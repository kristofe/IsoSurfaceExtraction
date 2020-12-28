import os, sys
import glob
import platform
import math
import numpy as np
import torch

from cffi import FFI

class QuadraticMarchingCubes:
  def __init__(self):
    self.ffi = FFI()
    with open("Src/exported_routines.h") as header:
      header_str = header.read()
      cstr = ""
      ignore = False
      for line in header_str.splitlines():
        if(line.startswith("#if")):
          ignore = True

        if(ignore == False):
          cstr += line

        if(line.startswith("#end")):
          ignore = False

      self.ffi.cdef(cstr)

    correctWorkingDirectory = os.getcwd()
    libname_start = correctWorkingDirectory + "/build/libquadratic_iso"

    if(platform.system() == "Darwin"):
        if os.path.exists(libname_start + ".dylib"):
          libname = libname_start + ".dylib"
        else:
          libname = libname_start + "d.dylib"
    elif(platform.system() == "Windows"):
      libname = correctWorkingDirectory + "/build/Debug/quadratic_iso.dll"
    else:
      if os.path.exists(libname_start + ".so"):
        libname = libname_start + ".so"
      else:
        libname = libname_start + "d.so"

    self.isosurf = self.ffi.dlopen(libname)
    os.chdir(os.getcwd())

    print(self.isosurf)
  
  def run(self, isovalue, np_sdf_data, dim, path=None):
    #Allocate the maximum possible amount of memory used for buffers
    #passed into the C++ code.
    np_tris = np.zeros((dim*dim*dim,3), dtype=np.int32)
    np_verts = np.zeros((dim*dim*dim,3), dtype=np.float32)
    ffi_vert_count = self.ffi.new("int*")
    ffi_tri_count = self.ffi.new("int*")

    #Run the C++ code.
    self.isosurf.run_quadratic_mc(
                            self.ffi.cast("int", dim), 
                            self.ffi.cast("float*",np_sdf_data.ctypes.data),
                            self.ffi.cast("float", isovalue), 
                            self.ffi.cast("float*",np_verts.ctypes.data), ffi_vert_count, 
                            self.ffi.cast("int*",np_tris.ctypes.data), ffi_tri_count)

    #Trim off unused memory
    np_verts = np_verts[:ffi_vert_count[0],:]
    np_tris = np_tris[:ffi_tri_count[0],:]

    if path is not None:
      with open(path, "w") as f:
        f.write("# OBJ file\n")
        for i in range(ffi_vert_count[0]):
          f.write(f"v {np_verts[i,0]:3.4f} {np_verts[i,1]:3.4f} {np_verts[i,2]:3.4f}\n")
        for i in range(ffi_tri_count[0]):
          f.write(f"f {np_tris[i,0]+1:d} {np_tris[i,1]+1:d} {np_tris[i,2]+1:d}\n")

    return np_verts, np_tris



if __name__ == "__main__":
  import test_sdf
  dim = 32 
  sdf = test_sdf.TestSDF(dim)
  qmc = QuadraticMarchingCubes()
  qmc.run(0.0, sdf.data.numpy(), dim)