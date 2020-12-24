import os, sys
import glob
import platform
import math
import numpy as np
import torch

from cffi import FFI

ffi = FFI()

class FFIData:
  # Probably have to use torch.matmul. - Matrix Matrix multiply of Tensors
  #  If the first argument is 2-dimensional and the second argument is
  #  1-dimensional, the matrix-vector product is returned.
  def __init__(self, dim):
    self.dim = dim
    self.data = torch.FloatTensor(dim, dim, dim)
    self.grads = torch.ones(dim, dim, dim, 3)
    #self.verts = torch.rand(dim, dim, dim, 12)
    self.verts = torch.ones(dim, dim, dim, 12) * 0.5

    #FIXME: Make sure that you cleanup the memory allocated here
    #FIXME: Make sure that you cleanup the memory allocated here
    #FIXME: Make sure that you cleanup the memory allocated here
    self.ffi_data  = ffi.new(f"float[{dim*dim*dim}]")
    self.ffi_grads  = ffi.new(f"float[{dim*dim*dim*3}]")
    self.ffi_verts  = ffi.new(f"float[{dim*dim*dim*12}]")
    self.initialize()

  def copy_to_ffi(self):
    l = 0
    ig = 0
    iv = 0
    for i in range(self.dim):
      for j in range(self.dim):
        for k in range(self.dim):
          self.ffi_data[l] = self.data[k, j,i] # OpenGL orders arrays differently
          for ii in range(3):
            self.ffi_grads[ig] = self.grads[k, j,i, ii]
            ig = ig + 1
          for ii in range(12):
            self.ffi_verts[iv] = self.verts[k, j,i, ii]
            iv = iv + 1
          l = l + 1

  def calc_sphere(self, i, j, k):
    d2 = (self.dim-1)/2
    i = i - d2
    j = j - d2
    k = k - d2
    l = math.sqrt(i*i + j*j + k*k)
    d = l - d2
    return d
  
  def calc_torus(self, i, j, k):
    tube_radius = 0.15
    ring_radius = 0.25
    d2 = (self.dim-1)/2
    x = (i - d2)/self.dim
    y = (j - d2)/self.dim
    z = (k - d2)/self.dim
    qx = math.sqrt(x*x + z*z) - ring_radius
    lq = math.sqrt(qx*qx + y*y)
    return lq-tube_radius

  def initialize(self):
    eps = 0.00001
    # row versus column matrix. which is it?
    for i in range(self.dim):
      for j in range(self.dim):
        for k in range(self.dim):
          #d = self.calc_sphere(i,j,k)
          #gx = self.calc_sphere(i+eps,j,k) - d
          #gy = self.calc_sphere(i,j+eps,k) - d
          #gz = self.calc_sphere(i,j,k+eps) - d
          d = self.calc_torus(i,j,k)
          gx = self.calc_torus(i+eps,j,k) - d
          gy = self.calc_torus(i,j+eps,k) - d
          gz = self.calc_torus(i,j,k+eps) - d
          self.data[i,j,k] = d
          self.grads[i,j,k,0] = gx/eps
          self.grads[i,j,k,1] = gy/eps
          self.grads[i,j,k,2] = gz/eps

os.chdir(os.getcwd())

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

  ffi.cdef(cstr)

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

isosurf = ffi.dlopen(libname)
print(isosurf)
isosurf.test()

dim = 16 


#FIXME: Make sure that you cleanup the memory allocated here
#FIXME: Make sure that you cleanup the memory allocated here
#FIXME: Make sure that you cleanup the memory allocated here
ffi_tris  = ffi.new(f"int[{dim*dim*dim*3}]")
ffi_verts = ffi.new(f"float[{dim*dim*dim*12}]")
ffi_vert_count = ffi.new("int*")
ffi_tri_count = ffi.new("int*")

path = ffi.new("char[]","mc.obj".encode('ascii'))
isovalue = 0.0
sdf = FFIData(dim)
sdf.copy_to_ffi()
#isosurf.run_quadratic_mc(path, ffi.cast("int", dim), sdf.ffi_data, ffi.cast("float", isovalue))
isosurf.run_quadratic_mc(path, ffi.cast("int", dim), sdf.ffi_data, ffi.cast("float", isovalue), ffi_verts, ffi_vert_count, ffi_tris, ffi_tri_count)
print(f'FROM PYTHON: vert count {ffi_vert_count[0]}, tri count {ffi_tri_count[0]}')

np_verts = np.zeros((ffi_vert_count[0], 3), dtype=np.float32)
for i in range(ffi_vert_count[0]):
  np_verts[i,0] = ffi_verts[i*3 + 0]
  np_verts[i,1] = ffi_verts[i*3 + 1]
  np_verts[i,2] = ffi_verts[i*3 + 2]

np_tris = np.zeros((ffi_tri_count[0], 3), dtype=np.int32)
for i in range(ffi_tri_count[0]):
  np_tris[i,0] = ffi_tris[i*3 + 0]
  np_tris[i,1] = ffi_tris[i*3 + 1]
  np_tris[i,2] = ffi_tris[i*3 + 2]

with open("verify_mc.obj", "w") as f:
  f.write("# OBJ file\n")
  for i in range(ffi_vert_count[0]):
    f.write(f"v {np_verts[i,0]:3.4f} {np_verts[i,1]:3.4f} {np_verts[i,2]:3.4f}\n")
  for i in range(ffi_tri_count[0]):
    f.write(f"f {np_tris[i,0]+1:d} {np_tris[i,1]+1:d} {np_tris[i,2]+1:d}\n")

