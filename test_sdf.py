import os, sys
import math
import numpy as np
import torch

class TestSDF:
  def __init__(self, dim):
    self.dim = dim
    self.data = torch.FloatTensor(dim, dim, dim)
    self.grads = torch.ones(dim, dim, dim, 3)
    self.verts = torch.ones(dim, dim, dim, 12) * 0.5

    self.initialize()

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
