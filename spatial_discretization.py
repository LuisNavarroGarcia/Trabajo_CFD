import numpy as np
from scipy.sparse import csr_matrix


def spatial_discretization(mesh = None, fluid_prop = None, bc = None, u = None, w = None,
 t = None, type_storage = 1, diffusion_integrator = None, convection_integrator = None):
        
  A_dif, BC_dif, bc_dif = diffusion_integrator(mesh = mesh, fluid_prop = fluid_prop, bc = bc, u = u, w = w, t = t)
  A_conv, BC_conv, bc_conv = convection_integrator(mesh = mesh, fluid_prop = fluid_prop, bc = bc, u = u, w = w, t = t)

  if type_storage == 1:
      A = csr_matrix(A_dif + BC_dif + bc_dif)
      b = csr_matrix(bc_dif + bc_conv)

  else:
      A = A_dif + BC_dif + bc_dif
      b = bc_dif + bc_conv

  return A, b