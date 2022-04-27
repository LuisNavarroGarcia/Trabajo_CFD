import numpy as np
from scipy.sparse import csr_matrix


def spatial_discretization(malla, fluid_prop, cdc, u, w, t , type_storage, diffusion_integrator, convection_integrator):
        
  A_dif, BC_dif, bc_dif = diffusion_integrator(malla, fluid_prop, cdc, u, w, t)
  A_conv, BC_conv, bc_conv = convection_integrator(malla, fluid_prop, cdc, u, w, t)

  if type_storage == 1:
      A = csr_matrix(A_dif + BC_dif + bc_dif)
      b = csr_matrix(bc_dif + bc_conv)

  else:
      A = A_dif + BC_dif + bc_dif
      b = bc_dif + bc_conv

  return A, b