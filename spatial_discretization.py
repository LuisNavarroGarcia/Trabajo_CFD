import numpy as np


def spatial_discretization(mesh = None, fluid_prop = None, bc = None, u = None, w = None,
 t = None, diffusion_integrator = None, convection_integrator = None):
        
    '''
   function spatial_discretization : Function that determines the spatial discretization using Finite Volumes. To do this, matrices obtained from diffusion and
   convection integrators are summed. A, BC and bc matrices are calculated, for convection and diffusion.
   Them, A matrices and BC matrices are summed on one hand, and bc matrices on the other hand.
   
   PARAMETERS
   ----------
   bc_type : vector that constains the type of boundary condition for each frontier.
   bc_handler : vector of values at each boundary 'n'.
   
   BC : NxN matrix with convection terms associated to boundary conditions
   
   bc : Nx1 vector with convection terms associated to boundary conditions
   x : position in x axis [m]
   y : position in y axis [m]
   t : time [s]
   w: state vector [T]
   mesh : class with all the mesh properties:
      - Rn : coordinates of nodes forming the mesh
      - Cells : index of nodes forming each cell
      - Volumes : area (2D mesh) of each cell
      - Neighbours : neighboring cells for each cell and/or boundaries
      - Normals : external normal vector for each face of each cell
      - Areas : length of each face of each cell
      - Faces : coordinates of the center of each face for each cell
      
   fluid_prop : class containing the properties of the fluid.
      - k : thermal conductivity [W/m*K]
      - rho : density [kg/m^3]
      - cv : specific heat [J/kg*K]
      
   OUTPUT
   ----------
   A: N x N matrix that multiplicates temperature vector in cells.
   b: N x 1 matrix with all independent terms. 
   '''

    A_dif, BC_dif, bc_dif = diffusion_integrator(mesh = mesh, fluid_prop = fluid_prop, bc = bc, u = u, w = w, t = t)
    A_conv, BC_conv, bc_conv = convection_integrator(mesh = mesh, fluid_prop = fluid_prop, bc = bc, u = u, w = w, t = t)

    A = A_dif + BC_dif + A_conv + BC_conv
    b = bc_dif + bc_conv

    return A, b