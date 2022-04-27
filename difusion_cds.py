import numpy as np
import fluid_prop

def difusion_cds(mesh, fluid_prop, boundary_cond, t):

     num_cells = len(mesh.cells)
     K = np.zeros(num_cells)
     BC = np.zeros(num_cells)
     bc = np.zeros(num_cells, 1)



     for i in range(num_cells):

          for j in range(3):
               n_ij = np.array([mesh.normals[i, j, 0], mesh.normals[i, j, 1]])
               A_ij = mesh.areas[i, j]
               V_i = mesh.V[i]

               if mesh.neighbours[i, j] > 0:
               
                    neighbour_index = mesh.neighbours[i, j]
                    r_i = np.array(mesh.Rc[i, 0]), mesh.Rc[i, 1]
                    r_j = np.array(mesh.Rc[neighbour_index, 0], mesh.Rc[neighbour_index, 1])
                    d_ij = r_i - r_j
                    d_norm = np.linalg.norm(d_ij)
                    dn = np.dot(d_ij, n_ij)
                    cond = -A_ij / V_i * fluid_prop.k * (dn / (d_norm * d_norm)) \
                                   / (fluid_prop.cv * fluid_prop.rho)
                    K[i, neighbour_index] = cond
                    K[i, i] = K[i, i] - cond

               else:

                    r_i = np.array([mesh.Rc[i, 0], mesh.Rc[i, 1]])
                    r_face = np.array([mesh.faces[i, j, 0], mesh.faces[i, j, 1]])
                    n_ij = np.array([mesh.normals[i, j, 0], mesh.normals([i, j, 1])])
                    d_ic = r_i - r_face
                    d_norm = np.linalg.norm(d_ic)
                    dn = np.dot(d_ic, n_ij)
                    frontier = np.absolute(mesh.neighbours[i, j])
                    cond = -A_ij / V_i * fluid_prop.k * (dn / (d_norm * d_norm)) \
                                   / (fluid_prop.cv * fluid_prop.rho)  
                    [BC_i, bc_i] = '''boundary_cond.type{frontier, 0}........'''

                    BC[i, i] = BC_i
                    bc[i] = bc_i

     return K, BC, bc
          


