import numpy as np
from boundary_conditions import neumann_convection, dirichlet_convection

def conv_upwind_order1(mesh = None, fluid_prop = None, bc = None, u = None, w = None, t = None):

    N = len(mesh.cells)
    C_matrix = np.zeros((N, N)) # matriz NxN
    bc_vec = np.zeros((N,1)) # matriz Nx1
    BC = np.zeros((N,N)) # matriz NxN


    for i in range(N):
        V_i = mesh.V[i] # volumen celda i
        Centroids = [mesh.Rc[i, 0], mesh.Rc[i, 1]] # Matriz 1x2 con las ...
        #coordenadas del centroide de la celda #se sobrescribe cada vez
        v = u([Centroids[0]], [Centroids[1]], t) #campo de velocidad dependiente del centroide y el tiempo
        for j in range(3):
            k = mesh.neighbours[i, j] # indice de la celda vecina de la cara j en la celda i
            A_ij = mesh.areas[i, j] # area cara j, celda i
            n_ij = [mesh.normals[i, j, 0], mesh.normals[i, j, 1]] #normal exterior de la cara j de la celda i
            vn = np.dot(v, n_ij)
            conv = -A_ij * vn * fluid_prop.rho * fluid_prop.cv #rho*cv??????????????????????  
            
            if k >= 0:  #si al lado hay celda
                if vn < 0:     #si la celda vecina esta aguas arriba
                    V_k = mesh.V[k]
                    C_matrix[i, k] = conv/V_i
                    C_matrix[k, k] = C_matrix[k, k]-conv/V_k   
            else:  #si al lado hay frontera
                centroid_face = [mesh.faces[i, j, 0], mesh.faces[i, j, 1]]    #centroides de las caras
                if bc.bc_type[np.absolute(k)-1] == 1:
                    [BC_i, bc_i, _, _] = neumann_convection(
                        bc_handler = bc.bc_handler[np.absolute(k)-1] , BC = BC,
                    bc = bc_vec , conv = conv, iteration = i, x = centroid_face[0], y = centroid_face[1], t = t,
                    mesh = mesh, fluid_prop = fluid_prop
                    )
                else:
                    [BC_i, bc_i, _, _] = dirichlet_convection(
                        bc_handler = bc.bc_handler[np.absolute(k)-1] , BC = BC,
                    bc = bc_vec , conv = conv, iteration = i, x = centroid_face[0], y = centroid_face[1], t = t,
                    mesh = mesh, fluid_prop = fluid_prop
                    )
    return C_matrix, BC_i, bc_i
                
                
                
    #            frontera = np.absolute(k) #frontera??
    #            r_cara = [Mesh.faces[i,j,1], Mesh.faces[i,j,2]]
    #            conv = -A_ij*vn/V_i #*rho*cv????????????????
    #            [BC_i, bc_i]=cdc.tipo{frontera,2}(cdc.handler{frontera]},...
    #                BC,bc,conv,i,r_cara(1),r_cara(2),t)
    #            if vn<0
    #                bc[i]=bc_i
    #            else
    #                BC[i,i]=BC_i

                