import numpy as np
from fluid_prop import rho, cv
from boundary_conditions import neumann_convection, dirichlet_convection

def conv_upwind_order1(mesh, cdc_tipo, u, t):

    N = len(mesh.cells)
    C_matrix = np.zeros((N, N)) # matriz NxN
    bc = np.zeros((N,1)) # matriz Nx1
    BC = np.zeros((N,N)) # matriz NxN


    for i in range(N):

        V_i = mesh.V[i] # volumen celda i
        Centroids = [mesh.Rc[i, 1], mesh.Rc[i, 2]] # Matriz 1x2 con las ...
        #coordenadas del centroide de la celda #se sobrescribe cada vez
        v = u[Centroids[1], Centroids [2], t] #campo de velocidad dependiente del centroide y el tiempo

        for j in range(3):
            k = mesh.neighbours[i, j] # indice de la celda vecina de la cara j en la celda i
            A_ij = mesh.areas[i, j] # area cara j, celda i
            n_ij = [mesh.normals[i, j, 1], mesh.normals[i, j, 2]] #normal exterior de la cara j de la celda i
            vn = np.dot(v, n_ij)
            conv=-A_ij * vn * rho * cv #rho*cv??????????????????????   

            if k > 0:  #si al lado hay celda
                if vn < 0:
                         #si la celda vecina esta aguas arriba
                    V_k = mesh.V[k]
                    C_matrix[i,k] = conv/V_i
                    C_matrix[k,k] = C_matrix[k,k] - conv/V_k  

            else:  #si al lado hay frontera
                if cdc_tipo == 1:
                    [BC_i, bc_i] = neumann_convection()
                else:
                    [BC_i, bc_i] = dirichlet_convection()
                
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

            