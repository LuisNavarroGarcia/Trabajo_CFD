import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

class ErrorPlot():
    
    def __init__(self): 
            self.data_point_x = np.array([0])
            self.data_point_y = np.array([0])
            self.error_array = np.array([0])
        
        
    def __call__(self, point, error_value, labelplot, it):
        
        it = it + 1

        if it == 1:
            plt.plot(1, 1, 'bo-', 1, error_value, 'go-')
            self.data_point_x = np.array([1])
            self.data_point_y = np.array([1])
            self.error_array = np.array([error_value])
            # Addition of tittle, labels, ...
            plt.grid()
            plt.title('Representacion de dos funciones')
            plt.xlabel('Iterations')
            plt.ylabel(labelplot)
            plt.plot()
            
        if it >= 2:
            self.data_point_x = np.append(self.data_point_x, it)
            self.data_point_y = np.append(self.data_point_y, point)
            self.error_array = np.append(self.error_array, error_value)
            plt.plot(self.data_point_x, self.data_point_y, 'bo-', self.data_point_x, self.error_array, 'go-')
            # Addition of tittle, labels, ...
            plt.grid()
            plt.title('Representacion de dos funciones')
            plt.xlabel('Iterations')
            plt.ylabel(labelplot)
            plt.plot()
            
        return self.data_point_x, self.data_point_y, it


def contour_plot(w, mesh, num_map):
    """
    function contour_plot : plots the solution in a contour plot 
    The solution is interpolated in the nodes of the mesh. 
    
    PARAMETERS
    ----------
    
    w : state value evaluated in time t. It is a VECTOR COLUM.
    
    mesh : struc with all the mesh porperties.
    
    num_map : to define the colours
        1 indicates VIRIDIS
        2 ndicates INFERNO
        
    INTERNAL PARAMETERS
    ----------
    
    x_nodes : COLUM VECTOR with the x component of the mesh nodes
    
    y_nodes : COLUMN VECTOR with the y component of the mesh nodes
    
    asp_ratio : aspect ratio of the plot. In base of TAMA?? or the domain. 
        It is a SCALAR VALUE. 
        
    rec : reconstructed temperatures in the nodes, 
        It is a COLUMN VECTOR. 
        
    newpoitns : new points created when scaling the mesh 
        It is a SCALAR VALUE. 
        
    Z : temperature matrix of the mesh. 
        It is a N*2 MATRIX. 
        
    maps : available maps. Possible options: VERIDI and INFERNO.
    """
    
    ##The data is read 
    conectivity_list = mesh.cells
    #print('convectivity_list:','\n', conectivity_list)
    x_nodes = mesh.Rn[:, 0] * 1000; #multiply the x position of the nodes by 1000
    y_nodes = mesh.Rn[:, 1] * 1000; #multiply the y position of the nodes by 1000
    asp_ratio = (max(x_nodes) - min(x_nodes)) / (max(y_nodes) - min(y_nodes))
    #print('xnodes: ','\n', x_nodes)
    #print('y_nodes: ','\n', y_nodes)
    #print('asp_ratio: ','\n', asp_ratio)
    
    ##A recostruction of the nodes is applied
    x_size = np.size(x_nodes)
    #print('x_size:', '\n', x_size)
    rec = np.zeros(x_size)
    
    for i in range(x_size):
        n = np.where(conectivity_list == i)[0] 
        # print('n: ', n)
        summatory = 0
        cont = 0
        
        for j in range(len(n)):
            summatory = summatory + w[n[j]]
            cont = cont + 1
        try:
            rec[i] = summatory / cont
        except Exception: 
            rec[i] = summatory
            
    # print('cont: ','\n', cont)
    # print('rec: ','\n', rec)
    ##Interpolation of the mesh solution 
    newpoints = 501
    xv = np.linspace(min(x_nodes), max(x_nodes), newpoints)
    #print('xv', '\n', xv)
    yv = np.linspace(min(y_nodes), max(y_nodes), newpoints)
    #print('tv', '\n', yv)
    X, Y = np.meshgrid(xv, yv)
    #print('X', '\n', X)
    #print('Y', '\n', Y)
    #print('X size:', '\n', len(X), '\n', len(X[0]))
    #print('Y size:', '\n', len(Y), '\n', len(Y[0]))
    
    
    Z = griddata((x_nodes, y_nodes), np.transpose(rec), (X, Y), method = 'linear')
    #print('Z', '\n', Z)
    #print('Y size:', '\n', np.size(Z))
    
    #Selection of the colour map type
    colourmap = ['binary','viridis', 'inferno', 'plasma', 'magma', 'cividis']
    
    #Selection of the lettering 
    plt.rcParams.update({'font.family':'fantasy'})
    plt.rcdefaults()
    fig,ax=plt.subplots(1,1)
    cp = ax.contourf(X, Y, Z, cmap=colourmap[num_map])
    fig.colorbar(cp) # Add a colorbar to a plot
    ax.set_title('Filled Contours Plot')
    ax.set_xlabel('x (cm)')
    ax.set_ylabel('y (cm)')
    plt.show()

            
        