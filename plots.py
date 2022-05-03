import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

class ErrorPlot():
    
    def __init__(self): 
        self.data_point_it_for_max = np.array([])
        self.data_point_it_for_mean = np.array([])
        self.data_point_y_max = np.array([])
        self.data_point_y_mean = np.array([])
        self.error_array_max = np.array([])
        self.error_array_mean = np.array([])
        
        
    def __call__(self, point, error_value, it, plot_type):
        
        if plot_type == 0: #Maximum type error
            self.data_point_it_for_max = np.append(self.data_point_it_for_max, it)
            self.data_point_y_max = np.append(self.data_point_y_max, point)
            self.error_array_max = np.append(self.error_array_max, error_value)
            
            plt.plot(self.data_point_it_for_max, self.data_point_y_max, 'bo')
            plt.plot(self.data_point_it_for_max, self.error_array_max, 'go')
            
            # Addition of tittle, labels, ...
            plt.title('Convergence of the error')
            plt.xlabel('Iterations')
            plt.ylabel('Error [-]')
            plt.legend(['Max error evolution', 'Max error threshold'])
            plt.plot()
            

        if plot_type == 1: #Mean type error
            self.data_point_it_for_mean = np.append(self.data_point_it_for_mean, it)
            self.data_point_y_mean = np.append(self.data_point_y_mean, point)
            self.error_array_mean = np.append(self.error_array_mean, error_value)

           
            plt.plot(self.data_point_it_for_mean, self.data_point_y_mean, 'ro')
            plt.plot(self.data_point_it_for_mean, self.error_array_mean, 'yo')
            
            # Addition of tittle, labels, ...
            plt.grid()
            plt.title('Convergence of the error')
            plt.xlabel('Iterations')
            plt.ylabel('Error [-]')
            plt.legend(['Mean error evolution', 'Mean error threshold'])
            plt.plot()

        plt.legend(['Max error evolution', 'Max error threshold', 'Mean error evolution', 'Mean error threshold'])
            
        return self.data_point_it_for_max, self.data_point_it_for_mean, self.data_point_y_max, self.data_point_y_mean, self.error_array_max, self.error_array_mean, it

def contour_plot(w, mesh, num_map, num_interp):
    """
    function contour_plot : plots the solution in a contour plot 
    The solution is interpolated in the nodes of the mesh. 
    
    PARAMETERS
    ----------
    
    w : state value evaluated in time t. It is a VECTOR COLUM.
    
    mesh : struc with all the mesh porperties.
    
    num_map : variable to choose the colours of the contour plot
        (0) indicates 'binary'
        (1) indicates 'viridis'
        (2) indicates 'inferno'
        (3) indicates 'plasma'
        (4) indicates 'magma'
        (5) indicates 'cviridis'

    num_interp : variable to choose the interpolation method used for the data
        (0) indicates 'nearest' interpolation
        (1) indicates 'linear' interpolation
        (2) indicates 'cubic' interpolation
        It is recommended to set the cubic interpolation, as it is the most precise one
        


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
    
    #Selection of the interpolation method
    interp_method = ['nearest', 'linear', 'cubic'] 
    
    Z = griddata((x_nodes, y_nodes), np.transpose(rec), (X, Y), method = interp_method[num_interp])
    #print('Z', '\n', Z)
    #print('Y size:', '\n', np.size(Z))
    
    #Selection of the colour map type
    colourmap = ['binary','viridis', 'inferno', 'plasma', 'magma', 'cividis']
    #print('Z', '\n', Z)
    #print('Y size:', '\n', np.size(Z))
    
    #Selection of the lettering 
    plt.rcParams.update({'font.family':'fantasy'})

    #Figure creation
    levels_contour = 500
    plt.rcdefaults()
    fig,ax=plt.subplots(1, 1, figsize=(asp_ratio*5*20, 5))

    cp = ax.contourf(X, Y, Z, levels_contour, cmap=colourmap[num_map])
    fig.colorbar(cp) # Add a colorbar to a plot
    ax.set_title('Temperature Contour Plot')
    ax.set_xlabel('x (cm)')
    ax.set_ylabel('y (cm)')
    plt.show()