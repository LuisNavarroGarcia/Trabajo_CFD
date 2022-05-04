from matplotlib import colors, lines
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

class ErrorPlot():
    
    def __init__(self): 
        self.data_point_it_for_max = np.array([])
        self.data_point_it_for_mean = np.array([])
        self.data_point_y_max = np.array([])
        self.data_point_y_mean = np.array([])

        plt.ion()
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        
        
    def __call__(self, point, error_value, it, plot_type, act_plot):
        
        if plot_type == 0: #Maximum type error
            self.data_point_it_for_max = np.append(self.data_point_it_for_max, it)
            self.data_point_y_max = np.append(self.data_point_y_max, point)

            plt.axhline(y = error_value, color = 'g', linestyle = 'dashed', label = 'Max error threshold')
            
            self.line1 = plt.plot(self.data_point_it_for_max, self.data_point_y_max, 'b-', label = 'Max error evolution') 
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()

            # Addition of tittle, labels, ...
            plt.title('Convergence of the error')
            plt.xlabel('Iterations')
            plt.ylabel('Error [-]')
            plt.plot()
            

        if plot_type == 1: #Mean type error
            self.data_point_it_for_mean = np.append(self.data_point_it_for_mean, it)
            self.data_point_y_mean = np.append(self.data_point_y_mean, point)

            plt.axhline(y = error_value, color = 'y', linestyle = 'dashed', label = 'Mean error threshold')
            
            self.line2, = plt.plot(self.data_point_it_for_mean, self.data_point_y_mean, 'r-', label = 'Mean error evolution') 
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()

            # Addition of tittle, labels, ...
            plt.grid()
            plt.title('Convergence of the error')
            plt.xlabel('Iterations')
            plt.ylabel('Error [-]')
            plt.plot()
               
        if act_plot[1] == 1 and act_plot[2] == 1:
            plt.legend(['Max error threshold', 'Max error evolution', 'Mean error threshold', 'Mean error evolution'])
        elif act_plot[1] == 1 and act_plot[2] == 0:
            plt.legend(['Max error threshold', 'Max error evolution'])
        elif act_plot[1] == 0 and act_plot[2] == 1:
            plt.legend(['Mean error threshold', 'Mean error evolution'])
        
        return self.data_point_it_for_max, self.data_point_it_for_mean, self.data_point_y_max, self.data_point_y_mean, it




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
    x_nodes = mesh.Rn[:, 0] * 1000; #multiply the x position of the nodes by 1000
    y_nodes = mesh.Rn[:, 1] * 1000; #multiply the y position of the nodes by 1000
    asp_ratio = (max(x_nodes) - min(x_nodes)) / (max(y_nodes) - min(y_nodes))

    ##A recostruction of the nodes is applied
    x_size = np.size(x_nodes)
    rec = np.zeros(x_size)
    
    for i in range(x_size):
        n = np.where(conectivity_list == i+1)[0] 
        summatory = 0
        cont = 0
        
        for j in range(len(n)):
            summatory = summatory + w[n[j]]
            cont = cont + 1
        try:
            rec[i] = summatory / cont
        except Exception: 
            rec[i] = summatory
            
    ##Interpolation of the mesh solution 
    newpoints = 501
    xv = np.linspace(min(x_nodes), max(x_nodes), newpoints)
    yv = np.linspace(min(y_nodes), max(y_nodes), newpoints)
    X, Y = np.meshgrid(xv, yv)
    
    #Selection of the interpolation method
    interp_method = ['nearest', 'linear', 'cubic'] 
    
    Z = griddata((x_nodes, y_nodes), np.transpose(rec), (X, Y), method = interp_method[num_interp])
    
    #Selection of the colour map type
    colourmap = ['binary','viridis', 'inferno', 'plasma', 'magma', 'cividis']
    
    #Selection of the lettering 
    plt.rcParams.update({'font.family':'fantasy'})

    #Figure creation
    levels_contour = 500
    plt.rcdefaults()
    fig,ax=plt.subplots(1, 1, figsize=(asp_ratio*5*20, 5))

    cp = ax.contourf(X, Y, Z, levels_contour, cmap=colourmap[num_map])
    fig.colorbar(cp) # Add a colorbar to a plot
    ax.set_title('Temperature Contour Plot')
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    plt.show()