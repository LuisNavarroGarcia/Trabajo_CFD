from matplotlib import colors, lines
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

class ErrorPlot():
    
    """class ErrorPlot : class to plot the evolution of the mean or maximum 
    error with the iterations, to check their convergence. 
    It will only appear if the activation plot option is activated. 

    PARAMETERS
    ----------

    point : variable indicating the error calculated. It will be added as 
    a new point in the plot, whether as a mean or as a maximum error, 
    depending on the plot_type variable. 
    
    error_value : it indicates the threshold for the error. It can be refered to 
    the mean error or to the maximum error, depending on the plot_type variable. 
     
    it : variable containing the iteration that will be added to the plot. 
    
    plot_type : variable to define if the point to be addend in the plot 
    is refered to the maximum error evolution, or the mean eror one. 
    (0) If the point added refers to maximum error
    (1) If the point added refers to mean error

    
    act_plot : Variable defining if the plot is activated (and will 
    appear during the simulation) or not. 
    This option is set before runing the simulation by the user. 
    (0) If the plot is not activated. 
    (1) If the plot is not activated. 

    INTERNAL PARAMETERS
    -------------------

    self.fig : variable to create the figre. 
    It gathers all the figure plot settings

    self.ax : variable generated to define the axis of the plot 
    
    OUTPUTS
    ----------

    self.data_point_it_for_max : vector array. It gathers each time the iteration 
    (only the multiple ones of the sampling frequency, for the cases when maximum error 
    is activated) for which a new point in the figure will appear. 

    self.data_point_it_for_mean : vector array. It gathers each time the iteration 
    (only the multiple ones of the sampling frequency, for the cases when mean error 
    is activated) for which a new point in the figure will appear. 

    self.data_point_y_max : it gathers for each iteration the corresponfing 
    mean error, to be drawn if the plot. 

    self.data_point_y_mean : it gathers for each iteration the corresponfing 
    maximum error, to be drawn if the plot. 
    """

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




def contour_plot(w, mesh, colourmap, interp_method):
    """
    function contour_plot : plots the solution, interpolating the 
    results from each cell to get a smooth temperature contour.  
    The solution is interpolated using the final result (last iteration of
    the state vector, w).
    
    PARAMETERS
    ----------
    
    w : state value evaluated in time t. It is a VECTOR COLUM.
    
    mesh : struc with all the mesh porperties.
    
    colourmap : varianles used to choose the colours of the contour plot
        Some examples to set in this variable: 'binary', 'viridis'
        'inferno', 'plasma', 'magma', 'cviridis', etc.
        
    interp_method : variable used to choose the interpolation method 
        of the data from w. 
        Some examples to set in this variable: 'nearest' , 'linear', 'cubic'.
        It is recommended to set the cubic interpolation, as it is the most precise one.  


    INTERNAL PARAMETERS
    ----------
    
    x_nodes : COLUM VECTOR. It gathers the x component of the mesh nodes
    
    y_nodes : COLUMN VECTOR. It gathers the the y component of the mesh nodes
    
    asp_ratio : aspect ratio of the plot. It is used to help in the scaling of the 
    contour figure. Variation of the width over the hight of the figure. 
        
    rec : reconstructed temperatures in the nodes. 
        It is a COLUMN VECTOR. 
        
    newpoitns : new points created when scaling the mesh 
        It is a SCALAR VALUE. 
        
    Z : temperature matrix of the mesh obtained through the iterpolation. 
        It is a N*2 MATRIX, corresponding to the 2 dimensions of the domain
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
    
    Z = griddata((x_nodes, y_nodes), np.transpose(rec), (X, Y), method = interp_method)
    
    #Selection of the lettering 
    plt.rcParams.update({'font.family':'fantasy'})

    #Figure creation
    levels_contour = 500
    plt.rcdefaults()
    fig,ax=plt.subplots(1, 1, figsize=(asp_ratio*5*20, 5))

    cp = ax.contourf(X, Y, Z, levels_contour, cmap=colourmap)
    fig.colorbar(cp) # Add a colorbar to a plot
    ax.set_title('Temperature Contour Plot')
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    fig.tight_layout()
    plt.show()