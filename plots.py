import numpy as np
import matplotlib.pyplot as plt

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
            plt.show()
            
        if it >= 2:
            self.data_point_x = np.append(self.data_point_x, point[0])
            self.data_point_y = np.append(self.data_point_y, point[1])
            self.error_array = np.append(self.error_array, error_value)
            plt.plot(self.data_point_x, self.data_point_y, 'bo-', self.data_point_x, self.error_array, 'go-')
            # Addition of tittle, labels, ...
            plt.grid()
            plt.title('Representacion de dos funciones')
            plt.xlabel('Iterations')
            plt.ylabel(labelplot)
            plt.show()
            
        return self.data_point_x, self.data_point_y, it

            
        