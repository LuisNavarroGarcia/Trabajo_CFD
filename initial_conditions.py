import numpy as np

class InitCond():

    def __init__(self, t0 = 0, u = None):

        """
        class InitCond : stores the initial conditions of the problem.
        
        PROPERTIES
        ----------
        - t0 : initial time [s]
        - T : initial temperature array (as a function of position)
        - u : initial velocity array
        
        """

        self.t0 = t0
        self.T = lambda x, y : np.full((1, len(x)), 500)
        try:
            self.u = lambda x, y : u(x, y, self.t0)
        except AttributeError:
            print('No Velocity Field entered')
            raise