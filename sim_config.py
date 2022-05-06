class SimConfig():

    def __init__(self, courant, t_final):

        """
        class SimConfig : stores several properties
        
        PROPERTIES
        ----------
        - courant : courant constant
        - tfinal : simulation time [s]
        
        """

        self.courant = courant
        self.tfinal = t_final