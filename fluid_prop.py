class FluidProp():
    '''
    class FluidProp : stores some properties of the fluid used to simulate
        - k : thermal conductivity [W/m*K]
        - rho : density [kg/m^3]
        - cv : specific heat [J/kg*K]
    '''
    def __init__(self, rho = None, cv = None, k = None):
        self.rho = rho
        self.cv = cv
        self.k = k