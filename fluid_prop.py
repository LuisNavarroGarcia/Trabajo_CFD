k = 0.025 # Thermal Conductivity [W/(m·K)]
cv = 717.5  # Specific Heat [J/(kg·K)]
rho = 1.225 # Density [kg/m^3]

class FluidProp():

    def __init__(self, rho = rho, cv = cv, k = k):
        self.rho = rho
        self.cv = cv
        self.k = k