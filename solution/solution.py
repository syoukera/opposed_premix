import base_array

class BaseSolution():
    '''Base solution class for simulation'''

    def __init__(self):
        self.R_array = base_array.DensityArray(5)
        self.P_array = base_array.PressureArray(5)
        self.T_array = base_array.TemperatureArray(5)
        self.V_array = base_array.VelocityArray(5)

