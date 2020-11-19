import base_array

class BaseSolution():
    '''Base solution class for simulation'''

    def __init__(self, name=None):
        self.name = name

        self.R_array = base_array.DensityArray(parent=self, n=5)
        self.P_array = base_array.PressureArray(parent=self, n=5)
        self.T_array = base_array.TemperatureArray(parent=self, n=5)
        self.V_array = base_array.VelocityArray(parent=self, n=5)

    def solve(self):
        '''Solve probrems using TDMA and SIMPLE method'''


    def get_properties(self):
        '''Get state variables of simulation'''
