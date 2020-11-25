import numpy as np
import pandas as pd
import base_array

class BaseSolution():
    '''Base solution class for simulation'''

    def __init__(self, name=None):
        self.name = name

        self.length = 10 # cm
        self.num_grid = 2001

        # grid
        self.dy = self.length/(self.num_grid - 1)
        self.y = np.linspace(0, self.length, self.num_grid) # cm

        # Path to csv file exported using CHEMKIN-PRO
        self.path_ck = 'data/export_OpposedDiffusioon_CH4_GRI.csv'
        self.df_ck = pd.read_csv(self.path_ck)

        self.R_array = base_array.DensityArray(parent=self)
        self.P_array = base_array.PressureArray(parent=self)
        self.T_array = base_array.TemperatureArray(parent=self)
        self.V_array = base_array.AxialVelocityArray(parent=self)
        self.G_array = base_array.RadialVelocityArray(parent=self)

    def solve(self):
        '''Solve probrems using TDMA and SIMPLE method'''


    def get_properties(self):
        '''Get state variables of simulation'''
