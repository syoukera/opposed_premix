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

        # calculation time
        self.start_time = 0.0
        self.end_time   = 1e-4
        self.dt         = 1e-6

        # calclation step
        self.total_step = int(self.end_time/self.dt)

        # Path to csv file exported using CHEMKIN-PRO
        self.path_ck = 'data/export_OpposedDiffusioon_CH4_GRI.csv'
        self.df_ck = pd.read_csv(self.path_ck)

        # Declear array of state variables
        self.R_array = base_array.DensityArray(parent=self)
        self.P_array = base_array.PressureArray(parent=self)
        self.T_array = base_array.TemperatureArray(parent=self)
        self.V_array = base_array.AxialVelocityArray(parent=self)
        self.G_array = base_array.RadialVelocityArray(parent=self)

        # Declear array of parameters
        self.mu_array = base_array.ParameterArray(parent=self, name='Mixture_viscosity (g/cm-sec)')
        self.cp_array = base_array.ParameterArray(parent=self, name='Specific_heat_Cp (erg/g-K)')
        self.lm_array = base_array.ParameterArray(parent=self, name='Mixture_thermal_conductivity (erg/cm-K-sec)')
        self.TPG_array = base_array.ParameterArray(parent=self, name='Pressure_gradient (g/cm3-s2)')

    def solve(self):
        '''Solve probrems using TDMA and SIMPLE method'''
        
        # Inititalise calclation
        self.time = self.start_time
        # Initialize function
        self.interporate_arrays()
        self.average_arrays()
        self.initialize_arrays()

        for n_step in range(self.total_step):
            self.time = n_step*self.dt
            self.time_step()

    def time_step(self):
        '''Progress a time step'''

        self.G_array.solve()

    def initialize_arrays(self):
        '''Initialize solved arrays'''
        
        self.G_array.initialize()

    def interporate_arrays(self):
        '''Interpolate required arrays'''

        # State variables array
        self.V_array.interpolate()
        self.R_array.interpolate()

        # Parameter array
        self.mu_array.interpolate()
        self.TPG_array.interpolate()

    def average_arrays(self):
        '''Average operation for required arrays'''

        self.V_array.average_variables()
        self.mu_array.average_variables()

