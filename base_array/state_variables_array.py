from .base_array import *
import cantera as ct
import numpy as np
import pickle

class StateVariablesArray(BaseArray):
    '''Variable array for state variables'''

    def __init__(self, parent, var=None):
        super().__init__(parent, var)
        self.coef_a = np.zeros(self.num_grid)
        self.coef_b = np.zeros(self.num_grid)
        self.coef_c = np.zeros(self.num_grid)
        self.coef_d = np.zeros(self.num_grid)
        
    def solve_TDMA(self):
        '''Solve linear equations using TDMA method'''
        # preparation    
        P = np.zeros(self.num_grid)
        Q = np.zeros(self.num_grid)
        
        # step 1
        # a, b, c, d = calc_coef(0)
        P[0] = self.coef_b[0]/self.coef_a[0]
        Q[0] = self.coef_d[0]/self.coef_a[0]
        
        # step 2
        for i in range(1, self.num_grid):
            # a, b, c, d = calc_coef(i)
            P[i] = self.coef_b[i]/(self.coef_a[i] - self.coef_c[i]*P[i-1])
            Q[i] = (self.coef_d[i] + self.coef_c[i]*Q[i-1])/(self.coef_a[i] - self.coef_c[i]*P[i-1])
            
        # step 3
        self.variable_array[self.num_grid-1] = Q[self.num_grid-1]
        
        # step 4
        for i in range(self.num_grid-2, -1, -1):
            self.variable_array[i] = P[i]*self.variable_array[i+1] + Q[i]

    def solve(self):
        '''High level function to Solve variable array'''
        self.calc_coef()
        self.solve_TDMA()

class DensityArray(StateVariablesArray):
    '''Variable array for density'''
    
    def __init__(self, parent, var=None):
        super().__init__(parent, var)
        self.name = 'Density (g/cm3)'
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0
