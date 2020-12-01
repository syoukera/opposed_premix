from .base_array import *
from .state_variables_array import *
import cantera as ct
import numpy as np

class AxialVelocityArray(StateVariablesArray):
    '''Variable array for axial velocity'''

    def __init__(self, parent, var=None):
        super().__init__(parent, var)
        self.name = 'Axial_velocity (cm/sec)'

    def initialize(self):
        '''Initialize variable array for simulation'''
        self.variable_array = np.linspace(100, -100, self.num_grid)
    
    def average_variables(self):
        '''
        Average variables at intermediate grid point
        for velocity in stagard grid
        'u' means Upper case
        '''

        self.variable_array_u = self.variable_array.copy()
        self.variable_array_u[1:] += self.variable_array[:-1].copy()
        self.variable_array_u[1:] /= 2
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''

        # Read variables from parent soluction
        mu    = self.parent_solution.mu.variable_array 
        mu_s  = self.parent_solution.mu.variable_array_s 
        R_s   = self.parent_solution.R.variable_array_s
        V     = self.parent_solution.V.variable_array
        V_old = self.parent_solution.V_old
        P     = self.parent_solution.P.variable_array
        G     = self.parent_solution.G.variable_array

        for p in range(self.num_grid):
        
            # Upper boundary condition
            if p == 0:
                self.coef_a[p] = 1.0
                self.coef_b[p] = 0.0
                self.coef_c[p] = 0.0
                self.coef_d[p] = 100
                continue
                
            # Lower boundary condition
            if p == self.num_grid-1:
                self.coef_a[p] = 1.0
                self.coef_b[p] = 0.0
                self.coef_c[p] = 0.0
                self.coef_d[p] = -100
                continue

            
            # Inner points of numerical grid    
            C_E = 4*mu[p+1]/(3*R_s[p]*self.dy**2)
            C_P = 4*mu[p]/(3*R_s[p]*self.dy**2)

            self.coef_a[p] = 1/self.dt + C_E + C_P + np.abs(V[p]/self.dy)
            self.coef_b[p] = C_E + np.max(0, -V[p]/self.dy)
            self.coef_c[p] = C_P + np.max(0, +V[p]/self.dy)
            self.coef_d[p] = V_old[p]/self.dt - (P[p+1] - P[p])/(R_s[p]/self.dy) \
                           - C_E*G[p+1]*self.dy + C_P*G[p]*self.dy \
                           + 2*mu_s[p]*(G[p+1] - G[p])/(R_s[p]*self.dy)
