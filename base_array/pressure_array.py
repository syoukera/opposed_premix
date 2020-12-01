from .base_array import *
from .state_variables_array import *
import cantera as ct
import numpy as np
import scipy.interpolate as interp

class PressureArray(StateVariablesArray):
    '''Variable array for pressure'''

    def __init__(self, parent, var=None):
        super().__init__(parent, var)
        self.name = 'Pressure (g/cm/s2)'

    def initialize(self):
        '''Initialize variable array for simulation'''
        self.variable_array = 1013250

    def interpolate(self):
        '''
        Interpolate and assign variables from other value arrays
        Converted to cgs-unit [g/cm/s2] from [atm]
        '''

        df_ck = self.parent_solution.df_ck

        dis = df_ck['Distance (cm)'].to_numpy()
        phi = df_ck['Pressure (atm)'].to_numpy()
        phi_cgs = phi*1013250

        f = interp.interp1d(dis, phi_cgs, kind="cubic")
        self.variable_array = f(self.y)
    

    def calc_coef(self):
        '''Calculate coefficients for TDMA'''

        # Read variables from parent soluction
        R = self.parent_solution.R.variable_array
        R_s = self.parent_solution.R.variable_array_s
        R_old = self.parent_solution.R_old
        G = self.parent_solution.G.variable_array
        V = self.parent_solution.V.variable_array

        '''Calculation of d'''
        coef_a_V = self.parent_solution.V.coef_a
        d = 1/np.multiply(R_s, coef_a_V)*self.dy

        for p in range(self.num_grid):
        
            # Upper boundary condition
            if p == 0:
                self.coef_a[p] = 1.0
                self.coef_b[p] = 0.0
                self.coef_c[p] = 0.0
                self.coef_d[p] = 1013250
                continue
                
            # Lower boundary condition
            if p == self.num_grid-1:
                self.coef_a[p] = 1.0
                self.coef_b[p] = 0.0
                self.coef_c[p] = 0.0
                self.coef_d[p] = 1013250
                continue
            
            # Inner points of numerical grid    
            self.coef_a[p] = R_s[p]*d[p]/self.dy + R_s[p-1]*d[p-1]/self.dy
            self.coef_b[p] = R_s[p]*d[p]/self.dy
            self.coef_c[p] = R_s[p-1]*d[p-1]/self.dy
            self.coef_d[p] = - (R[p] - R_old[p])/self.dy - 2*R[p]*G[p] \
                             - (R_s[p]*V[p] - R_s[p-1]*V[p-1])/self.dy