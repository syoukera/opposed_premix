from .base_array import *
from .state_variables_array import *
import cantera as ct
import numpy as np

class RadialVelocityArray(StateVariablesArray):
    '''Variable array for radial velocity'''

    def __init__(self, parent, var=None):
        super().__init__(parent, var)
        self.name = 'Radial_velocity (1/s)'

    def initialize(self):
        '''Initialize variable array for simulation'''
        center = int(self.num_grid/2) + 1
        self.variable_array = np.zeros(self.num_grid)
        self.variable_array[:center] = np.linspace(0, 80, center)
        self.variable_array[-center:] = np.linspace(80, 0, center)
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''

        # Read variables from parent soluction
        mu_s = self.parent_solution.mu.variable_array_s
        R    = self.parent_solution.R.variable_array
        V_u  = self.parent_solution.V.variable_array_u
        TPG  = self.parent_solution.TPG.variable_array[0]


        for p in range(self.num_grid):
        
            # Upper boundary condition
            if p == 0:
                self.coef_a[p] = 1.0
                self.coef_b[p] = 0.0
                self.coef_c[p] = 0.0
                self.coef_d[p] = 0.0
                continue
                
            # Lower boundary condition
            if p == self.num_grid-1:
                self.coef_a[p] = 1.0
                self.coef_b[p] = 0.0
                self.coef_c[p] = 0.0
                self.coef_d[p] = 0.0
                continue
            
            # Inner points of numerical grid    
            C_p = mu_s[p]*self.dt/(R[p]*self.dy**2)    
            C_w = mu_s[p-1]*self.dt/(R[p]*self.dy**2)

            self.coef_a[p] = 1 + C_p + C_w + self.variable_array[p]*self.dt
            self.coef_b[p] = C_p
            self.coef_c[p] = C_w
            
            # term v*dG/dy
            method = 'upwind'
            if method == 'center':
                D = V_u[p]*self.dt/(2*self.dy)
                self.coef_b[p] += - D
                self.coef_c[p] += + D
            elif method == 'upwind':
                D = V_u[p]*self.dt/self.dy
                if V_u[p] > 0:
                    self.coef_a[p] += D
                    self.coef_c[p] += D
                else:
                    self.coef_a[p] += - D
                    self.coef_b[p] += - D
            else:
                sys.exit('Unknow method is choosen')        
                
            self.coef_d[p] = self.variable_array[p] - TPG*self.dt/R[p]