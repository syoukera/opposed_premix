from .base_array import *
import numpy as np
import scipy.interpolate as interp

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

class TemperatureArray(StateVariablesArray):
    '''Variable array for temperature'''

    def __init__(self, parent, var=None):
        super().__init__(parent, var)
        self.kind = 'Temperature (K)'
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0

class DensityArray(StateVariablesArray):
    '''Variable array for density'''
    
    def __init__(self, parent, var=None):
        super().__init__(parent, var)
        self.kind = 'Density (g/cm3)'
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0

class AxialVelocityArray(StateVariablesArray):
    '''Variable array for axial velocity'''

    def __init__(self, parent, var=None):
        super().__init__(parent, var)
        self.kind = 'Axial_velocity (cm/sec)'
    
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
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0

class RadialVelocityArray(StateVariablesArray):
    '''Variable array for radial velocity'''

    def __init__(self, parent, var=None):
        super().__init__(parent, var)
        self.kind = 'Radial_velocity (1/s)'

    def initialize(self):
        '''Initialize variable array for simulation'''
        center = int(self.num_grid/2) + 1
        self.variable_array = np.zeros(self.num_grid)
        self.variable_array[:center] = np.linspace(0, 80, center)
        self.variable_array[-center:] = np.linspace(80, 0, center)
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''

        for p in range(self.num_grid):
        
            # Upper boundary condition
            if p == 0:
                self.coef_a = 1.0
                self.coef_b = 0.0
                self.coef_c = 0.0
                self.coef_d = 0.0
                continue
                
            # Lower boundary condition
            if p == self.num_grid-1:
                self.coef_a = 1.0
                self.coef_b = 0.0
                self.coef_c = 0.0
                self.coef_d = 0.0
                continue
            
            # Inner points of numerical grid    
            C_p = mu_s[p]*dt/(R[p]*dy**2)    
            C_w = mu_s[p-1]*dt/(R[p]*dy**2)

            self.coef_a = 1 + C_p + C_w + self.variable_array[p]*dt
            self.coef_b = C_p
            self.coef_c = C_w
            
            # term v*dG/dy
            method = 'upwind'
            if method == 'center':
                self.coef_D = V_l[p]*dt/(2*dy)
                self.coef_b += - D
                self.coef_c += + D
            elif method == 'upwind':
                D = V_l[p]*dt/dy
                if V_l[p] > 0:
                    self.coef_a += D
                    self.coef_c += D
                else:
                    self.coef_a += - D
                    self.coef_b += - D
            else:
                sys.exit('Unknow method is choosen')        
                
            self.coef_d = self.variable_array[p] - TPG*dt/R[p]
    
        # self.coef_a = np.ones(self.num_grid) * 1.0
        # self.coef_b = np.ones(self.num_grid) * 0.0
        # self.coef_c = np.ones(self.num_grid) * 0.0
        # self.coef_d = np.ones(self.num_grid) * 1.0

class PressureArray(StateVariablesArray):
    '''Variable array for pressure'''

    def __init__(self, parent, var=None):
        super().__init__(parent, var)
        self.kind = 'Pressure (atm)'
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0
