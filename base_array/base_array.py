import numpy as np
import scipy.interpolate as interp
import sys


class BaseArray():
    '''Base class for variables array'''
    
    def __init__(self, parent, var=None):
        self.kind = 'Base'

        # Assign parent solution
        self.parent_solution = parent

        # Read parent parameters 
        self.y = self.parent_solution.y
        self.num_grid = self.parent_solution.num_grid

        # Define my own variables
        self.variable_array = np.zeros(self.num_grid)
        
        if var is None:
            pass
        else:
            if vals.shape == self.variable_array.shape:
                self.variable_array = var
            else:
                sys.exit('Exception: Thrown vals has different shape.')

    def average_variables(self):
        '''
        Average variables at intermediate grid point
        for variable in stagard grid
        's' means Small case
        '''

        self.variable_array_s = self.variable_array.copy()
        self.variable_array_s[:-1] += self.variable_array[1:].copy()
        self.variable_array_s[:-1] /= 2
    
    def interpolate(self):
        '''Interpolate and assign variables from other value arrays'''

        df_ck = self.parent_solution.df_ck

        dis = df_ck['Distance (cm)'].to_numpy()
        phi = df_ck[self.kind].to_numpy()

        f = interp.interp1d(dis, phi, kind="cubic")
        self.variable_array = f(self.y)

        
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
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0

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
