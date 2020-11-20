import numpy as np
import scipy.interpolate as interp
import sys


class BaseArray():
    '''Base class for variables array'''
    
    def __init__(self, parent, n, var=None):
        self.kind = 'Base'

        self.parent_solution = parent

        self.num_grid = n
        self.variable_array = np.zeros(self.num_grid)
        self.y = self.parent_solution.y
        
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

    def __init__(self, parent, n, var=None):
        super().__init__(parent, n, var)
        self.coef_a = np.zeros(self.num_grid)
        self.coef_b = np.zeros(self.num_grid)
        self.coef_c = np.zeros(self.num_grid)
        self.coef_d = np.zeros(self.num_grid)
        
    def solve_TDMA(self):
        '''Solve TDMA for own variable_array'''
        pass

class TemperatureArray(StateVariablesArray):
    '''Variable array for temperature'''

    def __init__(self, parent, n, var=None):
        super().__init__(parent, n, var)
        self.kind = 'Temperature (K)'
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0

class DensityArray(StateVariablesArray):
    '''Variable array for density'''
    
    def __init__(self, parent, n, var=None):
        super().__init__(parent, n, var)
        self.kind = 'Density (g/cm3)'
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0

class VelocityArray(StateVariablesArray):
    '''Variable array for velocity'''

    def __init__(self, parent, n, var=None):
        super().__init__(parent, n, var)
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

class PressureArray(StateVariablesArray):
    '''Variable array for pressure'''

    def __init__(self, parent, n, var=None):
        super().__init__(parent, n, var)
        self.kind = 'Pressure (atm)'
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0
