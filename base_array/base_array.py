import numpy as np
import sys

class BaseArray():
    '''Base class for variables array'''
    
    def __init__(self, parent, n, var=None):
        self.parent_solution = parent

        self.num_grid = n
        self.variable_array = np.zeros(self.num_grid)
        
        if var is None:
            pass
        else:
            if vals.shape == self.variable_array.shape:
                self.variable_array = var
            else:
                sys.exit('Exception: Thrown vals has different shape.')

    def average_variables(self):
        '''Average variables for using in stagard grid'''
        pass                

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
    
    def interpolate(self):
        '''Interpolate and assign variables from other value arrays'''
        pass

class TemperatureArray(StateVariablesArray):
    '''Variable array for temperature'''
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0

    def get_density_array(self):
        pass

class DensityArray(StateVariablesArray):
    '''Variable array for density'''
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0

class VelocityArray(StateVariablesArray):
    '''Variable array for velocity'''
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0

class PressureArray(StateVariablesArray):
    '''Variable array for pressure'''
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0
