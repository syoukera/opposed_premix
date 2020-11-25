import numpy as np
import scipy.interpolate as interp
import sys


class BaseArray():
    '''Base class for variables array'''
    
    def __init__(self, parent, var=None):
        self.name = 'Base'

        # Assign parent solution
        self.parent_solution = parent

        # Read parent parameters 
        self.num_grid = self.parent_solution.num_grid
        self.dy = self.parent_solution.dy
        self.y = self.parent_solution.y
        self.dt = self.parent_solution.dt

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
        phi = df_ck[self.name].to_numpy()

        f = interp.interp1d(dis, phi, kind="cubic")
        self.variable_array = f(self.y)


class ParameterArray(BaseArray):
    '''Variable array for state variables'''

    def __init__(self, parent, name, var=None):
        super().__init__(parent, var)
        self.name = name
        self.interpolate()