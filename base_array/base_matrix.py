import numpy as np

class BaseMatrix():
    '''Base class for variables matrix'''
    
    def __init__(self, parent):
        self.name = 'Base'

        # Assign parent solution
        self.parent_solution = parent

        # Read parent parameters 
        self.num_grid = self.parent_solution.num_grid
        self.num_species = self.parent_solution.num_species
        self.dy = self.parent_solution.dy
        self.y = self.parent_solution.y
        self.dt = self.parent_solution.dt

        # Define my own variables
        self.variable_array = np.zeros((self.num_grid, self.num_species))
        
    def average_variables(self):
        '''
        Average variables at intermediate grid point
        for variable in stagard grid
        's' means Small case
        '''

        self.variable_array_s = self.variable_array.copy()
        self.variable_array_s[:-1, :] += self.variable_array[1:, :].copy()
        self.variable_array_s[:-1, :] /= 2

class ParameterMatrix(BaseMatrix):
    '''Variable matrix for parameter variables'''

    def __init__(self, parent, name):
        super().__init__(parent)
        self.name = name
