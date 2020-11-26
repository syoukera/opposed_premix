from .base_array import *
from .state_variables_array import *
import cantera as ct
import numpy as np
import pickle
import scipy.interpolate as interp

class MoleFractionArray(StateVariablesArray):
    '''Variable array for mole fraction of a species'''

    def __init__(self, parent, name, var=None):
        super().__init__(parent, var)
        self.name = name

    def interpolate(self):
        '''Interpolate and assign variables from other value arrays'''

        df_ck = self.parent_solution.df_ck

        dis = df_ck['Distance (cm)'].to_numpy()
        
        name = 'Mole_fraction_' + self.name + ' ()'
        phi = df_ck[name].to_numpy()

        f = interp.interp1d(dis, phi, kind="cubic")
        self.variable_array = f(self.y)

class MassFractionArray(StateVariablesArray):
    '''Variable array for mass fraction of a species'''

    def __init__(self, parent, name, var=None):
        super().__init__(parent, var)
        self.name = name
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0

class SpeciesList():
    '''Base class for list of mole fraction and mass fraction'''

    def __init__(self, parent, var=None):
        # Assign parent solution
        self.parent_solution = parent

        # Get species name of cti and csv
        # CTI means names from get from cantera chemistry set
        self.name_species_cti = self.parent_solution.name_species_cti

        # CK means names imported from CHEMKIN-PRO result
        with open('data/species_name_ck.txt', 'rb') as f:
            self.name_species_ck = pickle.load(f)

        # Defilen list of species
        self.list = []

    def interpolate_species_arrays(self):
        '''Interpolate each array in list'''
        
        for arr in self.list:
            if arr.name in self.name_species_ck:
                arr.interpolate()
            else:
                arr.variable_array = np.zeros(self.parent_solution.num_grid)

    def get_numpy_matrix(self):
        '''
        Retrun numpy matrix for cantera input
        shape: (num_species, num_grid)
        '''

        for i, arr in enumerate(self.list):
            if i == 0:
                mat = arr.variable_array.reshape((-1, 1))
            else:
                mat = np.concatenate((mat, arr.variable_array.reshape((-1, 1))), axis=1)

        return mat


class MoleFractionList(SpeciesList):
    '''List of Array for Mole Fractions'''

    def __init__(self, parent, var=None):
        super().__init__(parent, var)

        # Make list of MoleFractionArray
        for name in self.parent_solution.name_species_cti:
            arr = MoleFractionArray(self.parent_solution, name)
            self.list.append(arr)

class MassFractionList(SpeciesList):
    '''List of Array for Mass Fractions'''

    def __init__(self, parent, var=None):
        super().__init__(parent, var)
        
        # Make list of MassFractionArray
        for name in self.parent_solution.name_species_cti:
            arr = MassFractionArray(self.parent_solution, name)
            self.list.append(arr)