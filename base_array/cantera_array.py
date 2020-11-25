from .base_array import *
import cantera as ct
import numpy as np

class CanteraArray(BaseArray):
    '''Variable array for cantera variables'''

    def __init__(self, parent, var=None):
        super().__init__(parent, var)

        gas = ct.Solution('gri30.xml')
        states = ct.SolutionArray(gas, (num_grid))
        name_species_cti = states.species_names

        # CHMKINの名前を取得
        name_species_raw = df_ck.iloc[:, 9:56].columns
        name_species = []
        for n in name_species_raw:
            name_species.append(n.split('_')[2].split(' ')[0])

