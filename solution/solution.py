import cantera as ct
import numpy as np
import pandas as pd
from tqdm import tqdm
import base_array

class BaseSolution():
    '''Base solution class for simulation'''

    def __init__(self):
        self.length = 10 # cm
        self.num_grid = 2001

        # grid
        self.dy = self.length/(self.num_grid - 1)
        self.y = np.linspace(0, self.length, self.num_grid) # cm

        # calculation time
        self.start_time = 0.0
        self.end_time   = 1e-4
        self.dt         = 1e-6

        # calclation step
        self.total_step = int(self.end_time/self.dt)

        # Path to csv file exported using CHEMKIN-PRO
        self.path_ck = 'data/export_OpposedDiffusioon_CH4_GRI.csv'
        self.df_ck = pd.read_csv(self.path_ck)

        # Chemistry file in cantera
        self.gas = ct.Solution('gri30.xml')
        self.ct_array = ct.SolutionArray(self.gas, (self.num_grid))
        self.name_species_cti = self.ct_array.species_names
        self.num_species = len(self.name_species_cti)

        # Declear array of state variables
        self.R = base_array.DensityArray(parent=self)
        self.P = base_array.PressureArray(parent=self)
        self.T = base_array.TemperatureArray(parent=self)
        self.V = base_array.AxialVelocityArray(parent=self)
        self.G = base_array.RadialVelocityArray(parent=self)

        # Declear array of old varialbes
        self.R_old = np.zeros(self.num_grid)
        self.P_old = np.zeros(self.num_grid)
        self.T_old = np.zeros(self.num_grid)
        self.V_old = np.zeros(self.num_grid)
        self.G_old = np.zeros(self.num_grid)
        self.Y_old = np.zeros((self.num_grid, self.num_species))

        # Delcear array list of Chemical species
        self.X_list = base_array.MoleFractionList(parent=self)
        self.Y_list = base_array.MassFractionList(parent=self)  

        # Declear array of parameters
        self.TPG = base_array.ParameterArray(parent=self, name='Pressure_gradient (g/cm3-s2)')

        # Declear array of parameters obtained using cantera
        self.Mn   = self.ct_array.molecular_weights
        self.cp_m = base_array.ParameterArray(parent=self, name='Specific_heat_Cp (erg/g-K)')
        self.cp_i = base_array.ParameterMatrix(parent=self, name='Specific_heat_Cp (erg/g-K)')
        self.lm   = base_array.ParameterArray(parent=self, name='Mixture_thermal_conductivity (erg/cm-K-sec)')
        self.Df   = base_array.ParameterMatrix(parent=self, name='Multicomponent diffusion coefficients (cm2/s)')
        self.Dt   = base_array.ParameterMatrix(parent=self, name='Species thermal diffusion coefficients(g/cm/s)')
        self.mu   = base_array.ParameterArray(parent=self, name='Mixture_viscosity (g/cm-sec)')
        self.h    = base_array.ParameterMatrix(parent=self, name='Species partial molar enthalpies (erg/mol)')
        self.wdot = base_array.ParameterMatrix(parent=self, name='Net production rates for each species (mol/cm3/s)')

    def solve(self):
        '''Solve probrems using TDMA and SIMPLE method'''
        
        # Inititalise calclation
        self.time = self.start_time
        # Initialize function
        self.interporate_arrays()
        self.initialize_arrays()
        self.setup_cantera_array()
        self.average_arrays()

        for n_step in tqdm(range(self.total_step)):

            self.save_old_value()

            # Prepare parameters
            self.setup_cantera_array()
            self.average_arrays()

            # Time evolution
            self.time = n_step*self.dt
            self.time_step()

    def interporate_arrays(self):
        '''Interpolate required arrays'''

        # State variables array
        self.R.interpolate()
        self.V.interpolate()
        self.G.interpolate()
        self.T.interpolate()
        self.P.interpolate()
        self.X_list.interpolate()

        # Parameter array
        self.TPG.interpolate()
        
    def initialize_arrays(self):
        '''Initialize solved arrays'''
        
        # self.R.initialize()
        # self.V.initialize()
        # self.G.initialize()
        self.T.initialize()
        # self.P.initialize()
        # self.Y_list.initialize()

    def setup_cantera_array(self):
        '''First setup for cantera array and parameters'''

        self.assign_TPX_to_cantera()
        self.Y_list.assign_numpy_matrix(self.ct_array.Y)
        self.get_parameters_from_cantera()


    def assign_TPX_to_cantera(self):
        '''Assign TPX to cantera, unit: K, Pa, -'''

        self.ct_array.TPX = self.T.variable_array,      \
                            self.P.variable_array*1e-1, \
                            self.X_list.get_numpy_matrix()
        self.ct_array.transport_model='Multi'

    def get_parameters_from_cantera(self):
        '''Get thermdynamic and transport properties from cantera'''

        # Specific heat capacity at constant pressure [J/kg/K]. convert to [erg/g/K]
        self.cp_m.variable_array = self.ct_array.cp_mass * 1e4
        # Array of species partial molar specific heat capacities at constant pressure [J/kmol/K]. convert to [erg/mol/K]
        cp_i_mole = self.ct_array.partial_molar_cp * 1e4
        # Array of species partial mass specific heat capacities at constant pressure [erg/g/K].
        self.cp_i.variable_array = np.divide(cp_i_mole, self.Mn)
        # Thermal conductivity. [W/m/K]. convert to [erg/cm/K/s]
        self.lm.variable_array = self.ct_array.thermal_conductivity * 1e5
        # Multicomponent diffusion coefficients [m^2/s]. convert to [cm2/s]
        self.Df.variable_array = self.ct_array.mix_diff_coeffs_mass * 1e4
        # species thermal diffusion coefficients [kg/m/s]. convert to [g/cm/s]
        self.Dt.variable_array = self.ct_array.thermal_diff_coeffs * 1e4
        # Viscosity [Pa-s]. convert to [g/cm/s]
        self.mu.variable_array = self.ct_array.viscosity * 1e1
        # Array of species partial molar enthalpies [J/kmol]. convert to [erg/mol]
        self.h.variable_array = self.ct_array.partial_molar_enthalpies * 1e4
        # Net production rates for each species. [kmol/m^3/s] for bulk phases. convert to [mol/cm3/s]
        self.wdot.variable_array = self.ct_array.net_production_rates * 1e-3

    def average_arrays(self):
        '''Average operation for required arrays'''

        self.V.average_variables()
        self.R.average_variables()
        self.mu.average_variables()
        self.lm.average_variables()

    def save_old_value(self):
        '''Save old value of state variables'''
        self.R_old = self.R.variable_array
        self.P_old = self.P.variable_array
        self.T_old = self.T.variable_array
        self.V_old = self.V.variable_array
        self.G_old = self.G.variable_array
        self.Y_old = self.Y_list.get_numpy_matrix()

    def time_step(self):
        '''Progress a time step'''

        self.T.solve()

