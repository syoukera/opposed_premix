from .base_array import *
from .state_variables_array import *
import cantera as ct
import numpy as np

class TemperatureArray(StateVariablesArray):
    '''Variable array for temperature'''

    def __init__(self, parent, var=None):
        super().__init__(parent, var)
        self.name = 'Temperature (K)'

    def initialize(self):
        '''Initialize array of temperature'''
        
        self.variable_array = np.ones(self.num_grid) * 298
        self.variable_array[int(self.num_grid/4):-int(self.num_grid/4)] = 2000

    def calc_j_diff(self):
        '''Calculate diffusion flux'''
        y  = self.parent_solution.y
        R  = self.parent_solution.R.variable_array
        T  = self.parent_solution.T.variable_array
        X  = self.parent_solution.X_list.get_numpy_matrix()
        Y  = self.parent_solution.Y_list.get_numpy_matrix()
        Df = self.parent_solution.Df.variable_array
        Dt = self.parent_solution.Dt.variable_array


        dXdy = np.gradient(X, y, axis=0)
        dTdy = np.gradient(T, y, axis=0)

        dXdy_mass = np.multiply(dXdy, np.divide(Y, X))
        j_f = np.multiply(np.multiply(R.reshape((-1, 1)), Df), dXdy_mass)
        j_t = np.multiply(Dt, np.divide(dTdy, T).reshape((-1, 1)))
        self.j_diff = - j_f - j_t
        
        # Cnvert nan to 0.0
        self.j_diff[np.isnan(self.j_diff)] = 0.0

    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''
        lm_s  = self.parent_solution.lm.variable_array_s
        R     = self.parent_solution.R.variable_array
        cp_m  = self.parent_solution.cp_m.variable_array
        cp_i  = self.parent_solution.cp_i.variable_array
        V_u   = self.parent_solution.V.variable_array_u
        T_old = self.parent_solution.T_old
        P     = self.parent_solution.P.variable_array
        P_old = self.parent_solution.P_old
        h     = self.parent_solution.h.variable_array
        wdot  = self.parent_solution.wdot.variable_array

        self.calc_j_diff()

        for p in range(self.num_grid):
        
            # Upper boundary condition
            if p == 0:
                self.coef_a[p] = 1.0
                self.coef_b[p] = 0.0
                self.coef_c[p] = 0.0
                self.coef_d[p] = 298
                continue
                
            # Lower boundary condition
            if p == self.num_grid-1:
                self.coef_a[p] = 1.0
                self.coef_b[p] = 0.0
                self.coef_c[p] = 0.0
                self.coef_d[p] = 298
                continue

            # Inner points of numerical grid    
            C_p = lm_s[p]/(R[p]*cp_m[p]*self.dy**2)
            C_w = lm_s[p-1]/(R[p]*cp_m[p]*self.dy**2)
            q_p = 1/(R[p]*cp_m[p])*np.dot(cp_i[p, :], self.j_diff[p, :])
            
            self.coef_a[p] = 1/self.dt + C_p + C_w + np.abs(V_u[p]/self.dy) + np.abs(q_p/self.dy)
            self.coef_b[p] = C_p + max(0, -V_u[p]/self.dy) + max(0, -q_p/self.dy)
            self.coef_c[p] = C_w + max(0, +V_u[p]/self.dy) + max(0, +q_p/self.dy)
            self.coef_d[p] = T_old[p]/self.dt + (P[p] - P_old[p])/(R[p]*self.dt) - np.dot(h[p, :], wdot[p, :])/(R[p]*cp_m[p])
