from .base_array import *
import cantera as ct
import numpy as np
import pickle
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

    def solve(self):
        '''High level function to Solve variable array'''
        self.calc_coef()
        self.solve_TDMA()

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

class DensityArray(StateVariablesArray):
    '''Variable array for density'''
    
    def __init__(self, parent, var=None):
        super().__init__(parent, var)
        self.name = 'Density (g/cm3)'
    
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
        self.name = 'Axial_velocity (cm/sec)'
    
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
        self.name = 'Radial_velocity (1/s)'

    def initialize(self):
        '''Initialize variable array for simulation'''
        center = int(self.num_grid/2) + 1
        self.variable_array = np.zeros(self.num_grid)
        self.variable_array[:center] = np.linspace(0, 80, center)
        self.variable_array[-center:] = np.linspace(80, 0, center)
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''

        # Read variables from parent soluction
        mu_s = self.parent_solution.mu.variable_array_s
        R    = self.parent_solution.R.variable_array
        V_u  = self.parent_solution.V.variable_array_u
        TPG  = self.parent_solution.TPG.variable_array[0]


        for p in range(self.num_grid):
        
            # Upper boundary condition
            if p == 0:
                self.coef_a[p] = 1.0
                self.coef_b[p] = 0.0
                self.coef_c[p] = 0.0
                self.coef_d[p] = 0.0
                continue
                
            # Lower boundary condition
            if p == self.num_grid-1:
                self.coef_a[p] = 1.0
                self.coef_b[p] = 0.0
                self.coef_c[p] = 0.0
                self.coef_d[p] = 0.0
                continue
            
            # Inner points of numerical grid    
            C_p = mu_s[p]*self.dt/(R[p]*self.dy**2)    
            C_w = mu_s[p-1]*self.dt/(R[p]*self.dy**2)

            self.coef_a[p] = 1 + C_p + C_w + self.variable_array[p]*self.dt
            self.coef_b[p] = C_p
            self.coef_c[p] = C_w
            
            # term v*dG/dy
            method = 'upwind'
            if method == 'center':
                D = V_u[p]*self.dt/(2*self.dy)
                self.coef_b[p] += - D
                self.coef_c[p] += + D
            elif method == 'upwind':
                D = V_u[p]*self.dt/self.dy
                if V_u[p] > 0:
                    self.coef_a[p] += D
                    self.coef_c[p] += D
                else:
                    self.coef_a[p] += - D
                    self.coef_b[p] += - D
            else:
                sys.exit('Unknow method is choosen')        
                
            self.coef_d[p] = self.variable_array[p] - TPG*self.dt/R[p]

class PressureArray(StateVariablesArray):
    '''Variable array for pressure'''

    def __init__(self, parent, var=None):
        super().__init__(parent, var)
        self.name = 'Pressure (g/cm/s2)'

    def interpolate(self):
        '''
        Interpolate and assign variables from other value arrays
        Converted to cgs-unit [g/cm/s2] from [atm]
        '''

        df_ck = self.parent_solution.df_ck

        dis = df_ck['Distance (cm)'].to_numpy()
        phi = df_ck['Pressure (atm)'].to_numpy()
        phi_cgs = phi*1013250

        f = interp.interp1d(dis, phi_cgs, kind="cubic")
        self.variable_array = f(self.y)
    
    def calc_coef(self):
        '''Calculate coefficients for TDMA'''
        self.coef_a = np.ones(self.num_grid) * 1.0
        self.coef_b = np.ones(self.num_grid) * 0.0
        self.coef_c = np.ones(self.num_grid) * 0.0
        self.coef_d = np.ones(self.num_grid) * 1.0