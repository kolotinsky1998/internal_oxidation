import numpy as np
from converter import Converter
from diffusion import Poisson_Nernst
np.seterr(all=None, divide="warn", over="warn", under="warn", invalid="warn")
############################################################
############### Define system parameters####################
############################################################
T = 1900 #K
D_ref = 1
c_ref = 1.0e24 # 1/sm^3
eps_r = 9.8 # c-units
l_D = 9.4e-9 # m
L = 1.0e-6 # m
a = 1.0e-6 # m 
delta = 1.0e-9 # m
k_ob  = 1.0e3 
k_alb = 1.0e2
k_ehb = 1.0e3
k_sb  = 1.0e3
S = 4 #number of species
I = 3 #number of cells in x-direction
J = 10 #number of cells in y-direction
N = 150 #number of time iterations
## indexes
## 0 - oxygen vacancy, 1 - aluminium vacancy, 2 - electrons, 3 - holes ##
D = np.zeros(S)
D[0] = 0.01*D_ref
D[1] = 0.01*D_ref
D[2] = 1*D_ref
D[3] = 1*D_ref
## charges
z = np.zeros(S)
z[0] = 2
z[1] = -3
z[2] = -1
z[3] = 1
## initial concentrations
c_0 = np.zeros((I,J,S))
for i in range(I):
    for j in range(J):
        c_0[i][j][0] = 5.0e-3*c_ref
        c_0[i][j][1] = (10.0/3.0)*1.0e-3*c_ref
        c_0[i][j][2] = 0.1*c_ref
        c_0[i][j][3] = 0.1*c_ref
        
############################################################
############################################################
############################################################

converter = Converter(S, D, z, a, delta, L, T, eps_r, I, J, N,c_ref)
converter.PrintDimmensionalUnits()
poisson_nernst = Poisson_Nernst(converter, c_0)

c_image = np.zeros((N,S,J,I))
for n in range(N-1):
    poisson_nernst.UpdateConcentration()
    c_time = poisson_nernst.GetConcentration()
    for i in range(I):
        for j in range(J):
            for s in range(S):
                c_image[n][s][j][i] = c_time[i][j][s]
pot = poisson_nernst.CalculateAndPrintPotential()

x_=converter.sublattice_x
y_=converter.sublattice_y