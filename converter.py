#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
class Converter:
    def __init__(self, S, D, z, a, delta, L, T, eps_r, I, J, N, c_ref):
        #physical conditions
        self.T = T
        self.eps_r = eps_r
        self.c_ref = c_ref
        #Some physical fundamental quatities
        self.k_b = 1.3807e-23 # Дж/К
        self.e_0 = 1.6022e-19 #Кл
        self.eps_0 = 8.8542e-12 #Ф/м
        
        #Coefficient of convertation from lattice units to physical units
        #Phycial_Unit = Lattice_Unit * converer_unit
        #length
        self.converter_length = L
        #diffusion coefficient
        self.converter_D = np.amax(D)
        #potential
        self.converter_fi = self.k_b*self.T/self.e_0
        #time
        self.converter_t = (L**2)/np.amax(D)
        #spiecies' parameters
        self.S = S
        self.D = np.zeros(S) 
        for s in range (S):
            self.D[s] = D[s] / self.converter_D
        self.z = z
        
        #geometric properties
        self.a = a / self.converter_length
        self.delta = delta / self.converter_length
        self.L = L / self.converter_length
        self.l_d = ((self.eps_0*self.eps_r*self.k_b*self.T)/(self.e_0*self.e_0*self.c_ref))**0.5/self.converter_length
        
        #generating lattice
        self.I = I 
        self.J = J
        self.lattice_x = np.zeros(I+1)
        self.lattice_y = np.zeros(J+1)
        for i in range(I+1):
            self.lattice_x[i] = (self.a*0.5*float(i))/(float(I))
        
        for j in range(J+1):
            self.lattice_y[j] = (self.L*float(j))/(float(J))
         
        #generating sublattice
        self.sublattice_x = np.zeros(I)
        self.sublattice_y = np.zeros(J)
        for i in range(I):
            self.sublattice_x[i] = 0.5*(self.lattice_x[i+1] + self.lattice_x[i])
        
        for j in range(J):
            self.sublattice_y[j] = 0.5*(self.lattice_y[j+1] + self.lattice_y[j])
            
        #generating coordinate steps for sublattice - i direction
        self.dx = np.zeros(I-1)
        for i in range (I-1):
            self.dx[i]=0.5*(self.lattice_x[i+2] + self.lattice_x[i+1]) - \
            0.5*(self.lattice_x[i+1] + self.lattice_x[i])
            
        #generating coordinate steps for sublattice - j direction
        self.dy = np.zeros(J-1)
        for j in range (J-1):
            self.dy[j] = 0.5*(self.lattice_y[j+2] + self.lattice_y[j+1]) - \
            0.5*(self.lattice_y[j+1] + self.lattice_y[j])
        
        #generating coordinate steps for lattice - i direction
        self.deltax = np.zeros(I)
        for i in range (I):
            self.deltax[i] = self.lattice_x[i+1] - self.lattice_x[i]
            
        #generating coordinate steps for lattice - j direction
        self.deltay = np.zeros(J)
        for j in range (J):
            self.deltay[j] = self.lattice_y[j+1] - self.lattice_y[j]
            
        #generating time loop array
        self.N = N
        self.t = np.zeros(self.N)
        self.dt = min((np.amin(self.dy)**2),(np.amin(self.dx)**2)) / np.amax(self.D)
        self.tmax = self.N * self.dt
        for n in range(self.N):
            self.t[n] = self.dt*n
            
    #convert lattice potential to physical one 
    def convert_fi(self,fi):
        return fi*self.converter_fi
    
    #return sublattice in metres i direction
    def real_sublattice_x(self):
        real_x = np.zeros(self.I)
        for i in range(self.I):
            real_x[i] = self.sublattice_x[i]*self.converter_length
        return real_x
    
    #return sublattice in metres j direction
    def real_sublattice_y(self):
        real_y = np.zeros(self.J)
        for j in range(self.J):
            real_y[j] = self.sublattice_y[j]*self.converter_length
        return real_y
    #return current calculation time in seconds
    def convert_lattice_time(self, time):
        return time*self.converter_t
    
    def PrintDimmensionalUnits(self):
        print("################# Dimmensional units #####################")
        print("Thickness of the scale: L = ", self.L)
        print("Shielding lenght l_D = ", self.l_d)
        print("Length of the hexagonal side: a = ", self.a)
        print("Thickness of the grain boundary: delta = ", self.delta)
        print("Diffusion coefficients for all species: O = ", self.D[0], ", Al = ", self.D[1], ", e = ", self.D[2], ", h = ", self.D[3])
        print("Charges for all species: O = ", self.z[0], ", Al = ", self.z[1], ", e = ", self.z[2], ", h = ", self.z[3])
        print("Integration time step: dt = ", self.dt)
        print("Simulation time: t_max = ", self.tmax)
        
        
        

if __name__ == "__main__":
    S = 1
    D = np.array(1)
    z = np.array(1)
    a = 1
    delta = 1
    L = 2
    T = 1
    eps_r = 2
    I = 3
    J = 3
    N = 2
    converter = Converter(S, D, z, a, delta, L, T, eps_r, I, J, N)
    