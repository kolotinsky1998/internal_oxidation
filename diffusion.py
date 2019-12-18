import scipy
from scipy import optimize
from scipy.optimize import newton_krylov
import warnings
import numpy as np
import math 
import math as m
from converter import Converter

warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")

class Poisson_Nernst:
    def __init__(self, converter,c_0):
        """Constructor"""
        self.converter = converter
        self.S = converter.S
        self.I = converter.I
        self.J = converter.J
        self.z = converter.z
        self.dx = converter.dx
        self.dy = converter.dy
        self.deltax = converter.deltax
        self.deltay = converter.deltay
        self.D = converter.D
        self.a = converter.a
        self.l_d = converter.l_d
        self.t = converter.t
        self.delta = converter.delta
        self.c_ref = converter.c_ref
        self.c_0 = np.zeros(self.S*self.I*self.J)
        self.V_TRI = (3.0**0.5)*(self.a**2)*self.delta/8.0
        for i in range (self.I):
            for j in range(self.J):
                for s in range(self.S):
                    k = i + j*self.I + s*self.I*self.J
                    self.c_0[k] = c_0[i][j][s] / self.c_ref
        self.time = 0
        
        ######################################################
        #####  Calculation Green function for potential   ####
        ######################################################
        self.potential_matrix = np.zeros((self.I,self.J,self.I,self.J)) # Green function
        self.L = converter.L
        self.a = converter.a
        self.a_romb = self.a * m.sqrt(3.0) 
        
        rcut = 5.0 # cat radious (real-space)
        kcut = 65.0 # cat radious (k-space)
        kunit = 4 * m.pi / m.sqrt(3.0) /  self.a_romb # unit vector of k-lattice (module)
        alpha = 15.0 # cut parameter
        for i in range(self.I):
            for j in range(self.J): 
                xij = converter.sublattice_x[i]
                yij = converter.sublattice_y[j]
                zij = 0
                for k in range(self.I):
                    for l in range(self.J):
                        contrib = 0.0
                        q = self.deltax[k]*self.deltay[l]*self.delta #unit concentration cell's charge
                        xkl = converter.sublattice_x[k]
                        ykl = converter.sublattice_y[l]
                        zkl = 0
                        x_kl_im, y_kl_im, z_kl_im = self.coordimage(xkl,ykl,zkl)
                        # real part of interaction inside cut sphere
                        for n1 in range(-int(rcut/ self.a_romb), +int(rcut/ self.a_romb) + 1):
                            for n2 in range(-int(rcut/ self.a_romb), +int(rcut/ self.a_romb) + 1):
                                Rijkln = []
                                for ind in range(len(x_kl_im)):
                                    Rijkln.append( self.distance(xij,x_kl_im[ind] +                             (n1+n2)* self.a_romb*m.cos(30.0*m.pi/180.0),yij,y_kl_im[ind],zij,z_kl_im[ind] + (n1-n2)* self.a_romb*m.sin(30*m.pi/180) ) )
                                for rad in Rijkln:
                                    if rad < rcut:
                                        if rad == 0.0:
                                            contrib += 0.0
                                        else:
                                            contrib += self.potreal(q, rad, alpha)

                        # inverse part of interaction inside cut sphere
                        for k1 in range(-int(kcut/kunit), +int(kcut/kunit)+1):
                            for k2 in range(-int(kcut/kunit), +int(kcut/kunit)+1):
                                if k1 == 0 and k2 == 0:
                                    for ind in range(len(x_kl_im)):
                                        mnozh = 2.0 * m.sqrt(m.pi) /  self.a_romb**2 / m.sin(60.0*m.pi/180.0) * q
                                        distancey = y_kl_im[ind] - yij
                                        potinv = mnozh * ( 1.0 / alpha * m.exp(-alpha**2*distancey**2) + m.sqrt(m.pi)*distancey*m.erf(alpha*distancey) )
                                        contrib -= potinv
                                else:
                                    potinv = 0.0
                                    for ind in range(len(x_kl_im)):
                                        kk = kunit**2*( k1**2 + k2**2 - 2*k1*k2*m.cos(60.0*m.pi/180.0) )
                                        kmod = m.sqrt(kk)
                                        if kmod < kcut:
                                            distancey = y_kl_im[ind] - yij
                                            mnozh1 = m.pi /  self.a_romb**2 / m.sin(60*m.pi/180) * q / kmod
                                            mnozh2 = m.exp(kmod*distancey)*( 1 - m.erf(kmod/2.0/alpha + alpha*distancey) ) + m.exp(-kmod*distancey)*(1-m.erf(kmod/2.0/alpha - alpha*distancey))
                                            mnozh = mnozh1*mnozh2
                                            distancexm = x_kl_im[ind] - xij
                                            distancezm = z_kl_im[ind] - zij
                                            potinv += mnozh * m.cos( (k1+k2)*kunit*m.cos(60.0*m.pi/180.0)*distancexm + (k1-k2)*kunit*m.sin(60.0*m.pi/180.0)*distancezm )
                                    contrib += potinv
                        if i == k and j == l:
                            contrib -= 2.0 * alpha / m.sqrt(m.pi) * q
                        self.potential_matrix[i][j][k][l] = contrib                     
        #####################################################
        #####     Calculation of chemical constants     #####
        #####################################################
        self.kob = 1e3
        self.kalb = 1e2
        self.ksb = 1e3
        self.kehb = 1e3
        self.Ko = np.zeros((self.I,self.J))
        self.Kal = np.zeros((self.I,self.J))
        self.Ks = np.zeros((self.I,self.J))
        self.Keh = np.zeros((self.I,self.J))
        self.Plo = 1
        self.Pref = 1
        self.Phi = 1e1
        k = np.zeros(self.S, dtype=np.int)
        for j in range(self.J):
            for i in range(self.I):
                for s in range(self.S):
                    k[s] =  i + j*self.I + s*self.I*self.J
                if j == 0 or j== self.J-1:
                    self.Ko[i][j] = self.c_0[k[3]]**2.0 / self.Pref**(0.5) / self.c_0[k[0]] 
                    self.Kal[i][j] = self.c_0[k[1]]**(2.0/3.0) / self.c_0[k[2]]**2.0 / self.Pref**(0.5) 
                    self.Ks[i][j] = self.c_0[k[0]]*self.c_0[k[1]]**(2.0/3.0)
                    self.Keh[i][j] = self.c_0[k[2]]*self.c_0[k[3]]
                else:
                    self.Ko[i][j] = 0.0
                    self.Kal[i][j] = 0.0
                    self.Ks[i][j] = 0.0
                    self.Keh[i][j] = self.c_0[k[2]]*self.c_0[k[3]]
                    
    def coordimage(self,x,y,z):
        ximage = []
        yimage = []
        zimage = []
        ximage.append(x)
        yimage.append(y)
        zimage.append(z)
        ximage.append(-x)
        yimage.append(y)
        zimage.append(z)
        ximage.append(-0.5*self.a-(0.5*self.a-x)*m.cos(60.0*m.pi/180.0))
        yimage.append(y)
        zimage.append((0.5*self.a - x)*m.sin(60.0*m.pi/180.0))
        ximage.append(-0.5*self.a-(0.5*self.a-x)*m.cos(60.0*m.pi/180.0))
        yimage.append(y)
        zimage.append(-(0.5*self.a - x)*m.sin(60.0*m.pi/180.0))
        ximage.append(0.5*self.a+(0.5*self.a-x)*m.cos(60.0*m.pi/180.0))
        yimage.append(y)
        zimage.append((0.5*self.a-x)*m.sin(60.0*m.pi/180.0))
        ximage.append(0.5*self.a+(0.5*self.a-x)*m.cos(60.0*m.pi/180.0))
        yimage.append(y)
        zimage.append(-(0.5*self.a-x)*m.sin(60.0*m.pi/180.0))
        return ximage, yimage, zimage
    
    def potreal(self, charge, dist, alpha):
        return charge / dist * (1 - m.erf(alpha * dist))
    
    def distance(self, x1,x2,y1,y2,z1,z2):
        dx = x2-x1
        dy = y2-y1
        dz = z2-z1
        return (dx**2 + dy**2 + dz**2) ** 0.5
                    
    def GetConcentration(self):
        c = np.zeros((self.I, self.J, self.S))
        for i in range(self.I):
            for j in range(self.J):
                for s in range(self.S):
                    k = i + j*self.I + s*self.I*self.J
                    c[i][j][s] = self.c_0[k] 
        return c
    
    #return current time in seconds
    def Current_Phys_Time(self):
        return self.converter.convert_lattice_time(self.t[self.time])
    
    def UpdateConcentration(self):
        if self.time % 10 == 0:
            print("Current simulation time: t = ", self.t[self.time])
        self.c_0 = newton_krylov(self.system, self.c_0)
        self.UpdateTime()
        
    def UpdateTime(self):
        self.time = self.time + 1
        
    def B(self, x_):
        if (math.exp(x_)-1.0 !=0):
            return x_/(math.exp(x_)-1.0)
        else:
            return 1.
    
    def potential(self, c):
        pot = np.zeros((self.I,self.J))
        
        for i in range(self.I):
            for j in range(self.J):
                for k in range(self.I):
                    for l in range(self.J):
                        full_concentration = 0.0
                        for s in range (self.S):
                            point = k + l*self.I + s*self.I*self.J
                            full_concentration = full_concentration + c[point]*self.z[s]
                        pot[i][j] += self.potential_matrix[i][j][k][l]*(full_concentration)/(self.l_d*self.l_d)
        
        '''for i in range(self.I):
            for j in range(self.J):
                pot[i][j] = -0.6 + float(j)*1.2/float(self.J)'''
        return pot
    
    def CalculateAndPrintPotential(self):
        pot = np.zeros((self.I,self.J))
        for i in range(self.I):
            for j in range(self.J):
                for k in range(self.I):
                    for l in range(self.J):
                        full_concentration = 0
                        for s in range (self.S):
                            point = k + l*self.I + s*self.I*self.J
                            full_concentration = full_concentration + self.c_0[point]*self.z[s]
                        pot[i][j] += self.potential_matrix[i][j][k][l]*(full_concentration)/(self.l_d*self.l_d)
        return pot
        
    def reaction(self, c):
        reaction_array = np.zeros((self.I,self.J,self.S))
        
        k = np.zeros(self.S, dtype=np.int)
        for j in range(self.J):
            for i in range(self.I):
                for s in range(self.S):
                    k[s] =  i + j*self.I + s*self.I*self.J
                if j == 0:
                    Ro = self.kob * (self.Ko[i][j]*self.Phi**0.5*c[k[0]] - c[k[3]]**2.0)
                    Ral = self.kalb * (self.Kal[i][j]*self.Phi**0.5*c[k[2]]**2.0 - c[k[1]]**(2.0/3.0))
                    Rs = -self.ksb * (c[k[1]]**(2.0/3.0)*c[k[0]] - self.Ks[i][j])
                    Reh = -self.kehb * (c[k[2]]*c[k[3]] - self.Keh[i][j])
                    reaction_array[i][j][0] = - Ro + Rs
                    reaction_array[i][j][1] = Ral + (2.0/3.0)*Rs
                    reaction_array[i][j][2] = - 3*Ral + Reh
                    reaction_array[i][j][3] = 2*Ro + Reh
                elif j == self.J-1:
                    Ro = self.kob * (self.Ko[i][j]*self.Plo**0.5*c[k[0]] - c[k[3]]**2.0)
                    Ral = self.kalb * (self.Kal[i][j]*self.Plo**0.5*c[k[2]]**2.0 - c[k[1]]**(2.0/3.0))
                    Rs = -self.ksb * (c[k[1]]**(2.0/3.0)*c[k[0]] - self.Ks[i][j])
                    Reh = -self.kehb * (c[k[2]]*c[k[3]] - self.Keh[i][j])
                    reaction_array[i][j][0] = - Ro + Rs
                    reaction_array[i][j][1] = Ral + (2.0/3.0)*Rs
                    reaction_array[i][j][2] = - 3*Ral + Reh
                    reaction_array[i][j][3] = 2*Ro + Reh
                else:
                    Reh = -self.kehb * (c[k[2]]*c[k[3]] - self.Keh[i][j])
                    reaction_array[i][j][0] = 0.0
                    reaction_array[i][j][1] = 0.0
                    reaction_array[i][j][2] = Reh
                    reaction_array[i][j][3] = Reh
        return reaction_array

    
    
    def system(self, c):

        fi = self.potential(c)
        R = self.reaction(c)
        S = self.S
        I = self.I
        J = self.J
        z = self.z
        f = np.zeros(S*I*J)
        dt = self.t[self.time+1] - self.t[self.time]
        D = self.D
        dx = self.dx 
        dy = self.dy
        deltax = self.deltax
        deltay = self.deltay
        a = self.a
        L = self.L
        for s in range(S):
            for j in range(1, J-1):
                for i in range(1, I-1):
                
                    k = i + j*I + s*I*J
                
                    fi_ = fi[i+1][j] - fi[i][j]
                    J1 = -D[s] * (c[k+1] * self.B(-z[s]*fi_) - c[k] * self.B(z[s]*fi_)) / dx[i]
                    
                    fi_ = fi[i][j] - fi[i-1][j]
                    J2 = -D[s] * (c[k] * self.B(-z[s]*fi_) - c[k-1] * self.B(z[s]*fi_)) / dx[i-1]
                    
                    fi_ = fi[i][j+1] - fi[i][j]
                    J3 = -D[s] * (c[k+I] * self.B(-z[s]*fi_) - c[k] * self.B(z[s]*fi_)) / dy[j]
                
                    fi_ = fi[i][j] - fi[i][j-1]
                    J4 = -D[s] * (c[k] * self.B(-z[s]*fi_) - c[k-I] * self.B(z[s]*fi_)) / dy[j-1]
                
                    f[k] = c[k] - self.c_0[k] + \
                    2.0 * dt * (J1 - J2) / (dx[i] + dx[i-1]) + \
                    2.0 * dt * (J3 - J4) / (dy[j] + dy[j-1]) - \
                    dt * R[i][j][s]
                
        for s in range(S):
            J0 = 0
            j = 0
            for i in range(I):
                k = i + j*I + s*I*J
                fi_ = fi[i][j+1] - fi[i][j]
                J0 = J0 - D[s] * deltax[i] * (c[k+I] * self.B(-z[s]*fi_) - c[k] * self.B(z[s]*fi_)) / dy[j]
            for i in range(I):
                j = 0
                k = i + j*I + s*I*J
                f[k] = c[k] - self.c_0[k] - (0.5*self.delta*J0/self.V_TRI + R[i][j][s])*dt 
            J0 = 0 
            j = J-1
            for i in range(I):
                k = i + j*I + s*I*J
                fi_ = fi[i][j] - fi[i][j-1]
                J0 = J0 - D[s] * deltax[i] * (c[k] * self.B(-z[s]*fi_) - c[k-I] * self.B(z[s]*fi_)) / dy[j-1]
            for i in range(I):
                j = J-1
                k = i + j*I + s*I*J
                f[k] = c[k] - self.c_0[k] + (0.5*self.delta*J0/self.V_TRI + R[i][j][s])*dt 
        
            for j in range(1,J-1):
                i = I-1
                k = i + j*I + s*I*J
   
                fi_ = fi[i][j] - fi[i-1][j]
                J2 = -D[s] * (c[k] * self.B(-z[s]*fi_) - c[k-1] * self.B(z[s]*fi_)) / dx[i-1]
                    
                fi_ = fi[i][j+1] - fi[i][j]
                J3 = -D[s] * (c[k+I] * self.B(-z[s]*fi_) - c[k] * self.B(z[s]*fi_)) / dy[j]
                
                fi_ = fi[i][j] - fi[i][j-1]
                J4 = -D[s] * (c[k] * self.B(-z[s]*fi_) - c[k-I] * self.B(z[s]*fi_)) / dy[j-1]
                
                f[k] = c[k] - self.c_0[k] + \
                1.0 * dt * (0. - J2) / (dx[i-1]) + \
                2.0 * dt * (J3 - J4) / (dy[j] + dy[j-1]) - \
                dt * R[i][j][s]  
            i = 0
            for j in range(1,J-1):
                k = i + j*I + s*I*J
                
                fi_ = fi[i+1][j] - fi[i][j]
                J1 = -D[s] * (c[k+1] * self.B(-z[s]*fi_) - c[k] * self.B(z[s]*fi_)) / dx[i]
                    
                fi_ = fi[i][j+1] - fi[i][j]
                J3 = -D[s] * (c[k+I] * self.B(-z[s]*fi_) - c[k] * self.B(z[s]*fi_)) / dy[j]
                
                fi_ = fi[i][j] - fi[i][j-1]
                J4 = -D[s] * (c[k] * self.B(-z[s]*fi_) - c[k-I] * self.B(z[s]*fi_)) / dy[j-1]
                
                f[k] = c[k] - self.c_0[k] + \
                1.0 * dt * (J1 - 0.0) / (dx[i]) + \
                2.0 * dt * (J3 - J4) / (dy[j] + dy[j-1]) - \
                dt * R[i][j][s]
        return f

if __name__ == "__main__":
    S = 1
    D = np.zeros(1)
    z = np.zeros(1)
    a = 1
    delta = 1
    L = 2
    T = 1
    eps_r = 2
    I = 4
    J = 4
    N = 2
    converter = Converter(S, D, z, a, delta, L, T, eps_r, I, J, N)
    poisson_nernst = Poisson_Nernst(converter)
    for n in range(N-1):
        poisson_nernst.UpdateConcentration()
    
        