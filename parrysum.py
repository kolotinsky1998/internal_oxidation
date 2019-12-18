"""
Created on Thu Oct 10 09:16:02 2019

@author: vladislavnikolaev
"""

import math as m
import cmath as cm
import numpy as np

def distance(x1,x2,y1,y2,z1,z2):
    """ Возвращает расстояние между двумя точками """
    dx = x2-x1
    dy = y2-y1
    dz = z2-z1
    return (dx**2 + dy**2 + dz**2) ** 0.5

def potreal(charge, dist, alpha):
    """ Определяет действительную часть взаимодействия в сумме Перри """
    return charge / dist * (1 - m.erf(alpha * dist))

def coord(i,j,mesh): 
    """ Определяет координаты узла между двумя узлами сетки """
    return 0.5 * ( mesh[0][i] + mesh[0][i+1] ), 0.5 * ( mesh[1][j] + mesh[1][j+1] ), 0.0

def coordimage(x,y,z,L):
    ximage = []
    yimage = []
    zimage = []
    ximage.append(x)
    yimage.append(y)
    zimage.append(z)
    ximage.append(-x)
    yimage.append(y)
    zimage.append(z)
    ximage.append(-0.5*L-(0.5*L-x)*m.cos(60.0*m.pi/180.0))
    yimage.append(y)
    zimage.append((0.5*L - x)*m.sin(60.0*m.pi/180.0))
    ximage.append(-0.5*L-(0.5*L-x)*m.cos(60.0*m.pi/180.0))
    yimage.append(y)
    zimage.append(-(0.5*L - x)*m.sin(60.0*m.pi/180.0))
    ximage.append(0.5*L+(0.5*L-x)*m.cos(60.0*m.pi/180.0))
    yimage.append(y)
    zimage.append((0.5*L-x)*m.sin(60.0*m.pi/180.0))
    ximage.append(0.5*L+(0.5*L-x)*m.cos(60.0*m.pi/180.0))
    yimage.append(y)
    zimage.append(-(0.5*L-x)*m.sin(60.0*m.pi/180.0))
    return ximage, yimage, zimage


def pot_matrix(mesh, delta): # функция от сетки X*Y (x бежит по границе поверхности, y - по границе зерна), толщины межзеренной границы
    I = len(mesh[0]) - 1 # задаю размер массива по Х
    J = len(mesh[1]) - 1 # задаю размер массива по У
    
    potential_matrix = np.zeros((I,J,I,J)) # матрица, которую я собираюсь вернуть в конце функции, она будет задавать взаимодействие ячейки i,j с ячейкой k,l
    rcut, kcut = 5.0, 65.0 # радиусы обрезания для прямого и обратного пространств
    L = 1.0 # длина стороны шестиугольника a_hex, должна быть равна 1.0
    a = 2.0 * L * m.sin( 60 * m.pi / 180 ) # длина стороны ромба через L
    kunit = 4 * m.pi / m.sqrt(3.0) / a # единичный вектор обратной гексоганальной решетки (модуль)
    alpha = 15.0 # параметр обрезания, должен быть больше нуля, НЕ ДОЛЖЕН влиять на результат

    for i in range(I):
        for j in range(J): # рассматриваем произвольную ячейку с координатами (i,j)
            xij, yij, zij = coord(i,j,mesh)
            for k in range(I):
                for l in range(J): # рассчитываем коэффициент взаимодействия с ячейкой (k,l)

                    contrib = 0.0
                    q = (mesh[0][k+1] - mesh[0][k]) * (mesh[1][l+1] - mesh[1][l]) * delta # / 2.0 # заряд ячейки (k,l) с единичной концентрацией
                    xkl, ykl, zkl = coord(k,l,mesh)
                    x_kl_im, y_kl_im, z_kl_im = coordimage(xkl,ykl,zkl,L)
                    # рассчитываю реальную часть взаимодействия внутри радиуса обрезания rcut
                    for n1 in range(-int(rcut/a), +int(rcut/a) + 1):
                        for n2 in range(-int(rcut/a), +int(rcut/a) + 1):
                            Rijkln = []
                            for ind in range(len(x_kl_im)):
                                Rijkln.append( distance(xij,x_kl_im[ind] + (n1+n2)*a*m.cos(30.0*m.pi/180.0),yij,y_kl_im[ind],zij,z_kl_im[ind] + (n1-n2)*a*m.sin(30*m.pi/180) ) )
                            for rad in Rijkln:
                                if rad < rcut:
                                    if rad == 0.0:
                                        contrib += 0.0
                                    else:
                                        contrib += potreal(q, rad, alpha)

                    # рассчитываю взаимодействие в обратной решетке
                    for k1 in range(-int(kcut/kunit), +int(kcut/kunit)+1):
                        for k2 in range(-int(kcut/kunit), +int(kcut/kunit)+1):
                            if k1 == 0 and k2 == 0:
                                for ind in range(len(x_kl_im)):
                                    mnozh = 2.0 * m.sqrt(m.pi) / a**2 / m.sin(60.0*m.pi/180.0) * q
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
                                        mnozh1 = m.pi / a**2 / m.sin(60*m.pi/180) * q / kmod
                                        mnozh2 = m.exp(kmod*distancey)*( 1 - m.erf(kmod/2.0/alpha + alpha*distancey) ) + m.exp(-kmod*distancey)*(1-m.erf(kmod/2.0/alpha - alpha*distancey))
                                        mnozh = mnozh1*mnozh2
                                        distancexm = x_kl_im[ind] - xij
                                        distancezm = z_kl_im[ind] - zij
                                        potinv += mnozh * m.cos( (k1+k2)*kunit*m.cos(60.0*m.pi/180.0)*distancexm + (k1-k2)*kunit*m.sin(60.0*m.pi/180.0)*distancezm )
                                contrib += potinv
                    if i == k and j == l:
                        contrib -= 2.0 * alpha / m.sqrt(m.pi) * q
                    potential_matrix[i][j][k][l] = contrib

    return potential_matrix

mesh = np.array([np.linspace(0.0, 0.5, 2), np.linspace(0.0, 1.0, 3)])
print(mesh)
delta = 0.001

potmatr = pot_matrix(mesh,delta)
print(potmatr)  
              
