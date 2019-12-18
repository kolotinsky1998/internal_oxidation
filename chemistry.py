import math as m
import cmath as cm
import numpy as np

def equilibrium_constants(I,J,c0,Pref):
    Ko = np.zeros((I,J))
    Kal = np.zeros((I,J))
    Ks = np.zeros((I,J))
    Keh = np.zeros((I,J))
    # определяем константы реакций
    for j in range(J):
        for i in range(I):
            if j == 0:
                Ko[i][j] = c0[i][j][3]**2.0 / Pref**(0.5) / c0[i][j][0] # where is the concentration of neutral oxygen in this constant?
                Kal[i][j] = c0[i][j][1]**(2.0/3.0) / c0[i][j][2]**2.0 / Pref**(0.5) # the same here
                Ks[i][j] = c0[i][j][0]*c0[i][j][1]**(2.0/3.0)
                Keh[i][j] = c0[i][j][2]*c0[i][j][3]
            elif j == J-1:
                Ko[i][j] = c0[i][j][3]**2.0 / Pref**(0.5) / c0[i][j][0] # where is the concentration of neutral oxygen in this constant?
                Kal[i][j] = c0[i][j][1]**(2.0/3.0) / c0[i][j][2]**2.0 / Pref**(0.5) # the same here
                Ks[i][j] = c0[i][j][0]*c0[i][j][1]**(2.0/3.0)
                Keh[i][j] = c0[i][j][2]*c0[i][j][3]
            else:
                Ko[i][j] = 0.0
                Kal[i][j] = 0.0
                Ks[i][j] = 0.0
                Keh[i][j] = c0[i][j][2]*c0[i][j][3]
            
    return Ko, Kal, Ks, Keh #четыре двумерных массива констант равновесия в начальный момент времени для разных реакций, вызывается в начале кода, используется в ходе всего выполнения, эти массивы хранятся всех расчетов, неизменяются (потому что являются функциями начальных концентраций, начальных!)

def reactions(I,J,S,kob,kalb,ksb,kehb,Ko,Kal,Ks,Keh,c,Phi,Plo): # I,J,S - размеры массивов, затем 4 константы обратной реакции - ДАНЫ В ДИССЕРТАЦИИ НА СТРАНИЦЕ 49; дальше 4 массива констант равновесия, полученных в функции выше; и наконец, массив концентраций (трехмерный, в обратном порядке, ну ладно) и значения давлений сверху и снизу - это в диссертации.

    # индексы: 0 - вакансии кислорода, 1 - алюминия, 2 - электроны, 3 - дырки
    reaction_array = np.zeros((I,J,S))
    for j in range(J):
        for i in range(I):
            if j == 0:
                Ro = kob * (Ko[i][j]*Phi**0.5*c[i][j][0] - c[i][j][3]**2.0)
                Ral = kalb * (Kal[i][j]*Phi**0.5*c[i][j][2]**2 - c[i][j][1]**(2.0/3.0))
                Rs = -ksb * (c[i][j][1]**(2.0/3.0)*c[i][j][0] - Ks[i][j])
                Reh = -kehb * (c[i][j][2]*c[i][j][3] - Keh[i][j])
                reaction_array[i][j][0] = - Ro + Rs
                reaction_array[i][j][1] = Ral + (2.0/3.0)*Rs
                reaction_array[i][j][2] = - 3*Ral + Reh
                reaction_array[i][j][3] = 2*Ro + Reh
            elif j == J-1:
                Ro = kob * (Ko[i][j]*Plo**0.5*c[i][j][0] - c[i][j][3]**2)
                Ral = kalb * (Kal[i][j]*Plo**0.5*c[i][j][2]**2 - c[i][j][1]**(2.0/3.0))
                Rs = -ksb * (c[i][j][1]**(2.0/3.0)*c[i][j][0] - Ks[i][j])
                Reh = -kehb * (c[i][j][2]*c[i][j][3] - Keh[i][j])
                reaction_array[i][j][0] = - Ro + Rs
                reaction_array[i][j][1] = Ral + (2.0/3.0)*Rs
                reaction_array[i][j][2] = - 3*Ral + Reh
                reaction_array[i][j][3] = 2*Ro + Reh
            else:
                Reh = -kehb * (c[i][j][2]*c[i][j][3] - Keh[i][j])
                reaction_array[i][j][0] = 0.0
                reaction_array[i][j][1] = 0.0
                reaction_array[i][j][2] = Reh
                reaction_array[i][j][3] = Reh

    return reaction_array

# это проверка

I,J,S = 4,5,4
c0 = np.zeros((I,J,S))
for i in range(I):
    for j in range(J):
        for s in range(S):
            if s == 0:
                c0[i,j,s]=0.005
            elif s == 1:
                c0[i,j,s]=(10.0/3.0)*0.001
            else:
                c0[i,j,s]=0.1

Phi = 1e5
Plo = 1.0
Pref = Plo

Ko, Kal, Ks, Keh = equilibrium_constants(I,J,c0,Pref)
print (Ko)
kob = 1e3
kalb = 1e2
ksb = kehb = 1e3
c=c0
reaction_array = reactions(I,J,S,kob,kalb,ksb,kehb,Ko,Kal,Ks,Keh,c,Phi,Plo)
print(reaction_array[:,:,1])

