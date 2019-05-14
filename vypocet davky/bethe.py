# -*- coding: utf-8 -*-
"""
Created on Fri May  3 10:17:21 2019

@author: michal.sestak
"""

from math import sqrt, log, pi
from scipy.constants import N_A, epsilon_0, e
from scipy.interpolate import interp1d, CubicSpline
from scipy.integrate import quad
import numpy as np
import pandas as pd
from konstanty import *
import time
import matplotlib.pyplot as plt


#Define some natural constants

mp = 1.672621637e-27        #kg
me = 9.1093821499999992e-31 #kg
ma = 6.64465620e-27         #kg, mass of the alfa particle
qe = e                      #C
eps0 = epsilon_0            #F/m
c = 299792458               #m/s
#rho_vzduch = 0.00120479        #g/cm^3


#and nnow some problem specific ones

Ekin = 5.49e6*qe               #6 MeV in J
Ekin_0 = Ekin
m0   = ma                   #mass of projectile
Z1   = 2                    #nuclear charge of projectile (alfa)
EB   = 85.7*qe              #85.7 eV in J, i.e., Mean Excitation Energy of Dry Air
#ne   = 3.3456e29            #electron density of water in m^-3

E_alfa=np.array([5.489, 6.002, 7.689]) #v MeV

#def dosah_alfa_vzduch(E):
#    '''
#    Input:
#        E(float): kineticka energie [MeV]
#    Output:
#        R(float): dosah ve vzduchu [mm]
#    '''
#    R=3.18*E**(3/2)
#    return R
#
#dosah=dosah_alfa_vzduch(E_alfa)*0.1

data=pd.read_csv('ASTAR vzduch.txt', skiprows = 8, sep="\t")
#databaze=pd.DataFrame([data.loc[:,'E[MeV]'], data.loc[:,'S_total[MeV cm2/g]'], data.loc[:,'CSDARange[g/cm2]']]).T

x = np.asarray(data.loc[:,'E[MeV]'])
S = np.asarray(data.loc[:,'S_total[MeV cm2/g]'])*rho_vzduch
R = np.asarray(data.loc[:,'CSDARange[g/cm2]'])/rho_vzduch

S_spline = CubicSpline(x, S,extrapolate=True)
R_spline = CubicSpline(x, R)

xs = np.logspace(-3, 3, num=1000)
#xs=[round(x,3) for x in xs]
#pd.Series(xs).to_csv('file.txt',index=False)

Ss = S_spline(xs)
Rs = R_spline(xs)
plt.loglog(xs,Ss)
plt.grid()
plt.show()

#def find_nearest(Ekin,y,x=x):
#    '''
#    Input:
#        Ekin(float): energie alfa castice, [MeV]
#    Output:
#        dEdx(float): celkova brzdna schopnost, [MeV/cm]
#    '''
#    idx = np.abs((Ekin - x)).argmin()
#    if Ekin<x[idx]:
#        x = x[idx-1:idx+1]
#        y = y[idx-1:idx+1]
#        value = interp1d(x, y, kind='linear')(Ekin)
#    elif Ekin>x[idx]:
#        x = x[idx:idx+2]
#        y = y[idx:idx+2]
#        value = interp1d(x, y, kind='linear')(Ekin)
#    elif Ekin==x[idx]:
#        value = y[idx]
#    return value
#
#def find_nearest_dEdx(Ekin, y=S):
#    return find_nearest(Ekin, y)*rho_vzduch
#
#def find_nearest_R(Ekin, y=R):
#    return find_nearest(Ekin, y)/rho_vzduch

dosah_alfa=np.array([R_spline(E) for E in E_alfa])


#def usla_draha(Ekin, fce=find_nearest_dEdx, dx=1e-1):
#    '''
#    Input:
#        Ekin(float): pocatecni kineticka energie v MeV
#    Output:
#        x(float): delka projite drahy v cm
#    '''
#    x=0         #position in cm
#    dE=0     #energy loss in MeV
##    print('Ekin = '+str(round(Ekin, 3))+' MeV, s = '+str(s)+' cm;')
#    head = 'Delka drahy [cm], E_kin [MeV], S [MeV/cm]'
#    print(head)
#    while Ekin > 0:
#        
#        dE = fce(Ekin)*dx     #units J/m*dx
#        x = x+dx
#        Ekin = Ekin - dE
##        if dE < 0:
##            break
#        print(str(x) + ', ' + str(Ekin) + ', ' + str(dE/dx))
#    return x
#
#s=usla_draha(6)

def zbyla_energie(s, Ekin, fce=S_spline, dx=1e-1):
    '''
    Input:
        s(float): draha alfa castice ve vzduchu nez dorazi k cipu [cm]
    Output:
        E_zbyla(float): zbyla energie alfa castice po jejim pruchodu vzduchu po draze dlouhe s
    '''
    x=0         #position in cm
    dE=0     #energy loss in MeV
#    print('Ekin = '+str(Ekin)+'MeV')
#    print('delka drahy = '+str(s)+' cm')    
    while True:
        x = x+dx
        if x > s:
            break
        dE = fce(Ekin)*dx     #units MeV/cm*dx
        Ekin = Ekin - dE
        if dE < 0:
            print('dE vyslo zaporne! dE = '+str(dE))
            break
        if Ekin < 0:
            print('Ekin vysla zaporne! Ekin = '+str(Ekin))
            break
#    print('------------------------------------------')
#    print('Delka drahy [cm], E_kin [MeV], S [MeV/cm]')
#    print(str(x-dx) + ', ' + str(Ekin) + ', ' + str(dE/dx)+'\n')
    E_zbyla = Ekin
    return E_zbyla

def geometrie_alfa(s):
    '''
    Input:
        s zadavat v cm
    Output:
        podil alfa castic ktere leti smerem na cip (delano pres prostorovy uhel)
    '''
    f=1/2*(1-s/np.sqrt(s**2+r_cip**2))
    return f

def vypocet_alfa():
    '''
    Output:
        I_E_list(ndarray): integral vsech zbylych energii alf jdouci z koule o polomeru rovnu dosahu te dane
                            alfa castice; beru, ze v kazdem elementu objemu vznika jedna alfa, tj. objemova 
                            aktivita je 
    '''
    I_E_list=[]
    for i, Ekin in enumerate(E_alfa):
        R_max = dosah_alfa[i] #[cm]
        fce = lambda s: 4*np.pi*s**2*zbyla_energie(s, Ekin=Ekin)*geometrie_alfa(s)
        I_E = quad(fce, r_pouzdro, R_max+r_pouzdro) #integral vsech energii od povrchu pouzdra do dosahu
        I_E_list.append(I_E[0])
    return np.array(I_E_list)

I_E=vypocet_alfa()
V_slupka=4/3*np.pi*((dosah_alfa+r_pouzdro)**3-r_pouzdro**3)
E_pouzdro=I_E/V_slupka
#zbyla_energie(dosah_alfa[0],E_alfa[0])
#zbyla_energie(dosah_alfa[1],E_alfa[1])
#zbyla_energie(dosah_alfa[2],E_alfa[2])