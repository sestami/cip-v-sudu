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
ma = 6.64465620e-27         #kg,the mass of the alfa particle
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


#------------------------------------------------------------------------------------------------
#ALFA
E0_alfa=np.array([5.489, 6.002, 7.689]) #v MeV
data_alfa=pd.read_csv('ASTAR vzduch.txt', skiprows = 8, sep="\t")

E_alfa_data = np.asarray(data_alfa.loc[:,'E[MeV]'])
S_alfa_data = np.asarray(data_alfa.loc[:,'S_total[MeV cm2/g]'])*rho_vzduch
R_alfa_data = np.asarray(data_alfa.loc[:,'CSDARange[g/cm2]'])/rho_vzduch

S_alfa_spline = CubicSpline(E_alfa_data, S_alfa_data,extrapolate=True)
R_alfa_spline = CubicSpline(E_alfa_data, R_alfa_data)

def S_alfa(E):
    if E>0:
        return S_alfa_spline(E)
    elif E<=0:
        return 0

Es = np.logspace(-3, 3, num=1000)
Es = np.linspace(0, 10**3, num=10000)
Ss = S_alfa_spline(Es)
Rs = R_alfa_spline(Es)
#plt.plot(Es,Ss)
#plt.grid()
#plt.show()
dosah_alfa=np.array([R_alfa_spline(E) for E in E0_alfa])

#------------------------------------------------------------------------------------------------
#BETA
E0_beta=np.array([0.219, 0.638]) #v MeV
data_beta=pd.read_csv('ESTAR vzduch.txt', skiprows = 8, sep="\t")
data_beta_dosah=pd.read_csv('ESTAR vzduch dosah.txt', skiprows = 8, sep="\t")

E_beta_data = np.asarray(data_beta.loc[:,'E[MeV]'])
S_beta_data = np.asarray(data_beta.loc[:,'S_total[MeV cm2/g]'])*rho_vzduch
R_beta_data = np.asarray(data_beta_dosah.loc[:,'CSDARange[g/cm2]'])/rho_vzduch

S_beta_spline = CubicSpline(E_beta_data, S_beta_data,extrapolate=True)
R_beta_spline = CubicSpline(data_beta_dosah.loc[:,'E[MeV]'], R_beta_data)

def S_beta(E):
    if E>0:
        return S_beta_spline(E)
    elif E<=0:
        return 0

dosah_beta=np.array([R_beta_spline(E) for E in E0_beta])
#------------------------------------------------------------------------------------------------

def zbyla_energie(s, Ekin, fce=S_alfa_spline, dx=1e-1):
    '''
    PRO JAKOUKOLIV NABITOU CASTICI
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
            print('castice prosla celou drahu!')
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

def geometrie_nabiteCastice(s):
    '''
    PRO JAKOUKOLIV NABITOU CASTICI
    Input:
        s: vzdalenost od cipu, v niz alfa castice vznikla [cm]
    Output:
        podil alfa castic ktere leti smerem na cip (delano pres prostorovy uhel)
    '''
    f=1/2*(1-s/np.sqrt(s**2+r_cip**2))
    return f
#------------------------------------------------------------------------------------------------
def vypocet(E0=E0_alfa, dosah=dosah_alfa):
    '''
    Output:
        I_E_list(ndarray): I_E je stredni energie, ktera zbyde alfa castici o dane
                           energii po dojiti k cipu, prenasobena objemem koule,
                           z niz emitovane alfa castice mohou k cipu dojit 
                           (tj. o polomeru rovnem dosahu te dane alfa castice ve vzduchu)
    '''
    I_E_list=[]
    for i, Ekin0 in enumerate(E0):
        R_max = dosah[i] #[cm]
        fce = lambda s: 4*np.pi*s**2*zbyla_energie(s, Ekin0, fce=S_alfa)*geometrie_nabiteCastice(s)
        I_E = quad(fce, r_pouzdro, R_max+r_pouzdro) #integral vsech energii od povrchu pouzdra do dosahu
        I_E_list.append(I_E[0])
    return np.array(I_E_list)
#------------------------------------------------------------------------------------------------

I_E_beta=vypocet(E0=E0_beta, dosah=dosah_beta)
I_E_alfa=vypocet()

#TO DO: pozor, beta ma vetsi dosah nez jsou rozmery sudu!!!!
#       I_E_beta vychazi zaporne!!! WTF

#I_E_alfa=vypocet_alfa()
#V_slupka=4/3*np.pi*((dosah_alfa+r_pouzdro)**3-r_pouzdro**3)
#E_pouzdro=I_E/V_slupka
#zbyla_energie(dosah_alfa[0],E_alfa[0])
#zbyla_energie(dosah_alfa[1],E_alfa[1])
#zbyla_energie(dosah_alfa[2],E_alfa[2])