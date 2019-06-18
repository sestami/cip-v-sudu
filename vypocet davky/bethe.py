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
import logging


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
E0_beta=np.array([0.220, 0.639]) #v MeV
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

def faktor_pruniku(s):
    if s<=polomer_sudu:
        return 1
    elif s**2>polomer_sudu**2+(vyska_sudu/2)**2: #uvazovana kulova slupka je uplne mimo sud
        print("Rozmer vetsi nez rozmery sudu, zadny prispevek.")
        return 0
    elif s>polomer_sudu and s<=vyska_sudu/2: #uvazovana kulova slupka je castecne mimo sud
        h=s-polomer_sudu
        faktor_pruniku=1-h/s
        print("polomer slupky je vetsi nez polomer sudu, faktor pruniku = "+str(faktor_pruniku))
        return faktor_pruniku
    elif s>polomer_sudu and s>vyska_sudu/2:
        #DORESIT!!!!
        h1=s-polomer_sudu
        h2=s-vyska_sudu/2
        faktor_pruniku=1-h1/s-h2/s
        print("polomer slupky je vetsi nez polomer sudu a nez polovina vysky sudu, faktor pruniku = "+str(faktor_pruniku))
        if faktor_pruniku<0:
            logging.error("faktoru pruniku je mensi nez nula!!!")
            return False
        return faktor_pruniku
    else:
        logging.error("nejaka chyba nebo nedotazenost ve funkci 'faktor_pruniku'!")
        return False
    
def zbyla_energie(s, Ekin, fce, dx=1e-1):
    '''
    PRO JAKOUKOLIV NABITOU CASTICI
    Input:
        s(float): draha alfa castice ve vzduchu nez dorazi k cipu [cm]
        Ekin(float): energie dane castice [MeV]
    Output:
        E_zbyla(float): zbyla energie alfa castice po jejim pruchodu vzduchu po draze dlouhe s [MeV]
    '''
    print()
#    f=open('bethe_output.txt','w+')
#    f.write('x[cm]  Ekin[MeV]  S[MeV/cm]\n')
    x=0         #position in cm
    dE=0     #energy loss in MeV
#    print('Ekin = '+str(Ekin)+'MeV')
#    print('delka drahy = '+str(s)+' cm')
    faktorPruniku=faktor_pruniku(s)
    if faktorPruniku==0:
        return 0
    while True:
        x = x+dx
        if x > s:
            print('Castice prosla celou drahu!')
            break
        dE = fce(Ekin)*dx     #units MeV/cm*dx
        Ekin = Ekin - dE
        if Ekin < 0:
            print('Ekin vysla zaporne! Ekin = '+str(Ekin)+' Nastavuji Ekin = 0')
            Ekin=0
            break
        if dE < 0:
            print('dE vyslo zaporne! dE = '+str(dE))
            break
#        f.write(str(x)+'  '+str(Ekin)+'  '+str(dE/dx)+'\n')
#    print('------------------------------------------')
#    print('Delka drahy [cm], E_kin [MeV], S [MeV/cm]')
#    print(str(x-dx) + ', ' + str(Ekin) + ', ' + str(dE/dx)+'\n')
#    f.close()
    E_zbyla = Ekin
    return faktorPruniku*E_zbyla

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
def vypocet(E0=E0_alfa, dosah=dosah_alfa, dx=1e-1):
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
        #!!!!! (fce brzdne schopnosti se zadava manualne sem! lepsi implementace hazela chybu)
        fce = lambda s: 4*np.pi*s**2*zbyla_energie(s, Ekin0, S_beta, dx=dx)*geometrie_nabiteCastice(s)
        #!!!!!
        I_E = quad(fce, r_pouzdro, R_max+r_pouzdro) #integral vsech energii od povrchu pouzdra do dosahu
        I_E_list.append(I_E)
    return np.array(I_E_list)
#------------------------------------------------------------------------------------------------

I_E_beta=vypocet(E0=E0_beta, dosah=dosah_beta)
#I_E_alfa=vypocet()
I_E_beta_zaznam=np.array([[1.72585526e+00, 2.76434966e-04], [7.05350163e+00, 1.85648118e-02]])
E_zbyla=zbyla_energie(dosah_beta[1], E0_beta[1], fce=S_beta)

#TO DO: pozor, beta ma vetsi dosah nez jsou rozmery sudu!!!! (zbyva udelat pripad r_sf>h/2, viz radek 113)
#       vypocet I_E_beta trva hrozne dlouho
#       POZOR!! zdali se pouziva S_alfa nebo S_beta se rozhoduje uvnitr fce vypocet, na radku 148!!! (jinak mi to
#            hazi chybu "Kernel died, restarting")

#I_E_alfa=vypocet_alfa()
#V_slupka=4/3*np.pi*((dosah_alfa+r_pouzdro)**3-r_pouzdro**3)
#E_pouzdro=I_E/V_slupka
#zbyla_energie(dosah_alfa[0],E_alfa[0])
#zbyla_energie(dosah_alfa[1],E_alfa[1])
#zbyla_energie(dosah_alfa[2],E_alfa[2])