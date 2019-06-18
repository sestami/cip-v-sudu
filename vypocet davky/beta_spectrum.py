# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:06:00 2019

@author: michal.sestak
"""
from math import sqrt, log, pi
#from scipy.constants import N_A, epsilon_0, e
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
import numpy as np
import pandas as pd
#from konstanty import *
#import time
import matplotlib.pyplot as plt
from random import uniform

plt.close('all')
#ENERGIE V MeV!!!!!!!!
df_Pb=pd.read_csv('beta_214Pb.txt', sep='\t', decimal=',')
df_Bi=pd.read_csv('beta_214Bi.txt', sep='\t', decimal=',')

spline_Pb = CubicSpline(df_Pb.iloc[:,0], df_Pb.iloc[:,1])
spline_Bi = CubicSpline(df_Bi.iloc[:,0], df_Bi.iloc[:,1])
def spektrum(spline, E):
    if spline==spline_Pb:
       if E<=1:
           return spline(E)
       elif E>1:
           return 0
    if spline==spline_Bi:
        if E<=3.2:
            return spline(E)
        elif E>3.2:
            return 0

def comulative(spline, E_interest):
    f=lambda E: spektrum(spline, E)
    return quad(f, 0, E_interest)

def graf_spektrum(df, spline):
    fig,ax=plt.subplots()
    E=np.linspace(0, df.iloc[-1,0],num=100)
    hustota=np.array([spektrum(spline, el) for el in E])
    #komulativni=np.array([comulative(Pb_spline, el) for el in E])[:,0]

    ax.plot(E, hustota)
    ax.grid()
    ax.set_xlabel('$E$ [MeV]')
    ax.set_ylabel('zastoupení [-]')
    


def rejection_method(spektrum,M,df,N=10**6):
    #N...number of generated random numbers
    #a,b,M...define the area of generating random numbers for
    #        rejection
    a=0 #spektrum vzdy zacina nulou
    b=df.iloc[-1,0]
    random_numbers=[]
    for i in np.arange(N):
        while True:
            x=uniform(a,b)
            y=uniform(0,M)
            if y<spektrum(x):
                break
        random_numbers.append(x)
    return np.array(random_numbers)

maximum_Pb=0.155
maximum_Bi=0.165

def graf_nahodnaCisla(nahodnaCisla):
    fig, ax=plt.subplots()
    ax.hist(nahodnaCisla, rwidth=0.9, bins=20, label='generated betas')
    ax.set_xlabel('$E$ [MeV]')
    ax.set_ylabel('pocet [-]')
    ax.legend()
    ax.grid()
    
def maximalni_chyba(sigma, N=10**6, epsilon=0.01):
    return sigma/np.sqrt(N*epsilon)

nahodnaCisla_Pb=rejection_method(lambda E: spektrum(spline_Pb, E), maximum_Pb, df_Pb)
#graf_spektrum(df_Pb, spline_Pb)
#graf_nahodnaCisla(nahodnaCisla_Pb)

nahodnaCisla_Bi=rejection_method(lambda E: spektrum(spline_Bi, E), maximum_Bi, df_Bi)
#graf_spektrum(df_Bi, spline_Bi)
#graf_nahodnaCisla(nahodnaCisla_Bi)

E_Pb_mean=nahodnaCisla_Pb.mean()
E_Bi_mean=nahodnaCisla_Bi.mean()

E_Pb_sigma=nahodnaCisla_Pb.std()
E_Bi_sigma=nahodnaCisla_Bi.std()

E_Pb_error=maximalni_chyba(E_Pb_sigma)
E_Bi_error=maximalni_chyba(E_Bi_sigma)


#zaznamy
E_Pb_stredni=0.21915398688762827 #v MeV
E_Bi_stredni=0.6387920731310772  #v MeV