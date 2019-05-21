# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:06:00 2019

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

df_Pb=pd.read_csv('beta_214Pb.txt', sep='\t', decimal=',')

Pb_spline = CubicSpline(df_Pb.iloc[:,0], df_Pb.iloc[:,1])

def Pb_comulative(E):
    return quad(Pb_spline,0, E)

E=np.linspace(0, df_Pb.iloc[-1,0],num=100)
plt.plot(E, Pb_spline(E))
plt.grid()
#moc to nefunguje