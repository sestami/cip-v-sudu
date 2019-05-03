"""
Bethe Bloch equation and calculation of energy loss along trajectory

author: Bob Wimmer
date: April 27, 2011
email: wimmer@physik.uni-kiel.de

"""

from math import sqrt, log, pi
from scipy.constants import N_A, epsilon_0, e
import numpy as np
from scipy.integrate import quad
#Define some natural constants

mp = 1.672621637e-27        #kg
me = 9.1093821499999992e-31 #kg
ma = 6.64465620e-27         #kg, mass of the alfa particle
qe = e                      #C
eps0 = epsilon_0            #F/m
c = 299792458               #m/s
rho_vzduch = 1.20479        #kg/m^3


#and nnow some problem specific ones

Ekin = 5.49e6*qe               #6 MeV in J
Ekin_0 = Ekin
m0   = ma                   #mass of projectile
Z1   = 2                    #nuclear charge of projectile (alfa)
EB   = 85.7*qe              #85.7 eV in J, i.e., Mean Excitation Energy of Dry Air
#ne   = 3.3456e29            #electron density of water in m^-3

E_alfa=np.array([5.490, 6.002, 7.689]) #v MeV

def dosah_alfa_vzduch(E):
    '''
    Input:
        E(float): kineticka energie [MeV]
    Output:
        R(float): dosah ve vzduchu [mm]
    '''
    R=3.18*E**(3/2)
    return R

dosah_alfa=dosah_alfa_vzduch(E_alfa)

def ne():
    '''
    electron density in the material
    '''
    atomic_number_list = np.array([6, 7, 8, 18])
    fraction_by_weight = np.array([0.000124, 0.755267, 0.231781, 0.012827])
    Z_effective        = sum(atomic_number_list*fraction_by_weight)
    return N_A*Z_effective*rho_vzduch/28.97 #brano z https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=104, https://en.wikipedia.org/wiki/Molar_mass
    
ne = ne()

#relativisticky pristup
def beta(v):
    return v/c

def gamma(v):               #Lorentz gamma
    return 1./sqrt(1-v*v/c/c)

def v_of_Ekin_m0_relativisticky(Ekin, m0): #invert kinetic energy, E_kin, for speed, v.
    b2 = 1.-1./(1.+Ekin/m0/c/c)**2
    return sqrt(b2)*c

def dEdx_relativisticky(Ekin,Z1=Z1,m0=ma,EB=EB,ne=ne): #Bethe-Bloch equation
    '''
    [J/m]
    '''
    v = v_of_Ekin_m0_relativisticky(Ekin, m0)
    b2 = beta(v)**2
    C = Z1**2*qe**4/4/pi/eps0**2/me
    ln_term = log(2.*me*v**2/EB)
    return C/v**2*ne*(ln_term  - log(1.-b2) - b2)

#nerelativisticky pristup
def v_of_Ekin_m0(Ekin, m0): #invert kinetic energy, E_kin, for speed, v.
    return sqrt(2*Ekin/m0)

def dEdx(Ekin,Z1=Z1,m0=ma,EB=EB,ne=ne): #Bethe-Bloch equation
    '''
    [J/m]
    '''
    v = v_of_Ekin_m0(Ekin, m0)
    ln_term = log(2*me*v**2/EB)
    return Z1**2*qe**4/4/pi/eps0**2/me/v**2*ne*ln_term

#aproximujici pristup
def dEdx_aprox(Ekin,Z1=Z1,m0=ma,EB=EB,ne=ne): #Bethe-Bloch equation
    '''
    [J/m]
    '''
    v = v_of_Ekin_m0(Ekin, m0)
    ln_term = log(2*me*v**2/EB)
    return Z1/me/v**2*rho_vzduch*1/2*ln_term

#Energy loss in first layer
print(str(dEdx_relativisticky(Ekin)/(qe*1e9)) + ' in MeV/mm')
print(str(dEdx_relativisticky(Ekin)/(qe*1e8)) + ' in MeV/cm')

#initialize position, energy loss, and dx
def zbyla_energie(s, Ekin, fce=dEdx_relativisticky, dx=1e-3):
    '''
    Input:
        s(float): draha alfa castice ve vzduchu nez dorazi k cipu [m]
    Output:
        E_zbyla(float): zbyla energie alfa castice po jejim pruchodu vzduchu po draze dlouhe s
    '''
    x=0         #position in m
    dE = 0     #energy loss
    print('Ekin = '+str(round(Ekin/(qe*1e6), 3))+' MeV, s = '+str(s)+' m;')
    head = 'Delka drahy [m], E_kin [MeV], S [MeV/cm]'
    print(head)
    while x < s:
        dE = fce(Ekin)*dx     #units J/m*dx
        x = x+dx
        Ekin = Ekin - dE
        if dE < 0:
            break
    string = str(x) + ', ' + str(Ekin/(qe*1e6)) + ', ' + str(dE/(qe*1.e8)/dx) + '\n'
    print(string)
    E_zbyla = Ekin
    return E_zbyla

def vypocet():
    I_E_list=[]
    for Ekin in E_alfa:
        R = dosah_alfa_vzduch(Ekin)*10**(-3) #[m]
        Ekin = Ekin*10**6*qe #[J]
        fce = lambda s: zbyla_energie(s, Ekin=Ekin)
        I_E = quad(fce, 0.01, R) #integral vsech energii
        I_E_list.append(I_E[0])
    return np.array(I_E_list)

I_E=vypocet()/qe/1e6 