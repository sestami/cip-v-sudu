# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 13:47:59 2019

@author: michal.sestak

TO DO:
    - vztahnout vysledne davky k objemove aktivite a k dobe trvani expozice (v jednotkach za den, resp. za hodinu)
    - davka od dcer posazenych na cipu
    - GAMA: zahrnout soucinitele absorpce
DONE:
    - vztahnuti davek (od alfa, gama) k dobe trvani expozice [den^(-1)]
APROXIMACE:
    ALFY
    ...dopsat
    GAMY
    1) zanedbani zeslabeni v pouzdru cipu
    2) beru cip a pouzdro jako koule danych polomeru
    3) integruji v cylindrickych souradnicich v mezich r\in [r_pouzdro,polomer_sudu] a z\in [r_pouzdro,polomer_sudu],
    tj. zanedbavam prispevky od r\in [0,r_pouzdro] (valec na pouzdrem nahoru a dolu) a od z\in [0,r_pouzdro] (prunik
    valce o polomeru r_pouzdro a vysce taktez r_pouzdro (tento valec predstavuje pouzdro cipu) a valce predstavujici sud);
     TO DO: udpravit rozmery valce na rozmery pouzdra (nebo mozna ne? kvuli konzistentnosti, aby se integrovalo jenom
     tam, kde se uvazuje vznik gam)
    
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, dblquad

#ROZMERY
polomer_sudu=54/2 #v cm
vyska_sudu=83 #v cm

V_cip=(6.727)*(7.217)*0.15*10**(-3) #objem cipu v cm^3
V_sud=np.pi*polomer_sudu**2*vyska_sudu #objem sudu v cm^3

kryti_ze_strany_1=6.51*10**(-1) #tloustka pouzdra na stranach cipu (odpovidajici rozmer je sirka) v cm
kryti_ze_strany_2=6.55*10**(-1) #tloustka pouzdra na stranach cipu (odpovidajici rozmer je delka) v cm
kryti_vertikalne=0.69*10**(-1) #tloustka pouzdra nahore a dole cipu v cm

rho_cip=2.330 #hustota cipu (kremiku) v g/cm^3
m_cip=V_cip*rho_cip*10**(-3) #hmotnost cipu v kg

#hustoty jednotlivych prvku ze kterych je pouzdro v g/cm^3; [O, Si, C, N]
rho_pouzdro_list=np.array([1.33151E-03,2.330,2,1.16528E-03])

#atomove koncentrace jednotlivych prvku v pouzdru; [O, Si, C, N]
c_pouzdro_list=np.array([61.33,23.80,10.36,4.5])/100

#vahove koncentrace jednotlivych prvku v pouzdru; [O, Si, C, N]
w_pouzdro_list=np.array([53.41,36.39,6.77,3.43])/100

#efektivni hustota pouzdra
rho_pouzdro=sum(c_pouzdro_list*rho_pouzdro_list)

#vstupni jednorazova davka
A_in=10**4

#objemová aktivita v sudu
a_metry=A_in/V_sud*10**6 #v Bq/m^3
a=A_in/V_sud #v Bq/cm^3

#doba trvani expozice
T_tri_mesice=3*30*24*60*60 # tri mesice, v s
T_jeden_den=24*60*60 # jeden den, v s

#premenove konstanty
l0=np.log(2)/(3.81*24*60*60)
l1=np.log(2)/(3.1*60)
l2=np.log(2)/(26.8*60)
l3=np.log(2)/(19.9*60)

#rovnovazny faktor
F=0.1

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#DAVKA OD ALF
#sdelena energie cipu od alf z jedne premeny
E_alfa_1=1/2*(5490+6002+7687) #v keV 

#aktivita není konstantni; tento pripad je relevantni, ten spravny
def D_alfa(t):
    '''
    TO DO:
        Braggova krivka pro alfy (dosah alf)
        pridat prispevky dcer posazenych na pouzdru (a vypocitat podil alf ktere do cipu doleti)
    Output:
        davka od alf pri promenne aktivite, v Gy
    '''
    E_alfa=E_alfa_1*a*V_cip*(1-np.exp(-l0*t))/l0*1.6*10**(-16)
    return E_alfa/m_cip

DAlfa=D_alfa(T_jeden_den) #v Gy/(den)

def kontrola_linearity_davkoveho_integralu():
    casy=np.array([i*24*60*60 for i in np.arange(1,3*32,5)])
    D=D_alfa(casy)
    plt.plot(casy/60/60/24,D/a_metry*casy,'x',label='davkovy integral vztahnuty na (bq/m^3)/den')
    plt.xlabel('[dny]')
    plt.ylabel(r'$\left[\frac{Gy\cdot den}{Bq/m^3}\right]$')
    plt.grid()
    plt.legend()
    plt.show()

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#DAVKA OD GAMY

#korekce na geometrii, cip bereme jako kouli o polomeru 3 mm, L je vzdalenost od cipu
r_cip=0.3 # aproximace rozmeru cipu jako kouli, toto je jeji polomer v cm
r_pouzdro=1 #v cm
def geometrie_gama(r,z,my):
    '''
    Input:
        r,z zadavat v cm
        my je lin. soucinitel zeslabeni (lze zadavat jako (hmotnostni soucinitel zeslabeni*hustota))
    '''
    L=np.sqrt(r**2+z**2)
    f=np.pi*(1-L/np.sqrt(L**2+r_cip**2))*np.exp(-my*L)*r
    return f

rho_air=1 #v g/cm^3
#prvni slozka je energie gamy (v keV), druha slozka je soucinitel zeslabeni ve vzduchu (linearni=hmotnostni), treti slozka je vytezek
gama_511=[511,8.712E-02,0.08] #od Rn
gama_352=[352,9.800E-02,0.38]
gama_300=[300,1.067E-01,0.27]
gama_609=[609,8.055E-02,0.46]
gama_1180=[1180,5.687E-02,0.21]
gama_1764=[1764,4.80E-02,0.15]
gama_2204=[2204,4.447E-02,0.05]

E_list=np.array([gama_511[0],gama_352[0],gama_300[0],gama_609[0],gama_1180[0],gama_1764[0],gama_2204[0]])
mu_list=np.array([gama_511[1],gama_352[1],gama_300[1],gama_609[1],gama_1180[1],gama_1764[1],gama_2204[1]])
Y_list=np.array([gama_511[2],gama_352[2],gama_300[2],gama_609[2],gama_1180[2],gama_1764[2],gama_2204[2]])

absCoeff_list=np.array([2.971E-02,2.968E-02,2.932E-02,2.951E-02,2.7E-02,2.445E-02,2.3E-02])

#integral pres objem sudu, kdyby v sudu byl 1 Bq/m^3
def I_vzduch(mu): # I jako integral
    f_mu = lambda r, z: geometrie_gama(r,z,mu)
    I_h=dblquad(f_mu, r_pouzdro, polomer_sudu, lambda x: r_pouzdro, lambda x: vyska_sudu/2)
    return a*2*I_h[0]

#I_list=np.array([I_vzduch(mu) for mu in mu_list]) #aktivita v cipu pro jednotkovou objemovou aktivitu
#A_list=np.array([a*Y*I_list[i] if i==0 else F*a*Y*I_list[i] for i,Y in enumerate(Y_list)]) #a musi byt v Bq/cm^3 !!!; zahrnuje obj. aktivitu, vytezek, integral pres objem sudu

# ZANEDBAVAM ZESLABENI V POUZDRU!!!
# A_sum=np.sum(A_list)
#E_jednotliveGamy=E_list*A_list

# PRISPEVEK OD STEN SUDU
S_sud=2*np.pi*polomer_sudu*(vyska_sudu+polomer_sudu) #v cm^2
S_pouzdro=2*((6.554+6.727+6.815)*(6.515+7.217+6.518)+(6.554+6.727+6.815)*(0.69+0.15)+(6.515+7.217+6.518)*(0.69+0.15))*10**(-2) #v cm^2

#geometrie_gama_plá=lambda z,my: 0.5*(1-np.sqrt(polomer_sudu**2+z**2)/np.sqrt(polomer_sudu**2+z**2+r**2))*np.exp(-my*np.sqrt(polomer_sudu**2+z**2))
def I_sud(mu):
    '''
    Pozor, zde se musi obj. aktivity prepocitavat na plosne, pocitaji se tady totiz plosne integraly ()
    
    1 -> prispevky od plaste
    2 -> prispevky od podstav
    '''
    f_1_mu = lambda z: geometrie_gama(polomer_sudu,z,mu)
    I_1_h=quad(f_1_mu, 0, vyska_sudu/2) #pocitam pouze prispevky od horni poloviny sudu
    I_1=2*np.array(I_1_h) #zahrnuje plosnou aktivitu na plasti
    
    f_2_mu = lambda r: geometrie_gama(r,vyska_sudu/2,mu)
    I_2_h=quad(f_2_mu, 0, polomer_sudu)
    I_2=2*np.array(I_2_h) #dve podstavy; zahrnuje plosnou aktivitu na podstavach
    
    I=a*V_sud/S_sud*(I_1+I_2)
    return S_sud/(S_sud+S_pouzdro)*I[0]

def I_pouzdro(mu):
    '''
    je to takove od oka...
    '''
    return S_pouzdro/(S_pouzdro+S_sud)*(a*V_sud)*1/4

#prispevek od radonu
I_Rn=a*I_vzduch(mu_list[0])
A_Rn=Y_list[0]*I_Rn
A_Rn_deponovane=A_Rn*(1-np.exp(-absCoeff_list[0]*rho_cip*2*r_cip))
E_Rn=A_Rn_deponovane*E_list[0]

def gama_zpracovani(I_fce, F):
    '''
    Output:
        energie deponovana od dcer v cipu v keV
    '''
    I_dcery=np.array([I_fce(mu) for mu in mu_list[1:]]) #aktivita v cipu pro jednotkovou objemovou aktivitu
    A_dcery=np.array([F*Y*I_dcery[i] for i,Y in enumerate(Y_list[1:])]) #a musi byt v Bq/cm^3 !!!; zahrnuje vytezek, integral pres objem sudu
    A_deponovane=[A_dcery[0]*(1-np.exp(-absCoeff*rho_cip*2*r_cip)) for i,absCoeff in enumerate(absCoeff_list[1:])]
    #ZANEDBAVAM ZESLABENI V POUZDRU!!!
    #A_sum=np.sum(A_list)
    return A_deponovane*E_list[1:]

def D_gama(t):
    '''
    TO DO:
        zeslabeni v pouzdru?
    DONE:
        koeficient absorpce v cipu
    Predpoklady:
        aktivita není konstantni
        vsechny gamy dosle do cipu se zde deponuji
    Input:
        doba trvani pro urceni casoveho integralu davky
    Output:
        davka od gam pri promenne aktivite, v J/kg
    '''
    E_vzduch=gama_zpracovani(I_vzduch,F)
    E_sud=gama_zpracovani(I_sud,1-F)
    E_pouzdro=gama_zpracovani(I_pouzdro,1-F)
    
    E_celkove=(E_Rn+np.sum(E_vzduch+E_sud+E_pouzdro))*1.6*10**(-16)*(1-np.exp(-l0*t))/l0 #zahrnuti casoveho integralu
    return E_celkove/m_cip

DGama=D_gama(T_jeden_den)

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#DAVKA OD BETY

#vstupni udaje
E_endpoints=np.array([672,729,1542,3272])
E_middle=1/2*E_endpoints #overit
Y_bety=np.array([50,43,35,18.2])/100

rho_vzduch=0.0012 # v g/cm3
#dalsi hustoty jsou: rho_cip, rho_pouzdro