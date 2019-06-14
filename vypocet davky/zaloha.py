# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 11:17:38 2019

@author: michal.sestak
"""

#hustoty jednotlivych prvku ze kterych je pouzdro
rho_O=1.33151E-03
rho_Si=2.330
rho_C=2
rho_N=1.16528E-03

#atomove koncentrace jednotlivych prvku v pouzdru
c_O=61.33/100
c_Si=23.80/100
c_C=10.36/100
c_N=4.5/100

#vahove koncentrace jednotlivych prvku v pouzdru
w_O=53.41/100
w_Si=36.39/100
w_C=6.77/100
w_N=3.43/100

#davka od alf, bere se jen objem cipu, je uvazovano, ze aktivita je po celou dobu konstantni
#def D_alfa_konst(t):
#    '''
#    Output:
#        davka od alf pri konstantni aktivite, v Gy
#    NEPOUZIVAT!!!
#    '''
#    E_alfa=E_alfa_1*a*V_cip*t*1.6*10**(-16)
#    return E_alfa/m_cip

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
#def dosah_alfa_libovolnaLatka(E,rho,A_r):
#    '''
#    Input:
#        rho(float): hustota latky v g/cm^3
#        A_r(float): relativni atomova hmotnost
#    Output:
#        R_X(float): dosah v dane latce, urcen s relativni chybou +-15 %
#    '''
#    return 0.3*dosah_alfa_vzduch(E)/rho*np.sqrt(A_r)
#
##Puvodni pristup
#E_alfa_1=1/2*(5.490+6.002+7.687) #v MeV
#def D_alfa_puvodni(t):
#    '''
#    TO DO:
#        Braggova krivka pro alfy (dosah alf)
#        pridat prispevky dcer posazenych na pouzdru (a vypocitat podil alf ktere do cipu doleti)
#    Output:
#        davka od alf pri promenne aktivite, v Gy
#    '''
#    E_alfa=E_alfa_1*a*V_cip*(1-np.exp(-l0*t))/l0*1.6*10**(-13)
#    return E_alfa/m_cip



#def vyvoj_aktivity_a_integralu_aktivity():
#    '''
#    - kdyz se to da jenom v ramci jednoho dne, tak to linearni vubec neni
#    - v ramci tri mesicu ano
#    '''
#    casy=np.array([i*24*60*60 for i in np.arange(1,10*3.82,2)])
##    D=D_alfa(casy)
#    plt.figure()
#    plt.plot(casy/60/60/24,a*V_sud*(1-np.exp(-l0*casy))/l0,'x',label='integral aktivity v sudu')
#    plt.xlabel('[dny]')
#    plt.ylabel(r'$\int_0^T A dt$ [Bq$\cdot$s]')
#    plt.grid()
#    plt.legend()
#    plt.show()
#    
#    plt.figure()
#    plt.plot(casy/60/60/24,a*V_sud*np.exp(-l0*casy),'x',label='vyvoj aktivity v sudu')
#    plt.xlabel('[dny]')
#    plt.ylabel(r'$A$ [Bq]')
#    plt.grid()
#    plt.legend()
#    plt.show()
    

#aproximace funkce pro vypocteni geometrickeho integralu
#L=np.linspace(3,400,num=1000)
#
#plt.plot(L,korekce_geometrie(L),label='$korekce$')
#plt.plot(L,1/L**(1.7),label='$1/L^{1,7}$')
#plt.plot(L,korekce_geometrie(L)-1/L**(1.7),label='$korekce-1/L^{1.7}$')
#plt.grid()
#plt.xlabel('L')
#plt.ylabel('korekce na geometrii (bere se v uvahu pouze uzky svazek)')
#plt.legend()
#plt.show()

#def geometrie_gama_plast(z,my):
#    '''
#    Input:
#        z zadavat v cm
#        my je lin. soucinitel zeslabeni (lze zadavat jako (hmotnostni soucinitel zeslabeni*hustota))
#    '''
#    L=np.sqrt(polomer_sudu**2+z**2)
#    f=np.pi*polomer_sudu*(1-L/np.sqrt(L**2+r_cip**2))*np.exp(-my*L)
#    return f
#
#def geometrie_gama_podstavy(r,my):
#    '''
#    Input:
#        r zadavat v cm
#        my je lin. soucinitel zeslabeni (lze zadavat jako (hmotnostni soucinitel zeslabeni*hustota))
#    '''
#    L=np.sqrt(r**2+(vyska_sudu/2)**2)
#    f=np.pi*(1-L/np.sqrt(L**2+r_cip**2))*np.exp(-my*L)
#    return f


#I_list=np.array([I_vzduch(mu) for mu in mu_list]) #aktivita v cipu pro jednotkovou objemovou aktivitu
#A_list=np.array([a*Y*I_list[i] if i==0 else F*a*Y*I_list[i] for i,Y in enumerate(Y_list)]) #a musi byt v Bq/cm^3 !!!; zahrnuje obj. aktivitu, vytezek, integral pres objem sudu

# ZANEDBAVAM ZESLABENI V POUZDRU!!!
# A_sum=np.sum(A_list)
#E_jednotliveGamy=E_list*A_list

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#PRISPEVEK OD BETY

#def dosah_beta_1(E_max,rho):
#    '''
#    plati pouze pro hlinik; pro vzduch pouzij NIST XCOM
#    '''
#    R=0.412*E_max**(1.265-0.0954*(np.log(E_max)))
#    return R/rho
#
#def dosah_beta_2(E_max,rho):
#    '''
#    plati pouze pro hlinik
#    '''
#    R=0.543*E_max-0.160
#    return R/rho

#def geometrie_beta(r,z,mu_beta):
#    '''
#    - Zatim se nelisi geometrie_beta, ale predpokladam, ze se lisit bude
#    - tak 
#    Input:
#        r,z zadavat v cm
#        my_beta je lin. soucinitel zeslabeni (lze zadavat jako (hmotnostni soucinitel zeslabeni*hustota))
#    '''
#    L=np.sqrt(r**2+z**2)
#    f=np.pi*(1-L/np.sqrt(L**2+r_cip**2))*np.exp(-mu_beta*L)*r
#    return f




bethe.py
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

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#DAVKA OD BETY
'''Jak na to:
    - nejprve vypocitat dosahy z empirickych vzorcu -> NOPE, dosahy urciny z tabulek NIST XCOM; radove odpovidaji dosahum
                                                       urcenym z techto empirickych vzorcu (ktere plati pro hlinik)
    - budu to delat pres soucinitele absorpce
    - budu uvazovat, ze beta zareni dosle do cipu ma energii rovnou jedne ctvrtine jeho E_max (to by melo davku znacne
      nadhodnocovat)
 
'''

#vstupni udaje
E_endpoints=np.array([672,729,1520,3272]) #v keV
E_middle=1/2*E_endpoints #overit
Y_list_beta=np.array([49,42,35,18.2])/100 #resit toto je nesmysl !!!

def soucinitel_absorpce(E_max, rho):
    return 22*E_max**(-4/3)/rho

#hmotnostni soucinitele absorpce
mu_list_beta=np.array([soucinitel_absorpce(E_max,rho_vzduch) for E_max in E_endpoints])


def beta_zpracovani(I_fce, F):
    '''
    Output:
        energie deponovana od beta zareni dcer v cipu v keV
    Predpoklady:
        dva zpusoby:
            1) od kazde castice se deponuje energie rovna jedne desetine maximalni energie daneho beta zareni
            2) jako u gamy, energie deponovana od jedne absorbovane castice se bere jedna desetina max. 
               energie dane bety
    '''
    I_dcery_po_vzduchu=np.array([I_fce(mu*rho_vzduch) for mu in mu_list_beta]) #aktivita dosla na povrch pouzdra
    I_dcery=np.array([I_0*np.exp(-mu_list_beta[i]*rho_pouzdro) for i,I_0 in enumerate(I_dcery_po_vzduchu)]) #aktivita v cipu pro jednotkovou objemovou aktivitu
    A_dcery=np.array([F*Y*I_dcery[i] for i,Y in enumerate(Y_list_beta)]) #a musi byt v Bq/cm^3 !!!; zahrnuje vytezek, integral pres objem sudu
    #1) to, co se dostalo do cipu, se plne deponuje
    A_deponovane=A_dcery 
    #2) viz v popisu funkce
#    A_deponovane=np.array([A_dcery[i]*(1-np.exp(-absCoeff*rho_cip*2*r_cip)) for i,absCoeff in enumerate(mu_list_beta)])
    return A_deponovane*1/10*E_endpoints

def D_beta(t):
    '''
    TO DO:
        
    DONE:
        
    Predpoklady:
        viz fce beta_zpracovani
    Input:
        doba trvani pro urceni casoveho integralu davky
    Output:
        davka od gam pri promenne aktivite, v J/kg
    '''
    E_vzduch=beta_zpracovani(I_vzduch,F)
    E_sud=beta_zpracovani(I_sud,1-F)
    E_pouzdro=beta_zpracovani(I_pouzdro,1-F)
    
    E_celkove=np.sum(E_vzduch+E_sud+E_pouzdro)*1.6*10**(-16)*(1-np.exp(-l0*t))/l0 #zahrnuti casoveho integralu
    return E_celkove/m_cip