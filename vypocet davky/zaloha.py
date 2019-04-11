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
