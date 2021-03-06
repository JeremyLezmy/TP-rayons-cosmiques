# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 16:32:25 2019

@author: jerem
"""

from pylab import*
from numpy import polynomial as P
import numpy as np
import matplotlib.pyplot as plt
#############Lecture des fichiers#############,r"D:\Téléchargement\tp2019-1",r"D:\Téléchargement\tp2019-2",r"D:\Téléchargement\tp2019-3"
liste_fichiers = [r"D:\Téléchargement\tp2019pb2"]
liste_size = [2001]
tt_masses = []
delta_mu=[]
tab_dE1=[]
beta_gamma=[]
vtest=[]
for rr in range(0, len(liste_fichiers)):
    nom_fichier = liste_fichiers[rr]
    nb_event = liste_size[rr]

    fd=open(nom_fichier, 'r')
    all_lignes = fd.readlines()
    c = 2.99792456*10**8


    nb = 0
    for l in range(0, int(len(all_lignes)/nb_event)):
         #initialisation
         lignes = all_lignes[l*nb_event:l*nb_event+nb_event]
         temps = []
         col1 = []
         col2 = []
         col3 = []
         col4 = []
         for i in range(0,nb_event):
             temps.append(i*2.*10.**-10)
             tab = lignes[i].split('\t')
             col1.append(float(tab[1]))
             col2.append(float(tab[2]))
             col3.append(float(tab[3]))
             col4.append(float(tab[4]))
    
         tab_mini=[]
         tab_largeur = []
         tt=[]
         vv=[]
         ############Calcul du temps moyen#################
         for col in [col1, col2, col3]:
             indice_mini=0
             for j in range(1,len(col)):
                 if col[j]<col[indice_mini]:
                     indice_mini=j
             moyenne = 0.
             if nom_fichier==r"D:\Téléchargement\tp2019-4":
                 taille=len(col)-2650
             else:
                 taille = len(col)-650                                  #####ATTENTION MODIF
             for j in range(0, taille):
                 moyenne=moyenne+col[j]
             moyenne =moyenne/taille
             v = min(indice_mini, abs(indice_mini-len(col)-1))-1
             
             t_m =[]
             v_m =[]
             
             for k in range(0, len(col)):
                 t_m.append(temps[k])
                 v_m.append(col[k]-moyenne)
             vtest.append(v_m)        
             
             tt.append(t_m[v_m.index(min(v_m))])
             vv.append(min(v_m))
            
             ##############Calcul des vitesses############
         if tt[0]<tt[1] and tt[1]<tt[2] :
             beta1 = (22.75*10**(-2))/(1.*c*(tt[1]-tt[0]))
             beta2 = (5.5*10**(-2) )/(1.*c*(tt[2]-tt[1]-(17.5*10**(-2)/(beta1*c))))                   ###ATTENTION MODIF
             beta3 = beta1
        
        
         
         #############Condition pour avoir des evenements coherents##########
         
             if beta2>0.  and beta2<beta3 and beta1<1.:
                 #########Calcul des gamma##########
                 gamma1 = 1./np.sqrt((1.-beta1**2.))
                 gamma2 = 1./np.sqrt((1.-beta2**2.))
                 mec2 = 0.511 #Mev
                 Pi = np.pi
                 dx = 1.5#cm
                 Z = 82.
                 A = 207.2 #g/mol
                 z = -1.
                 M= 105.66
                 I = 823.*10**(-6)
                 K= 0.307075
                 #############calcul de Tmax et de la perte d'energie avec bethebloch############
    
                 T_max=(2.*mec2*(beta1**2)*(gamma1**2))/(1.+(2*gamma1*mec2/M)+(mec2/M)**2)
                 P1=(K*(z**2)*Z)/(A*(beta1**2))
                
                 dE1 =11.35*dx*P1*(0.5*np.log(abs((2*mec2*beta1**2*gamma1**2*T_max)/I**2))-beta1**2)
                 
                 dE1=6.24*10**12*dx*10**(-2)*(4*np.pi/(9.11*10**(-31)*9*10**16))*(Z*11.35*10**3*6.022*10**(23)/(A*10**(-3)*beta1**2))*((((1.6*10**-19)**2)/(4*np.pi*8.854*10**-12))**2)*(np.log(abs((2*mec2*beta1**2)/(I*(1-beta1**2))))-beta1**2)
                 #•print(dE1)
                 beta_gamma.append(beta1*gamma1)
                
                 
                 ###########Calcul des masses###########
                 mc2 = (dE1)/(gamma1-gamma2)
                 print(l, " Masse : ", mc2)
                 if mc2>0:
                     tt_masses.append(mc2)
                 ############Calcul des erreurs############
                 delta_x1=0.2
                 delta_x2=0.2
                 delta_t1=1*10**-10
                 delta_t2=1*10**-10
                 delta_dx = 0.2
                 d1 = 22.75
                 d2 = 17.5
                 delta_beta1=(beta1*((delta_x1/d1)+(2*delta_t1/(tt[1]-tt[0]))))
                 delta_beta2=(beta2*((delta_x2/d2)+(2*delta_t2/(tt[2]-tt[1]-(11.75*10**(-2)/(beta1*c))))))       #######ATTENTION MODIF
    
                 delta_gamma1=(gamma1**3*beta1*delta_beta1)
                
                 delta_gamma2=(gamma2**3*beta2*delta_beta2)
                
                 f1=2.*delta_beta1/beta1
                 f2=2.*delta_gamma1/gamma1
                 f3=2.*mec2*delta_gamma1/M
                 f4=1.+2.*gamma1*(mec2/M)+(mec2/M)**2
                 delta_T_max = (T_max*(f1+f2-(f3/f4)))
                
                 f6=K*z**2*Z/(A*beta1**3)
                 f7=(np.log(2.*mec2*beta1**2*gamma1**2*T_max/I**2)+1.)*delta_beta1
                 f8=K*z**2*Z*delta_gamma1/(A*beta1**2*gamma1)
                 f9=K*z**2*Z*delta_T_max/(2.*A*beta1**2*T_max)
                
                 #delta_dE1 = f6*f7+f8+f9
                 delta_dE1 = dE1 * (delta_dx/dx) + 2*(delta_beta1/beta1) * abs(-dE1 + 6.24*10**12*dx*10**(-2)*(4*np.pi/(9.11*10**(-31)*9*10**16))*(Z*11.35*10**3*6.022*10**(23)/(A*10**(-3)*beta1**2))*((((1.6*10**-19)**2)/(4*np.pi*8.854*10**-12))**2)*((gamma1/beta1)**2 -1))
                 
                 #print("beta1=",beta1,"beta2=",beta2,"gamma1=",gamma1,"gamma2=",gamma2,"dE1=",dE1,"T_max=",T_max,"masse=",mc2)
                
                
                 #delta_m = mc2 *(delta_dE1/dE1 + 10**-1*delta_dx/dx+(delta_gamma1+delta_gamma2)/(gamma1-gamma2))
                 delta_m = delta_dE1/(gamma1 -gamma2) + dE1*((delta_gamma1+delta_gamma2)/(gamma1-gamma2)**2)
                 delta_mu.append(delta_m)
                 #print("delta_beta1=",delta_beta1,"delta_beta2=",delta_beta2,"delta_gamma1=",delta_gamma1,"delta_gamma2=",delta_gamma2,"delta_dE1=",delta_dE1,"delta_T_max=",delta_T_max,"delta_m=",delta_m)
                 nb = nb+ 1
    
                
                
