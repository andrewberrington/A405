
# coding: utf-8

# In[2]:

import numpy as np
import pandas as pd
import h5py
from pprint import pformat
from matplotlib import pyplot as plt
from scipy import optimize


# In[17]:

import a405thermo.thermlib as therm
reload(therm)

def rootfinder_thetaes(temp_guess2, pressure2, theta_es):
    '''Rootfinder for theta_es
    inputs: temp_guess (K), pressure (Pa)
    outputs: difference between temp_guess and temp'''
    
    R_v = 461.5 # J kg^-1 K^-1
    T_0 = 273.15 # K
    C_pv = 1870 # J kg^-1 K^-1
    l_v0 = 2.501e6 # J kg^-1
    e_s0 = 6.11e2 # Pa
    Tp = 273.16 # K
    C_l = 4187 # J kg^-1 K^-1
    C_pd = 1004 # J kg^-1 K^-1
    P_0 = 1e5 # Pa
    R_d = 287.05 # J kg^-1 K^-1
    
    phi_0 = l_v0/Tp
    
    l_v = (C_pv-C_l)*(temp_guess2-T_0)+l_v0 # l_v function of temp
    
    phi_l = C_l*np.log(temp_guess2/Tp)
    
    e_s = e_s0*np.exp((C_pv*np.log(temp_guess2/Tp)+phi_0-C_l*np.log(temp_guess2/Tp)-(l_v/temp_guess2))/R_v) #e_s function of temp
    
    r_s2 = ((e_s*0.622)/(pressure2-e_s))
    
    c_p = (C_pd+(r_s2*C_l))
    
    theta = temp_guess2*(P_0/pressure2)**(R_d/c_p)
    
    return theta_es - theta*((l_v*r_s)/(c_p*temp_guess2))

def calc_thetaes(pressure2,theta_es):
    '''Calculation for theta_es using rootfinder
    inputs: pressure (Pa)
    outputs: theta_es (K)'''
    
    temp2 = optimize.zeros.brentq(rootfinder_thetaes, 1, 350, args=(pressure2,theta_es))
    return temp2


# In[18]:

pressarray2 = np.arange(1000,390,-10)*100 # for Pa
theta_es = 273.15 # K
thetemps2 = []

for the_pressure in pressarray2:
    thetemps2.append(calc_thetaes(the_pressure, theta_es))
    
thetemps2 = np.array(thetemps2)-273.15
thepress2 = pressarray2/100

print(thetemps2)

