# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:29:19 2019

@author: WangPeide
"""
from math import log10, sqrt
# Input parameters 
R_m = 360   # tensile strength
R_p02 = 220 # yield strength
R_z = 125    # surface roughness
alph_k = 1  # stress concentration factor
n = 1       # notch sensitivity caused by stress gradient influence and localized plastic deformation at the notch base
mechanical_method = 1 # 1 for normative, 2 for test
MAT = 1     # 1-cast iron 2-cast steel
R = -1      # stress ratio
j = 3       # quality level for component
j_0 = 1     # constant for material and test method
S_pu = 2/3  # survival probability
t = 120     # wall thickness
gamma_M = 1.265 # partial safety factor for material
sign = 0    # thickness correction. 0-No 1-Yes
# calculation of SN curve
# sigm_b
if mechanical_method == 1:
    sigm_b = 1.06*R_m
    simg_s = 1.106*R_p02
else:
    sigm_b = R_m
    sigm_s = R_p02
# surface roughness factor
F_o = 1-0.22*(log10(R_z)**0.64)*log10(sigm_b)+0.45*(log10(R_z)**0.53)
# notch factor
belt_k = alph_k/n
# total influence factor fok
F_ok = sqrt(belt_k**2-1+1/(F_o**2))

# fatigue strength of specimen
if MAT == 1:
    sigm_w = 0.27*sigm_b +100
    M = 0.00035*sigm_b +0.08
else:
    sigm_w = 0.27*sigm_b +85
    M = 0.00035*sigm_b +0.05

# fatigue shtrngth of component
sigm_wk = sigm_w/F_ok

# slopes of SN curve m1 and m2
m1 = 5.5/(F_ok**2)+6
m2 = 2*m1-1

# factor for influence of mean stress
u = 1/(M+1)*sigm_wk/sigm_b
a = (1+R)/(1-R)*sigm_wk/sigm_b
p = (1/(M+1)-1+u**2)/(u**2-u)
if a==0:
    Fm = 1
else:
    if p <= 1:
        Fm = -1*(1+p*a)/(2*a**2*(1-p))+sqrt(1/(1-p)/a**2+((1+p*a)/2/a**2/(1-p))**2)
    else:
        Fm = -1*(1+p*a)/(2*a**2*(1-p))-sqrt(1/(1-p)/a**2+((1+p*a)/2/a**2/(1-p))**2)
        
# stress amplitude at knee of SN curve
sigm_A = sigm_wk*Fm

# number of load cycles at knee of SN curve
N_D = 10**(6.8-3.6*(1/m1))

# upgrading factors
# quelity level
S_d = 0.85**(j-j_0) 
# thickness-dependent tensile value Rm
S_t = (t/25)**((-0.15)*sign)
# total upgrading factor
S = S_pu*S_d*S_t

# upgraded stress range at knee of SN curve
delt_sigm_A_s = 2*sigm_A*S/gamma_M
sigm_d = sigm_A*S/gamma_M
# upper limit of fatigue life line
delt_sigm_1 = R_p02*(1-R)/gamma_M
sigm_1 = R_p02*(1-R)/gamma_M
# number of load cycles at upper fatigue limit
delt_N_1 = N_D*(delt_sigm_A_s/delt_sigm_1)**m1
N_1 = N_D*(2*sigm_d/sigm_1)**m1