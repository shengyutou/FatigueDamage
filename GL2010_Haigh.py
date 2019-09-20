# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 17:22:55 2019

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