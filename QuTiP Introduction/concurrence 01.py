# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 13:25:35 2021

@author: SIAVASH
"""

from matplotlib import pyplot as plt
from qutip import *
import numpy as np
from qutip.qip.gates import *
from numpy import e, real, sort, sqrt


#Quantum States

s1 = basis(2,0)
s2 = basis(2,1)
s_p = (s1 + s2)/sqrt(2)

psi1 = tensor(s1,s1)
psi2 = (tensor(s1,s1) + tensor(s2,s2))/sqrt(2)

rho1 = ket2dm(psi1)
rho2 = ket2dm(psi2)

rho3 = (tensor((s1*s1.dag()),(s1*s1.dag())) + tensor(s_p*s_p.dag(), s2*s2.dag()))/2

#Calculating Von Neumman Entropy
V1 = entropy_vn (rho1)
V2 = entropy_vn (rho2)
V3 = entropy_vn (rho3)

#Calculating Concurrence
C1 = concurrence (rho1)
C2 = concurrence (rho2)
C3 = concurrence (rho3)

#Calculating Linear entropy
L1 = entropy_linear (rho1)
L2 = entropy_linear (rho2)
L3 = entropy_linear (rho3)

