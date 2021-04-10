# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 12:40:16 2021

@author: SIAVASH
"""

from matplotlib import pyplot as plt
from qutip import *
import numpy as np
from qutip.qip.gates import *
from numpy import e, real, sort, sqrt

D = []
CC = []
LL = []
VV = []

#Basis
u1 = basis(2,0)
u2 = basis(2,1)
u_p = (u1 + u2)/np.sqrt(2)
#intial density matrix
s = np.linspace(0,1,100)
p = np.ndarray.tolist(s)

#p = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
for px in p: 

    rho_0 = px*tensor((u1*u1.dag()),(u1*u1.dag())) + (1 - px)*tensor(u_p*u_p.dag(), u2*u2.dag())
    
    #von Neumman Measurement
    gamma1 = u1 * u1.dag()
    gamma2 = u2 * u2.dag()
    
    p1 = tensor(gamma1,qeye(2)) * rho_0 * tensor(gamma1.dag(),qeye(2))
    p2 = tensor(gamma2,qeye(2)) * rho_0 * tensor(gamma2.dag(),qeye(2))
    p11 = np.trace(p1)
    p22 = np.trace(p2) 
    
    rho_B1 = p1.ptrace(1) / p11
    rho_B2 = p2.ptrace(1) / p22
    
    rho_A = rho_0.ptrace(0)
    rho_B = rho_0.ptrace(1)
    
    #Eigen Values
    E_A = rho_A.eigenenergies()
    E_B = rho_B.eigenenergies()
    E_AB = rho_0.eigenenergies()
    E_B1 = rho_B1.eigenenergies()
    E_B2 = rho_B2.eigenenergies()
    
    
    #Calculating Vomn Neumman Entropy
    S_rho_AB = entropy_vn (rho_0)
    S_rho_A = entropy_vn (rho_A)
    S_rho_B = entropy_vn (rho_B)
    S_rho_B1 = entropy_vn (rho_B1)
    S_rho_B2 = entropy_vn (rho_B2)
    
    # Quantum Discord Calculations
    p1_min = min(E_B1)* S_rho_B1
    p2_min = min(E_B2)* S_rho_B2
    p1_max = max(E_B1)* S_rho_B1
    p2_max = max(E_B2)* S_rho_B2
    
    #S_B = p1_min * S_rho_B1 + p2_min * S_rho_B2
    Discord = S_rho_A - S_rho_AB + p22*S_rho_B2 + p11*S_rho_B1
    
    #Calculating Concurrence
    C0 = concurrence(rho_0)
    M0 = entropy_linear(rho_0)
    
    D.append(Discord)
    CC.append(C0)
    LL.append(M0)
    VV.append(S_rho_AB)
    
#plt.plot(p,D,'b-', label = 'Discord')
#plt.plot(p,CC,'k.', label = 'Concurrence' )
#plt.plot(p,VV,'r-', label = 'linear entropy')

# Plotting Figs
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, constrained_layout=True, sharey=True)
l1, = ax1.plot(p, D)
l2, = ax1.plot(p,CC)
#ax1.set_title('Discord and Concurrence')
ax1.set_xlabel('p')
ax1.set_ylabel('D and C')
fig.legend((l1, l2), ('Discord', 'Concurrence'), 'upper left')

ax2.plot(p,VV)
ax2.set_xlabel('p')
ax2.set_title('von Neumman entropy')

ax3.plot(p,LL)
ax3.set_xlabel('p')
ax3.set_title('Linear entropy')

#fig.suptitle('QI', fontsize=16)
