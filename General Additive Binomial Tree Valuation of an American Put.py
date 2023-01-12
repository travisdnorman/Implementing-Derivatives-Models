import math
import numpy as np
from matplotlib import pyplot as plt

K = 100 #strike
T = 1
S = 100
r = 0.06
N = 3
sig = 0.2
dh = 0.00 #dividend rate
tau = 2/3 #time of dividend
Hd = 72 #down barrier price

dT = T/N
nu =r - 0.5*sig**2
dxu = math.sqrt(sig**2*dT+(nu*dT)**2)
dxd = -dxu
pu = 1/2 + 1/2*(nu*dT/dxu)
pd = 1-pu

disc = math.exp(-r*dT)
dpu = disc*pu
dpd = disc*pd
edxud = math.exp(dxu-dxd)
edxd = math.exp(dxd)

St = np.zeros(N+3)
C = np.zeros(N+3)

St[0] = S*(1-dh)*math.exp((N+2)*dxd)
for j in range(1,N+3):
    St[j] = St[j-1]*edxud

#add Asset Values to plot    
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
plt.scatter(np.full(N+3,N),St)
for v in St:
    ax1.text(N,v+2,round(v,2),ha="center")

for j in range(N+3):
    if St[j] > Hd:
        C[j] = max(0,K-St[j])
    else:
        C[j] = 0

#add Option values to plot    
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
plt.scatter(np.full(N+3,N),C)
for v in C:
    ax2.text(N,v+2,round(v,4),ha="center")
    
for i in reversed(range(N)):
    if i*dT < tau and (i+1)*dT >= tau:
        adj = 1-dh
    else:
        adj = 1
    for j in range(i+3):
        St[j] = St[j]/(adj*edxd)
        if St[j] > Hd:
            C[j] = dpd*C[j]+dpu*C[j+1]
            C[j] = max(C[j], K-St[j])
        else:
            C[j] = 0
        
    plt.figure(1)
    plt.scatter(np.full(i+3,i),St[:i+3])
    for v in St[:i+3]:
        ax1.text(i,v+2,round(v,2),ha="center")
    
    plt.figure(2)
    plt.scatter(np.full(i+3,i),C[:i+3])
    for v in C[:i+3]:
        ax2.text(i,v+2,round(v,4),ha="center")
    
plt.figure(1)
plt.show()

plt.figure(2)
plt.show()

print("American Put = ", C[1])

#Sensitivity Calculations
Delta = (C[2]-C[0])/(St[2]-St[0])
print("Delta = ", Delta)
Gamma = (((C[2]-C[1])/(St[2]-St[1]))-((C[1]-C[0])/(St[1]-St[0])))/(1/2*(St[2]-St[0]))
print("Gamma = ", Gamma)
