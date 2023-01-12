import math
import numpy as np
from matplotlib import pyplot as plt

K = 100 #strike
T = 1
S = 100
r = 0.06
N = 3
sig = 0.20

#set Trigeorgis coedfficients
dT = T/N
nu = r-0.5*sig**2
dxu = math.sqrt(sig**2*dT+(nu*dT)**2)
dxd = -dxu
pu = 1/2 + 1/2*(nu*dT/dxu)
pd = 1-pu

disc =math.exp(-r*dT)

St = np.zeros(N+1)
C = np.zeros(N+1)

St[0] = S*math.exp(N*dxd)
for j in range(1,N+1):
    St[j] = St[j-1]*math.exp(dxu-dxd)

for j in range(N+1):
    C[j] = max(0,St[j]-K)

plt.figure(1)
plt.scatter(np.full(N+1, N),St)

plt.figure(2)
plt.scatter(np.full(N+1, N),C)
    
for i in reversed(range(N)):
    for j in range(i+1):
        C[j] = disc*(pu*C[j+1]+pd*C[j])
        
    St[0] = S*math.exp(i*dxd)
    for j in range(1,i+1):
        St[j] = St[j-1]*math.exp(dxu-dxd)
    plt.figure(1)
    plt.scatter(np.full(i+1, i),St[:i+1])
    
    plt.figure(2)
    plt.scatter(np.full(i+1, i),C[:i+1])

plt.figure(1)    
plt.show()

plt.figure(2)    
plt.show()
        
print("European Call Value = ", C[0])