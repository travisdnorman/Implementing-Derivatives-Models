from math import exp, sqrt
import numpy as np

K = 100 #strike
T = 1
S = 100
r = 0.06
N = 3
sig = 0.20
div = 0.03
dx = 0.2

dt = T/N
nu = r-div-0.5*sig**2
edx = exp(dx)
pu = 0.5*dt*((sig/dx)**2 + nu/dx)
pm = 1 - dt*(sig/dx)**2-r*dt
pd = 0.5*dt*((sig/dx)**2-nu/dx)

assert pu > 0, f"Pu greater than 0 required, got: {pu}"
assert pm > 0, f"Pm greater than 0 required, got: {pm}"
assert pd > 0, f"Pd greater than 0 required, got: {pd}"
assert dx >= sig*sqrt(3*dt), f"Stability and Convergence condition not satisfied: dx={dx}, sig={sig}, dt={dt}"

St = np.zeros(N*2+1)
C = np.zeros((2,N*2+1))

St[0] = S*exp(-N*dx)

for j in range(1,N*2+1):
    St[j] = St[j-1]*edx

for j in range(N*2+1):
    C[0,j] = max(0,K-St[j])
    
for i in range(N-1,-1,-1):
    for j in range(1,N*2):
        C[1,j] = pu*C[0,j+1]+pm*C[0,j]+pd*C[0,j-1]
        
    C[1,0] = C[1,1] + (St[1]-St[0])
    C[1,N*2] = C[1,N*2-1]
    
    #apply early exercise condition
    for j in range(N*2+1):
        C[0,j] = max(C[1,j], K-St[j])
        
print("American Put = " + str(C[0,N]))