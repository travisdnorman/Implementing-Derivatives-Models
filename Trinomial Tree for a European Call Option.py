from math import exp
import numpy as np
from matplotlib import pyplot as plt

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
pu = 0.5*((sig**2*dt+nu**2*dt**2)/dx**2 + nu*dt/dx)
pm = 1 - (sig**2*dt+nu**2*dt**2)/dx**2
pd = 0.5*((sig**2*dt+nu**2*dt**2)/dx**2-nu*dt/dx)
disc = exp(-r*dt)

St = np.zeros(N*2+1)
C = np.zeros((N+1,N*2+1))

St[0] = S*exp(-N*dx)

for j in range(1,N*2+1):
    St[j] = St[j-1]*edx

for j in range(N*2+1):
    C[N][j] = max(0,St[j]-K)
    
for i in range(N-1,-1,-1):
    for j in range((N-i),i+N+1):
        C[i,j] = disc*(pu*C[i+1,j+1]+pm*C[i+1,j]+pd*C[i+1,j-1])
        
print("European Call = " + str(C[0,N]))