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
pu = -0.5*dt*((sig/dx)**2 + nu/dx)
pm = 1 + dt*(sig/dx)**2+r*dt
pd = -0.5*dt*((sig/dx)**2-nu/dx)

pmp = np.zeros(2*N+1)
pp = np.zeros(2*N+1)

def solve_implicit_tridiagonal_system(C,pu,pm,pd,lamba_L,lamba_U):
    #substitute boundary condition at j = 0 into j = 1
    pmp[1] = pm+pd
    pp[1] = C[0,1] + pd*lamba_L
    
    #eliminate upper diagonal
    for j in range(2,2*N):
        pmp[j] = pm - pu*pd/pmp[j-1]
        pp[j] = C[0,j] - pp[j-1]*pd/pmp[j-1]
        
    #use boundary condition at j=2N and equation at j=2N-1
    C[1,2*N] = (pp[2*N-1] + pmp[2*N-1]*lamba_U)/(pu + pmp[2*N-1])
    C[1,2*N-1] = C[1,2*N] - lamba_U
    
    #back substitution
    for j in range(2*N-2,-1,-1):
        C[1,j] = (pp[j] - pu*C[1,j+1])/pmp[j]

St = np.zeros(N*2+1)
C = np.zeros((2,N*2+1))

St[0] = S*exp(-N*dx)

for j in range(1,N*2+1):
    St[j] = St[j-1]*edx

for j in range(N*2+1):
    C[0,j] = max(0,K-St[j])
    
#compute derivative boundary condition
lamba_L = -1*(St[1]-St[0])
lamba_U = 0
    
for i in range(N-1,-1,-1):
    solve_implicit_tridiagonal_system(C,pu,pm,pd,lamba_L,lamba_U)
    
    #apply early exercise condition
    for j in range(N*2+1):
        C[0,j] = max(C[1,j], K-St[j])
        
print("American Put = " + str(C[0,N]))