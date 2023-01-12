import math
import numpy as np

K = 100 #strike
T = 1
S = 100
r = 0.06
N = 3
u = 1.1
d = 1/u
dT = T/N
p = (math.exp(r*dT)-d)/(u-d)
disc = math.exp(-r*dT)

St = np.zeros(N+1)
C = np.zeros(N+1)

for j in range(N+1):
    St[j] = S*(d**(N-j)*u**j)

for j in range(N+1):
    C[j] = max(0,St[j]-K)
    
for i in reversed(range(N)):
    for j in range(i+1):
        C[j] = disc*(p*C[j+1]+(1-p)*C[j])