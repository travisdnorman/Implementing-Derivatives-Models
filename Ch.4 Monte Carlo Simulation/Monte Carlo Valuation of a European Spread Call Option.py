from math import exp, sqrt
import numpy as np

K = 1 #strike
T = 1
S1 = 100
S2 = 110
r = 0.06
N = 1
M = 100
sig1 = 0.20
sig2 = 0.30
rho = 0.50
div1 = 0.03
div2 = 0.04

dt = T/N
nu1dt = (r-div1-0.5*sig1**2)*dt
nu2dt = (r-div2-0.5*sig2**2)*dt
sig1sdt = sig1*sqrt(dt)
sig2sdt = sig2*sqrt(dt)
srho = sqrt(1-rho**2)

sum_CT = 0
sum_CT2 = 0

for j in range(M):
    St1 = S1
    St2 = S2
    cv1 = 0
    cv2 = 0
    
    for i in range(1,N+1):
        z1 = np.random.normal()
        epsilon2 = np.random.normal()
        z2 = rho*z1+srho*epsilon2
        St1 = St1*exp(nu1dt+sig1sdt*z1)
        St2 = St2*exp(nu2dt+sig2sdt*z2)
        
    CT = max(0, St1-St2-K)
    sum_CT += CT
    sum_CT2 += CT*CT
    
call_value = sum_CT/M*exp(-r*T)
SD = sqrt((sum_CT2 - sum_CT*sum_CT/M)*exp(-2*r*T)/(M-1))
SE = SD/sqrt(M)

print(f"Call Value = {call_value}")
print(f"Standard Dev of Call Value = {SD}")
print(f"Standard Error of Call Value = {SE}")