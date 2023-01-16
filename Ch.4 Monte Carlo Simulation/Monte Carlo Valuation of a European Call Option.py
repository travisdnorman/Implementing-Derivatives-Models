from math import exp, sqrt, log
import numpy as np

K = 100 #strike
T = 1
S = 100
r = 0.06
N = 10
M = 100
sig = 0.20
div = 0.03

dt = T/N
nudt = (r-div-0.5*sig**2)*dt
sigsdt = sig*sqrt(dt)
lnS = log(S)

sum_CT = 0
sum_CT2 = 0

for j in range(M):
    lnSt = lnS
    for i in range(1,N+1):
        epsilon = np.random.normal()
        lnSt = lnSt + nudt + sigsdt*epsilon
        
    ST = exp(lnSt)
    CT = max(0, ST-K)
    sum_CT += CT
    sum_CT2 += CT*CT
    
call_value = sum_CT/M*exp(-r*T)
SD = sqrt((sum_CT2 - sum_CT*sum_CT/M)*exp(-2*r*T)/(M-1))
SE = SD/sqrt(M)

print(f"Call Value = {call_value}")
print(f"Standard Dev of Call Value = {SD}")
print(f"Standard Error of Call Value = {SE}")