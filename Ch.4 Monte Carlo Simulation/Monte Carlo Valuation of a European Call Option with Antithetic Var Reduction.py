from math import exp, sqrt, log
import numpy as np

K = 100 #strike
T = 1
S = 100
r = 0.06
N = 1
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
    lnSt1 = lnS
    lnSt2 = lnS
    
    for i in range(1,N+1):
        epsilon = np.random.normal()
        lnSt1 = lnSt1 + nudt + sigsdt*epsilon
        lnSt2 = lnSt2 + nudt + sigsdt*-epsilon
        
    ST1 = exp(lnSt1)
    ST2 = exp(lnSt2)
    CT = 0.5*(max(0, ST1-K) + max(0, ST2-K))
    sum_CT += CT
    sum_CT2 += CT*CT
    
call_value = sum_CT/M*exp(-r*T)
SD = sqrt((sum_CT2 - sum_CT*sum_CT/M)*exp(-2*r*T)/(M-1))
SE = SD/sqrt(M)

print(f"Call Value = {call_value}")
print(f"Standard Dev of Call Value = {SD}")
print(f"Standard Error of Call Value = {SE}")