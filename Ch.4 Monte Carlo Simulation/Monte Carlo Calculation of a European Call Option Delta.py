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

dt = T
nudt = (r-div-0.5*sig**2)*dt
sigsdt = sig*sqrt(dt)
lnS = log(S)

sum_CT = 0
sum_CT2 = 0

for j in range(M):
    epsilon = np.random.normal()
    e = exp(nudt + sigsdt*epsilon)
    
    ST = S*e
    
    if (ST > K):
        CT = e
    else:
        CT = 0
        
    sum_CT += CT
    sum_CT2 += CT*CT
    
delta_value = sum_CT/M*exp(-r*T)
SD = sqrt((sum_CT2 - sum_CT*sum_CT/M)*exp(-2*r*T)/(M-1))
SE = SD/sqrt(M)

print(f"Delta Value = {delta_value}")
print(f"Standard Dev of Delta Value = {SD}")
print(f"Standard Error of Delta Value = {SE}")