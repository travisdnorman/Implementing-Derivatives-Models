from math import exp, sqrt, log
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
rho12 = 0.8
rho13 = 0.5
rho14 = -0.6
div1 = 0.03
div2 = 0.04
alpha1 = 0.5
alpha2 = 0.4
Vbar1 = 0.2
Vbar2 = 0.3
xi1 = .05
xi2 = .02


dt = T/N
sdt = sqrt(dt)
alpha1dt = alpha1*dt
alpha2dt = alpha2*dt
xi1sdt = xi1*sqrt(dt)
xi2sdt = xi2*sqrt(dt)
lnS1 = log(S1)
lnS2 = log(S2)
srho12 = sqrt(1-rho12**2)
srho13 = sqrt(1-rho13**2)
srho14 = sqrt(1-rho14**2)

sum_CT = 0
sum_CT2 = 0

for j in range(M):
    lnSt1 = lnS1
    lnSt2 = lnS2
    cv1 = 0
    cv2 = 0
    Vt1 = Vbar1
    Vt2 = Vbar2
    
    for i in range(1,N+1):
        epsilon = [np.random.normal() for _ in range(4)]
        z1 = epsilon[0]
        z2 = rho12*z1+srho12*epsilon[1]
        z3 = rho13*z1+srho13*epsilon[2]
        z4 = rho14*z1+srho14*epsilon[3]
        
        Vt1 = Vt1 + alpha1dt*(Vbar1-Vt1) + xi1sdt*sqrt(Vt1)*z3
        Vt2 = Vt2 + alpha2dt*(Vbar2-Vt2) + xi2sdt*sqrt(Vt2)*z4
        
        lnSt1 = lnSt1 + (r-div1-0.5*Vt1)*dt + sqrt(Vt1)*sdt*z1
        lnSt2 = lnSt2 + (r-div2-0.5*Vt2)*dt + sqrt(Vt2)*sdt*z2
        
    St1 = exp(lnSt1)
    St2 = exp(lnSt2)
        
    CT = max(0, St1-St2-K)
    sum_CT += CT
    sum_CT2 += CT*CT
    
call_value = sum_CT/M*exp(-r*T)
SD = sqrt((sum_CT2 - sum_CT*sum_CT/M)*exp(-2*r*T)/(M-1))
SE = SD/sqrt(M)

print(f"Call Value = {call_value}")
print(f"Standard Dev of Call Value = {SD}")
print(f"Standard Error of Call Value = {SE}")