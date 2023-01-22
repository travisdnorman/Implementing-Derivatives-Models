from math import exp, sqrt, log
import numpy as np
import scipy.stats as scistat

def d1(S, K, r, sigma, T):
    return (np.log(S/K) + (r+sigma*sigma/2)*T)/(sigma*sqrt(T))

def d2(S, K, r, sigma, T):
    return d1(S, K, r, sigma, T) - sigma*sqrt(T)

class Call:        
    def Price(S, K, r, sigma, T):
        return np.maximum(S - K, 0) if T==0 else S*scistat.norm.cdf(d1(S, K, r, sigma, T)) - K*np.exp(-r*T)*scistat.norm.cdf(d2(S, K, r, sigma, T))

    def Delta(S, K, r, sigma, T):
        return scistat.norm.cdf(d1(S, K, r, sigma, T))

    def Gamma(S, K, r, sigma, T):
        return scistat.norm.pdf(d1(S, K, r, sigma, T))/(S*sigma*sqrt(T))

    def Vega(S, K, r, sigma, T):
        return S*scistat.norm.pdf(d1(S, K, r, sigma, T))*sqrt(T)

    def Theta(S, K, r, sigma, T):
        aux1 = -S*scistat.norm.pdf(d1(S, K, r, sigma, T))*sigma/(2*sqrt(T))
        aux2 = -r*K*np.exp(-r*T)*scistat.norm.cdf(d2(S, K, r, sigma, T))
        return aux1+aux2

    def Rho(S, K, r, sigma, T):
        return K*T*np.exp(-r*T)*scistat.norm.cdf(d2(S, K, r, sigma, T))
    
K = 100 #strike
T = 1
S = 100
r = 0.06
N = 52
M = 1000
sig = 0.20
div = 0.03
alpha = 0.5
Vbar = 0.10
xi = .02

beta1 = -0.88
beta2 = -0.42
beta3 = -0.0003


dt = T/N
sig2 = sig**2
sdt = sqrt(dt)
alphadt = alpha*dt
xisdt = xi*sqrt(dt)
erddt = exp((r-div)*dt)
egam1 = exp(2*(r-div)*dt)
egam2 = -2*erddt+1
eveg1 = exp(-alphadt)
eveg2 = Vbar - Vbar*eveg1

sum_CT = 0
sum_CT2 = 0

for j in range(M):
    St1 = S
    St2 = S
    Vt = sig2
    maxSt1 = St1
    maxSt2 = St2
    cv1 = 0
    cv2 = 0
    cv3 = 0
    
    for i in range(1,N+1):
        t = (i-1)*dt
        delta1 = Call.Delta(St1,K,r,sig,T-t)
        delta2 = Call.Delta(St2,K,r,sig,T-t)
        gamma1 = Call.Gamma(St1,K,r,sig,T-t)
        gamma2 = Call.Gamma(St2,K,r,sig,T-t)
        vega1 = Call.Vega(St1,K,r,sig,T-t)
        vega2 = Call.Vega(St2,K,r,sig,T-t)
        
        epsilon = np.random.normal()
        Vtn = Vt + alphadt*(Vbar-Vt) + xisdt*sqrt(Vt)*epsilon
        
        epsilon = np.random.normal()
        Stn1 = St1*exp((r-div-0.5*Vt)*dt + sqrt(Vt)*sdt*epsilon)
        Stn2 = St2*exp((r-div-0.5*Vt)*dt + sqrt(Vt)*sdt*(-epsilon))
        
        cv1 = cv1 + delta1*(Stn1-St1*erddt) + delta2*(Stn2-St2*erddt)
        cv2 = cv2 + gamma1*((Stn1-St1)**2-St1**2*(egam1*exp(Vt*dt)+egam2)) \
                    + gamma2*((Stn2-St2)**2-St2**2*(egam1*exp(Vt*dt)+egam2))
        cv3 = cv3 + vega1*((Vtn-Vt)-(Vt*eveg1+eveg2-Vt)) \
                    + vega2*((Vtn-Vt)-(Vt*eveg1+eveg2-Vt))
                    
        Vt = Vtn
        St1 = Stn1
        St2 = Stn2
        
        maxSt1 = max(maxSt1, St1)
        maxSt2 = max(maxSt2, St2)
        
        
    CT = 0.5*(max(0,maxSt1-K) + max(0,maxSt2-K) + beta1*cv1 + beta2*cv2 + beta3*cv3)
    sum_CT += CT
    sum_CT2 += CT*CT
    
call_value = sum_CT/M*exp(-r*T)
SD = sqrt((sum_CT2 - sum_CT*sum_CT/M)*exp(-2*r*T)/(M-1))
SE = SD/sqrt(M)

print(f"Call Value = {call_value}")
print(f"Standard Dev of Call Value = {SD}")
print(f"Standard Error of Call Value = {SE}")