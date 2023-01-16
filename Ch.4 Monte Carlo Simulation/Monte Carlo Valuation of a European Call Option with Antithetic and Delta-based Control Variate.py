from math import exp, sqrt, log
import numpy as np
import scipy.stats as scistat

K = 100 #strike
T = 1
S = 100
r = 0.06
N = 10
M = 100
sig = 0.20
div = 0.03

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

dt = T/N
nudt = (r-div-0.5*sig**2)*dt
sigsdt = sig*sqrt(dt)
erddt = exp((r-div)*dt)

beta1 = -1

sum_CT = 0
sum_CT2 = 0

for j in range(M):
    St1 = S
    St2 = S
    cv1 = 0
    cv2 = 0
    
    for i in range(1,N+1):
        t = (i-1)*dt
        delta1 = Call.Delta(St1,K,r,sig,T-t)
        delta2 = Call.Delta(St2,K,r,sig,T-t)
        epsilon = np.random.normal()
        Stn1 = St1*exp(nudt+sigsdt*epsilon)
        Stn2 = St2*exp(nudt+sigsdt*-epsilon)
        cv1 = cv1 + delta1*(Stn1-St1*erddt)
        cv2 = cv2 + delta2*(Stn2-St2*erddt)
        St1 = Stn1
        St2 = Stn2
        
    CT = 0.5*(max(0, St1-K) + beta1*cv1 + max(0,St2-K) + beta1*cv2)
    sum_CT += CT
    sum_CT2 += CT*CT
    
call_value = sum_CT/M*exp(-r*T)
SD = sqrt((sum_CT2 - sum_CT*sum_CT/M)*exp(-2*r*T)/(M-1))
SE = SD/sqrt(M)

print(f"Call Value = {call_value}")
print(f"Standard Dev of Call Value = {SD}")
print(f"Standard Error of Call Value = {SE}")