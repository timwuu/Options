import math

# Abramowiz Stegun Approximation
# HASTINGS.  MAX ERROR = .000001
def normalcdf(X):
	T=1/(1+.2316419*math.fabs(X))
	D= 0.3989422804*math.exp(-X*X/2.0)   # 1/sqrt(2*pi)
	Prob=D*T*(.31938153+T*(-.356563782+T*(1.781477937+T*(-1.821255978+T*1.330274429))))
	if (X>0.0):
		Prob=1-Prob

	return Prob

def d1(S,K,r,sigma,T):
    """
    * Calc d1(S,K,r,sigma,T).
    * @param {number} S - the stock price.
    * @param {number} K - the strike price.
    * @param {number} r - the interest rate.
    * @param {number} sigma - the volatility.
    * @param {number} T - the duration.
    * @return d1.
    * @customfunction """
    return (math.log(S/K)+(r+sigma*sigma/2.0)*T)/(sigma*math.sqrt(T))

def DELTACALL(S,K,r,sigma,T):
    """  Calc delta of CALL options.
    @param {number} S - the stock price.
    @param {number} K - the strike price.
    @param {number} r - the interest rate.
    * @param {number} sigma - the volatility.
    * @param {number} T - the duration.
    * @return delta of CALL options.
    * @customfunction """
    
    return normalcdf((math.log(S/K)+(r+sigma*sigma/2.0)*T)/(sigma*math.sqrt(T)))

def CALL(S,K,r,sigma,T):
    _d1= d1(S,K,r,sigma,T)
    _d2= _d1-sigma*math.sqrt(T)
  
    return S*normalcdf(_d1) - K*normalcdf(_d2)*math.exp(-r*T) 


print( normalcdf(0.0))

S = 665.85
K = 732
r = 0.0527
sigma = 0.1864 * 0.9
T = 4./52.

print( CALL(S, K, r, sigma, T))