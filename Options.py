import math
import data
import sys

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

def PUT(S,K,r,sigma,T):
    _d1= d1(S,K,r,sigma,T)
    _d2= _d1-sigma*math.sqrt(T)

    return S*(normalcdf(_d1)-1.0) + K*math.exp(-r*T)*(1.0-normalcdf(_d2))


#strategy call credit spread

def strategy_call_cs( dataset, strike_a, strike_b, vol_discount=1.0, leverage=1.0):

    prv_dat = None
    prv_tmp = (None,None,None,None)

    payoff = 0.0
    ret = 1.0

    for dat in dataset:
        S = dat[2]
        K1 = round(S * strike_a)
        K2 = round(S * strike_b)
        r = dat[3]
        sigma = dat[1] * vol_discount
        T = 4./52.

        CALL1= CALL(S, K1, r, sigma, T)
        CALL2= CALL(S, K2, r, sigma, T)

        credit = CALL1-CALL2

        base = (K1+K2)/2.0

        if prv_dat is not None:
            _K1 = prv_tmp[0]
            _K2 = prv_tmp[1]
            _credit = prv_tmp[2]
            _base = prv_tmp[3]

            payoff = _credit
            if S > _K1:
                payoff += _K1 - min( S, _K2)

            ret *= (1.0 + payoff/_base*leverage)

        prv_dat = dat
        prv_tmp = ( K1, K2, credit, base)

        #print( prv_tmp, payoff, ret)

    return ret


strike_a = float(sys.argv[1])
strike_b = float(sys.argv[2])
vol_discount = float(sys.argv[3])
leverage = float(sys.argv[4])

ret1 = strategy_call_cs( data.dat, strike_a, strike_b, vol_discount, leverage)

print( ret1)