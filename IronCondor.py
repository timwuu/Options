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


#strategy iron condor

def strategy_iron_condor( dataset, strike_a, strike_b, strike_c, strike_d, vol_discount=1.0, leverage=1.0):

    prv_dat = None
    prv_tmp = (None,None,None,None,None,None,None)

    payoff1 = 0.0
    payoff2 = 0.0
    ret = 1.0

    for dat in dataset:
        S = dat[2]
        K1 = round(S * strike_a)
        K2 = round(S * strike_b)
        K3 = round(S * strike_c)
        K4 = round(S * strike_d)
        r = dat[3]
        sigma = dat[1] * vol_discount
        T = 4./52.

        PUT1= PUT(S, K1, r, sigma, T)
        PUT2= PUT(S, K2, r, sigma, T)
        CALL1= CALL(S, K3, r, sigma, T)
        CALL2= CALL(S, K4, r, sigma, T)

        credit1 = PUT2-PUT1
        credit2 = CALL1-CALL2

        base = S  #(K1+K2)/2.0

        if prv_dat is not None:
            _K1 = prv_tmp[0]
            _K2 = prv_tmp[1]
            _K3 = prv_tmp[2]
            _K4 = prv_tmp[3]
            _credit1 = prv_tmp[4]
            _credit2 = prv_tmp[5]
            _base = prv_tmp[6]

            payoff1 = _credit1
            if S < _K2:
                payoff1 -= _K2 - max( S, _K1)

            payoff2 = _credit2
            if S > _K3:
                payoff2 += _K3 - min( S, _K4)

            ret *= (1.0 + (payoff1+payoff2)/_base*leverage)

        prv_dat = dat
        prv_tmp = ( K1, K2, K3, K4, credit1, credit2, base)

        #print( prv_tmp, payoff, ret)

    return ret


strike_a = float(sys.argv[1])
strike_b = float(sys.argv[2])
strike_c = float(sys.argv[3])
strike_d = float(sys.argv[4])
vol_discount = float(sys.argv[5])
leverage = float(sys.argv[6])

ret1 = strategy_iron_condor( data.dat, strike_a, strike_b, strike_c, strike_d, vol_discount, leverage)

print( ret1)