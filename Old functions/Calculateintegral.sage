import pdb
def errorterm(N,a1,a2,b1,b2,k1,k2):
    a = (-1)^k1*i^k2*factorial(k1+2)*(k2+2)/(8*pi^2*N^2)*zeta_hat(-b1/N,k1+2,k1,0,N)*zeta_hat(a2/N,2,1,0,N)*zeta_hat(-b2/N,-k2+1,k2,0,N)
    if k1 == 0 and mod(b1,N)==0:
        b = i^k2*(k2+2)/(2*N^2)*zeta_hat(-a1/N,0,0,N)*zeta_hat(a2/N,1,1,0,N)*zeta_hat(-b2/N,-k2,k2+1,0,N)
    else:
        b = 0
    if mod(k2,2)==1:
        c = (k1+2)*(k2+2)/(2*N^2)*(2*pi)^(k1+2*k2+2)*zeta_hat(a1/N,-k2,0,0,N)*hurwitz_zeta(-k2-1,-a2/N)*
    return a+b
    

def zeta_hat(x,s, sign=1,eps=0, N=1):
#    Input: Real number x, positive integer s, sign 1 or -1
#    Output: zeta_hat(x,s)+(-1)^sign zeta_hat(-x,s)
    pdb.set_trace()
    if (-1)^sign == -1 and s ==1:
        return 2*i*pi * (1/2-x+x.floor())    
    elif (-1)^sign == (-1)^s:
        return -bernoulli_polynomial(x-x.floor(),s)*(2*i*pi)^s/factorial(s)
    else:
        zetaN = CyclotomicField(N).gen()
        tmp1 = sum([CC(zetaN)^(N*x*u)*hurwitz_zeta(s,u/N) for u in range(1,N)])
        tmp2 = sum([CC(zetaN)^(-N*x*u)*hurwitz_zeta(s,-u/N) for u in range(1,N)])
        return N^(-s)*(tmp1 + (-1)^sign * tmp2)
