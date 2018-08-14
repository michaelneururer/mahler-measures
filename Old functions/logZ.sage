R.<q> = PowerSeriesRing(QQ)
phi = DirichletGroup(8)[2]
def logZ(t,N=100):
    s = 0
    for n in range(1,N):
        for m in srange(1,N/n):
            s+= (-1)^n * QQ(phi(m)) * n * t^(m*n)
    return s

def fun(t,M=12):
    return t - sum([2^k * t^(2^k) for k in range(1,M)])

def testfun(t,N, symbolic = False):
    s = -1
    for n in range(1,N):
        for m in srange(1,N/n):
            if symbolic:
                s+= 2*phi(n)*n*(2*q^(2*m*n)-q^(m*n))
            else:
                s+= 2*phi(n)*n*(2*exp(-2*2*pi*t*m*n)-exp(-2*pi*t*m*n))
    return s
