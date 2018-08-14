B = ModularForms(Gamma1(5),2).basis()
f = B[2]
f2 = f.qexp(100)
Pow.<q> = PowerSeriesRing(QQ)
def phi15_qexp(prec):
    prod = Pow(1)
    for n in range(prec):
        if n%5 == 1 or n%5 == 4:
            prod *= (sum([q^(n*l) for l in range(ceil(prec/n))])+O(q^prec))^5

        if n%5 ==2 or n%5 == 3:
            prod *= ((1-q^n)+O(q^prec))^5
    return prod
p = q^(-1)*phi15_qexp(100)
#Lets find (h(z) - h(0))(h(z) - h(1/2))(h(z) h(2/5))
#We have h(2/5) = 0, h(1/2) = 5*zeta5^3 + 5*zeta5^2 + 8, h(0) = - 5*zeta5^3 - 5*zeta5^2 + 3
#We see h(1/2)+h(0) = 11, h(1/2)*h(0) = -1. Hence (h(z)-h(1/2))(h(z)-h(0)) = h(z)^2-11h(z)-1
pre_fac = p*(p^2-11*p-1)
def cusp_rep(c):
    #Finds the cusp 0,2/5,1/2 or infty to which c is equivalent to
    for x in Gamma1(5).cusps():
        if Cusp(c).is_gamma1_equiv(x,5)[0]:
            return x

def canonical_f(k,m):
    #k even and non-negative
    if m < -k:
        raise ValueError, 'We must have m \geq -k'
    if k == 0:
        L = [1]
        for i in srange(1,m+1):
            tmp = L[0] * p^(i) #This is a form starting with q^(-i)
            for j in srange(-i+1,1):
                tmp -= tmp[j] * L[k-j]
            L.append(tmp)
        return L
    else:
        L = [f2^(k/2)]
        for i in srange(-k+1,m+1):
            tmp = L[0] * p^(k+i) #This is a form starting with q^(-i)
            for j in srange(-i+1,k+1):
                tmp -= tmp[j] * L[k-j]
            L.append(tmp)
        return L
def canonical_g(k,m):
    #k even and positive
    #Recursively finds all g_{k,m} from g_{k,-k+3} to f_{k,m}
    if m < -k+3:
        raise ValueError, 'We must have m \geq -k+3'
    L = [f2^(k/2)*pre_fac]
    for i in srange(-k+4,m+1):
        tmp = L[0] * p^(k-3+i) #This is a form starting with q^(-i)
        for j in srange(-i+1,k-2):
            tmp -= tmp[j] * L[k-3-j]
        L.append(tmp)
    return L

Rq.<q> = LaurentSeriesRing(QQ);

def Up(series, p):
    qvar = (series.parent()).gen();
    upper = floor((series.prec()-0.1)/p) + 1;
    lower = min(ceil(series.valuation()/p), 0);
    ans = O(qvar^upper);
    if p == 1:
        return series;
    else:
        for k in range(lower, upper):
            if series[k*p] != 0:
                ans = ans + series[k*p]*qvar^(k);
        return ans ;
