import pdb

from sage.modular.cusps import Cusp
from sage.modular.arithgroup.congroup_gamma0 import Gamma0_class
from sage.rings.infinity import Infinity
from sage.rings.complex_field import ComplexField, ComplexNumber
from sage.rings.rational_field import RationalField
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.big_oh import O
from sage.functions.transcendental import hurwitz_zeta
from sage.functions.generalized import sgn
from sage.functions.trig import cot
from sage.functions.other import factorial, ceil, floor
from sage.combinat.combinat import bernoulli_polynomial as ber_pol
from sage.matrix.constructor import Matrix
from sage.arith.misc import gcd, divisors
from sl2_operations import shift_cusp

Gamma0 = Gamma0_class
QQ = RationalField()
CC = ComplexField()

def hurwitz_hat(d, N, l, zetaN = None):
    r"""
    Computes the function zeta_hat(d/N, l) from Brunaults paper "Regulateurs
    modulaires via la methode de Rogers--Zudilin" as an element of Q(zetaN)
    at a non-positive integer l. We use the formula
    sum_{x in Z/NZ} zetaN**(xu) zeta(x/N, s) = N**s zeta_hat(u/N, s).
    INPUT:
    - d - int
    - N - int, positive int
    - l - int, non-positive integer
    OUTPUT:
    - an element of the ring Q(zetaN), equal to zeta_hat(d/N, l)
    """
    if zetaN == None:
        zetaN = CyclotomicField(N).gen()
    d = d%N
    k = 1-l
    if k < 1:
        raise ValueError, 'k has to be positive'
    if k == 1:
        s = -QQ(1)/QQ(2) # This is the term in the sum corresponding to x = 0
        s += sum([zetaN**(x*d)*(QQ(1)/QQ(2)-QQ(d)/QQ(N)+floor(QQ(d)/QQ(N)))
                    for x in range(1,N)])
    else:
        s = sum([-zetaN**(x*d)*(ber_pol(QQ(x)/QQ(N)-floor(QQ(x)/QQ(N)),k))/k
                    for x in range(N)])
    return N**(k-1) * s


def eis_E(cv, dv, N, k, Q=None, param_level=1, prec=10):
    r"""
    Computes the coefficient of the Eisenstein series for $\Gamma(N)$.
    Not intended to be called by user.
    INPUT:
    - cv - int, the first coordinate of the vector determining the \Gamma(N)
      Eisenstein series
    - dv - int, the second coordinate of the vector determining the \Gamma(N)
      Eisenstein series
    - N - int, the level of the Eisenstein series to be computed
    - k - int, the weight of the Eisenstein seriess to be computed
    - Q - power series ring, the ring containing the q-expansion to be computed
    - param_level - int, the parameter of the returned series will be
      q_{param_level}
    - prec - int, the precision.  The series in q_{param_level} will be truncated
      after prec coefficients
    OUTPUT:
    - an element of the ring Q, which is the Fourier expansion of the Eisenstein
      series
    """
    if Q == None:
        Q = PowerSeriesRing(CyclotomicField(N),'q')
    R = Q.base_ring()
    zetaN = R.zeta(N)
    q = Q.gen()

    cv = cv%N
    dv = dv%N

    #if dv == 0 and cv == 0 and k == 2:
    #   raise ValueError("E_2 is not a modular form")

    if k == 1:
        if cv == 0 and dv == 0:
            raise ValueError("that shouldn't have happened...")
        elif cv == 0 and dv != 0:
            s = QQ(1)/QQ(2) * (1 + zetaN**dv)/(1 - zetaN**dv)
        elif cv != 0:
            s = QQ(1)/QQ(2) - QQ(cv)/QQ(N) + floor(QQ(cv)/QQ(N))
    elif k > 1:
        if cv == 0:
            s = hurwitz_hat(QQ(dv), QQ(N), 1-k, zetaN)
        else:
            s = 0
    for n1 in xrange(1, prec): # this is n/m in DS
        for n2 in xrange(1, prec/n1 + 1): # this is m in DS
            if Mod(n1, N) == Mod(cv, N):
                s += n2**(k-1) * zetaN**(dv*n2) * q**(n1*n2)
            if Mod(n1, N) == Mod(-cv, N):
                s += (-1)**k * n2**(k-1) * zetaN**(-dv*n2) * q**(n1*n2)
    return s+O(q**floor(prec))

def eis_G(a, b, N, k, Q=None, t=1, prec = 20):
    if Q == None:
        Q = PowerSeriesRing(QQ,'q')
    R = Q.base_ring()
    q = Q.gen()
    a = ZZ(a%N)
    b = ZZ(b%N)
    s = 0

    if k==1:
        if a==0 and not b==0:
            s = QQ(1)/QQ(2)-QQ(b%N)/QQ(N)
        elif b==0 and not a==0:
            s = QQ(1)/QQ(2) - QQ(a%N)/QQ(N)
    elif k>1:
        if b==0:
            s = -N**(k-1)*ber_pol(QQ(a%N)/QQ(N),k)/QQ(k)

    #If a == 0 or b ==0 the loop has to start at 1
    starta, startb = 0, 0
    if a == 0:
        starta = 1
    if b == 0:
        startb = 1

    for v in srange(starta, (prec/t + a)/N):
        for w in srange(startb,(prec/t/abs((-a+v*N)) + b)/N+1):
            s += q**(t*(a+v*N)*(b+w*N)) * (a+v*N)**(k - 1)
            if (-a+v*N)>0 and (-b+w*N)>0:
                s += (-1)**k * q**(t*(-a+v*N)*(-b+w*N)) * (-a+v*N)**(k - 1)
    return s+O(q**floor(prec))

def eis_F(cv, dv, N, k, Q=None, prec=10, t=1):
    """
    Computes the coefficient of the Eisenstein series for $\Gamma(N)$.
    Not indented to be called by user.
    INPUT:
    - cv - int, the first coordinate of the vector determining the \Gamma(N)
      Eisenstein series
    - dv - int, the second coordinate of the vector determining the \Gamma(N)
      Eisenstein series
    - N - int, the level of the Eisenstein series to be computed
    - k - int, the weight of the Eisenstein seriess to be computed
    - Q - power series ring, the ring containing the q-expansion to be computed
    - param_level - int, the parameter of the returned series will be
      q_{param_level}
    - prec - int, the precision.  The series in q_{param_level} will be truncated
      after prec coefficients
    OUTPUT:
    - an element of the ring Q, which is the Fourier expansion of the Eisenstein
      series
    """
    if Q == None:
        Q = PowerSeriesRing(CyclotomicField(N),'q{}'.format(N))
    R = Q.base_ring()
    zetaN = R.zeta(N)
    q = Q.gen()
    s = 0
    if k == 1:
        if cv%N == 0 and dv%N != 0:
            s = QQ(1)/QQ(2) * (1 + zetaN**dv)/(1 - zetaN**dv)
        elif cv%N != 0:
            s = QQ(1)/QQ(2) - QQ(cv)/QQ(N) + floor(QQ(cv)/QQ(N))
    elif k > 1:
        s = - ber_pol(QQ(cv)/QQ(N) - floor(QQ(cv)/QQ(N)),k)/QQ(k)
    for n1 in xrange(1, ceil(prec/QQ(t))): # this is n/m in DS
        for n2 in xrange(1, ceil(prec/QQ(t)/QQ(n1)) + 1): # this is m in DS
            if Mod(n1, N) == Mod(cv, N):
                s += N**(1-k)*n1**(k-1) * zetaN**(dv*n2) * q**(t*n1*n2)
            if Mod(n1, N) == Mod(-cv, N):
                s += (-1)**k * N**(1-k) * n1**(k-1) * zetaN**(-dv*n2) * q**(t*n1*n2)
    return s+O(q**floor(prec))

def eis_H(a, b, N, k, Q=None, t=1, prec = 10):
    if Q == None:
        Q = PowerSeriesRing(CyclotomicField(N),'q')
    R = Q.base_ring()
    zetaN = R.zeta(N)
    q = Q.gen()
    a = ZZ(a%N)
    b = ZZ(b%N)
    s = 0

    if k==1:
        if a==0 and not b==0:
            s = -QQ(1)/QQ(2) * (1+zetaN**b)/(1-zetaN**b)
        elif b==0 and not a==0:
            s = -QQ(1)/QQ(2) * (1+zetaN**a)/(1-zetaN**a)
        elif a != 0 and b != 0:
            s = -QQ(1)/QQ(2) * ((1+zetaN**a)/(1-zetaN**a)+(1+zetaN**b)/(1-zetaN**b))
    elif k>1:
        s = hurwitz_hat(-b,N,1-k,zetaN)

    for m in srange(1, prec/t):
        for n in srange(1,prec/t/m+1):
            s += (zetaN**(-a*m-b*n)+(-1)**k*zetaN**(a*m+b*n))*n**(k-1)*q**(m*n)
    return s+O(q**floor(prec))

def eis_phipsi(phi, psi, k, prec=10, t=1,cmplx = False):
    r"""
    Return Fourier expansion of Eisenstein series at the cusp oo.

    INPUT:

    - ``phi`` -- Dirichlet character.
    - ``psi`` -- Dirichlet character.
    - ``k`` -- integer, the weight of the Eistenstein series.
    - ``prec`` -- integer (default: 10).
    - ``t`` -- integer (default: 1).

    OUTPUT:

    The Fourier expansion of the Eisenstein series $E_k^{\phi,\psi, t}$ (as
    defined by [Diamond-Shurman]).

    EXAMPLES:
    sage: phi = DirichletGroup(3)[1]
    sage: psi = DirichletGroup(5)[1]
    sage: E = eisenstein_series_at_inf(phi, psi, 4)
    """
    N1, N2 = phi.level(), psi.level()
    N = N1*N2
    #The Fourier expansion of the Eisenstein series at infinity is in the field Q(zeta_Ncyc)
    Ncyc = lcm([euler_phi(N1), euler_phi(N2)])
    if cmplx == True:
        CC = ComplexField(53)
        pi = ComplexField().pi()
        I = ComplexField().gen()
        R = CC
        zeta = CC(exp(2*pi*I/Ncyc))
    else:
        R = CyclotomicField(Ncyc)
        zeta = R.zeta(Ncyc)
        phi, psi = phi.base_extend(R), psi.base_extend(R)
    Q = PowerSeriesRing(R, 'q')
    q = Q.gen()
    s = O(q**prec)

    #Weight 2 with trivial characters is calculated separately
    if k==2 and phi.conductor()==1 and psi.conductor()==1:
        if t==1:
            raise TypeError('E_2 is not a modular form.')
        s = 1/24*(t-1)
        for m in srange(1,prec):
            for n in srange(1,prec/m+1):
                s += n * (q**(m*n)-t*q**(m*n*t))
        return s+O(q**prec)

    if psi.level()==1 and k==1:
        s -= phi.bernoulli(k)/ k
    elif phi.level()==1:
        s -= psi.bernoulli(k)/ k

    for m in srange(1, prec/t):
        for n in srange(1,prec/t/m+1):
            s += 2* phi(m) * psi(n) * n**(k-1) * q**(m*n*t)
    return s+O(q**prec)
