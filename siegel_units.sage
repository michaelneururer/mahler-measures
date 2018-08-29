load('sl2_operations.py')

def is_gamma_N_equiv(x,y,N):
    #Checks if two cusps a and b are equivalent under Gamma(N)
    a,b,c,d = x.numerator(), x.denominator(), y.numerator(), y.denominator()
    if a%N == c%N and b%N == d%N:
        return True
    elif a%N == (-c)%N and b%N == (-d)%N:
        return True
    else:
        return False

def order_of_siegel_unit(a,b,N, cusp = Cusp(Infinity)):
    r"""
    Calculates the absolute order of a Siegel unit on Gamma(N) at a cusp.
    """
    #Find matrix that maps infinity to cusp
    n,d = cusp.numerator(), cusp.denominator()
    gamma = complete_column(n,d)
    #Use that up to a scalar g_(a,b) o gamma = g_((a,b)gamma)
    new_a, new_b = (Matrix([[a,b]]) * gamma).rows()[0]
    x = QQ(new_a)/QQ(N) - floor(QQ(new_a)/QQ(N))
    return 1/2 * bernoulli_polynomial(x,2)

def divisor_of_siegel_unit(a,b,group):
    r"""
        INPUT:
            a - int.
            b - int.
            group - either a congruence group or a natural number N, which stands for the group Gamma(N).
        OUTPUT:
            Returns divisor of the Siegel unit g_{a,b} with respect to the given group.
    """
    if group in ZZ:
        group = Gamma(group)
    N = group.level()
    C = group.cusps()
    v = vector(QQ,[0 for c in C])
    if Mod(a,N)==0 and Mod(b,N)==0:
        return v
    v = [order_of_siegel_unit(a,b,N, cusp = c) for c in C]
    return vector(v)
def siegel_unit_divisor_matrix(group):
    if group in ZZ:
        group = Gamma(group)
    N = group.level()
    C = group.cusps()
    mat = []
    for a in range(N):
        for b in range(N):
            mat.append(divisor_of_siegel_unit(a,b,group))
    return Matrix(QQ,mat)

def siegel_quotients_trivial(N):
    #Returns the SUQ's with trivial divisor on Gamma(N)
    mat = siegel_unit_divisor_matrix(Gamma(N))
    d = mat.denominator()
    matZZ = Matrix(ZZ,mat.denominator()*mat)
    B = matZZ.left_kernel().basis()
    res = []
    for b in B:
        tmp = []
        for a in range(N):
            tmp.append(b[a*N:(a+1)*N])
        res.append(Matrix(tmp))
    return res

def find_siegel_unit_rep(div, group=None):
    r"""
        INPUT:
            div - a divisor of the cusps of the congruence group, the cusps should be ordered just as in group.cusps(). Note that we take orders with respect to the local coordinate corresponding to q at infinity, i.e., a simple zero at a cusp c of width d corresponds to the term 1/d * (c) in the divisor.
            group - a congruence subgroup of the form GammaH or Gamma(N)
        OUTPUT: Exponents of the Siegel units needed to obtain this divisor

        EXAMPLES:
            #Gamma1(6):
            sage: Gamma1(6).cusps()
            [0, 1/3, 1/2, Infinity]
            sage: Gamma1(6).cusp_width(0), Gamma1(6).cusp_width(1/3)
            (6,2)
            sage: a = [-1/6,1/2,0,0] #This corresponds to the function s which has a simple zero at 1/3 and a pole at 0
            sage: find_siegel_unit_rep(a, group = Gamma1(6))
            [ 0 -4  4  0  0  0]
            [ 0  0  0  0  0  0]
            [ 0  0  0  0  0  0]
            [ 0  0  0  0  0  0]
            [ 0  0  0  0  0  0]
            [ 0  0  0  0  0  0]
            #This corresponds to the modular function (g_{0,2}/g_{0,1})^4.

            #Gamma1(8):
            sage: Gamma1(8).cusps()
            [0, 1/4, 1/3, 3/8, 1/2, Infinity]
            sage: Gamma1(8).cusp_width(0), Gamma1(8).cusp_width(1/3)
            (8,8)
            sage: a = [-1/8, 0, 1/8, 0, 0, 0] #This corresponds to the function Z which has a simple zero at 1/3 and a pole at 0
            sage: find_siegel_unit_rep(a, group = Gamma1(8))
            [ 0 -2  0  2  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            #This corresponds to the modular function (g_{0,3}/g_{0,1})^2.

    """

    if group in ZZ:
        group = Gamma(group)
    N = group.level()
    C = group.cusps()
    C_GammaN = Gamma(N).cusps()
    #We pull the divisor back to a divisor for Gamma(N)
    if group == Gamma(N):
        div_GammaN = div
    else:
        div_GammaN = zero_vector(QQ,len(C_GammaN))
        for i in range(len(C_GammaN)):
            for j in range(len(C)):
                if C_GammaN[i].is_gamma_h_equiv(C[j],group)[0]:
                    div_GammaN[i] += div[j]
    mat = siegel_unit_divisor_matrix(Gamma(N))
    d = mat.denominator()
    matZZ = Matrix(ZZ,d*mat)
    H,U = matZZ.hermite_form(transformation=True)
    res_tmp = H.solve_left(d*vector(div_GammaN))
    res = res_tmp * U
    output_matrix = []
    for a in range(N):
        output_matrix.append(res[a*N:(a+1)*N])
    return Matrix(output_matrix)

def siegel_unit_argument(a,b,N,cusp = Cusp(Infinity)):
    r"""
    Determines the argument of a Siegel unit at a cusp.
    """
    a = a%N
    b = b%N
    c,d = cusp.numerator(),cusp.denominator()
    if d == 0:
        if a == 0:
            return pi*P1(QQ(b)/QQ(N))
        else:
            return 0
    return pi*bernoulli_polynomial(QQ(a)/QQ(N), 2)*QQ(c)/QQ(d) \
           - 2* pi*sum([P1(QQ(1)/QQ(d)*(k-QQ(a)/QQ(N))) \
           *P1(QQ(c)/QQ(d)*(k-QQ(a)/QQ(N))-QQ(b)/QQ(N)) \
           for k in range(1,d+1)]) + (a==0)*pi*P1(QQ(b)/QQ(N))

def evaluate_siegel_unit_numerical(a,b,N,cusp = Cusp(Infinity), prec_tau = 10, prec_prod = 100):
    r"""
    Evaluates a quotient of Siegel units at a cusp numerically. The input is in the form of an NxN matrix where the (i,j)-coefficient corresponds to the exponent of g_{i,j}. Can also be called with cusp an element of the upper half plane.

        EXAMPLES:
        sage: M8 = 2*Eij(0,3,8) - 2*Eij(0,1,8)
        sage: for c in Gamma1(8).cusps():
         ...:    print c, CC(evaluate_siegel_unit_quotient_numerical(M8,c,prec_prod=5000,prec_tau=10))*(-I)
        0 2577.97166150864 - 9.09494701772928e-13*I
        1/4 -1.00000000000012 + 3.55271367880050e-15*I
        1/3 -0.000274501076652590 + 0.000274802562797830*I
        3/8 0.171572875253790 + 1.79457422433884e-15*I
        1/2 1.00000042624885 + 4.26249035068160e-7*I
        Infinity 5.82842712474619
        sage: M14_8 = Eij(0,2,8)*6 + 2*Eij(0,4,8) - 4* Eij(0,1,8) -4*Eij(0,3,8)
        sage: for c in [Cusp(Infinity),0,1/4,1/2]:
        ....:     print c, evaluate_siegel_unit_quotient_numerical(M14_8,c,prec_prod=10000,prec_tau=10)
        Infinity -2.00000000000000 - 2.89441106094618e-16*I
        0 -644.993012352625 + 7.16260561905370e-14*I
        1/4 4.12632005004012e-27 - 1.05826060025810e-39*I
        1/2 -1.00000000000000 - 1.07081010725096e-13*I
    """
    a,b = a%N, b%N
    if cusp.parent() == QQ or cusp == Cusp(Infinity):
        c,d = cusp.numerator(), cusp.denominator()
        gamma = complete_column(c,d)
        tau2 = prec_tau * I
        tau = moebius_transform(gamma,tau2)
    else:
        tau = cusp
    M = (QQ(1)/QQ(2)*bernoulli_polynomial(QQ(a)/QQ(N),2)).denominator()
    qMN = CC(exp(2*pi*I*tau/(M*N))) #This is to avoid taking fractional roots of q
    zetaN = CC(exp(2*pi*I/N))
    res = qMN**(M*N*QQ(1)/QQ(2) * bernoulli_polynomial(QQ(a)/QQ(N),2))
    res *= (1-qMN**(M*a)*zetaN**b)
    for n in range(1,prec_prod):
        res *= (1-qMN**(M*(n*N+a))*zetaN**(b))*(1-qMN**(M*(n*N-a))*zetaN**(-b))
    return res

def siegel_unit_zeta(a,b,N,gamma):
    r"""
    Determines the argument of a root of unity zeta such that g_{a,b}(\gamma\tau) = \zeta g_{(a,b)\gamma}(\tau).
    """
    a, b = a%N, b%N
    a2, b2 = vector(QQ, [a,b]) * gamma
    a2, b2 = a2%N, b2%N
    corr = 0
    if gamma[1][0] < 0:
        gamma = -gamma
        a2,b2 = -a2%N,-b2%N
        #Correction added to the argument to amount for g_{0,b} = -zetaN**(-b)*g_{0,-b}
        corr = 2*b2/QQ(N)+1 if a2==0 else 0
    c,e,d,f = gamma.list()
    if d == 0:
        return pi*e*bernoulli_polynomial(QQ(a)/QQ(N),2)
    s = QQ(c)/QQ(d)*bernoulli_polynomial(QQ(a)/QQ(N),2) + QQ(f)/QQ(d)*bernoulli_polynomial(QQ(a2)/QQ(N),2)
    s+= (a==0)*P1(QQ(b)/QQ(N)) - (a2==0)*P1(QQ(b2)/QQ(N))
    s+=- 2*sum([P1(QQ(1)/QQ(d)*(k-QQ(a)/QQ(N)))*P1(QQ(c)/QQ(d)*(k-QQ(a)/QQ(N))-QQ(b)/QQ(N)) for k in range(1,d+1)])
    s+= corr
    while True:
        if 0 <= s < 2:
            break
        if s >= 2:
            s = s-2
        if s < 0:
            s += 2
    s = pi*s
    return s

def P1(x):
    fracx = x - floor(x)
    if fracx == 0:
        return 0
    else:
        return fracx-QQ(1)/QQ(2)

class SiegelUnitQuotient(SageObject):

    def __init__(self,M,group=None):
        r"""
        Create a quotient of Siegel units.
        """
        self.M = M
        self.N = M.dimensions()[0]
        if group == None:
            self.group = Gamma(self.N)
        else:
            self.group = group
    def __repr__(self):
        print 'Siegel unit quotient for group {} represented by matrix'.format(self.group)
        return str(self.M)

    def level(self):
        return self.N
    def __neg__(self):
        return (-1)*self
    def __mul__(self,other):
        #At the moment the group associated to a product of Siegel unit quotients is set to Gamma(N) in case they are Siegel unit quotients unless self.group contains other.group or the other way around.
        if isinstance(other,SiegelUnitQuotient):
            if self.N==other.N:
                if (self.group).is_subgroup(other.group):
                    G = other.group
                elif (other.group).is_subgroup(self.group):
                    G = self.group
                else:
                    G = Gamma(self.N)
                return SiegelUnitQuotient(self.M + other.M, G)
        if other in ZZ:
            return SiegelUnitQuotient(other*self.M,self.group)
        raise ValueError

    __rmul__ = __mul__

    def __div__(self,other):
        if isinstance(other,SiegelUnitQuotient):
            return self * ((-1)*other)
        raise ValueError
    __rdiv__ = __div__

    def __pow__(self,other):
        if other in ZZ:
            return self * other

    def q_expansion(self, prec=10):
        #Returns the q_expansion without the first factor q^(1/QQ(2)*B_2(a/QQ(self.N))) with relative precision prec.
        K = CyclotomicField(self.N, 'zeta{}'.format(self.N))
        zetaN = K.gen()
        Pow = PowerSeriesRing(K,'q{}'.format(self.N))
        qN = Pow.gen()
        res = 1 + O(qN^prec)
        for i in range(self.N):
            for j in range(self.N):
                if self.M[i][j]>0:
                    res*= (1-qN^i*zetaN^j)**self.M[i][j]
                    res*= prod([((1-qN^(self.N*n+i)*zetaN^j)*(1-qN^(self.N*n-i)*zetaN^(-j)))**(self.M[i][j]) for n in range(1,ceil(prec/QQ(self.N))+1)])
                if self.M[i][j]<0:
                    if i>0:
                        res *= (sum([qN^(i*l)*zetaN^(j*l) for l in range(ceil(prec/QQ(i))+1)]))**(-self.M[i][j])
                    else:
                        res *= (1-zetaN^j)**(self.M[i][j])
                    for n in range(1,ceil(prec/QQ(self.N))+1):
                        a = sum([qN^((self.N*n+i)*l)*zetaN^(j*l) for l in range(ceil(prec/QQ(self.N*n))+1)])
                        b = sum([qN^((self.N*n-i)*l)*zetaN^(-j*l) for l in range(ceil(prec/QQ(self.N*n-i)) )])
                        res *= ((a*b)**(-self.M[i][j])+ O(qN^prec))
        return res
    qexp = q_expansion

    def log_q_expansion(self,prec=10):
        #Returns a list containing the tau-term plus the constant term and the q_expansion of log(self)
        K = CyclotomicField(self.N, 'zeta{}'.format(self.N))
        zetaN = K.gen()
        Pow = PowerSeriesRing(K,'q{}'.format(self.N))
        qN = Pow.gen()
        res = O(qN^(prec))
        tau_term= 0
        const_term = 0
        for i in range(self.N):
            for j in range(self.N):
                if self.M[i][j]!=0:
                    tau_term += pi*I*bernoulli_polynomial(i/QQ(self.N),2)
                    if i==0:
                         const_term += log(1-exp(2*pi*I*j/QQ(self.N)))
                    for m in range(1,prec):
                        for n in range(1, ceil(prec/QQ(m))):
                            if n%self.N==i:
                                res += -self.M[i][j]*zetaN^(j*m)/QQ(m)*qN^(m*n)
                            if n%self.N==-i:
                                res += -self.M[i][j]*zetaN^(-j*m)/QQ(m)*qN^(m*n)
        return const_term,tau_term,res
    log_qexp = log_q_expansion

    def log_abs_q_expansion(self,x,prec=10):
        #Returns the expansion in qy = e^(-2pi n y) of log|self(x+iy)|
        Pow = PowerSeriesRing(RR, 'qy')
        qy = Pow.gen()
        res = O(qy^prec)
        const_term,tau_term,log_exp = self.log_q_expansion(prec)
        res = CC(const_term+x*tau_term).real_part()
        for j in range(prec):
            res+=(CC(log_exp[j])*CC(exp(2*pi*I*x/self.N))).real_part()*qy^j
        return res


    def log(self,tau,prec=10):
        #Returns log(abs(self))(tau).
        c, Q = self.log_q_expansion(prec)
        Pow = PowerSeriesRing(CyclotomicField(self.N),'qN')
        qN = Pow.gen()
        Q = Pow(Q)
        res = CC(c.subs(tau=tau))
        res += (Q.polynomial().subs(qN=CC(exp(2*pi*I/self.N*tau))))
        return res
    def divisor(self):
        N = self.N
        return sum([self.M[i][j]*divisor_of_siegel_unit(i,j,group=self.group) for i in range(N) for j in range(N)])

    def evaluate(self,tau = Cusp(Infinity),numerical=False,prec_tau=100,prec_prod=100):
        r"""
        Evaluates a quotient of Siegel units at a cusp. The input is in the form of an NxN matrix where the (i,j)-coefficient corresponds to the exponent of g_{i,j}.
            INPUT:
                M - Matrix representing a Siegel unit quotient.
                tau - An element of Q or Cusp(Infinity). If numerical == True, then tau can be an element of the upper half plane
                numerical - If set to True then the infinite product defining the Siegel unit will be approximated.
                prec_tau - Only relevant if numerical == True. If gamma moves infinity to cusp, then we evaluate the infinite product at gamma*(prec_tau*I)
                prec_prod - Only relevant if numerical == True. Number of factors of the infinite product used for approximation.
            OUTPUT: (arg, r, order)
                (arg,r) is  value of the Siegel unit in polar coordinates and order is the order in the case of a zero or a pole. (arg,r) in the case of a pole can be ignored.
                In the future the function should output leading terms in case of a pole or zero.
            TODO:
                At the moment the output is a bit inconsistent. When one evaluates at Cusp(Infinity) one gets an argument of 0 and a value in a cyclotomic field, while evaluating at other cusps might produce a non-zero argument and always a real number as value. Both outputs are correct but I will at some point make it more consistent.
            EXAMPLES:
                #Gamma1(8):
                sage:  M8 = -2*Eij(0,1,8) + 2*Eij(0,3,8)
                sage: for c in Gamma1(8).cusps():
                ....:     print c, evaluate_siegel_unit_quotient(M8,c)
                0 Pole of order 1/8
                (1/2*pi, 0)
                1/4 (3/2*pi, 1)
                1/3 Zero of order 1/8
                (7/6*pi, 0)
                3/8 (1/2*pi, 0.171572875253810)
                1/2 (1/2*pi, 1)
                Infinity 2*zeta8^3 + 3*zeta8^2 + 2*zeta8

                #Gamma1(6):
                sage: a = [-1/6,1/2,0,0]
                sage: M= find_siegel_unit_rep(a, 6, group = Gamma1(6))
                sage: for c in Gamma1(6).cusps():
                ....:     print c,evaluate_siegel_unit_quotient(M,c)
                0 Pole of order 1/6
                (2/3*pi, 0)
                1/3 Zero of order 1/2
                (4/3*pi, 0)
                1/2 (2/3*pi, 1)
                Infinity 9*zeta6 - 9

                #Gamma(2,6):
                sage: M = Eij(3,0,6) + Eij(3,1,6) - Eij(3,2,6) - Eij(3,3,6)
                sage: for c in [Cusp(Infinity),0,1/5,1/4,1/3,2/3]:
                ....:     print c, evaluate_siegel_unit_quotient(M,c)
                Infinity 1
                0 Zero of order 1/6
                (0, 0)
                1/5 Pole of order 1/6
                (-8/5*pi, 0)
                1/4 (-pi, 1)
                1/3 (-pi, 1/2)
                2/3 (-pi, 2)

                #Gamma(4): #These are all correct as can be checked with evaluate_siegel_unit_quotient_numerical.
                sage: M4 = Eij(0,2,4) - 4*Eij(0,1,4)
                sage: for c in Gamma(4).cusps():
                ....:     print c,evaluate_siegel_unit_quotient(M4,c)
                0 (pi, 1)
                1/2 Zero of order 1/4
                (3/4*pi, 0)
                1 (1/2*pi, 1)
                2 (0, 1)
                3 (-1/2*pi, 1)
                Infinity Pole of order 1/4
                0

                #Gamma(4):
                sage: M = Eij(1,0,4) - Eij(0,1,4)
                sage: for c in Gamma(4).cusps():
                ....:     print c, evaluate_siegel_unit_quotient_numerical(M,c)
                0 2.78360368557622e-26 + 2.78360368557622e-26*I
                1/2 1.22304233836406e-8 - 5.65516717797777e-9*I
                1 0.980596860843385 - 0.196035117679344*I
                2 2.92654249242635e8 - 1.21221359122951e8*I
                3 0.832004450972997 - 0.554768928849868*I
                Infinity 1.90997040924671e25 + 1.90997040924671e25*I
                sage: for c in Gamma(4).cusps():
                ....:     print c, evaluate_siegel_unit_quotient(M,c)

                0 Zero of order 3/32
                (1/4*pi, 0)
                1/2 Zero of order 1/32
                (-5/32*pi, 0)
                1 (-1/16*pi, 1) #Correct numerically
                2 Pole of order 1/32
                (-1/8*pi, 0)
                3 (-3/16*pi, 1) #Not correct numerically. The argument should be 9/16*pi*I, not -3/16(pi*I).
                Infinity Pole of order 3/32
        """
        M = self.M
        N = self.N
        if not tau==Cusp(Infinity) and not tau in QQ:
            numerical = True
        if not numerical:
            cusp = tau
            if cusp == Cusp(Infinity):
                K = CyclotomicField(N)
                zetaN = K.gen()
                order_at_inf = 0
                value_at_inf = 1
                for i in range(N):
                    for j in range(N):
                        if M[i][j] != 0:
                            order_at_inf += 1/2*M[i][j]*bernoulli_polynomial(QQ(i)/QQ(N), 2)
                            if i%N == 0:
                                value_at_inf *= (1-zetaN**j)**(M[i][j])
                if order_at_inf > 0:
                    return 0,0,order_at_inf
                elif order_at_inf < 0:
                    return 0,Infinity,order_at_inf
                else:
                    return 0,value_at_inf,0
            else:
                a = cusp.numerator()
                b = cusp.denominator()
                argument = 0
                gamma = complete_column(a,b)
                M2 = Matrix([[0 for i in range(N)] for j in range(N)])
                for i in range(N):
                    for j in range(N):
                        if M[i][j] != 0:
                            i2, j2 = vector([i,j]) * gamma
                            M2 += M[i][j] * Eij(i2%N,j2%N,N)
                            argument += M[i][j] * siegel_unit_argument(i,j,N,cusp)
                return argument,abs(SiegelUnitQuotient(M2).evaluate(Cusp(Infinity))[1]),SiegelUnitQuotient(M2).evaluate(Cusp(Infinity))[2]
        if numerical:
            res = 1
            for i in range(N):
                for j in range(N):
                    if M[i][j] != 0:
                        res *= evaluate_siegel_unit_numerical(i,j,N,tau,prec_tau,prec_prod)**M[i][j]
            return res
    eval = evaluate

    def apply_matrix(self, gamma):
        r"""
        Apply a GL_2(Z/NZ) matrix to a Siegel unit quotient given as a matrix.
        INPUT:
            M - NxN-matrix representing a Siegel unit quotient, we denote by g_M.
            gamma - GL_2(Z/NZ) matrix.

        OUTPUT:
            M2 - The new Siegel unit quotient.
            arg - The argument of a root of unity zeta, such that g_M \circ \gamma = zeta * g_M2.

        EXAMPLES:
            #This is the calculation needed to confirm Lemma 8.2:
            sage: M2 = 2*Eij(6,3,8)-2*Eij(2,1,8) #Corresponds to the Siegel unit quotient (g_{6,3}/g_{2,1})^2
            sage: gamma = Matrix([[1,0],[-2,1]])
            sage: siegel_unit_quotient_zeta(M2,gamma)
            [ 0 -2  0  2  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0], -1/2*pi
            #This shows that (g_{6,3}/g_{2,1})^2 \circ gamma = -i (g_{0,3}/g_{0,1})^2

        """
        N = self.N
        M = self.M
        M2 = Matrix([[0 for i in range(N)] for j in range(N)])
        argument = 0
        for i in range(N):
            for j in range(N):
                if M[i][j] != 0:
                    i2, j2 = vector(QQ,[i,j]) * gamma
                    M2 += M[i][j] * Eij(i2%N,j2%N,N)
                    argument += M[i][j] * siegel_unit_zeta(i,j,N,gamma)
        return SiegelUnitQuotient(M2,Gamma(N)), argument

def testing(G):
    N = G.level()
    B = [Matrix([[0,1],[-1,0]]), Matrix([[2,3],[1,2]]) , Matrix([[1,4],[1,5]])]
    maximum = 0
    max_data = []
    tau = I
    for g in B:
        for a in range(N):
            for b in range(N):
                if a==0 and b==0:
                    continue
                print g,a,b,tau
                prec_prod = 100
                a2,b2 = vector(QQ,[a,b])*g
                while True:
                    x = evaluate_siegel_unit_numerical(a,b,N,moebius_transform(g,tau),prec_prod=prec_prod)
                    y = evaluate_siegel_unit_numerical(a2,b2,N,tau,prec_prod=prec_prod)
                    res1 = x/y
                    if abs(abs(res1)-1) > 0.00001:
                        prec_prod *= 5
                    else:
                        break
                print 'A', res1
                print a,b,N,g
                res2 = CC(exp(I*siegel_unit_zeta(a,b,N,g)))
                print 'B', res2
                diff = abs(res1 - res2)
                print 'Difference', diff
                if diff > maximum:
                    maximum = diff
                    max_data = [g,a,b,tau]

    return maximum, max_data
