load('sl2_operations.py')
load('eisenstein_series.py')

def beta(a,b, N):
    out = [[0 for i in range(N)] for j in range(N)]
    for i in range(N):
        for j in range(N):
            for x in range(N):
                for y in range(N):
                    out[(i-x)%N][(j-y)%N] += a[i][j]*b[x][y]
    return Matrix(out)

def mu(a,b,N):
    out = [[0 for i in range(N)] for j in range(N)]
    for i in range(N):
        for j in range(N):
            for x in range(N):
                for y in range(N):
                    out[(i+x)%N][(j+y)%N] += a[i][j]*b[x][y]
    return Matrix(out)

def tensorprod(D1,D2, N):
    #For two divisors it returns the tensor product as a dictionary that has keys (a,b,c,d) with values e, where e is the coefficient at (a,b)x(c,d)
    tensordiv = {}
    for i in range(N):
        for j in range(N):
            for x in range(N):
                for y in range(N):
                    if not D1[i][j]*D2[x][y]==0:
                        try:
                            tensordiv[(i,j,x,y)] += D1[i][j]*D2[x][y]
                        except KeyError:
                            tensordiv[(i,j,x,y)] = D1[i][j]*D2[x][y]
    return tensordiv

class EisensteinSymbol(SageObject):
    #Class to work with the Eisenstein symbol Eis^{k1,k2}(phi).

    def __init__(self, k1, k2, N, phi):
        r"""
        Input:
            k1 - int
            k2 - int
            N - int, the level
            phi - Dictionary that associates to (i,j,x,y) a coefficient
        """
        self.k1 = k1
        self.k2 = k2
        self.N = N
        self.phi = phi

    def __repr__(self):
        return 'Eisenstein symbol of weight ({},{}), level {}, and \n phi = {}'.format(self.k1,self.k2, self.N, self.phi)

    def latex(self):
        s = ''
        for X in self.phi.keys():
            try:
                #Omit the +-sign if self.phi[X] < 0
                if self.phi[X] > 0:
                    s += '+'
            except:
                s += '+'
            if not (self.phi[X] == 1 or self.phi[X] == -1):
                s += latex(self.phi[X])
            s += ' '
            s += '({}, {}, {}, {})'.format(X[0],X[1],X[2],X[3])
        return s

    def __eq__(self, other):
        self.clean_phi()
        other.clean_phi()
        if isinstance(other, EisensteinSymbol):
            if self.k1==other.k1 and self.k2==other.k2 and self.N==other.N and self.phi==other.phi:
                return True
        else:
            return False

    def __add__(self, other):
        if isinstance(other, EisensteinSymbol):
            if not [self.k1,self.k2,self.N] == [other.k1,other.k2,other.N]:
                raise ValueError, 'Can only add Eisenstein symbols of the same weight and level'
            out = {}
            for X in self.phi.keys():
                out[X] = self.phi[X]
            for X in other.phi.keys():
                try:
                    out[X] += other.phi[X]
                except KeyError:
                    out[X] = other.phi[X]

            return EisensteinSymbol(self.k1,self.k2,self.N,out)

        raise ValueError

    def __sub__(self,other):
        return self + (-1)*other

    def __mul__(self, scalar):
        out = {}
        for X in self.phi.keys():
            out[X] = scalar * self.phi[X]
        return EisensteinSymbol(self.k1,self.k2,self.N,out)

    __rmul__ = __mul__

    def __div__(self,scalar):
        return self.__mul__(1/scalar)

    __rdiv__ = __div__

    def clean_phi(self):
        #Removes unnecessary terms such as Eis^(0,1)((u,v),(0,0)) = 0
        out = {}
        for x in self.phi.keys():
            if self.k1%2 and (x[0]%self.N == -x[0]%self.N and x[1]%self.N==-x[1]%self.N):
                    continue
            if self.k2%2 and (x[2]%self.N == -x[2]%self.N and x[3]%self.N==-x[3]%self.N):
                continue
            a = 1 if x[0]%self.N == 0 else 0
            b = 3 if x[2]%self.N == 0 else 2
            if x[a]%self.N > self.N/2 and x[b]%self.N > self.N/2:
                try:
                    out[(-x[0]%self.N, -x[1]%self.N, -x[2]%self.N, -x[3]%self.N)] += (-1)**(self.k1+self.k2)*self.phi[x]
                except KeyError:
                    out[(-x[0]%self.N, -x[1]%self.N, -x[2]%self.N, -x[3]%self.N)] = (-1)**(self.k1+self.k2)*self.phi[x]
            elif x[a]%self.N > self.N/2:
                try:
                    out[(-x[0]%self.N, -x[1]%self.N, x[2]%self.N, x[3]%self.N)] += (-1)**(self.k1)*self.phi[x]
                except KeyError:
                    out[(-x[0]%self.N, -x[1]%self.N, x[2]%self.N, x[3]%self.N)] = (-1)**(self.k1)*self.phi[x]
            elif x[b]%self.N > self.N/2:
                try:
                    out[(x[0]%self.N, x[1]%self.N, -x[2]%self.N, -x[3]%self.N)] += (-1)**(self.k2)*self.phi[x]
                except KeyError:
                    out[(x[0]%self.N, x[1]%self.N, -x[2]%self.N, -x[3]%self.N)] = (-1)**(self.k2)*self.phi[x]
            else:
                try:
                    out[x] += self.phi[x]
                except KeyError:
                    out[x] = self.phi[x]
        for x in out.keys():
            if out[x] == 0:
                out.pop(x)
        self.phi = out
    def flip(self):
        #For a single summand E_D^{k1,k2}(u_1,u_2) it returns E_D^{k2,k1}(u_2,u_1). For a general Eisenstein symbol we extend linearly.
        out = {}
        for X in self.phi.keys():
            out[(X[2],X[3],X[0],X[1])] = self.phi[X]
        return EisensteinSymbol(self.k2,self.k1,self.N,out)
    def apply_matrix(self, g):
        #Applies a matrix g to a product divisor
        if gcd(det(g),self.N) > 1:
            raise ValueError, 'g has to be in GL_2(Z/NZ)'
        out = {}
        for X in self.phi.keys():
            i,j = vector([X[0],X[1]]) * g
            x,y = vector([X[2],X[3]]) * g
            i,j,x,y = i%self.N, j%self.N, x%self.N, y%self.N
            out[(i,j,x,y)] = self.phi[X]
        return EisensteinSymbol(self.k1,self.k2,self.N,out)

    def abs_conv(self):
        #Tests if the Eisenstein symbol of a divisor D converges absolutely when integrated over X{0,\infty}.
        if not (self.k1 == 0 and self.k2 == 1):
            raise NotImplementedError, 'At the moment we only deal with k1 == 0, k2 == 1'
        s1 = 0
        s2 = 0
        zetaN = CC(exp(2*pi*I/self.N))
        for X in self.phi.keys():
            s1 += self.phi[X] * bernoulli_polynomial((X[0]%self.N)/self.N, 2) * bernoulli_polynomial((X[2]%self.N)/self.N,3)
            if X[0] == 0:
                s2 += self.phi[X] * log(abs(1-zetaN**(X[1]))) * bernoulli_polynomial((X[2]%self.N)/self.N,3)
        return s1,s2

    def residues(self, alpha):
        #See Proposition 7.4
        resX, resY = 0,0
        zetaN = CC(exp(2*pi*I/self.N))
        for X in self.phi.keys():
            if X[0] == 0:
                resX += self.phi[X] * -4*pi**2/self.N**2 * log(abs(1-zetaN**X[1])) * bernoulli_polynomial((X[2]%self.N)/self.N,3) * alpha**2
                resY += self.phi[X] * -8*pi**2/self.N**2 * log(abs(1-zetaN**X[1])) * bernoulli_polynomial((X[2]%self.N)/self.N,3) * alpha
            if X[2] == 0:
                resX += self.phi[X] * 3*pi*I/self.N**2 * bernoulli_polynomial((X[0]%self.N)/self.N,2) * (hurwitz_hat_numerical(-X[3],self.N,2) - hurwitz_hat_numerical(X[3],self.N,2)) * alpha
        return CC(resX), CC(resY)

    def chasles_defect(self, alpha, beta, gamma):
        #Compute the integrals of eta along the cycles X{alpha,gamma}-X{alpha,beta}-X{beta,gamma} and Y{alpha,gamma}-Y{alpha,beta}-Y{beta,gamma}, where alpha, beta, gamma are cusps
        if alpha == beta or beta == gamma or alpha == gamma:
            return CC(0),CC(0)
        intX, intY = 0, 0
        #Compute the contribution at the cusp alpha
        a = alpha.numerator()
        c = alpha.denominator()
        g = complete_column(a,c)
        g_inv = g^(-1)
        b = g[0,1]
        d = g[1,0]
        z_beta = moebius_transform(g^(-1),beta)
        z_gamma = moebius_transform(g^(-1),gamma)
        eta_g = self.apply_matrix(g)
        res_beta = eta_g.residues(QQ(z_beta))
        res_gamma = eta_g.residues(QQ(z_gamma))
        intX_eta_g = res_beta[0] - res_gamma[0]
        intY_eta_g = res_beta[1] - res_gamma[1]
        intX += a*intX_eta_g + b*intY_eta_g
        intY += c*intX_eta_g + d*intY_eta_g
        #Compute the contribution at the cusp beta
        a = beta.numerator()
        c = beta.denominator()
        g = complete_column(a,c)
        g_inv = g^(-1)
        b = g[0,1]
        d = g[1,0]
        z_gamma = moebius_transform(g^(-1),gamma)
        z_alpha = moebius_transform(g^(-1),alpha)
        eta_g = self.apply_matrix(g)
        res_gamma = eta_g.residues(QQ(z_gamma))
        res_alpha = eta_g.residues(QQ(z_alpha))
        intX_eta_g = res_gamma[0] - res_alpha[0]
        intY_eta_g = res_gamma[1] - res_alpha[1]
        intX += a*intX_eta_g + b*intY_eta_g
        intY += c*intX_eta_g + d*intY_eta_g
        #Compute the contribution at the cusp gamma
        a = gamma.numerator()
        c = gamma.denominator()
        g = complete_column(a,c)
        g_inv = g^(-1)
        b = g[0,1]
        d = g[1,0]
        z_alpha = moebius_transform(g^(-1), alpha)
        z_beta = moebius_transform(g^(-1),beta)
        eta_g = self.apply_matrix(g)
        res_alpha = eta_g.residues(QQ(z_alpha))
        res_beta = eta_g.residues(QQ(z_beta))
        intX_eta_g = res_alpha[0] - res_beta[0]
        intY_eta_g = res_alpha[1] - res_beta[1]
        intX += a*intX_eta_g + b*intY_eta_g
        intY += c*intX_eta_g + d*intY_eta_g
        return intX, intY

    def bad_term(self, g=Matrix([[1,0],[0,1]])):
        #This is a temporary function that returns the "bad term" of an Eisenstein symbol. I.e., the terms that produce quasi-modular forms, where we are not sure if the main theorem applies.
        symbol = self.apply_matrix(g)
        bad_terms = {}
        for i in range(symbol.N):
            for j in range(symbol.N):
                #Check if the sum of the coefficients of (u1,j,i,0) is 0.
                check = sum([symbol.phi[X] for X in symbol.phi.keys() if X[1] == j and X[2] == i and X[3] == 0])
                if not check == 0:
                    for h in range(symbol.N):
                        if (h,j,i,0) in symbol.phi.keys():
                            bad_terms[(h,j,i,0)] = symbol.phi[(h,j,i,0)]
        return EisensteinSymbol(self.k1,self.k2,self.N,bad_terms)

    def integrate(self, prec = 20, short_form = 1,check_modular = False):
        #Returns the modular form that appears in the L-function when one integrates the Eisenstein symbol over g^* X{0,\infty}
        #If short_form is greater than 1 it is assumed that the resulting modular form is in q**short_form and we return f(tau/short_form).
        #If check_modular == True, then
        if self.k2 == 1 and check_modular == True:
            for i in range(self.N):
                for j in range(self.N):
                    #Check if the sum of the coefficients of (u1,j,i,0) is 0.
                    check = sum([self.phi[X] for X in self.phi.keys() if X[1] == j and X[2] == i and X[3] == 0])
                    if not check == 0:
                        print 'Might be quasi-modular, since the coefficients of ( * , {}, {}, 0) dont add up to 0'.format(j,i)
        s = 0
        for X in self.phi.keys():
            s += self.phi[X] * brunault_int([X[0],X[1]],[X[2],X[3]], identity_matrix(2), self.N, self.k1, self.k2, prec = prec)
        if short_form > 1:
            for i in range(len(s)):
                if s[i] != 0 and not short_form.divides(i):
                    raise ValueError, 'The q^(1/{})-expansion resulting from the integral is not a q^{}-expansion'.format(self.N, short_form/self.N)
            return vector([s[short_form*i] for i in srange(prec/short_form)])
        return s
    def A(self, prec = 20, short_form = 1):
        s = 0
        Pow = PowerSeriesRing(QQ,'q')
        q = Pow.gen()
        for X in self.phi.keys():
            a1,b1,a2,b2 = X[0],X[1],X[2],X[3]
            G1 = eis_G(b2,a1,self.N,self.k2+1, Pow, prec = prec) + eis_G(b2,-a1,self.N,self.k2+1,Pow,prec=prec)
            G2 = eis_G(b1,-a2,self.N,self.k1+1, Pow, prec = prec) - eis_G(b1,a2,self.N,self.k1+1,Pow,prec=prec)
            s += self.phi[X] * G1 * G2
        if short_form > 1:
            for i in range(len(s)):
                if s[i] != 0 and not short_form.divides(i):
                    raise ValueError, 'The q^(1/{})-expansion resulting from the integral is not a q^{}-expansion'.format(self.N, short_form/self.N)
            return vector([s[short_form*i] for i in srange(prec/short_form)])
        return 1/2*s+O(q^prec)

    def I_terms(self, prec = 20):
        #Returns the modular forms that appear in the L-functions of I1, I2, I3, on page 1144 of 'Regulateurs...'
        #The one for I1 is a weight k1+k2+2 quasi-modular. The L-value is taken at s = k1+k2+2
        #The one for I2 is a weight k1+1 modular form. The L-value is taken at s = k1+k2+2
        #The one for I3 is a weight k2+1 quasi-modular form. The L-value is taken at s = k2+1
        #WARNING: The result doesn't seem quite correct yet, perhaps it is correct up to a factor.
        I1,I2,I3 = 0,0,0
        K = CyclotomicField(lcm(4,self.N),'zeta{}'.format(lcm(4,self.N)))
        zeta4N = K.zeta(lcm(4,self.N))
        I = zeta4N**N
        Pow.<q> = PowerSeriesRing(K)
        for X in self.phi.keys():
            f = eis_H(X[3],X[0],self.N,self.k2+1,Pow,t=1,prec=prec) + eis_H(X[3],-X[0],self.N,self.k2+1,Pow,t=1,prec=prec)
            WN2f = eis_G(X[3],X[0],self.N,self.k2+1,Pow,t=1,prec=prec) + eis_G(X[3],-X[0],self.N,self.k2+1,Pow,t=1,prec=prec)
            WN2g = I**(self.k1+1)/QQ(self.N)*(eis_H(X[1],-X[2],self.N,self.k1+1,Pow,t=1,prec=prec) - eis_H(X[1],X[2],self.N,self.k1+1,Pow,t=1,prec=prec))
            g = eis_G(X[1],-X[2],self.N,self.k1+1,Pow,t=1,prec=prec) - eis_G(X[1],X[2],self.N,self.k1+1,Pow,t=1,prec=prec)
            I1 += self.phi[X] * f * WN2g + O(q**prec)
            I2 += -f[0] * self.phi[X] *  WN2g
            I3 += -g[0] * self.phi[X] * f
            #Now add (-1)^k times the I^(k2,k1)(u_2,u_2)
            f2 = eis_H(X[1],X[2],self.N,self.k1+1,Pow,t=1,prec=prec) + eis_H(X[1],-X[2],self.N,self.k1+1,Pow,t=1,prec=prec)
            WN2g2 = I**(self.k2+1)/QQ(self.N)*(eis_H(X[3],-X[0],self.N,self.k2+1,Pow,t=1,prec=prec) - eis_H(X[3],X[0],self.N,self.k2+1,Pow,t=1,prec=prec))
            g2 = eis_G(X[3],-X[1],self.N,self.k2+1,Pow,t=1,prec=prec) - eis_G(X[3],X[1],self.N,self.k2+1,Pow,t=1,prec=prec)
            I1 += (-1)^(self.k1+self.k2+1)*self.phi[X] * f2 * WN2g2 + O(q**prec)
            I2 += -(-1)^(self.k1+self.k2+1)*f2[0] * self.phi[X] *  WN2g2
            I3 += -(-1)^(self.k1+self.k2+1)*g2[0] * self.phi[X] * f2
        return I1,I2,I3

    def pre_factor(self):
        #Returns the Prefactor of the L-function in Theoreme 1.1 [Brunault-Regulateurs]
        k1,k2,N = self.k1,self.k2,self.N
        return (k1+2)*(k2+2)/(2*N^(k1+k2+2))* (2*pi)^(k1+k2+1) * I^(k1-k2+1)

def brunault_int(u,v,g, N, k1=0,k2=1,prec = 53):
    #u,v are vectors in (Z/N)^2 and g is in GL_2(Z/N)
    #returns - up to a scalar factor - the modular form that appears in the L-function in the integral of Eis^{k1,k2}(u,v) over xi(g). Equivalently the integral of Eis^{k1,k2}(gu,gv) over X{0,\infty}
    unew = [g[0][0]*u[0]+g[1][0]*u[1], g[0][1]*u[0]+g[1][1]*u[1]]
    vnew = [g[0][0]*v[0]+g[1][0]*v[1], g[0][1]*v[0]+g[1][1]*v[1]]
    unew = map(lambda x: ZZ(x%N), unew)
    vnew = map(lambda x: ZZ(x%N), vnew)

    if k1 == 1 and unew[1] == 0:
        raise NotImplementedError, 'k1 == 1 and b1 == 0 is not implemented'
    Pow = PowerSeriesRing(QQ,'q')
    q = Pow.gen()
    G1 = eis_G(vnew[1],unew[0],N, k2+1, Pow, prec = prec)
    G2 = eis_G(unew[1],-vnew[0],N, k1+1, Pow, prec = prec)
    G3 = eis_G(vnew[1],-unew[0],N, k2+1, Pow, prec = prec)
    G4 = eis_G(unew[1],vnew[0],N, k1+1, Pow, prec = prec)
    GG1 = G1*G2 + O(q^prec)
    GG2 = G3*G4 + O(q^prec)
    return vector(GG1.padded_list()) - vector(GG2.padded_list())

######################################

def lin_dep(L,k):
    A = pari('lindep({},{})'.format(L,k))
    return A, sum([A[i]*L[i] for i in range(len(L))])

def evaluateDirichletLFunction(chi, m):
    """
    Evaluates the Dirichlet L-function of chi at m<=0 using the formula in
    Zagier's book
    """
    if (type(m)!=Integer or m>0):
        raise ValueError('evaluation of Dirichlet L-functions only valid '\
                'when argument is integer <= 0')
    n = -m
    N = chi.level()
    s = 0
    for i in range(1, N+1):
        s += chi(i)*bernoulli_polynomial(i/N, n+1)
    return (-N^n/(n+1))*s

def hurwitz_hat_numerical(d, N, s):
    r"""
    Computes the function zeta_hat(d/N, l) from Brunaults paper "Regulateurs
    modulaires via la methode de Rogers--Zudilin" as an element of CC. We use the formula
    sum_{x in Z/NZ} zetaN**(xu) zeta(x/N, s) = N**s zeta_hat(u/N, s).
    INPUT:
    - d - int
    - N - int, positive int
    - l - int, non-positive integer
    OUTPUT:
    - an element of the ring Q(zetaN), equal to zeta_hat(d/N, l)
    """
    zetaN = exp(2*pi*I/QQ(N))
    d = d%N
    result = 0
    for x in range(1,N+1):
        result += CC(zetaN**(x*d)*hurwitz_zeta(s,QQ(x)/QQ(N)))
    return N**(-s) * result


def basis_matrices(M):
    #Input: Space of modular symbols
    #Output: matrices corresponding to the bottom rows given in the Manin-symbol representation of the basis
    #Warning: Only in weight 3 all the Manin symbols automatically come with the polynomial X. For higher weights this is not the case!
    B = M.basis()
    N = M.level()
    #The following is a list of the bottom rows of all the elements of B. It's surprisingly hard to extract them from a modular symbol.
    C = map(lambda x: x.manin_symbol_rep()[0][1].tuple()[1:3], B)
    matrices = []
    for c in C:
        tmp = [c[0],c[1]]
        if c[0] == 0:
            tmp[0] += N
        if c[1] == 0:
            tmp[1] += N
        matrices.append(complete_row(tmp[0],tmp[1],2))
    return matrices

def lin_dep(L,k):
    K = pari('lindep({},{})'.format(L,k))
    return K, sum([K[i] * L[i] for i in range(len(L))])
