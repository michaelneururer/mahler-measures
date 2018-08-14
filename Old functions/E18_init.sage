#Initialise data for E18_init

#Divisors are in Matrix form the (i,j)-th entry is the coefficient of (i,j)
DivX = [[-1,0,-1,0,1,0,1,0]]+7*[[0,0,0,0,0,0,0,0]]
DivY = [[-1 ,0 ,1 ,0 ,1 ,0 ,-1,0 ]]+7*[[0,0,0,0,0,0,0,0]]

Divg = list(-2*vector((-1 ,0 ,-1 ,1 ,0 ,0 ,0, 1)) + vector((-1,0,-1,0,1,0,1,0)) + vector((-1 ,0 ,1 ,0 ,1 ,0 ,-1,0)))
#g = (xy/(x-1)^2)
Divg = [Divg] + 7*[[0,0,0,0,0,0,0,0]]

DivZ = [[0,-1,0,1,0,0,0,0]]+7*[[0,0,0,0,0,0,0,0]]
DivZm1 = [[0,-1,0,1,2,0,-2,0]]+7*[[0,0,0,0,0,0,0,0]]

#t = (Z-1)^2/Z
Divt = [[0,-1,0,1,4,0,-4,0]]+7*[[0,0,0,0,0,0,0,0]]

ProdD = (2**6/3)*EisensteinSymbol(0,1,8,tensorprod(DivZ, mu(DivY,DivY,8),8));
#The next bit is explained in the paper, E = eta'
m1 = Matrix([[0,-1],[1,-2]]); m2 = Matrix([[-1,0],[-2,-1]]); m3 = Matrix([[0,-1],[1,0]])
E8 = ProdD.apply_matrix(m1) + 2*ProdD.apply_matrix(m2) + (-1)* ProdD.apply_matrix(m3)
E8.clean_phi()

print E8

#Sym8 = ModularSymbols(Gamma1(8),3);Sym8.basis()
#matrices8 = basis_matrices(Sym8)

#forms = []
#for m in matrices:
#    print m
#    forms.append(E.integrate(m, prec = 350))

#forms_short = [vector([g[8*i]/8 for i in srange(len(g)/8)]) for g in forms];
#for g in forms_short:
#    print Pow(list(g[0:8]))
print 'Integrating E over X{0,\infty}'
res = E8.integrate(prec = 44*8, short_form = 8)
print 'Result:', res
print 'Setting up matrix with coefficients of modular forms of level 8'
f = Newforms(Gamma1(8),3)[0]
chi0, chi1,chi2= DirichletGroup(1)[0], DirichletGroup(4)[1], DirichletGroup(8)[3]
f1 = eisenstein_series_at_inf(chi0,chi1,3, prec = 44).padded_list();
f2 = eisenstein_series_at_inf(chi0,chi1,3, t=2, prec = 44).padded_list();
f3 = eisenstein_series_at_inf(chi0,chi2,3, prec = 44).padded_list();
f4 = eisenstein_series_at_inf(chi1,chi0,3, prec = 44).padded_list();
f5 = eisenstein_series_at_inf(chi1,chi0,3, t=2, prec = 44).padded_list();
f6 = eisenstein_series_at_inf(chi2,chi0,3, prec = 44).padded_list();
L8 = [f1,f2,f3,f4,f5,f6, [0]+f.coefficients(43)]; L8 = Matrix(L8)

print 'Resulting modular form as a vector in the basis'
x8 = L8.T \vector(res)
print 'Mahler measure is the L-value at s=0  of the modular form', -1/(4*pi^2)*E8.pre_factor()*x8
