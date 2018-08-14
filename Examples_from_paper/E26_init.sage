Divx = Eij(0,0,6) - Eij(0,2,6) - Eij(0,3,6) + Eij(0,5,6)
Divy = Eij(0,0,6) + Eij(0,1,6) - Eij(0,3,6) - Eij(0,4,6)

Div_t = (Eij(3,0,6)+Eij(3,1,6)-Eij(3,2,6)-Eij(3,3,6))

bet = mu(Divy, Divy, 6)

ProdD = EisensteinSymbol(0,1,6,tensorprod(Div_t,bet,6));
#The factor 3 = N/2 comes from converting s to an Eisenstein symbol and the factor 2 = N/3 comes from converting bet to an Eisenstein symbol
ProdD *= 3 * 2;
ProdD.clean_phi()

g = Matrix([[2,1],[3,2]]); sigma = Matrix([[0,-1],[1,0]])
E1 = 3*ProdD.apply_matrix(g)
print 'Test absolute convergence of first term', E1.abs_conv()
E2 = 2*ProdD.apply_matrix(g*sigma)
print 'Test absolute convergence of second term', E2.abs_conv()
E26 = E1 - E2
E26.clean_phi()


#Sym6 = ModularSymbols(Gamma1(6),3);
#matrices6 = basis_matrices(Sym6)

print E26
print 'Test of absolute convergence when we integrate E over X{0,\infty}:', E26.abs_conv()
print 'Integrating E over X{0,\infty}'
res = E26.integrate(prec = 20*6, short_form = 3)
print 'Result:', res
print 'Setting up matrix with coefficients of modular forms of level 6'

chi0, chi3, chi4= DirichletGroup(1)[0], DirichletGroup(3)[1], DirichletGroup(4)[1]
e1 = eis_phipsi(chi0,chi3,3, prec = 40).padded_list()
e2 = eis_phipsi(chi0,chi3,3, t=2, prec = 40).padded_list()
e3 = eis_phipsi(chi0,chi3,3, t=4, prec = 40).padded_list()
e4 = eis_phipsi(chi0,chi4,3, t=1, prec = 40).padded_list()
e5 = eis_phipsi(chi0,chi4,3, t=3, prec = 40).padded_list()
f1 = eis_phipsi(chi3,chi0,3, prec = 40).padded_list()
f2 = eis_phipsi(chi3,chi0,3, t=2, prec = 40).padded_list()
f3 = eis_phipsi(chi3,chi0,3, t=4, prec = 40).padded_list()
f4 = eis_phipsi(chi4,chi0,3, t=1, prec = 40).padded_list()
f5 = eis_phipsi(chi4,chi0,3, t=3, prec = 40).padded_list()
N = Newforms(Gamma1(12),3,names = 'a')
f = N[0]
L12 = [e1,e2,e3,e4,e5,f1,f2,f3,f4,f5] + [[0]+f.coefficients(39)]
L12 = Matrix(L12)

print 'Resulting modular form as a vector in the basis'
x26 = L12.T \vector(res)
print 'Mahler measure is the L-value at s=0 of the modular form'
print -1/(4*pi^2)*E26.pre_factor()*x26
