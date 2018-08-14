#NEW NOTATION FOR DIVISORS:
DivX = Matrix([[ 1, 0, -1, -1, 0, 1 ]]+5*[[0,0,0,0,0,0]])
DivY = Matrix([[ 1, 1, 0, -1, -1, 0 ]]+5*[[0,0,0,0,0,0]])

Div_s = 4*Matrix([[0,-1,1,0,0,0]]+5*[[0,0,0,0,0,0]])

bet = mu(DivY, DivY, 6)

ProdD = EisensteinSymbol(0,1,6,tensorprod(Div_s,bet,6));
ProdD *= 2*3
ProdD.clean_phi()

m1 = Matrix([[0,-1],[1,-2]]); m2 = Matrix([[-1,0],[-2,-1]]); m3 = Matrix([[0,-1],[1,0]])
E6 = ProdD.apply_matrix(m1) + 2*ProdD.apply_matrix(m2) + (-1)* ProdD.apply_matrix(m3)
E6.clean_phi()

Sym6 = ModularSymbols(Gamma1(6),3);
matrices6 = basis_matrices(Sym6)

print E6
print 'Test of absolute convergence when we integrate E over X{0,\infty}:', E6.abs_conv()
print 'Integrating E over X{0,\infty}'
res = E6.integrate(prec = 20*6, short_form = 6)
print 'Result:', res
print 'Setting up matrix with coefficients of modular forms of level 6'

#Better:
chi0, chi1= DirichletGroup(1)[0], DirichletGroup(3)[1]
f1 = eis_phipsi(chi0,chi1,3, prec = 20).padded_list()
f2 = eis_phipsi(chi0,chi1,3, t=2, prec = 20).padded_list()
f3 = eis_phipsi(chi1,chi0,3, prec = 20).padded_list()
f4 = eis_phipsi(chi1,chi0,3, t=2, prec = 20).padded_list()
L6 = [f1,f2,f3,f4]
L6 = Matrix(L6)

print 'Resulting modular form as a vector in the basis'
x6 = L6.T \vector(res)
print 'Mahler measure is the L-value at s=0  of the modular form', -1/(4*pi^2)*E6.pre_factor()*x6
