#Divisors are in Matrix form the (i,j)-th entry is the coefficient of (i,j)
DivX = [[0,0,-2,0],[0,0,0,0],[0,0,2,0],[0,0,0,0]]
DivY = [[0,0,0,0],[-2,0,0,0],[0,0,0,0],[2,0,0,0]]

DivZ1 = [[0,-2,0,0],[0,0,0,0],[0,2,0,0],[0,0,0,0]]
DivZ2 = [[0,-4,1,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]

bet = mu(DivX,DivY,4)

ProdD1 = EisensteinSymbol(0,1,4,tensorprod(DivZ1, bet, 4));
ProdD2 = EisensteinSymbol(0,1,4,tensorprod(DivZ2, bet, 4));

ProdD1*= 2*4/3
ProdD2*= 2*4/3
#The factor 2 = N/2 comes from converting DivZ to an Eisenstein symbol and the factor 4/3 = N/3 comes from converting bet to an Eisenstein symbol

#The next bit is explained in the paper, E = eta'
m1 = Matrix([[1,0],[1,1]]); m2 = Matrix([[2,1],[1,1]]);
sigma = Matrix([[0,-1],[1,0]]);
m3 = Matrix([[1,0],[-1,1]]); m4 = Matrix([[2,-1],[-1,1]]);

X02_1 = ProdD1.apply_matrix(m1) + 2*ProdD1.apply_matrix(m2) + (-1)*ProdD1.apply_matrix(m2*sigma)
X02_2 = ProdD2.apply_matrix(m1) + 2*ProdD2.apply_matrix(m2) + (-1)*ProdD2.apply_matrix(m2*sigma)

X0m2_1= ProdD1.apply_matrix(m3) + 2*ProdD1.apply_matrix(m4) + ProdD1.apply_matrix(m4*sigma)
X0m2_2= ProdD2.apply_matrix(m3) + 2*ProdD2.apply_matrix(m4) + ProdD2.apply_matrix(m4*sigma)

Y02_1 = ProdD1.apply_matrix(m1) + (-1)*ProdD1.apply_matrix(m1*sigma) + ProdD1.apply_matrix(m2) + (-1)*ProdD1.apply_matrix(m2*sigma)
Y02_2 = ProdD2.apply_matrix(m1) + (-1)*ProdD2.apply_matrix(m1*sigma) + ProdD2.apply_matrix(m2) + (-1)*ProdD2.apply_matrix(m2*sigma)

Y0m2_1= (-1)*ProdD1.apply_matrix(m3) + (-1)*ProdD1.apply_matrix(m3*sigma) + (-1)*ProdD1.apply_matrix(m4) + (-1)*ProdD1.apply_matrix(m4*sigma)
Y0m2_2= (-1)*ProdD2.apply_matrix(m3) + (-1)*ProdD2.apply_matrix(m3*sigma) + (-1)*ProdD2.apply_matrix(m4) + (-1)*ProdD2.apply_matrix(m4*sigma)

r"""
m1 = Matrix([[1,0],[1,1]]); m2 = Matrix([[2,1],[1,1]]);
sigma = Matrix([[0,-1],[1,0]]);
eps = Matrix([[-1,0],[0,1]]);
E4_1 = ProdD.apply_matrix(m1) + ProdD.apply_matrix(m2) + (-1)*ProdD.apply_matrix(m1*sigma) + (-1)*ProdD.apply_matrix(m2*sigma)
E4_2 = (-1)*ProdD.apply_matrix(eps*m3) + (-1)*ProdD.apply_matrix(eps*m4) + (-1)*ProdD.apply_matrix(eps*m3*sigma) + (-1)*ProdD.apply_matrix(eps*m4*sigma)
"""

#E4 = E4_1 + (-1)*E4_2
#E4.clean_phi()

#print E4

int_X02_1 = X02_1.integrate(prec = 100)
int_X02_2 = X02_2.integrate(prec = 100)
int_X0m2_1 = X0m2_1.integrate(prec = 100)
int_X0m2_2 = X0m2_2.integrate(prec = 100)
int_Y02_1 = Y02_1.integrate(prec = 100)
int_Y02_2 = Y02_2.integrate(prec = 100)
int_Y0m2_1 = Y0m2_1.integrate(prec = 100)
int_Y0m2_2 = Y0m2_2.integrate(prec = 100)
#I1,I2,I3 = E4.I_terms(prec = 100)
#print 'Result:', res
print 'Setting up matrix with coefficients of modular forms of level 16'
K = CyclotomicField(4)
M = ModularForms(Gamma1(16),3).basis()
#MatM = Matrix(K,[m.coefficients(srange(100)) for m in M]
MatM = Matrix(QQ,[m.coefficients(srange(100)) for m in M])

print 'Integrating eta over X{0,2}'
try:
    x_1 = MatM.T \vector(int_X02_1)
    print 'Mahler measure is the L-value at s=0  of the modular form', -1/(4*pi^2)*X02_1.pre_factor()*x_1
except ValueError:
    print 'res is not a modular form'
try:
    x_2 = MatM.T \vector(int_X02_2)
    print 'Mahler measure is the L-value at s=0  of the modular form', -1/(4*pi^2)*X02_2.pre_factor()*x_2
except ValueError:
    print 'res is not a modular form'

print 'Integrating eta over X{0,-2}'
try:
    x_1 = MatM.T \vector(int_X0m2_1)
    print 'Mahler measure is the L-value at s=0  of the modular form', -1/(4*pi^2)*X0m2_1.pre_factor()*x_1
except ValueError:
    print 'res is not a modular form'
try:
    x_2 = MatM.T \vector(int_X0m2_2)
    print 'Mahler measure is the L-value at s=0  of the modular form', -1/(4*pi^2)*X0m2_2.pre_factor()*x_2
except ValueError:
    print 'res is not a modular form'

print 'Integrating eta over Y{0,2}'
try:
    x_1 = MatM.T \vector(int_Y02_1)
    print 'Mahler measure is the L-value at s=0  of the modular form', -1/(4*pi^2)*Y02_1.pre_factor()*x_1
except ValueError:
    print 'res is not a modular form'
try:
    x_2 = MatM.T \vector(int_Y02_2)
    print 'Mahler measure is the L-value at s=0  of the modular form', -1/(4*pi^2)*Y02_2.pre_factor()*x_2
except ValueError:
    print 'res is not a modular form'

print 'Integrating eta over Y{0,-2}'
try:
    x_1 = MatM.T \vector(int_Y0m2_1)
    print 'Mahler measure is the L-value at s=0  of the modular form', -1/(4*pi^2)*Y0m2_1.pre_factor()*x_1
except ValueError:
    print 'res is not a modular form'
try:
    x_2 = MatM.T \vector(int_Y0m2_2)
    print 'Mahler measure is the L-value at s=0  of the modular form', -1/(4*pi^2)*Y0m2_2.pre_factor()*x_2
except ValueError:
    print 'res is not a modular form'

print 'Determining the quasi-modular form'
print 'Determining the weight 3 part'

E = EisensteinForms(Gamma1(16),1).basis()
MatE = Matrix([m.coefficients(srange(100)) for m in E])
MatDE = copy(MatE)

# Derivate the Eisenstein series of weight 1
for i in range(MatDE.nrows()):
    for j in range(MatDE.ncols()):
        MatDE[i,j]*=j

#Mat=block_matrix(K,[[MatM],[MatDE]])
Mat=block_matrix(QQ,[[MatM],[MatDE]])

x02_1 = Mat.T \vector(int_X02_1)
x02_2 = Mat.T \vector(int_X02_2)
x0m2_1 = Mat.T \vector(int_X0m2_1)
x0m2_2 = Mat.T \vector(int_X0m2_2)
y02_1 = Mat.T \vector(int_Y02_1)
y02_2 = Mat.T \vector(int_Y02_2)
y0m2_1 = Mat.T \vector(int_Y0m2_1)
y0m2_2 = Mat.T \vector(int_Y0m2_2)

#A-terms
a02_2 = Mat.T \ vector(X02_2.A(prec=100).padded_list())
a0m2_2 = Mat.T \vector(X0m2_2.A(prec=100).padded_list())

# Weight 3 part of the quasi-modular form
F3_x02_1 = sum(x02_1[i]*MatM[i] for i in range(MatM.nrows()))
F3_x02_2 = sum(x02_2[i]*MatM[i] for i in range(MatM.nrows()))
F3_x0m2_1 = sum(x0m2_1[i]*MatM[i] for i in range(MatM.nrows()))
F3_x0m2_2 = sum(x0m2_2[i]*MatM[i] for i in range(MatM.nrows()))
F3_y02_1 = sum(y02_1[i]*MatM[i] for i in range(MatM.nrows()))
F3_y02_2 = sum(y02_2[i]*MatM[i] for i in range(MatM.nrows()))
F3_y0m2_1 = sum(y0m2_1[i]*MatM[i] for i in range(MatM.nrows()))
F3_y0m2_2 = sum(y0m2_2[i]*MatM[i] for i in range(MatM.nrows()))

#A-terms
F3_a02_2 = sum([a02_2[i]*MatM[i] for i in range(MatM.nrows())])
F3_a0m2_2 = sum([a0m2_2[i]*MatM[i] for i in range(MatM.nrows())])

N = Newforms(Gamma1(16),3,names='a')
f = N[0]
f8 = Newforms(Gamma1(8),3)[0]

chi0 = DirichletGroup(1)[0]
chi1,chi2, chi3, chi4 = DirichletGroup(4)[1], DirichletGroup(8)[3], DirichletGroup(16)[3], DirichletGroup(16)[7] #Odd characters
a1 = eis_phipsi(chi0,chi1,3, prec = 100).padded_list();
a2 = eis_phipsi(chi0,chi1,3, t=2, prec = 100).padded_list();
a3 = eis_phipsi(chi0,chi1,3, t=4, prec = 100).padded_list();
b1 = eis_phipsi(chi0,chi2,3, prec = 100).padded_list();
b2 = eis_phipsi(chi0,chi2,3, t=2, prec = 100).padded_list();
c1 = eis_phipsi(chi0,chi3,3, prec = 100).padded_list();
d1 = eis_phipsi(chi0,chi4,3, prec = 100).padded_list();
A1 = eis_phipsi(chi1,chi0,3, prec = 100).padded_list();
A2 = eis_phipsi(chi1,chi0,3, t=2, prec = 100).padded_list();
A3 = eis_phipsi(chi1,chi0,3, t=4, prec = 100).padded_list();
B1 = eis_phipsi(chi2,chi0,3, prec = 100).padded_list();
B2 = eis_phipsi(chi2,chi0,3, t=2, prec = 100).padded_list();
C1 = eis_phipsi(chi3,chi0,3, prec = 100).padded_list();
D1 = eis_phipsi(chi4,chi0,3, prec = 100).padded_list();
# The modular form f8(z)
f8_1 = [0]+f8.coefficients(100-1)
# The modular form f8(2z)
f8_2 = [0]*100
for i in range(1,50):
    f8_2[2*i] = f8_1[i]
# The modular form f
f_coeffs = [0]+f.coefficients(100-1)
L16 = [a1,a2,a3,b1,b2,c1,d1,A1,A2,A3,B1,B2,C1,D1,f8_1,f8_2,f_coeffs];
L16 = Matrix(L16)

print 'Resulting modular form as a vector in the basis'
v_x02_1 = L16.T \vector(F3_x02_1)
v_x02_2 = L16.T \vector(F3_x02_2)
v_x0m2_1 = L16.T \vector(F3_x0m2_1)
v_x0m2_2 = L16.T \vector(F3_x0m2_2)
v_y02_1 = L16.T \vector(F3_y02_1)
v_y02_2 = L16.T \vector(F3_y02_2)
v_y0m2_1 = L16.T \vector(F3_y0m2_1)
v_y0m2_2 = L16.T \vector(F3_y0m2_2)
print 'Integral over X{0,2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*X02_1.pre_factor()*v_x02_1
print 'Integral over X{0,2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*X02_2.pre_factor()*v_x02_2
print 'Integral over X{0,-2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*X0m2_1.pre_factor()*v_x0m2_1
print 'Integral over X{0,-2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*X0m2_2.pre_factor()*v_x0m2_2
print 'Integral over Y{0,2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*Y02_1.pre_factor()*v_y02_1
print 'Integral over Y{0,2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*Y02_2.pre_factor()*v_y02_2
print 'Integral over Y{0,-2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*Y0m2_1.pre_factor()*v_y0m2_1
print 'Integral over Y{0,-2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*Y0m2_2.pre_factor()*v_y0m2_2

print 'Determining the weight 1 part'

# Weight 1 part of the quasi-modular form
# This is an Eisenstein series of weight 1 for Gamma_1(16)
F1_x02_1 = sum(x02_1[MatM.nrows()+i]*MatE[i] for i in range(MatE.nrows()))
F1_x02_2 = sum(x02_2[MatM.nrows()+i]*MatE[i] for i in range(MatE.nrows()))
F1_x0m2_1 = sum(x0m2_1[MatM.nrows()+i]*MatE[i] for i in range(MatE.nrows()))
F1_x0m2_2 = sum(x0m2_2[MatM.nrows()+i]*MatE[i] for i in range(MatE.nrows()))
F1_y02_1 = sum(y02_1[MatM.nrows()+i]*MatE[i] for i in range(MatE.nrows()))
F1_y02_2 = sum(y02_2[MatM.nrows()+i]*MatE[i] for i in range(MatE.nrows()))
F1_y0m2_1 = sum(y0m2_1[MatM.nrows()+i]*MatE[i] for i in range(MatE.nrows()))
F1_y0m2_2 = sum(y0m2_2[MatM.nrows()+i]*MatE[i] for i in range(MatE.nrows()))

F1_a02_2 = sum([a02_2[MatM.nrows()+i]*MatE[i] for i in range(MatE.nrows())])
F1_a0m2_2 = sum([a0m2_2[MatM.nrows()+i]*MatE[i] for i in range(MatE.nrows())])

a1 = eis_phipsi(chi0,chi1,1, prec = 100).padded_list();
a2 = eis_phipsi(chi0,chi1,1, t=2, prec = 100).padded_list();
a3 = eis_phipsi(chi0,chi1,1, t=4, prec = 100).padded_list();
b1 = eis_phipsi(chi0,chi2,1, prec = 100).padded_list();
b2 = eis_phipsi(chi0,chi2,1, t=2, prec = 100).padded_list();
c1 = eis_phipsi(chi0,chi3,1, prec = 100).padded_list();
d1 = eis_phipsi(chi0,chi4,1, prec = 100).padded_list();
L16_2 = [a1,a2,a3,b1,b2,c1,d1];
L16_2 = Matrix(L16_2)

print 'Resulting modular form as a vector in the basis'
w_x02_1 = L16_2.T \vector(F1_x02_1)
w_x02_2 = L16_2.T \vector(F1_x02_2)
w_x0m2_1 = L16_2.T \vector(F1_x0m2_1)
w_x0m2_2 = L16_2.T \vector(F1_x0m2_2)
w_y02_1 = L16_2.T \vector(F1_y02_1)
w_y02_2 = L16_2.T \vector(F1_y02_2)
w_y0m2_1 = L16_2.T \vector(F1_y0m2_1)
w_y0m2_2 = L16_2.T \vector(F1_y0m2_2)
print 'Integral over X{0,2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*X02_1.pre_factor()*w_x02_1
print 'Integral over X{0,2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*X02_2.pre_factor()*w_x02_2
print 'Integral over X{0,-2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*X0m2_1.pre_factor()*w_x0m2_1
print 'Integral over X{0,-2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*X0m2_2.pre_factor()*w_x0m2_2
print 'Integral over Y{0,2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*Y02_1.pre_factor()*w_y02_1
print 'Integral over Y{0,2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*Y02_2.pre_factor()*w_y02_2
print 'Integral over Y{0,-2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*Y0m2_1.pre_factor()*w_y0m2_1
print 'Integral over Y{0,-2} is the L-value at s=0 of the modular form', -1/(4*pi^2)*Y0m2_2.pre_factor()*w_y0m2_2




# Using Eisenstein symbols of level 8

DivZ=Eij(1,5,8) + Eij(5,1,8) - Eij(1,1,8) - Eij(5,5,8)

bet8=Matrix(8,8)
bet8[2,4]=1
bet8[6,0]=1

ProdD = EisensteinSymbol(0,1,8,tensorprod(DivZ, bet8, 8));

ProdD.integrate(prec=200)
F=ProdD.A(prec=200)
vF=F.padded_list()

vG=100*[0]
for i in range(100):
    vG[i]=vF[2*i]

M = ModularForms(Gamma1(16),3).basis()
MatM = Matrix(QQ,[m.coefficients(srange(100)) for m in M])

try:
    x = MatM.T \vector(vG)
    print 'Mahler measure is the L-value at s=0  of the modular form', -1/(4*pi^2)*ProdD.pre_factor()*x
except ValueError:
    print 'res is not a modular form'

print 'Determining the modular form'

N = Newforms(Gamma1(16),3,names='a')
f = N[0]
f8 = Newforms(Gamma1(8),3)[0]

chi0 = DirichletGroup(1)[0]
chi1,chi2, chi3, chi4 = DirichletGroup(4)[1], DirichletGroup(8)[3], DirichletGroup(16)[3], DirichletGroup(16)[7] #Odd characters
a1 = eis_phipsi(chi0,chi1,3, prec = 100).padded_list();
a2 = eis_phipsi(chi0,chi1,3, t=2, prec = 100).padded_list();
a3 = eis_phipsi(chi0,chi1,3, t=4, prec = 100).padded_list();
b1 = eis_phipsi(chi0,chi2,3, prec = 100).padded_list();
b2 = eis_phipsi(chi0,chi2,3, t=2, prec = 100).padded_list();
c1 = eis_phipsi(chi0,chi3,3, prec = 100).padded_list();
d1 = eis_phipsi(chi0,chi4,3, prec = 100).padded_list();
A1 = eis_phipsi(chi1,chi0,3, prec = 100).padded_list();
A2 = eis_phipsi(chi1,chi0,3, t=2, prec = 100).padded_list();
A3 = eis_phipsi(chi1,chi0,3, t=4, prec = 100).padded_list();
B1 = eis_phipsi(chi2,chi0,3, prec = 100).padded_list();
B2 = eis_phipsi(chi2,chi0,3, t=2, prec = 100).padded_list();
C1 = eis_phipsi(chi3,chi0,3, prec = 100).padded_list();
D1 = eis_phipsi(chi4,chi0,3, prec = 100).padded_list();
# The modular form f8(z)
f8_1 = [0]+f8.coefficients(100-1)
# The modular form f8(2z)
f8_2 = [0]*100
for i in range(1,50):
    f8_2[2*i] = f8_1[i]
# The modular form f
f_coeffs = [0]+f.coefficients(100-1)
L16 = [a1,a2,a3,b1,b2,c1,d1,A1,A2,A3,B1,B2,C1,D1,f8_1,f8_2,f_coeffs];
L16 = Matrix(L16)

print 'Resulting modular form as a vector in the basis'
v3 = L16.T \vector(vG)
print 'Integral over X{0,oo} is the L-value at s=0 of the modular form', -1/(4*pi^2)*ProdD.pre_factor()*v3


#-----> We get:
print 'Integral over X{0,oo} is the L-value at s=0 of the modular form', (3/2048, -3/2048, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3/1024)

