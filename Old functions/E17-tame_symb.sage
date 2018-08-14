DivX = [ [1,0,0,0,0,1,-2], [0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]]
DivY = [ [1,-2,1,0,0,0,0], [0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]]
DivZ = [ [0,1,-3,3,-1,0,0], [0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]]
DivXm1 = [ [0,1,0,0,1,0,-2], [0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]]
DivYm1 = [ [0,-2,0,1,0,0,1], [0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]]
DivXZm1 = [ [0,0,1,1,0,0,-2], [0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]]
DivYZm1 = [ [0,-2,0,0,1,1,0], [0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]]

#Generators of principal horizontal divisors: u, u − 1, du − 1, v, v − 1,dv − 1
gens = map(lambda x: x[0],[DivX, DivXm1, DivXZm1, DivY, DivYm1, DivYZm1])
M = Matrix(gens)
#Divisors of X - A and Y-A
DivXmA = [-2,1,0,0,0,0,1]
DivYmA = [0,1,-2,1,0,0,0]

#Calculate tame symbols 
Pol.<d> = PolynomialRing(QQ)
#10000 stands for infinity... in the following code this will always be taken to the power of zero. It was just chosen as a number because 'Inf'^0 should be 1 but produces an error.
fcn_u = [0,1,1/d,1/d,1,0, 10000]
fcn_v = [0,10000,0,1,1/d,1/d,1]
def tame_symb(a,b, k,dexp1=0,dm1exp1=0,dexp2=0,dm1exp2=0):
    """
    finds the tame symbol at kA associated to two horizontal functions given as a vectors a = (a1, ..., a6) and b=(b1,...b6), where ai and bi correspond to the exponent of the i-th     generator above. dexp1 and dm1exp1 are the exponents of d and d-1 in the first input function and dexp2 and dexp2 for the second input function.
    """
    vx = sum([gens[i][k]*a[i] for i in range(6)])
    vy = sum([gens[i][k]*b[i] for i in range(6)])
    exps = [a[i]*vy-b[i]*vx for i in range(6)]
    prod = (-1)^(vx*vy)
    #print vx,vy,exps
    #Take care of problems regarding division by 0 or multiplication by oo.
    if k==0:
    #Problems come from u/v = 0/0 = -1
        prod *= (-1)^(exps[0])
        exps[0]=0
        exps[3]=0
    if k==1:
    #Problems come from v, v-1, vd-1, u-1=0. v(u-1)^2 = (v-1)(u-1)^2 = d-1 and (vd-1)(u-1)^2 = d (d-1)
    #Convert all v, v-1 and vd-1 to u-1's... in the end u-1 should have exponent 0
        prod *= (d-1)^(exps[3]+exps[4]+exps[5])*d^exps[5]
        exps[1] += -2*(exps[3]+exps[4]+exps[5])
        if exps[1] != 0:
            return 'exp of u-1 should be zero now!'
        exps[3]=0
        exps[4]=0
        exps[5]=0
    if k==2:
    #Problems come from v and ud-1. (ud-1)/v=(d-1)^(-2)
        prod *= (d-1)^(2*exps[2])
        exps[3] += exps[2]
        exps[2] = 0
        if exps[3] != 0:
            return 'error'
    if k==3:
    #Problems from v-1, ud-1. (ud-v)/(v-1) = (d-1)/d
        prod *= ((d-1)/d)^(exps[2])
        exps[4] += exps[2]
        exps[2] = 0
        if exps[4] != 0:
            return 'error'    
    if k==4:
        #Problems from u-1, vd-1. (u-1)/(vd-1) = d / (d-1)
        prod *= (d/(d-1))^(exps[1])
        exps[5] += exps[1]
        exps[1] = 0
        if exps[5] != 0:
            return 'error'  
    if k==5:
        #Problems from u, vd-1. u/(vd-1) = 1/(d-1)^2
        prod *= (d-1)^(-2*exps[0])
        exps[5] += exps[0]
        exps[0] = 0
        if exps[5] != 0:
            return 'error'         
           
    if k==6:
    #Problems from u, u-1, ud-1, v-1. u(v-1)^2 = (u-1)(v-1)^2 = d-1 and (ud-1)(v-1)^2 = d(d-1)
        prod *= (d-1)^(exps[0]+exps[1]+exps[2])*d^exps[2]
        exps[4] += -2*(exps[0]+exps[1]+exps[2])
        if exps[4] != 0:
            return 'exp of v-1 should be zero now!'
        exps[0]=0
        exps[1]=0
        exps[2]=0    
    prod *= fcn_u[k]^(exps[0]) * (fcn_u[k]-1)^(exps[1]) * (fcn_u[k]*d-1)^(exps[2]) * fcn_v[k]^(exps[3]) * (fcn_v[k]-1)^(exps[4]) * (fcn_v[k]*d-1)^(exps[5])
    prod *= d^(dexp1*vy - dexp2*vx) * (d-1)^(dm1exp1*vy - dm1exp2*vx)
    return prod

def product_rule(a,b, dexp1=0,dm1exp1=0,dexp2=0,dm1exp2=0):
    return prod([tame_symb(a,b,k, dexp1,dm1exp1,dexp2,dm1exp2) for k in range(7)])
    
def proj_eq(P,Q):
    for i in range(len(P)):
        if P[i] != 0:
            if Q[i] == 0:
                return False
            if [P[0]/P[i], P[1]/P[i], P[2]/P[i]] == [Q[0]/Q[i], Q[1]/Q[i], Q[2]/Q[i]]:
                return True
            else:
                return False
PPol.<u,v,d> = PolynomialRing(QQ)
PFrac = FractionField(PPol)
I  = PPol.ideal(d*(d-1)*u*v- (u*v-u-v)*(1+d*(u*v-u-v)))
PQuo=PPol.quotient_ring(I)
def convert_to_fcn(a, dexp=0,dm1exp=0):
    #converts exponents to a function
    f = d^dexp * (d-1)^dm1exp * u^a[0]*(u-1)^a[1]*(u*d-1)^a[2]*v^a[3]*(v-1)^a[4]*(v*d-1)^a[5]
    return f
def convert_to_fcns(a,b,dexp1=0,dm1exp1=0,dexp2=0,dm1exp2=0):
    f = d^dexp1 * (d-1)^dm1exp1 * u^a[0]*(u-1)^a[1]*(u*d-1)^a[2]*v^a[3]*(v-1)^a[4]*(v*d-1)^a[5]
    g = d^dexp2 * (d-1)^dm1exp2 * u^b[0]*(u-1)^b[1]*(u*d-1)^b[2]*v^b[3]*(v-1)^b[4]*(v*d-1)^b[5]
    return [f,g]
#Check Steinberg
def check_Steinberg(F,G,dexp1=0,dm1exp1=0,dexp2=0,dm1exp2=0):
    if PQuo((convert_to_fcn(F,dexp1,dm1exp1)+convert_to_fcn(G,dexp2,dm1exp2)-1).numerator())==0:
        return True
    else:
        return False
