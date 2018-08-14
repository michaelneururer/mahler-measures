#Horizontal divisors
DivX = (-1,0,-1,0,1,0,1,0)
DivY = (-1 ,0 ,1 ,0 ,1 ,0 ,-1,0 )
DivXm1 = (-1 ,0 ,-1 ,1 ,0 ,0 ,0, 1 )
DivYm1 = (-1 ,1 ,0 ,0 ,0 ,1 ,-1,0 )
DivXpY = (0 ,0 ,-1 ,0 ,2 ,0 ,-1,0 )
DivXpZY = (-1 ,1 ,-1 ,1 ,1 ,0 ,-1,0 )
DivXpZ = (-1 ,2 ,-1 ,0 ,0 ,0 ,0,0 )
DivYpZ = ( -1, 0, 0, 0, 0, 0, -1,2)
DivXp1oZ = (-1 ,0 ,-1 ,0 ,0 ,2 ,0,0 )
DivYp1oZ = (-1 ,0 ,0 ,2 ,0 ,0 ,-1,0 )
#Generating divisors: X,Y,X-1,Y-1,Y+Z, X+Y, X+ZY
gens = [DivX, DivY, DivXm1, DivYm1, DivYpZ, DivXpY, DivXpZY]
DivXmA = (0,-1,0,-1,0,1,0,1)
DivYmA = (0,-1,0,1,0,1,0,-1)

#Calculate tame symbols 
Pol.<Z> = PolynomialRing(QQ)
#10000 stands for infinity... in the following code this will always be taken to the power of zero. It was just chosen as a number because 'Inf'^0 should be 1 but produces an error.
fcn_X = [10000,-Z,10000,1,0,-1/Z,0,1]
fcn_Y = [10000,1,0,-1/Z,0,1,10000,-Z]
fcn_XpY = [-(Z-1)^2/Z, 1-Z, 10000,1-1/Z,0,1-1/Z,10000,1-Z]
def tame_symb(a,b, k,Zexp1=0,Zm1exp1=0,Zp1exp1=0,Zexp2=0,Zm1exp2=0,Zp1exp2=0):
    """
    finds the tame symbol at kA associated to two horizontal functions given as vectors a = (a1, ..., a7) and b=(b1,...b7), where ai and bi correspond to the exponent of the i-th     generator above. dexp1 and dm1exp1 are the exponents of d and d-1 in the first input function and dexp2 and dexp2 for the second input function.
    """
    vx = sum([gens[i][k]*a[i] for i in range(7)])
    vy = sum([gens[i][k]*b[i] for i in range(7)])
    exps = [a[i]*vy-b[i]*vx for i in range(7)]
    prod = (-1)^(vx*vy)
    #print vx,vy,exps
    #Take care of problems regarding division by 0 or multiplication by oo.
    if k==0:
    #Problems come from all functions except X+Y. X/Y = X/(Y-1) = X/(Y+Z) = -X/(X-1)=-1.
    #X/(X+YZ) = 1/(1-Z)
        prod *= (-1)^(exps[1]+exps[3]+exps[4])*(1-Z)^(exps[6])
        exps[0]=0
        exps[1]=0
        exps[2]=0
        exps[3]=0
        exps[4]=0
        exps[6]=0
    if k==1:
    #Problems come from Y-1,X+YZ. (Y-1)/(X+YZ) = 1/Z
        prod *= Z^(exps[6])
        exps[6]=0
        exps[3]=0
    if k==2:
    #Problems come from all functions except Y-1 and Y+Z. XY = -1
    #X/(X-1) = X/(X+Y) = X/(X+YZ) = 1
        prod *= (-1)^exps[1]
        exps[0]=0
        exps[1]=0
        exps[2]=0
        exps[5]=0
        exps[6]=0
    if k==3:
    #Problems from X-1, X+YZ. (X-1)*(X+YZ)=1
        exps[2] = 0
        exps[6] = 0   
    if k==4:
        #Problems from X,Y,X+Y,X+YZ. X/Y = -1, (X+Y)/X^2 = (Z-1)^2/Z
        #(X+YZ)/X = 1-Z
        prod *= (-1)^exps[1] * ((Z-1)^2/Z)^exps[5] * (1-Z)^exps[6]
        exps[0] = 0
        exps[1] = 0
        exps[5] = 0
        exps[6] = 0
    if k==5:
        #Problems from Y-1, so we can just omit it.
        exps[3] = 0                  
    if k==6:
        #Problems from all except X-1. XY = X(Y-1)=X(Y+Z)=X(X+Y) = -1 and X(X+YZ)= -Z
        prod *= (-1)^(exps[1]+exps[3]+exps[4]+exps[5])*(-Z)^exps[6]
        exps[0]=0
        exps[1]=0
        exps[3]=0
        exps[4]=0
        exps[5]=0
        exps[6]=0
    if k==7:
        #Problems from X-1, Y+Z. (Y+Z)/(X-1)^2 = -Z^2/(Z^2-1)
        prod*= (-Z^2/(Z^2-1))^exps[4]
        exps[4] = 0
        exps[2] = 0   
    prod *= fcn_X[k]^exps[0] * fcn_Y[k]^exps[1] * (fcn_X[k]-1)^exps[2] * (fcn_Y[k]-1)^exps[3] * (fcn_Y[k]+Z)^exps[4] * (fcn_XpY[k])^exps[5]*(fcn_X[k]+fcn_Y[k]*Z)^exps[6]
    prod *= Z^(Zexp1*vy - Zexp2*vx) * (Z-1)^(Zm1exp1*vy - Zm1exp2*vx) * (Z+1)^(Zp1exp1*vy - Zp1exp2*vx)
    return prod

def product_rule(a,b,Zexp1=0,Zm1exp1=0,Zp1exp1=0,Zexp2=0,Zm1exp2=0,Zp1exp2=0):
    return prod([tame_symb(a,b,k,Zexp1,Zm1exp1,Zp1exp1,Zexp2,Zm1exp2,Zp1exp2) for k in range(8)])
    
def proj_eq(P,Q):
    for i in range(len(P)):
        if P[i] != 0:
            if Q[i] == 0:
                return False
            if [P[0]/P[i], P[1]/P[i], P[2]/P[i]] == [Q[0]/Q[i], Q[1]/Q[i], Q[2]/Q[i]]:
                return True
            else:
                return False
                
PPol.<X,Y,Z> = PolynomialRing(QQ)
PFrac = FractionField(PPol)
I  = PPol.ideal(X^2*Y*Z + Y*Z + Y^2*Z*X + Z*X+Z^2*X*Y + X*Y -2)
PQuo=PPol.quotient_ring(I)
def convert_to_fcn(a, Zexp=0,Zm1exp=0,Zp1exp=0):
    #converts exponents to a function
    f = Z^Zexp * (Z-1)^Zm1exp * (Z+1)^Zp1exp * X^a[0]*Y^a[1]*(X-1)^a[2]*(Y-1)^a[3]*(Y+Z)^a[4]*(X+Y)^a[5]*(X+Y*Z)^a[6]
    return f
def convert_to_fcns(a,b,Zexp1=0,Zm1exp1=0,Zp1exp1=0,Zexp2=0,Zm1exp2=0,Zp1exp2=0):
    f = Z^Zexp1 * (Z-1)^Zm1exp1 * (Z+1)^Zp1exp1 * X^a[0]*Y^a[1]*(X-1)^a[2]*(Y-1)^a[3]*(Y+Z)^a[4]*(X+Y)^a[5]*(X+Y*Z)^a[6]
    g = Z^Zexp2 * (Z-1)^Zm1exp2 * (Z+1)^Zp1exp2 * X^b[0]*Y^b[1]*(X-1)^b[2]*(Y-1)^b[3]*(Y+Z)^b[4]*(X+Y)^b[5]*(X+Y*Z)^b[6]
    return [f,g]
#Check Steinberg
def check_Steinberg(F,G,Zexp1=0,Zm1exp1=0,Zp1exp1=0,Zexp2=0,Zm1exp2=0,Zp1exp2=0):
    if PQuo((convert_to_fcn(F,Zexp1,Zm1exp1,Zp1exp1)+convert_to_fcn(G,Zexp2,Zm1exp2,Zp1exp2)-1).numerator())==0:
        return True
    else:
        return False
