def paramXYZ(tau,z):
    w = pari([tau,1])
    X = pari( (w.ellwp(z)-w.ellwp(1/2))*(w.ellwp(1/8)-w.ellwp(1/4))^2/((w.ellwp(1/8)-w.ellwp(1/2))*(w.ellwp(z+1/4)-w.ellwp(1/4))*(w.ellwp(z)-w.ellwp(1/4))) )
    Y=pari((w.ellwp(-z)-w.ellwp(1/2))*(w.ellwp(1/8)-w.ellwp(1/4))^2/((w.ellwp(1/8)-w.ellwp(1/2))*(w.ellwp(-z+1/4)-w.ellwp(1/4))*(w.ellwp(-z)-w.ellwp(1/4))))
    Z = pari((w.ellwp(1/8)-w.ellwp(1/2))/(w.ellwp(1/4)-w.ellwp(1/2)))
    return X,Y,Z

def TestXYZ(tau,z):
#Test if the values of X,Y,Z lie on the variety X+1/X+Y*1/Y+Z+1/Z = 2
    X,Y,Z = paramXYZ(tau,z)
    return X+1/X+Y+1/Y+Z+1/Z

def Cycle(n,m,t1,t2):
#Shokurov-Cycle between -1/2 and 1/2
    tau = 1/2*exp(I*pi*(1-t1))
    z = t2*(n*tau+m)
    return [tau,z]

def ShokurovCycle(cusp1, cusp2, n,m, t1, t2):
#X,Y,Z evaluated on this half-circle. Perhaps I will add functionality for other cusps... at the moment
#only the half-circles between -1/2 and 1/2 work and the vertical lines to infinity.
    if cusp2 == 'inf':
        tau = cusp1 + (1/(1-t1)-1)*I
    else:
        tau = 1/2*exp(I*pi*(1-t1))
    z = t2*(n*tau+m)
    X,Y,Z = paramXYZ(tau,z)
    return X,Y,Z

def CycleValues(cusp1,cusp2,n,m,N, absolutevalues=False):
#Produces a list of ShokurovCycle values.
    list = []
    absolute = []
    for t1 in range(1,N):
        for t2 in range(1,N):
            try:
                s = ShokurovCycle(-1/2,1/2,n,m,t1/N,t2/N)
            except PariError:
                s = 'Too close to singularity'
            list.append([s,t1,t2])
            if absolutevalues:
                if len(s) == 3:
                    absolute.append([abs(s[0]),abs(s[1]),abs(s[2]),t1,t2])
                else:
                    absolute.append([s,t1,t2])
    if absolute:
        return list, absolute
    return list

def circ(t1):
    return 1/2*e^(I*pi*((1-t1)))

#All kinds of functions follow that I defined to check the limits I calculated.

def X(t1,a):
    tau = circ(t1); w = pari([1,tau])
    z = (1/2+a*t1)*tau
    X = pari( (w.ellwp(z)-w.ellwp(1/2))*(w.ellwp(1/8)-w.ellwp(1/4))^2/((w.ellwp(1/8)-w.ellwp(1/2))*(w.ellwp(z+1/4)-w.ellwp(1/4))*(w.ellwp(z)-w.ellwp(1/4))) )
    return CC(X)

def Y(t1,a):
    tau = circ(t1); w = pari([1,tau])
    z = (1/2+a*t1)*tau
    X = pari( (w.ellwp(z)-w.ellwp(1/2))*(w.ellwp(1/8)-w.ellwp(1/4))^2/((w.ellwp(1/8)-w.ellwp(1/2))*(w.ellwp(z-1/4)-w.ellwp(1/4))*(w.ellwp(z)-w.ellwp(1/4))) )
    return CC(X)

def Xtest(t1,a):
    tau = circ(t1); w = pari([1,(-tau-1)/(2*tau+1)])
    z = (1/2+a*t1)*tau
    X = pari(e^(2*pi*I/(8*tau+4))*(w.ellwp(z/(2*tau+1))-w.ellwp(1/4/(2*tau+1))))
    return CC(X)
def Xtest1(t1,a):
    tau = circ(t1); w = pari([1,(-tau-1)/(2*tau+1)])
    z = (1/2+a*t1)*tau
    X = pari(e^(2*pi*I/(8*tau+4))*(w.ellwp(z/(2*tau+1))+pi^2/3))
    return CC(X)
def Xtest2(t1,a):
    tau = circ(t1); w = pari([1,(-tau-1)/(2*tau+1)])
    z = (1/2+a*t1)*tau
    X = pari(e^(2*pi*I/(16*tau+8))*(w.ellwp(1/8/(2*tau+1))+pi^2/3))
    return CC(X)
def Xtest3(t1,a):
    tau = circ(t1); w = pari([1,(-tau-1)/(2*tau+1)])
    z = (1/2+a*t1)*tau
    X = pari((w.ellwp((z+1/4)/(2*tau+1))-w.ellwp(1/4/(2*tau+1))))
    return CC(X)
def XoverYconj(t1,a):
    tau = circ(t1); w = pari([1,(-tau-1)/(2*tau+1)])
    z = (1/2+a*t1)*tau
    XoY = pari( (w.ellwp(-z+1/4)-w.ellwp(1/4))/(w.ellwp(z+1/4)-w.ellwp(1/4)).conj())
    return CC(XoY)

def limit1(a):
    fa = f(a)
    return -4*pi^2/(e^(2*pi*I*fa)-2+e^(-2*pi*I*fa))
def limit2(a):
    fa = f(a)
    return -4*pi^2*(e^(2*pi*I*fa)-e^(-2*pi*I*fa))
def limit3(a):
    fa = f(a)
    return -4*pi^2*(-e^(2*pi*I*fa)-2-e^(-2*pi*I*fa))
def f(a):
    return CC(1/4 + a/(2-2*I*pi))

#PLOTTING

from sage.plot.colors import float_to_html
from sage.plot.colors import rainbow

def pointplot():
    list = CycleValues(-1/2,1/2,1,0,50,  False)
    list2 = [[l[0][0],l[0][1]] for l in list]
    arglist = [[RR(arg(l[0])),RR(arg(l[1]))] for l in list2]

    return arglist

#Before using nicerplot calculate list = CycleValues(-1/2,1/2,1,0,50,  False)

def nicerplot():
    colors = [float_to_html(1/50*i,0,0) for i in range(50)]
    bigtestlist = [[[arg(CC(l[0][0])),arg(CC(l[0][1]))] for l in list if l[1]==i] for i in range(1,50)]
    #Translate to fit the square
    for b in bigtestlist:
        for a in b:
            if a[0]<=0 and a[1]<=0:
                a[0]+= 2*pi
                a[1]+= 2*pi
            elif a[0]<0 and a[1]>=0:
                a[0]+= 2*pi
            elif a[0]>=0 and a[1]<=0:
                a[1]+=2*pi
    plotlist = scatter_plot(bigtestlist[0], facecolor = colors[0])
    plotlist += scatter_plot([[0,-1]], facecolor = colors[0])
    for i in range(1,49):
        plotlist += scatter_plot(bigtestlist[i], facecolor = colors[i])
        plotlist += scatter_plot([[2*pi/50*i,-1]], facecolor = colors[i])
    return plotlist

#Implicit function for psi given a fixed value of Z and a = (phi-pi)/(psi-pi)

def find_phi(a,u):
#Example: L = [(find_phi(a,-0.5), ((find_phi(a,-0.5)-pi)*a+pi)) for a in srange(-5.1,5.1,0.2)]
    return find_root(cos(phi)+cos((phi-pi)*a+pi) ==u, max(0,pi*(1-abs(1/a))),min(2*pi,pi*(1+abs(1/a))))
