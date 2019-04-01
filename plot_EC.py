import matplotlib.pyplot as plt
import numpy as np
import vectors

x_pts=[]
y_pts=[]
A = 0;
B = 0;
P = 0;

def addp2p(i,p,q):
    p = vectors.Point(15,6,0)
    print(checkPoint(p))

def checkPoint(p):
    print(p)
    for x in range(len(x_pts)):
        if (x_pts[x] == p.x and y_pts[x] == p.y):
            return (x_pts[x] == p.x and y_pts[x] == p.y)
    return False
    
# validates and adds two points
# returns resulting point
def add_ecp(p,q):
    if (checkPoint(p) and checkPoint(q)):
        λ = get_delta(p,q)
        Xr = get_Xr(λ, p.x, q.x)
        Yr = get_Yr(λ, p.x, Xr, q.y)
        r = vectors.Point(Xr,Yr,0)
        if(checkPoint(r)):
            #annotation = "({:d},{:d}) + ({:d},{:d}) = ({:d},{:d})".format(xp, yp, xq,yq,Xr,Yr)
            print(r.x)
            print(r.y)
            plt.plot(r.x, r.y, "or")
        else:
            print("Something went wrong....")
    else:
        return False

# Deterime the delta of two points.
def get_delta(p,q):
    if (p.x == q.x and p.y == q.y):
        #calculate same point delta using differiential calculus
        λ = ( (3* pow(p.x,2) + A) / (2*p.y))
    else:
        #calculate delta with rise over run
        num = (q.y - p.y)
        den = get_Inverse((q.x - p.x))
        λ = (num * den % P)
        #λ = ( (yq - yp / xq - xp) )
    return  λ

def get_Xr(λ, xp, xq):
    return ((pow(λ, 2) - xp -xq) % P)

def get_Yr(λ,xp,xr,yp):
    step1 = (xp - xr)
    step2 = (λ * step1)
    step3 = (step2 - yp)
    step4 = step3 % P
    return step4
    #return ((λ (xp - xr) -yp) % P)

def get_Inverse(x):
    for n in range(P):
        if n * x % P == 1:
            return n

#
#
def plot_ec(p,a,b ):
    #Set globals
    global A, B, P
    A = a
    B = b
    P = p
    for x in range(p):
        for y in range(p):
            LS=(y**2)%p
            RS=(x**3+a*x+b)%p
            if(LS==RS):
                x_pts.append(x)
                y_pts.append(y) 
    #print(x_pts)
    #print(y_pts)

    plt.plot(x_pts, y_pts, 'bs')
    plt.grid(True, axis='both')
    plt.grid(linestyle='solid')
    return plt

# plot ec: Ep(A,B) format and set global A and B for other functions
plot_ec(23,1,4)
#Print Globals
#print(P)
#print(A)
#print(B)
#add_ecp(8,8,13,11)   # Return same
addp2p(2,8,8)
#display the graphs
plt.show()