import matplotlib.pyplot as plt
import numpy as np

x_pts=[]
y_pts=[]
A = 0;
B = 0;
P = 0;

def checkPoint(px,py):
    for x in range(len(x_pts)):
        if (x_pts[x] == px and y_pts[x] == py):
            return (x_pts[x] == px and y_pts[x] == py)
    return False
    
# validates and adds two points
# returns resulting point
def add_ecp(xp,yp,xq,yq):
    if (checkPoint(xp,yp) and checkPoint(xq,yq)):
        λ = get_delta(xp,yp,xq,yq)
        Xr = get_Xr(λ, xp, xq)
        Yr = get_Yr(λ, xp, Xr, yp)
        if(checkPoint(Xr,Yr)):
            #annotation = "({:d},{:d}) + ({:d},{:d}) = ({:d},{:d})".format(xp, yp, xq,yq,Xr,Yr)
            print(Xr)
            print(Yr)
            plt.plot(Xr, Yr, "or")
        else:
            print("Something went wrong....")
    else:
        return False

# Deterime the delta of two points.
def get_delta(xp,yp,xq,yq):
    if (xp == xq and yp == yq):
        #calculate same point delta using differiential calculus
        λ = ( (3* pow(xp,2) + A) / (2*yp))
    else:
        #calculate delta with rise over run
        num = (yq - yp)
        den = get_Inverse((xq - xp))
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
add_ecp(8,8,13,11)   # Return same
#display the graphs
plt.show()