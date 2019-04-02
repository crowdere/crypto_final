import matplotlib.pyplot as plt
import numpy as np
import vectors
import random
import hashlib

x_pts=[]
y_pts=[]
A = 0
B = 0
P = 0

def dsadecrypt(hm,r,s,g,Q_):
    n = len(x_pts)
    #compute inverse of s mod n
    s_inv = 0
    for x in range(n):
        if x * s % n == 1:
            s_inv = x
    print('S_INV -> ',s_inv)
    #compute w
    w = s_inv % n
    #compute u
    u = (hm * w) % n
    #compute v
    v = (r * w) % n
    #compute X
    uG = addp2p(u,g)
    vQ = addp2p(v,Q_)
    result = add_ecp(uG,vQ)
    print(result.x % n)

def dsa():
    #generate hash of message
    m = 184921
    hm = hash(m) % 256
    print('H(M) ->',hm)
    #calculate n (number of points on plot)
    n = len(x_pts)
    print('N -> {0}'.format(n))
    #find random point on plot
    g = vectors.Point(3,10,0)
    #SENDER selects private key d
    d = 3
    print('D ->',d)
    #choose random K value
    k = random.randint(1,24)
    print('K ->',k)
    #compute dG (Q)
    Q_ = addp2p(d,g)
    print('Q -> ',Q_)
    #compute kG (P)
    P_ = addp2p(k,g)
    print('KG ->', P_)
    #compute r
    r = P_.x % n
    print('R ->', r)
    #compute s
    k_inv = 0
    #find inverse of k mod n
    for x in range(n):
        if x * k % n == 1:
            k_inv = x
    print('K_INV ->',k_inv)
    s = k_inv * (hm + (d*r)) % n
    print('S ->',s)
    dsadecrypt(hm,r,s,g,Q_)

def addp2p(i,p):
    q = p
    r = add_ecp(p,q)
    p_list = []
    p_list.append(p)
    print('{0} {1}'.format('1',p))
    for x in range(1,i):
        y = add_ecp(p_list[x-1],p_list[x-x])
        p_list.append(y)
        print(x+1,y)
    return p_list[i-1]

def checkPoint(p):
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
        Yr = get_Yr(λ, p.x, Xr, p.y)
        r = vectors.Point(Xr,Yr,0)
        if(checkPoint(r)):
            #annotation = "({:d},{:d}) + ({:d},{:d}) = ({:d},{:d})".format(xp, yp, xq,yq,Xr,Yr)
            #print(r.x)
            #print(r.y)
            plt.plot(r.x, r.y, "or")
            return r
        else:
            print("Something went wrong....")
    else:
        return False

# Deterime the delta of two points.
def get_delta(p,q):
    if (p.x == q.x and p.y == q.y):
        #calculate same point delta using differiential calculus
        num = (3 * pow(p.x, 2) + A) % P
        den = (2*p.y)
        r2 = get_Inverse(den)
        λ = num*(r2) % P
        #λ = ( (3* pow(p.x,2) + A) / (2*p.y))
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
    step1 = (int(round(xp) - int(round(xr))))
    step2 = (λ * step1)
    step3 = (step2 - int(round(yp)))
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
plot_ec(23,1,1)
#
#p = vectors.Point(3,3,0)
#q = vectors.Point(3,8,0)
#print(add_ecp(p,q))   # Return same
p = vectors.Point(0,10,0)
#addp2p(6,p)
dsa()
#display the graphs
plt.show()