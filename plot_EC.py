import matplotlib.pyplot as plt
import numpy as np
import vectors
import random
import hashlib

x_pts=[]
y_pts=[]
N = 0
A = 0
B = 0
P = 0
#-----
DSA_M = ""
DSA_R = 0
DSA_S = 0
DSA_HM = 0

def addp2p(i,p):
    q = p
    r = add_ecp(p,q)
    p_list = []
    p_list.append(p)
    #print('{0} {1}'.format('1',p))
    for x in range(1,i):
        y = add_ecp(p_list[x-1],p_list[x-x])
        p_list.append(y)
        #print(x+1,y)
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
        if get_Inverse(den) is None:
            return False
        else:
            r2 = get_Inverse(den)
            λ = num*(r2) % P
            return  λ
            #λ = ( (3* pow(p.x,2) + A) / (2*p.y))
    else:
        #calculate delta with rise over run
        num = (q.y - p.y)
        if get_Inverse((q.x - p.x)) is None:
            return False
        else:
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
    if x == 0:
        return None
    for n in range(P):
        if n * x % P == 1:
            return n

def get_InverseN(x, N):
    if x == 0:
        return None
    for n in range(N):
        if n * x % N == 1:
            return n
#
#
def plot_ec(p,a,b ):
    #Set globals
    global A, B, P, N
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
    N = len(x_pts)
    return plt
#=====================================
#Digital Signatures

#input values for private key then generate private key
def generateKeyPair(d, G):
    if (d <= (len(x_pts) -1) and d >= 1):
        Q = addp2p(d,G)
        return Q
    else:
        return None

def hashs(m):
    result = ""
    for char in m:
        result += str(ord(char))
    return int(result)

def sendMessage(k, G, m, d):
    if (k <= (len(x_pts) -1) and k >= 1):
        hm = hashs(m) % 256
        p_ = addp2p(k, G)
        r = 0
        while r == 0:
            r = p_.x % P
            k += 1
            p_ = addp2p(k, G)
        dr = d*r
        hmdr = hm + dr
        if get_InverseN(k,N) is None:
            return None
        else:
            s = ( get_InverseN(k, N) * (hmdr) )
            s2 = s % N
            #set Global
            global DSA_HM, DSA_M, DSA_R, DSA_S
            DSA_HM = hm
            DSA_M = m
            DSA_R = r
            DSA_S = s2
            return True
    else:
        return None


def verifyMessage(s, hm, r,Q,G):
    if get_InverseN(s,N) is None:
            return None
    else:
        w2 = get_InverseN(s, N) 
        w = w2 % N
        u = (hm * w ) % N
        v = (r * w) % N
        
        uG = addp2p(u, G)
        vQ = addp2p(v, Q)

        return add_ecp(uG, vQ)



#====================================
# plot ec: Ep(A,B) format and set global A and B for other functions
#plot_ec(23,1,1)

#Protocol for B to sign a hashed message, h(hm), with signature (r,s)to A
plot_ec(23,1,1)
print("=========================")
print('Valid Points on E{0}({1},{2})  are:'.format(P,A,B))
print(x_pts)
print(y_pts)
print("=========================")
#                 G.x G.y
G = vectors.Point(3,10,0)
print('G selected as: {0}'.format(G))
#                   d  G
d = 1
Q = generateKeyPair(d, G)
print('Public Key Q = {0} is a result of Private key (d) = {1} and Base (G)= {2}'.format(Q,d,G))
#               k G       m     d
if (sendMessage(1, G, "hello", d) == True):
    print("=========================")
    print('m,(r,s) = {0},({1},{2})'.format(DSA_M, DSA_R, DSA_S))
    print("=========================")

vm = verifyMessage(DSA_S, DSA_HM, DSA_R, Q, G)
print('A computes X = uG + vG = {}'.format(vm))
print('M=\"{}\" is valid because: x({}) mod n({}) = {} == r!'.format(DSA_M, vm.x, N, DSA_R))

plt.show()





#   
#works: 3,10 --- 9,16 -- 3,13  -- 
#no:    7,11 --- 4,0 -- 12,19 --1,7
#print('{0},({1},{2})'.format(DSA_M, DSA_R, DSA_S))
#1,3,

# for x in x_pts:
#     try:
#         G = vectors.Point(x_pts[x], y_pts[x], 0)
#         Q = generateKeyPair(1,G)
#         if (sendMessage(1, G, "TESCO IS THE BEST", 1) == True):
#             print(verifyMessage(DSA_S, DSA_HM, DSA_R, Q, G), " + ", G)
#     except:
#         continue