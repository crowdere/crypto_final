
import math
import time
import random

#Algorithm
#Input: N, e, d. 
#Output: p and q where pq=N.
#[Initialize] Set k←de−1.
#[Try a random g] Choose g at random from {2,…,N−1} and set t←k.
#[Next t] If t is divisible by 2, set t←t/2 and x←gtmodN. Otherwise go to step 2.
#[Finished?] If x>1 and y=gcd(x−1,N)>1 then set p←y and q←N/y, output (p,q) and terminate the algorithm. Otherwise go to step 3.

def prime_factors(n):
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    return factors
    

print(prime_factors(18721))
