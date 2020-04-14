#!/usr/bin/python3

from scipy.optimize import brentq

MAX_T = 30

def f(x,T):
    return ((1+x)**(T-1))*x - 1

a = 0 # negative
b = 1 # positive for T>1

print ("Critical values for the product hp")

for j in range(2, MAX_T + 1):
    print (j, brentq(f, a, b, float(j) ) )
