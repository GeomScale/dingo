#!/usr/bin/python3.6

import numpy as np
from src.fba import slow_fba 


m = 3
n = 6

A = np.zeros((2*n, n), dtype=np.float)
A[0:n] = np.eye(n)
A[n:] -=  np.eye(n,n, dtype=np.float)
print("\n This is the A matrix: \n")
print(A)

b = np.ones(2*n, dtype=np.float)
print("\n This is the vector b: \n")
print(b)

# Aeq = np.random.randint(-3, 3, size=(m, n))
Aeq = np.random.choice(np.arange(-3, 3), p=[0.05, 0.05, 0.3, 0.5, 0.1, 0.0], size=(m,n))
print("\n This is the Aeq matrix: \n")
print(Aeq)

beq = np.zeros(m)
print("\n This is the vector beq: \n")
print(beq)

obj_fun =  np.ones(n, dtype=np.float)

res = slow_fba(A, b, Aeq, beq, obj_fun)
print("optimal solution:")
print(res[0])
print("optimal value:")
print(res[1])

