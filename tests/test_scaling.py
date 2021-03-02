import sys
import numpy as np
from src.scaling import gmscale, apply_scaling 


A = np.random.choice(np.arange(-3, 3), p=[0.05, 0.05, 0.3, 0.5, 0.1, 0.0], size=(10,4))
print("\n This is the A matrix: \n")
print(A)

b = np.ones(10, dtype=np.float)
print("\n This is the vector b: \n")
print(b)

res = gmscale(A, 0, 0.99)
print("cscale:")
print(res[0])
print("rscale:")
print(res[1])

res = apply_scaling(A, b, res[0], res[1])

print("new A:")
print(res[0])
print("new b:")
print(res[1])

