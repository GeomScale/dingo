import numpy as np
from dingo import PolytopeSampler

# Definišemo kocku [-1,1]^3:
A = np.array([
    [ 1, 0, 0],
    [-1, 0, 0],
    [ 0, 1, 0],
    [ 0,-1, 0],
    [ 0, 0, 1],
    [ 0, 0,-1]
], dtype=float)
b = np.array([1, 1, 1, 1, 1, 1], dtype=float)

# Napravi sampler direktno iz (A, b)

# Generiši uzorke (na primer billiard_walk)
samples = PolytopeSampler.sample_from_polytope_no_multiphase(
        A, b, method = 'cdhr', n=1000, burn_in=100, thinning=1, variance=1.0, bias_vector=None, solver=None, ess=0
    )

print("Oblik matrice uzoraka:", samples.shape)
print("Min po osi:", samples.min(axis=1))
print("Max po osi:", samples.max(axis=1))

# Ako hoćeš plot
try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(samples[0,:], samples[1,:], samples[2,:], s=4, alpha=0.3)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title("Samples iz [-1,1]^3 kocke")
    plt.show()
except ImportError:
    pass
