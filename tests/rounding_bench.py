#!/usr/bin/python3

import os, sys, time, getopt
import numpy as np
import pickle
from PolyRound.api import PolyRoundApi
from PolyRound.static_classes.lp_utils import ChebyshevFinder
from PolyRound.settings import PolyRoundSettings
#import hopsy
import dingo
from dingo import MetabolicNetwork, PolytopeSampler

from dingo import set_default_solver
from volestipy import HPolytope
from scipy.linalg import eigh

def _ess_min(samples):
    """Minimum effective sample size across dimensions using Geyer's initial sequence estimator.

    samples: ndarray of shape (n_dims, n_samples)  (as returned by dingo samplers)
    Returns the minimum ESS over all dimensions.
    """
    # samples shape: (d, N) – transpose to (N, d)
    X = samples.T
    N, d = X.shape
    ess_vals = []
    for j in range(d):
        x = X[:, j] - X[:, j].mean()
        # Normalised autocorrelation via FFT
        n = len(x)
        fft_size = 1
        while fft_size < 2 * n:
            fft_size <<= 1
        f = np.fft.rfft(x, n=fft_size)
        acf_full = np.fft.irfft(f * np.conj(f))[:n]
        acf = acf_full / acf_full[0]   # normalise so lag-0 = 1
        # Initial positive sequence: sum pairs (lag 2k, 2k+1) while positive
        rho_sum = 0.0
        for k in range(1, n // 2):
            pair = acf[2 * k - 1] + acf[2 * k]
            if pair < 0:
                break
            rho_sum += pair
        eff_n = N / (1.0 + 2.0 * rho_sum)
        ess_vals.append(min(eff_n, N))
    return min(ess_vals)


def evaluate_rounding_quality(A, b, T_matrix=None, walk_method='billiard_walk', n_samples=None):
    """Compute three quality metrics for a rounded polytope {x : A x <= b}.

    Metrics
    -------
    1. **T_cond** – condition number of the rounding transformation T_matrix
       (σ_max / σ_min).  Close to 1 → near-isotropic transformation.
       Only computed when T_matrix is not None.
    2. **cov_ratio** – ratio of the largest to smallest eigenvalue of the
       empirical covariance of the samples drawn from the rounded polytope.
       Close to 1 → well-rounded (isotropic distribution).
    3. **ess_min** – minimum effective sample size (ESS) across all dimensions
       after a fixed number of billiard-walk steps.  Higher → better mixing.
    """
    d = A.shape[1]
    if n_samples is None:
        n_samples = max(200, 10 * d)

    burn_in  = int(5 * np.sqrt(d))
    thinning = max(1, burn_in)

    print(f"  [quality] sampling {n_samples} points in R^{d} (burn={burn_in}, thin={thinning}) …")
    samples = PolytopeSampler.sample_from_polytope_no_multiphase(
        A, b,
        method=walk_method,
        n=n_samples,
        burn_in=burn_in,
        thinning=thinning
    )
    # samples: (d, n_samples)

    # --- metric 1: T condition number ---
    if T_matrix is not None:
        sv = np.linalg.svd(T_matrix, compute_uv=False)
        sv = sv[sv > 1e-12]
        t_cond = sv.max() / sv.min() if sv.size > 0 else float('inf')
    else:
        t_cond = None

    # --- metric 2: covariance eigenvalue ratio ---
    X = samples.T                          # (n_samples, d)
    X_centered = X - X.mean(axis=0)
    cov = np.cov(X_centered, rowvar=False)
    eigvals = eigh(cov, eigvals_only=True) # sorted ascending
    cov_ratio = eigvals[-1] / eigvals[0] if eigvals[0] > 1e-14 else float('inf')

    # --- metric 3: minimum ESS ---
    ess = _ess_min(samples)

    return {
        "T_cond":    t_cond,    # None when no T_matrix supplied
        "cov_ratio": cov_ratio, # 1 = perfect isotropy
        "ess_min":   ess,       # higher = better mixing
    }


def test_rounding(rounding_method, transformed_polytope, name):
    if rounding_method == "PolyRound":
        start = time.time()
        rounded_polytope = PolyRoundApi.round_polytope(transformed_polytope)
        end   = time.time()
        A = rounded_polytope.A.to_numpy()
        b = rounded_polytope.b.to_numpy()
        print("Polytope derived from the " + name + " network, took "
              + str(end - start) + " sec to get rounded with PolyRound.")
        result = evaluate_rounding_quality(A, b, T_matrix=None)
    else:
        A = transformed_polytope.A.to_numpy()
        b = transformed_polytope.b.to_numpy()
        P = HPolytope(A, b)
        start = time.time()
        A_rounded, b_rounded, T_matrix, shift, round_value = P.rounding(rounding_method, None)
        end   = time.time()
        print("Polytope derived from the " + name + " network, took "
              + str(end - start) + " sec to get rounded with " + rounding_method + ".")
        result = evaluate_rounding_quality(A_rounded, b_rounded, T_matrix=T_matrix)

    # Print quality metrics
    print(f"  Quality metrics for {rounding_method}:")
    if result['T_cond'] is not None:
        print(f"    T condition number (σ_max/σ_min) : {result['T_cond']:.4g}  (1 = perfect)")
    print(f"    Cov eigenvalue ratio (λ_max/λ_min): {result['cov_ratio']:.4g}  (1 = isotropic)")
    print(f"    Min ESS across dimensions         : {result['ess_min']:.1f}")
    print()

def polyround_preprocess(model_path):

    print("Starting PolyRound preprocessing for " + model_path)

    name = model_path.split("/")[-1]

    # Import model and create Polytope object
    polytope = PolyRoundApi.sbml_to_polytope(model_path)
    print("Polyope for network " + name + " was built.")

    # Make a settings object for the polyround library - optional

    s = PolyRoundSettings()
    s.backend = 'gurobi'
    s.verbose = True
    s.check_lps = False  # Disable extra LP checks for better performance

    # Simplify the polytope
    start = time.time()
    simplified_polytope = PolyRoundApi.simplify_polytope(polytope, settings=s)
    end   = time.time()
    time_for_simplification = end - start
    print("Polytope derived from the " + name + " network, took " + str(time_for_simplification) + 
          " sec to get simplified.")

    # Polytope transformation
    start = time.time()
    transformed_polytope = PolyRoundApi.transform_polytope(simplified_polytope, settings=s)
    end   = time.time()
    time_for_transformation = end - start
    print("Polytope derived from the " + name + " network, took " + str(time_for_transformation) + 
          " sec to get transformed.")

    # Export simplified and transformed polytope as pickle file
    polytope_info = (
        transformed_polytope,
        name,
    )

    # Ensure output directory exists
    os.makedirs("simpl_transf_polytopes", exist_ok=True)
    
    with open(
        "simpl_transf_polytopes/polytope_" + name + ".pckl", "wb"
    ) as polyround_polytope_file:
        pickle.dump(polytope_info, polyround_polytope_file)

    test_rounding("PolyRound", transformed_polytope, name)
    test_rounding("isotropic_position", transformed_polytope, name)
    test_rounding("min_ellipsoid", transformed_polytope, name)
    test_rounding("john_position", transformed_polytope, name)
    test_rounding("log_barrier", transformed_polytope, name)
    test_rounding("volumetric_barrier", transformed_polytope, name)
    test_rounding("vaidya_barrier", transformed_polytope, name)

    
    # Export rounded polytope as pickle file
    #polytope_info = (
    #    rounded_polytope,
    #    name,
    #)

    #with open(
    #    "polyrounded_polytopes/polytope_" + name + ".pckl", "wb"
    #) as polyround_polytope_file:
    #    pickle.dump(polytope_info, polyround_polytope_file)

    # return rounded_polytope, name
    return


if __name__ == '__main__':
    # Wrap main execution in try/except to surface errors during runs
    try:
        # Enable Gurobi as the default solver for both PolyRound and dingo rounding methods
        set_default_solver("gurobi")
        
        if len(sys.argv) < 2:
            print("Usage: python tests/rounding_bench.py <model_filename>")
            sys.exit(1)

        current_directory = os.getcwd()
        network_name = sys.argv[1]
        dingo_directory = '/'.join(current_directory.split("/")[:-1])

        path_to_net = dingo_directory + "/ext_data/" + network_name
        print("path_to_net:", path_to_net)
        polyround_preprocess(path_to_net)
    except Exception as e:
        import traceback
        traceback.print_exc()
        sys.exit(1)