# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of the GeomScale project

import unittest
import os
import sys
import time
import warnings
import gc
import numpy as np

from dingo import MetabolicNetwork, PolytopeSampler
from dingo.pyoptinterface_based_impl import set_default_solver

# --- environment safety ---
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

# === Global parameters ===
N_SAMPLES = 50000
BURN_IN   = 200
NREFL     = None   # will default to ceil(sqrt(d)) inside generate_steady_states_sb_once

def _print_diag(diag: dict, elapsed_s: float, tag: str):
    secs_cpp = diag.get("seconds", None)
    secs_cpp_str = f"{secs_cpp:.6f}s" if isinstance(secs_cpp, float) else "None"
    print(
        f"\n[{tag}] Diagnostics:\n"
        f"  minESS  = {diag.get('minESS')}\n"
        f"  maxPSRF = {diag.get('maxPSRF')}\n"
        f"  N       = {diag.get('N')}\n"
        f"  phases  = {diag.get('phases')}\n"
        f"  seconds(C++) = {secs_cpp_str}\n"
        f"  elapsed(py)  = {elapsed_s:.6f}s\n"
    )

def _print_samples_summary(steady_states: np.ndarray, tag: str, n_dims_preview: int = 3):
    ss_min = np.nanmin(steady_states)
    ss_max = np.nanmax(steady_states)
    print(f"[{tag}] Steady states summary:")
    print(f"  shape = {steady_states.shape}, global min = {ss_min:.6g}, global max = {ss_max:.6g}")
    d, N = steady_states.shape
    dims = min(n_dims_preview, d)
    for i in range(dims):
        mean_i = np.nanmean(steady_states[i, :])
        std_i  = np.nanstd(steady_states[i, :])
        print(f"  dim {i}: mean = {mean_i:.6g}, std = {std_i:.6g}")
    print("")

def _run_bsb_once(model: MetabolicNetwork, tc: unittest.TestCase, tag: str):
    """Runs one Billiard Shake-and-Bake phase"""
    sampler = PolytopeSampler(model)
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    A, b, N, N_shift = sampler.get_polytope()
    d_eff = A.shape[1]
    n_rxns = N.shape[0]
    tc.assertGreater(d_eff, 0, "Effective dimension is zero")

    t0 = time.perf_counter()
    steady_states, diag = sampler.generate_steady_states_sb_once(
        n=N_SAMPLES, burn_in=BURN_IN, sampler="billiard_shake_and_bake", nreflections=NREFL
    )
    elapsed = time.perf_counter() - t0

    _print_diag(diag, elapsed, f"{tag} :: BILLIARD_SHAKE_AND_BAKE")
    _print_samples_summary(steady_states, f"{tag} :: BILLIARD_SHAKE_AND_BAKE")



    # sanity
    tc.assertEqual(steady_states.shape[0], n_rxns)
    tc.assertTrue(np.isfinite(diag["minESS"]))
    tc.assertGreater(diag["minESS"], 0)

    del steady_states, sampler, A, b, N, N_shift
    gc.collect()

class TestBilliardShakeAndBake(unittest.TestCase):
    def test_ecoli_core_json(self):
        input_file = os.path.join(os.getcwd(), "ext_data", "e_coli_core.json")
        model = MetabolicNetwork.from_json(input_file)
        _run_bsb_once(model, self, tag="e_coli_core.json")

    def test_ecoli_core_mat(self):
        input_file = os.path.join(os.getcwd(), "ext_data", "e_coli_core.mat")
        model = MetabolicNetwork.from_mat(input_file)
        _run_bsb_once(model, self, tag="e_coli_core.mat")

    def test_ecoli_core_sbml(self):
        input_file = os.path.join(os.getcwd(), "ext_data", "e_coli_core.xml")
        model = MetabolicNetwork.from_sbml(input_file)
        _run_bsb_once(model, self, tag="e_coli_core.xml")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        set_default_solver(sys.argv[1])
        sys.argv.pop(1)
    unittest.main()
