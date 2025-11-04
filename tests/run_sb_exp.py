#!/usr/bin/python3
import os, sys
import numpy as np
from dingo import PolytopeSampler
from time import process_time
import pickle


def sample_on_polyround_processed_polytope(p):
    name = os.path.basename(p)

    with open(p, "rb") as f:
        obj = pickle.load(f)
    polytope = obj[0]

    polyround_A = polytope.A.to_numpy()
    polyround_b = polytope.b.to_numpy()

    print(f"Dimensions = {polyround_A.shape[1]}, Ograničenja = {polyround_A.shape[0]}")
    start = process_time()

    steady_states, diag = PolytopeSampler.sample_from_polytope_sb_once(
        polyround_A,
        polyround_b,
        n=15000,
        burn_in=1000,
        sampler="sb",
        walk_len=48,
    )

    end = process_time()

    print(
        f"[{name}] minESS={diag['minESS']:.3f}  maxPSRF={diag['maxPSRF']:.3f}  "
        f"N={diag['N']}  phases={diag['phases']}  seconds={diag['seconds']}"
    )
    print(f"Total sampling time: {end - start:.2f} s")

    out_path = f"dingo_polyround_no_multiphase_{name}.pckl"
    with open(out_path, "wb") as f_out:
        pickle.dump({"samples": steady_states, "diagnostics": diag}, f_out)

    print(f"Saved results to {out_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: python run_bwr_exp.py <path_to_polytope_pickle>")
    file_name = sys.argv[1]
    sample_on_polyround_processed_polytope(file_name)
