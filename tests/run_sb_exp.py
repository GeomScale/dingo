# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2025 Vissarion Fisikopoulos
# Copyright (c) 2018-2025 Apostolos Chalkis
# Copyright (c) 2025-2025 Iva Janković

# Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

import os, sys, pickle, argparse
import numpy as np
from time import perf_counter
from volestipy import HPolytope

def _binary_search_k(P, S_all, L0, S_last, target_ess):
    """Find smallest k_rel in [1, S_last.shape[1]] where minESS >= target_ess."""
    lo, hi = 1, S_last.shape[1]
    best_k_rel, best_diag = hi, None

    def pref_diag(k_rel):
        return P.boundary_diag(S_all[:, : L0 + k_rel])

    while lo <= hi:
        mid = (lo + hi) // 2
        dmid = pref_diag(mid)
        if float(dmid.get("minESS", -1.0)) >= target_ess:
            best_k_rel, best_diag = mid, dmid
            hi = mid - 1
        else:
            lo = mid + 1

    if best_diag is None:
        best_diag = P.boundary_diag(S_all[:, : L0 + best_k_rel])
    return best_k_rel, best_diag

def sample_on_polyround_processed_polytope(p, target_ess, n_chunk, max_total,
                                           sampler, walk_len, burn_in_init):
    name = os.path.basename(p)

    with open(p, "rb") as f:
        obj = pickle.load(f)
    polytope = obj[0]

    A = np.ascontiguousarray(polytope.A.to_numpy(), dtype=np.float64)
    b = np.ascontiguousarray(polytope.b.to_numpy(), dtype=np.float64)
    assert A.flags['C_CONTIGUOUS'] and b.flags['C_CONTIGUOUS']

    d, m = A.shape[1], A.shape[0]
    m, n = A.shape
    d_meta = getattr(polytope, "dimension", None) or getattr(polytope, "d", None)
    print(f"m={m}, n={n}, d(from meta)={d_meta}")

    # nreflections = ceil(0.25 * d) for BSB, else 0
    if sampler.lower() == "bsb":
        nreflections = int(np.ceil(0.25 * d))
    else:
        nreflections = 0

    print(f"Dimensions = {d}, Constraints = {m}")

    P = HPolytope(A, b)

    S_parts, dur = [], []
    total = 0
    burn  = int(burn_in_init)

    while True:
        t0 = perf_counter()
        S_chunk = P.boundary_sample(
            sampler="sb",
            number_of_points=int(n_chunk),
            number_of_points_to_burn=int(burn),
            walk_len=int(walk_len),
            nreflections=int(nreflections),
        )
        t1 = perf_counter()

        dur.append(t1 - t0)
        S_parts.append(S_chunk)
        total += S_chunk.shape[1]
        burn = 0  # continue without burn-in

        S_all = np.concatenate(S_parts, axis=1)
        diag_all = P.boundary_diag(S_all)
        minESS_all = float(diag_all.get("minESS", np.nan))

        if minESS_all >= target_ess:
            # find minimal prefix in last chunk
            L0 = total - S_chunk.shape[1]
            k_rel, diag_at_k = _binary_search_k(P, S_all, L0, S_chunk, target_ess)
            k_abs = L0 + k_rel

            # sampling time = sum of past chunks + part of last chunk
            frac = k_rel / S_chunk.shape[1]
            sampling_seconds = sum(dur[:-1]) + frac * dur[-1]
            chunks_used = len(dur)

            print(
                f"[{name}] minESS={diag_at_k['minESS']:.3f}  "
                f"maxPSRF={diag_at_k.get('maxPSRF', np.nan):.3f}  "
                f"N={int(diag_at_k.get('N', k_abs))}  "
                f"chunks={chunks_used}  sampling_seconds={sampling_seconds:.3f}"
            )
            print(f"Total sampling time: {sampling_seconds:.2f} s")

            out_path = f"dingo_polyround_no_multiphase_{name}"
            with open(out_path, "wb") as f_out:
                pickle.dump(
                    {
                        "samples": S_all[:, :k_abs],
                        "diagnostics": {
                            **diag_at_k,
                            "N_at_threshold": int(k_abs),
                            "chunks": chunks_used,
                            "sampling_seconds": sampling_seconds,
                        },
                        "params": {
                            "sampler": sampler,
                            "walk_len": int(walk_len),
                            "nreflections": int(nreflections),
                            "target_ess": target_ess,
                            "n_chunk": int(n_chunk),
                            "max_total": int(max_total),
                            "burn_in_init": int(burn_in_init),
                        },
                    },
                    f_out,
                )
            print(f"Saved results to {out_path}")
            return

        if total >= max_total:
            sampling_seconds = sum(dur)
            chunks_used = len(dur)

            print(
                f"[{name}] minESS={diag_all.get('minESS', np.nan):.3f}  "
                f"maxPSRF={diag_all.get('maxPSRF', np.nan):.3f}  "
                f"N={int(diag_all.get('N', S_all.shape[1]))}  "
                f"chunks={chunks_used}  sampling_seconds={sampling_seconds:.3f}"
            )
            print(f"Total sampling time: {sampling_seconds:.2f} s  (target {target_ess} not reached)")

            out_path = f"dingo_polyround_no_multiphase_{name}"
            with open(out_path, "wb") as f_out:
                pickle.dump(
                    {
                        "samples": S_all,
                        "diagnostics": {
                            **diag_all,
                            "N_total_generated": int(S_all.shape[1]),
                            "chunks": chunks_used,
                            "sampling_seconds": sampling_seconds,
                        },
                        "params": {
                            "sampler": sampler,
                            "walk_len": int(walk_len),
                            "nreflections": int(nreflections),
                            "target_ess": target_ess,
                            "n_chunk": int(n_chunk),
                            "max_total": int(max_total),
                            "burn_in_init": int(burn_in_init),
                        },
                    },
                    f_out,
                )
            print(f"Saved results to {out_path}")
            return

def _make_parser():
    ap = argparse.ArgumentParser()
    ap.add_argument("polytope_pickle", type=str)
    ap.add_argument("--target-ess", type=float, default=1000.0)
    ap.add_argument("--n-chunk", type=int, default=5000)
    ap.add_argument("--max-total", type=int, default=200000)
    ap.add_argument("--sampler", type=str, default="bsb", choices=["sb", "bsb"])
    ap.add_argument("--walk-len", type=int, default=1)
    ap.add_argument("--burn-in-init", type=int, default=0)
    return ap

if __name__ == "__main__":
    parser = _make_parser()
    args = parser.parse_args()
    sample_on_polyround_processed_polytope(
        p=args.polytope_pickle,
        target_ess=args.target_ess,
        n_chunk=args.n_chunk,
        max_total=args.max_total,
        sampler=args.sampler,
        walk_len=args.walk_len,
        burn_in_init=args.burn_in_init,
    )
