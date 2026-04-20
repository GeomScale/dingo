from __future__ import annotations

import os
import csv
import sys
import json
import time
import pickle
import signal
import argparse
import subprocess
from typing import List, Tuple

import numpy as np

from PolytopeSampler import PolytopeSampler


CSV_HEADER = [
    "status",
    "error",
    "sampler",
    "model_id",
    "src_path",
    "samples_path",
    "m_constraints",
    "vars",
    "ess_target",
    "chunk_n",
    "burn_in_first",
    "walk_len",
    "nreflections",
    "max_calls",
    "calls",
    "k_star",
    "total_N",
    "N_diag",
    "minESS",
    "maxPSRF",
    "elapsed_sec",
    "loading_sec",
    "sampling_sec",
]


def append_row(csv_path: str, row: dict) -> None:
    os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)
    write_header = (not os.path.exists(csv_path)) or (os.path.getsize(csv_path) == 0)

    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_HEADER)
        if write_header:
            writer.writeheader()

        for k in CSV_HEADER:
            row.setdefault(k, "")
        writer.writerow(row)


def compute_nreflections(d: int, mode: str, fixed: int) -> int:
    if mode == "0":
        return 0
    if mode == "fixed":
        return int(fixed)
    return max(0, int(round(0.25 * d)))


def find_input_files(polytopes_dir: str, exts: Tuple[str, ...]) -> List[Tuple[str, str]]:
    items: List[Tuple[str, str]] = []
    for root, _, files in os.walk(polytopes_dir):
        for fn in sorted(files):
            if fn.lower().endswith(exts):
                path = os.path.join(root, fn)
                model_id = os.path.splitext(fn)[0]
                items.append((model_id, path))
    return items


def extract_Ab_from_object(obj):
    if isinstance(obj, dict):
        if "A" in obj and "b" in obj:
            A = np.asarray(obj["A"], dtype=np.float64)
            b = np.asarray(obj["b"], dtype=np.float64)
            name = str(obj.get("name", ""))
            return A, b, name
        if "P" in obj:
            P = obj["P"]
            A = np.asarray(P.A, dtype=np.float64)
            b = np.asarray(P.b, dtype=np.float64)
            name = str(obj.get("name", ""))
            return A, b, name

    if isinstance(obj, (tuple, list)):
        if len(obj) >= 2:
            if isinstance(obj[0], (np.ndarray, list)) and isinstance(obj[1], (np.ndarray, list)):
                A = np.asarray(obj[0], dtype=np.float64)
                b = np.asarray(obj[1], dtype=np.float64)
                name = str(obj[2]) if len(obj) >= 3 else ""
                return A, b, name

            P = obj[0]
            if hasattr(P, "A") and hasattr(P, "b"):
                A = np.asarray(P.A, dtype=np.float64)
                b = np.asarray(P.b, dtype=np.float64)
                name = str(obj[1]) if len(obj) >= 2 else ""
                return A, b, name

    if hasattr(obj, "A") and hasattr(obj, "b"):
        A = np.asarray(obj.A, dtype=np.float64)
        b = np.asarray(obj.b, dtype=np.float64)
        return A, b, ""

    raise ValueError("Unknown preprocessed file format: could not extract A and b.")


def load_preprocessed_Ab(src_path: str):
    low = src_path.lower()

    if low.endswith(".npz"):
        z = np.load(src_path, allow_pickle=True)
        if "A" not in z or "b" not in z:
            raise ValueError("NPZ file must contain keys 'A' and 'b'.")
        A = np.asarray(z["A"], dtype=np.float64)
        b = np.asarray(z["b"], dtype=np.float64)
        name = str(z["name"]) if "name" in z else ""
        return A, b, name

    with open(src_path, "rb") as f:
        obj = pickle.load(f)

    return extract_Ab_from_object(obj)


def save_samples(out_dir: str, model_id: str, sampler: str, S_star: np.ndarray) -> str:
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"S_star_{model_id}_{sampler}.npy")
    np.save(out_path, S_star)
    return out_path


def validate_Ab(A: np.ndarray, b: np.ndarray) -> None:
    if A.ndim != 2:
        raise ValueError(f"A must be 2D, got shape={A.shape}")
    if b.ndim != 1:
        raise ValueError(f"b must be 1D, got shape={b.shape}")
    if A.shape[0] != b.shape[0]:
        raise ValueError(f"Incompatible dimensions: A.shape={A.shape}, b.shape={b.shape}")
    if not np.isfinite(A).all():
        raise ValueError("A contains NaN or Inf.")
    if not np.isfinite(b).all():
        raise ValueError("b contains NaN or Inf.")


def sanitize_samples(S: np.ndarray, d: int) -> np.ndarray:
    S = np.asarray(S, dtype=np.float64)

    if S.ndim != 2 or S.size == 0:
        raise RuntimeError(f"Invalid sample matrix shape: {S.shape}")

    if S.shape[0] != d and S.shape[1] == d:
        S = S.T

    if S.shape[0] != d:
        raise RuntimeError(f"Expected samples of shape (d, N) with d={d}, got {S.shape}")

    if not np.isfinite(S).all():
        raise RuntimeError("Sample matrix contains NaN or Inf.")

    return np.asfortranarray(S)


def worker_run_one(args) -> dict:
    sampler_label = args.worker_sampler
    model_id = args.worker_model_id
    src_path = args.worker_src_path

    row = {
        "status": "fail",
        "error": "",
        "sampler": sampler_label,
        "model_id": model_id,
        "src_path": src_path,
        "samples_path": "",
        "ess_target": int(args.ess_target),
        "chunk_n": int(args.chunk_n),
        "burn_in_first": int(args.burn_in_first),
        "walk_len": int(args.walk_len),
        "max_calls": int(args.max_calls),
    }

    t0 = time.perf_counter()

    t_load0 = time.perf_counter()
    A, b, _ = load_preprocessed_Ab(src_path)
    A = np.ascontiguousarray(A, dtype=np.float64)
    b = np.ascontiguousarray(b, dtype=np.float64)
    validate_Ab(A, b)
    loading_sec = time.perf_counter() - t_load0

    d = int(A.shape[1])
    m = int(A.shape[0])
    nref = compute_nreflections(d, args.nref_mode, args.nref_fixed)

    row.update({
        "vars": d,
        "m_constraints": m,
        "nreflections": nref,
        "loading_sec": loading_sec,
    })

    t_samp0 = time.perf_counter()
    S_star, info = PolytopeSampler.boundary_sample_ess(
        A=A,
        b=b,
        ess_target=int(args.ess_target),
        chunk_n=int(args.chunk_n),
        burn_in_first=int(args.burn_in_first),
        sampler=sampler_label,
        walk_len=int(args.walk_len),
        nreflections=int(nref),
        max_calls=int(args.max_calls),
    )
    sampling_sec = time.perf_counter() - t_samp0
    elapsed = time.perf_counter() - t0

    S_star = sanitize_samples(S_star, d)

    if args.save_samples:
        samples_dir = os.path.join(args.out_dir, f"samples_{sampler_label}")
        samples_path = save_samples(samples_dir, model_id, sampler_label, S_star)
    else:
        samples_path = ""

    row.update({
        "status": "ok",
        "error": "",
        "samples_path": samples_path,
        "calls": int(info.get("calls", -1)),
        "k_star": int(info.get("k_star", -1)),
        "total_N": int(info.get("total_N", S_star.shape[1])),
        "N_diag": int(info.get("N", S_star.shape[1])),
        "minESS": float(info.get("minESS", float("nan"))),
        "maxPSRF": float(info.get("maxPSRF", float("nan"))),
        "elapsed_sec": elapsed,
        "sampling_sec": sampling_sec,
    })

    return row


def run_worker_mode(args) -> int:
    try:
        row = worker_run_one(args)
        print(json.dumps(row), flush=True)
        return 0
    except Exception as e:
        row = {
            "status": "fail",
            "error": repr(e),
            "sampler": args.worker_sampler,
            "model_id": args.worker_model_id,
            "src_path": args.worker_src_path,
            "samples_path": "",
            "m_constraints": -1,
            "vars": -1,
            "ess_target": int(args.ess_target),
            "chunk_n": int(args.chunk_n),
            "burn_in_first": int(args.burn_in_first),
            "walk_len": int(args.walk_len),
            "nreflections": -1,
            "max_calls": int(args.max_calls),
            "calls": -1,
            "k_star": -1,
            "total_N": -1,
            "N_diag": -1,
            "minESS": float("nan"),
            "maxPSRF": float("nan"),
            "elapsed_sec": float("nan"),
            "loading_sec": float("nan"),
            "sampling_sec": float("nan"),
        }
        print(json.dumps(row), flush=True)
        return 1


def run_job_in_subprocess(args, sampler_label: str, model_id: str, src_path: str) -> dict:
    cmd = [
        sys.executable,
        os.path.abspath(__file__),
        "--worker",
        "--worker_sampler", sampler_label,
        "--worker_model_id", model_id,
        "--worker_src_path", src_path,
        "--out_dir", args.out_dir,
        "--ess_target", str(args.ess_target),
        "--chunk_n", str(args.chunk_n),
        "--max_calls", str(args.max_calls),
        "--walk_len", str(args.walk_len),
        "--burn_in_first", str(args.burn_in_first),
        "--nref_mode", args.nref_mode,
        "--nref_fixed", str(args.nref_fixed),
    ]

    if args.save_samples:
        cmd.append("--save_samples")

    t0 = time.perf_counter()
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    elapsed = time.perf_counter() - t0

    stdout = (proc.stdout or "").strip()
    stderr = (proc.stderr or "").strip()

    if stdout:
        last_line = stdout.splitlines()[-1]
        try:
            row = json.loads(last_line)
            if "elapsed_sec" not in row or row["elapsed_sec"] in ("", None) or (
                isinstance(row["elapsed_sec"], float) and np.isnan(row["elapsed_sec"])
            ):
                row["elapsed_sec"] = elapsed
            return row
        except json.JSONDecodeError:
            pass

    if proc.returncode < 0:
        sig = -proc.returncode
        if sig == signal.SIGABRT:
            err = f"subprocess aborted with SIGABRT ({sig})"
        elif sig == signal.SIGSEGV:
            err = f"subprocess crashed with SIGSEGV ({sig})"
        else:
            err = f"subprocess terminated by signal {sig}"
    else:
        err = f"subprocess failed with returncode={proc.returncode}"

    if stderr:
        err += f" | stderr: {stderr[-1000:]}"
    elif stdout:
        err += f" | stdout: {stdout[-1000:]}"

    return {
        "status": "fail",
        "error": err,
        "sampler": sampler_label,
        "model_id": model_id,
        "src_path": src_path,
        "samples_path": "",
        "m_constraints": -1,
        "vars": -1,
        "ess_target": int(args.ess_target),
        "chunk_n": int(args.chunk_n),
        "burn_in_first": int(args.burn_in_first),
        "walk_len": int(args.walk_len),
        "nreflections": -1,
        "max_calls": int(args.max_calls),
        "calls": -1,
        "k_star": -1,
        "total_N": -1,
        "N_diag": -1,
        "minESS": float("nan"),
        "maxPSRF": float("nan"),
        "elapsed_sec": elapsed,
        "loading_sec": float("nan"),
        "sampling_sec": float("nan"),
    }


def build_arg_parser():
    ap = argparse.ArgumentParser()

    ap.add_argument("--polytopes_dir", help="Directory with preprocessed polytopes (.pckl or .npz)")
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--sampler", default="both", choices=["sb", "bsb", "both"])
    ap.add_argument("--ess_target", type=int, default=1000)
    ap.add_argument("--chunk_n", type=int, default=5000)
    ap.add_argument("--max_calls", type=int, default=1000)

    ap.add_argument("--walk_len", type=int, default=1)
    ap.add_argument("--burn_in_first", type=int, default=0)
    ap.add_argument("--nref_mode", default="0.25d", choices=["0.25d", "fixed", "0"])
    ap.add_argument("--nref_fixed", type=int, default=0)

    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--skip", type=int, default=0)
    ap.add_argument("--ext", default="pckl,npz")
    ap.add_argument("--save_samples", action="store_true")

    # internal worker mode
    ap.add_argument("--worker", action="store_true")
    ap.add_argument("--worker_sampler")
    ap.add_argument("--worker_model_id")
    ap.add_argument("--worker_src_path")

    return ap


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    if args.worker:
        raise SystemExit(run_worker_mode(args))

    if not args.polytopes_dir:
        raise RuntimeError("--polytopes_dir is required in normal mode.")

    os.makedirs(args.out_dir, exist_ok=True)

    exts = tuple("." + e.strip().lower().lstrip(".") for e in args.ext.split(",") if e.strip())
    items = find_input_files(args.polytopes_dir, exts)

    if args.skip > 0:
        items = items[int(args.skip):]
    if args.limit > 0:
        items = items[:int(args.limit)]

    if not items:
        raise RuntimeError(f"No files with extensions {exts} found in --polytopes_dir")

    samplers = ["sb", "bsb"] if args.sampler == "both" else [args.sampler]

    print(f"Found {len(items)} polytope files")
    print(f"Samplers: {samplers}")
    print(
        f"ESS={args.ess_target}, chunk_n={args.chunk_n}, walk_len={args.walk_len}, "
        f"burn_in_first={args.burn_in_first}, max_calls={args.max_calls}"
    )

    for sampler_label in samplers:
        csv_path = os.path.join(
            args.out_dir,
            f"results_preprocessed_boundary_ess_{sampler_label}_ESS{args.ess_target}_chunk{args.chunk_n}.csv"
        )

        print(f"\n=== sampler={sampler_label} ===")
        print(f"CSV -> {csv_path}")

        for idx, (model_id, src_path) in enumerate(items, start=1):
            print(f"[{idx}/{len(items)}] {model_id}")

            row = run_job_in_subprocess(
                args=args,
                sampler_label=sampler_label,
                model_id=model_id,
                src_path=src_path,
            )

            append_row(csv_path, row)

            print(
                f"  status={row['status']} | "
                f"minESS={row.get('minESS', '')} | "
                f"maxPSRF={row.get('maxPSRF', '')} | "
                f"error={row.get('error', '')[:120]}"
            )

        print(f"sampler={sampler_label}")

if __name__ == "__main__":
    main()