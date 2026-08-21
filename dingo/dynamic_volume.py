"""Wrapper for the supplied static relative-volume routine.

This module does not implement a second volume algorithm.  For every incoming
linear cut it calls MATLAB ``relative_volume`` from the supplied
MetabolicPolytopesCuttingPlanes source, adds that cut to the parent polytope,
and repeats.  No FBA objective or biomass floor is introduced.
"""
from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Iterable, Mapping

import numpy as np
from scipy.io import loadmat, savemat


_STATIC_OPTION_NAMES = {
    "r",
    "delta",
    "alpha",
    "nu",
    "N_utest",
    "window_size",
    "epsilon",
    "max_annealing_iter",
    "max_bisect_iter",
    "simdLen",
    "t_max_safety",
    "verb",
}
_WRAPPER_OPTION_NAMES = {"algorithm", "wrapper_timeout_sec"}
_SOURCE_DEFAULTS = {
    "r": 0.10,
    "delta": 0.05,
    "alpha": 0.20,
    "nu": 10,
    "N_utest": 320,
    "window_size": 5000,
    "epsilon": 0.10,
    "max_annealing_iter": 200,
    "max_bisect_iter": 100,
    "simdLen": 8,
    "t_max_safety": 1.01,
    "verb": 1,
}
_POSITIVE_INTEGER_OPTIONS = {
    "nu",
    "N_utest",
    "window_size",
    "max_annealing_iter",
    "max_bisect_iter",
    "simdLen",
}


@dataclass(frozen=True)
class HalfspaceCut:
    normal: np.ndarray
    threshold: float
    label: str = ""

    def __post_init__(self) -> None:
        normal = np.asarray(self.normal, dtype=float).reshape(-1)
        if normal.size == 0 or not np.all(np.isfinite(normal)):
            raise ValueError("cut normal must be a nonempty finite vector")
        if np.linalg.norm(normal) == 0:
            raise ValueError("cut normal must be nonzero")
        if not np.isfinite(self.threshold):
            raise ValueError("cut threshold must be finite")
        object.__setattr__(self, "normal", normal)


@dataclass(frozen=True)
class StaticVolumeResult:
    rho: float
    log_rho: float
    cut_ratios: np.ndarray
    cut_log_ratios: np.ndarray
    cut_n_phases: np.ndarray
    cut_total_steps: np.ndarray
    cut_converged: np.ndarray
    cut_labels: tuple[str, ...]
    algorithm: str = "static_relative_volume_sequential"
    objective_floor: bool = False


def cuts_from_bounds(
    parent_lb,
    parent_ub,
    child_lb,
    child_ub,
    *,
    atol: float = 1e-12,
) -> list[HalfspaceCut]:
    """Convert a nested bound update into sequential scalar halfspace cuts."""
    parent_lb = _vector(parent_lb, "parent_lb")
    parent_ub = _vector(parent_ub, "parent_ub")
    child_lb = _vector(child_lb, "child_lb")
    child_ub = _vector(child_ub, "child_ub")
    n = parent_lb.size
    if any(x.size != n for x in (parent_ub, child_lb, child_ub)):
        raise ValueError("parent and child bound vectors have different lengths")
    if atol < 0 or not np.isfinite(atol):
        raise ValueError("atol must be nonnegative and finite")
    if np.any(parent_lb > parent_ub) or np.any(child_lb > child_ub):
        raise ValueError("invalid parent or child bounds")
    if np.any(child_lb < parent_lb - atol) or np.any(child_ub > parent_ub + atol):
        raise ValueError("child bounds must define a subset of the parent bounds")
    cuts: list[HalfspaceCut] = []
    for j in range(n):
        if child_ub[j] < parent_ub[j] - atol:
            normal = np.zeros(n)
            normal[j] = 1.0
            cuts.append(HalfspaceCut(normal, float(child_ub[j]), f"ub[{j}]"))
        if child_lb[j] > parent_lb[j] + atol:
            normal = np.zeros(n)
            normal[j] = -1.0
            cuts.append(HalfspaceCut(normal, float(-child_lb[j]), f"lb[{j}]"))
    return cuts


def _vector(value, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=float).reshape(-1)
    if result.size == 0 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a nonempty finite vector")
    return result


def _matlab_quote(path: str | Path) -> str:
    return str(path).replace("'", "''")


def _string_from_mat(value) -> str:
    if isinstance(value, str):
        return value
    array = np.asarray(value)
    if array.dtype.kind in {"U", "S"}:
        return "".join(array.reshape(-1).astype(str).tolist()).strip()
    return str(value).strip()


class StaticVolumeUpdater:
    """Sequential wrapper around the unchanged supplied static estimator."""

    def __init__(
        self,
        backend: str = "batch",
        *,
        matlab_executable: str | None = None,
        matlab_source_dir: str | Path | None = None,
        polytope_sampler_matlab: str | Path | None = None,
    ) -> None:
        if backend != "batch":
            raise ValueError("only backend='batch' is supported")
        self.backend = backend
        self.matlab_executable = matlab_executable or os.environ.get(
            "DINGO_MATLAB_EXECUTABLE", "matlab"
        )
        self.matlab_source_dir = Path(
            matlab_source_dir
            or os.environ.get("DINGO_VOLUME_UPDATING_MATLAB_DIR", "")
            or Path(__file__).resolve().parent / "matlab" / "volume_updating"
        ).resolve()
        sampler = polytope_sampler_matlab or os.environ.get(
            "DINGO_POLYTOPE_SAMPLER_MATLAB"
        )
        self.polytope_sampler_matlab = Path(sampler).resolve() if sampler else None
        self._validate_installation()

    def _validate_installation(self) -> None:
        required = [
            "volume_update_batch.m",
            "relative_volume.m",
            "default_params.m",
            "ul_test.m",
            "new_window.m",
            "update_window.m",
            "init.m",
        ]
        missing = [name for name in required if not (self.matlab_source_dir / name).is_file()]
        if missing:
            raise FileNotFoundError(
                f"incomplete volume-updating MATLAB folder {self.matlab_source_dir}: {missing}"
            )
        executable = shutil.which(self.matlab_executable)
        if executable is None and not Path(self.matlab_executable).is_file():
            raise FileNotFoundError(
                f"MATLAB executable not found: {self.matlab_executable}. Set "
                "DINGO_MATLAB_EXECUTABLE if needed."
            )
        if self.polytope_sampler_matlab is not None:
            missing = [
                name
                for name in ("code", "bin")
                if not (self.polytope_sampler_matlab / name).is_dir()
            ]
            if missing:
                raise FileNotFoundError(
                    "PolytopeSamplerMatlab must contain code/ and bin/: "
                    f"{self.polytope_sampler_matlab} (missing {missing})"
                )

    @property
    def _uses_windows_matlab_from_wsl(self) -> bool:
        return os.name == "posix" and self.matlab_executable.lower().endswith(".exe")

    def _path_for_matlab(self, path: str | Path) -> str:
        path = str(Path(path).resolve())
        if not self._uses_windows_matlab_from_wsl:
            return path
        result = subprocess.run(
            ["wslpath", "-w", path],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()

    @staticmethod
    def _normalize_options(options: Mapping | None) -> tuple[dict, float]:
        options = dict(options or {})
        algorithm = options.pop("algorithm", "static")
        if algorithm not in {"static", "baseline", "relative_volume"}:
            raise ValueError(
                "algorithm must be 'static'/'baseline'. The wrapper does not "
                "provide a separate reuse estimator."
            )
        timeout = float(options.pop("wrapper_timeout_sec", 7200))
        unknown = set(options) - _STATIC_OPTION_NAMES
        if unknown:
            raise ValueError(
                "options are not part of the supplied static API: "
                + ", ".join(sorted(unknown))
            )
        for name, value in options.items():
            if isinstance(value, (bool, np.bool_)):
                raise ValueError(f"{name} must be numeric, not boolean")
            array = np.asarray(value)
            if array.size != 1:
                raise ValueError(f"{name} must be a numeric scalar")
            try:
                numeric = float(array.reshape(-1)[0])
            except (TypeError, ValueError) as error:
                raise ValueError(f"{name} must be a numeric scalar") from error
            if not np.isfinite(numeric):
                raise ValueError(f"{name} must be finite")
            if name in _POSITIVE_INTEGER_OPTIONS:
                if numeric <= 0 or not numeric.is_integer():
                    raise ValueError(f"{name} must be a positive integer")
                options[name] = int(numeric)
            else:
                options[name] = numeric
        if "r" in options and not 0 < options["r"] < 1:
            raise ValueError("r must be strictly between 0 and 1")
        if "delta" in options and options["delta"] <= 0:
            raise ValueError("delta must be positive")
        effective_r = options.get("r", _SOURCE_DEFAULTS["r"])
        effective_delta = options.get("delta", _SOURCE_DEFAULTS["delta"])
        if effective_r + effective_delta >= 1:
            raise ValueError("r + delta must be less than 1")
        if "alpha" in options and not 0 < options["alpha"] < 1:
            raise ValueError("alpha must be strictly between 0 and 1")
        if "epsilon" in options and options["epsilon"] <= 0:
            raise ValueError("epsilon must be positive")
        if "t_max_safety" in options and options["t_max_safety"] < 1:
            raise ValueError("t_max_safety must be at least 1")
        if "verb" in options and options["verb"] not in {0, 1, 2}:
            raise ValueError("verb must be 0, 1, or 2")
        if timeout <= 0 or not np.isfinite(timeout):
            raise ValueError("wrapper_timeout_sec must be positive and finite")
        return options, timeout

    @staticmethod
    def _normalize_cuts(cuts: Iterable[HalfspaceCut], dimension: int) -> list[HalfspaceCut]:
        result = []
        for cut in cuts:
            if not isinstance(cut, HalfspaceCut):
                raise TypeError("every cut must be a HalfspaceCut")
            if cut.normal.size != dimension:
                raise ValueError(
                    f"cut {cut.label!r} has dimension {cut.normal.size}, expected {dimension}"
                )
            result.append(cut)
        return result

    def _run_matlab(
        self,
        request_path: Path,
        output_path: Path,
        timeout: float,
    ) -> None:
        source = self._path_for_matlab(self.matlab_source_dir)
        request = self._path_for_matlab(request_path)
        output = self._path_for_matlab(output_path)
        sampler_setup = ""
        if self.polytope_sampler_matlab is not None:
            sampler = self._path_for_matlab(self.polytope_sampler_matlab)
            sampler_setup = (
                "setenv('DINGO_POLYTOPE_SAMPLER_MATLAB', "
                f"'{_matlab_quote(sampler)}'); "
            )
        command = (
            sampler_setup
            + f"addpath(genpath('{_matlab_quote(source)}')); "
            f"volume_update_batch('{_matlab_quote(request)}', "
            f"'{_matlab_quote(output)}');"
        )
        env = os.environ.copy()
        completed = subprocess.run(
            [self.matlab_executable, "-batch", command],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
        )
        if completed.returncode != 0:
            transcript = (completed.stdout + "\n" + completed.stderr).strip()
            raise RuntimeError(
                f"MATLAB volume update failed with code {completed.returncode}:\n"
                + transcript[-12000:]
            )
        if not output_path.is_file():
            raise RuntimeError("MATLAB exited successfully but produced no output MAT file")

    def estimate(
        self,
        S,
        lb,
        ub,
        cuts: Iterable[HalfspaceCut],
        *,
        options: Mapping | None = None,
    ) -> StaticVolumeResult:
        S = np.asarray(S, dtype=float)
        lb = _vector(lb, "lb")
        ub = _vector(ub, "ub")
        if S.ndim != 2 or S.shape[1] != lb.size or ub.size != lb.size:
            raise ValueError("S, lb and ub dimensions are inconsistent")
        if not np.all(np.isfinite(S)):
            raise ValueError("S must be finite")
        if np.any(lb > ub):
            raise ValueError("a lower bound exceeds its upper bound")
        cuts = self._normalize_cuts(cuts, lb.size)
        static_options, timeout = self._normalize_options(options)
        if not cuts:
            empty_float = np.empty(0, dtype=float)
            return StaticVolumeResult(
                1.0,
                0.0,
                empty_float,
                empty_float,
                np.empty(0, dtype=int),
                np.empty(0, dtype=int),
                np.empty(0, dtype=bool),
                (),
            )

        normals = np.vstack([cut.normal for cut in cuts])
        thresholds = np.asarray([cut.threshold for cut in cuts], dtype=float)
        with tempfile.TemporaryDirectory(prefix="volume_update_") as temp:
            temp_dir = Path(temp)
            request_path = temp_dir / "request.mat"
            output_path = temp_dir / "result.mat"
            savemat(
                request_path,
                {
                    "S": S,
                    "lb": lb.reshape(-1, 1),
                    "ub": ub.reshape(-1, 1),
                    "cut_normals": normals,
                    "cut_thresholds": thresholds.reshape(-1, 1),
                    "options": static_options,
                },
                do_compression=True,
            )
            self._run_matlab(request_path, output_path, timeout)
            result = loadmat(output_path, squeeze_me=True, struct_as_record=False)

        cut_ratios = np.atleast_1d(result["cut_ratios"]).astype(float)
        cut_logs = np.atleast_1d(result["cut_log_ratios"]).astype(float)
        cut_phases = np.atleast_1d(result["cut_n_phases"]).astype(int)
        cut_steps = np.atleast_1d(result["cut_total_steps"]).astype(int)
        cut_converged = np.atleast_1d(result["cut_converged"]).astype(bool)
        expected = len(cuts)
        arrays = [cut_ratios, cut_logs, cut_phases, cut_steps, cut_converged]
        if any(array.size != expected for array in arrays):
            raise RuntimeError("MATLAB result arrays do not match the number of cuts")
        if (
            not np.all(np.isfinite(cut_ratios))
            or np.any(cut_ratios < 0)
            or np.any(cut_ratios > 1 + 1e-6)
        ):
            raise RuntimeError("MATLAB returned an invalid per-cut ratio")
        expected_logs = np.log(np.maximum(cut_ratios, 1e-300))
        if not np.allclose(cut_logs, expected_logs, atol=1e-12, rtol=1e-10):
            raise RuntimeError("MATLAB per-cut ratios and log ratios are inconsistent")
        if np.any(cut_phases < 0) or np.any(cut_steps < 0):
            raise RuntimeError("MATLAB returned a negative phase or step count")
        rho = float(np.asarray(result["rho"]).squeeze())
        log_rho = float(np.asarray(result["log_rho"]).squeeze())
        if not np.isfinite(rho) or rho < 0 or rho > 1 + 1e-6:
            raise RuntimeError(f"invalid cumulative ratio returned by MATLAB: {rho}")
        if not np.isclose(log_rho, np.sum(cut_logs), atol=1e-10, rtol=1e-10):
            raise RuntimeError("MATLAB cumulative log ratio is inconsistent with cut ratios")
        algorithm = _string_from_mat(result.get("algorithm", "static_relative_volume_sequential"))
        if algorithm != "static_relative_volume_sequential":
            raise RuntimeError(f"unexpected MATLAB algorithm identifier: {algorithm!r}")
        objective_floor = bool(np.asarray(result.get("objective_floor", False)).squeeze())
        if objective_floor:
            raise RuntimeError("MATLAB result unexpectedly reports an objective floor")
        return StaticVolumeResult(
            rho=rho,
            log_rho=log_rho,
            cut_ratios=cut_ratios,
            cut_log_ratios=cut_logs,
            cut_n_phases=cut_phases,
            cut_total_steps=cut_steps,
            cut_converged=cut_converged,
            cut_labels=tuple(cut.label for cut in cuts),
            algorithm=algorithm,
            objective_floor=False,
        )


# Backward-compatible class name used by the earlier wrapper code.  The alias
# now points to the supplied static implementation; it does not imply a second
# dynamic or reuse algorithm.
MatlabVolumeUpdater = StaticVolumeUpdater


__all__ = [
    "StaticVolumeResult",
    "StaticVolumeUpdater",
    "HalfspaceCut",
    "MatlabVolumeUpdater",
    "cuts_from_bounds",
]
