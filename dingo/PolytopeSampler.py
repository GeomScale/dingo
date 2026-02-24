# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2024 Ke Shi

# Licensed under GNU LGPL.3, see LICENCE file

# Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

import numpy as np
import warnings
from typing import Optional, Tuple, Dict
import math
import time
from dingo.MetabolicNetwork import MetabolicNetwork
from dingo.utils import (
    map_samples_to_steady_states,
    get_matrices_of_low_dim_polytope,
    get_matrices_of_full_dim_polytope,
)

from dingo.pyoptinterface_based_impl import fba,fva,inner_ball,remove_redundant_facets

from volestipy import HPolytope


class PolytopeSampler:
    def __init__(self, metabol_net):

        if not isinstance(metabol_net, MetabolicNetwork):
            raise Exception("An unknown input object given for initialization.")

        self._metabolic_network = metabol_net
        self._A = None
        self._b = None
        self._N = None
        self._N_shift = None
        self._T = None
        self._T_shift = None
        self._parameters = {}
        self._parameters["nullspace_method"] = "sparseQR"
        self._parameters["opt_percentage"] = self.metabolic_network.parameters[
            "opt_percentage"
        ]
        self._parameters["distribution"] = "uniform"
        self._parameters["first_run_of_mmcs"] = True
        self._parameters["remove_redundant_facets"] = True

        self._parameters["tol"] = 1e-06
        self._parameters["solver"] = None
        self._last_hpoly = None

    def get_polytope(self):
        """A member function to derive the corresponding full dimensional polytope
        and a isometric linear transformation that maps the latter to the initial space.
        """

        if (
            self._A is None
            or self._b is None
            or self._N is None
            or self._N_shift is None
            or self._T is None
            or self._T_shift is None
        ):


            (
                max_flux_vector,
                max_objective,
            ) = self._metabolic_network.fba()

            if (
                self._parameters["remove_redundant_facets"]
            ):

                A, b, Aeq, beq = remove_redundant_facets(
                    self._metabolic_network.lb,
                    self._metabolic_network.ub,
                    self._metabolic_network.S,
                    self._metabolic_network.objective_function,
                    self._parameters["opt_percentage"],
                    self._parameters["solver"],
                )
            else:

                (
                    min_fluxes,
                    max_fluxes,
                    max_flux_vector,
                    max_objective,
                ) = self._metabolic_network.fva()

                A, b, Aeq, beq = get_matrices_of_low_dim_polytope(
                    self._metabolic_network.S,
                    self._metabolic_network.lb,
                    self._metabolic_network.ub,
                    min_fluxes,
                    max_fluxes,
                )

            if (
                A.shape[0] != b.size
                or A.shape[1] != Aeq.shape[1]
                or Aeq.shape[0] != beq.size
            ):
                raise Exception("Preprocess for full dimensional polytope failed.")

            A = np.vstack((A, -self._metabolic_network.objective_function))

            b = np.append(
                b,
                -np.floor(max_objective / self._parameters["tol"])
                * self._parameters["tol"]
                * self._parameters["opt_percentage"]
                / 100,
            )

            (
                self._A,
                self._b,
                self._N,
                self._N_shift,
            ) = get_matrices_of_full_dim_polytope(A, b, Aeq, beq)

            n = self._A.shape[1]
            self._T = np.eye(n)
            self._T_shift = np.zeros(n)

        return self._A, self._b, self._N, self._N_shift

    def generate_steady_states(
        self, ess=1000, psrf=False, parallel_mmcs=False, num_threads=1
    ):
        """A member function to sample steady states.

        Keyword arguments:
        ess -- the target effective sample size
        psrf -- a boolean flag to request PSRF smaller than 1.1 for all marginal fluxes
        parallel_mmcs -- a boolean flag to request the parallel mmcs
        num_threads -- the number of threads to use for parallel mmcs
        """

        self.get_polytope()

        P = HPolytope(self._A, self._b)

        self._A, self._b, Tr, Tr_shift, samples = P.mmcs(
            ess, psrf, parallel_mmcs, num_threads, self._parameters["solver"]
        )

        if self._parameters["first_run_of_mmcs"]:
            steady_states = map_samples_to_steady_states(
                samples, self._N, self._N_shift
            )
            self._parameters["first_run_of_mmcs"] = False
        else:
            steady_states = map_samples_to_steady_states(
                samples, self._N, self._N_shift, self._T, self._T_shift
            )

        self._T = np.dot(self._T, Tr)
        self._T_shift = np.add(self._T_shift, Tr_shift)

        return steady_states

    def generate_steady_states_no_multiphase(
        self, method = 'billiard_walk', n=1000, burn_in=0, thinning=1, variance=1.0, bias_vector=None, ess=1000
    ):
        """A member function to sample steady states.

        Keyword arguments:
        method -- An MCMC method to sample, i.e. {'billiard_walk', 'cdhr', 'rdhr', 'ball_walk', 'dikin_walk', 'john_walk', 'vaidya_walk', 'gaussian_hmc_walk', 'exponential_hmc_walk', 'hmc_leapfrog_gaussian', 'hmc_leapfrog_exponential}
        n -- the number of steady states to sample
        burn_in -- the number of points to burn before sampling
        thinning -- the walk length of the chain
        """

        self.get_polytope()

        P = HPolytope(self._A, self._b)

        if bias_vector is None:
            bias_vector = np.ones(self._A.shape[1], dtype=np.float64)
        else:
            bias_vector = bias_vector.astype('float64')

        samples = P.generate_samples(method.encode('utf-8'), n, burn_in, thinning, variance, bias_vector, self._parameters["solver"], ess)
        samples_T = samples.T

        steady_states = map_samples_to_steady_states(
                samples_T, self._N, self._N_shift
            )

        return steady_states

    def generate_steady_states_sb_once(
        self,n: int = 1000,burn_in: int = 0,sampler: str = "sb",nreflections: Optional[int] = None,walk_len: Optional[int] = None,
    ):
        """
        Single-call boundary sampler mapped to steady states.
        """
        self.get_polytope()
        if self._last_hpoly is None:
            self._last_hpoly = HPolytope(self._A, self._b)
        P = self._last_hpoly
        d = int(self._A.shape[1])

        S = P.boundary_sample(
            sampler=sampler,
            number_of_points=int(n),
            number_of_points_to_burn=int(burn_in),
            walk_len=int(walk_len),
            nreflections=int(nreflections),
        )

        if not np.isfinite(S).all():
            raise RuntimeError("boundary_sample returned NaN/Inf; check constraints or parameters.")

        return map_samples_to_steady_states(S, self._N, self._N_shift)

    @staticmethod
    def sample_from_polytope(
        A, b, ess=1000, psrf=False, parallel_mmcs=False, num_threads=1, solver=None
    ):
        """A static function to sample from a full dimensional polytope.

        Keyword arguments:
        A -- an mxn matrix that contains the normal vectors of the facets of the polytope row-wise
        b -- a m-dimensional vector, s.t. A*x <= b
        ess -- the target effective sample size
        psrf -- a boolean flag to request PSRF smaller than 1.1 for all marginal fluxes
        parallel_mmcs -- a boolean flag to request the parallel mmcs
        num_threads -- the number of threads to use for parallel mmcs
        """

        P = HPolytope(A, b)

        A, b, Tr, Tr_shift, samples = P.mmcs(
            ess, psrf, parallel_mmcs, num_threads, solver
        )


        return samples

    @staticmethod
    def sample_from_polytope_no_multiphase(
        A, b, method = 'billiard_walk', n=1000, burn_in=0, thinning=1, variance=1.0, bias_vector=None, solver=None, ess=1000
    ):
        """A static function to sample from a full dimensional polytope with an MCMC method.

        Keyword arguments:
        A -- an mxn matrix that contains the normal vectors of the facets of the polytope row-wise
        b -- a m-dimensional vector, s.t. A*x <= b
        method -- An MCMC method to sample, i.e. {'billiard_walk', 'cdhr', 'rdhr', 'ball_walk', 'dikin_walk', 'john_walk', 'vaidya_walk', 'gaussian_hmc_walk', 'exponential_hmc_walk', 'hmc_leapfrog_gaussian', 'hmc_leapfrog_exponential', 'shake_and_bake', 'billiard_shake_and_bake'}
        n -- the number of steady states to sample
        burn_in -- the number of points to burn before sampling
        thinning -- the walk length of the chain
        """
        if bias_vector is None:
            bias_vector = np.ones(A.shape[1], dtype=np.float64)
        else:
            bias_vector = bias_vector.astype('float64')

        P = HPolytope(A, b)

        samples = P.generate_samples(method.encode('utf-8'), n, burn_in, thinning, variance, bias_vector, solver, ess)

        samples_T = samples.T
        return samples_T

    @staticmethod
    def _parse_boundary_sampler_params(
        d: int,
        sampler: str,
        walk_len: Optional[int],
        nreflections: Optional[int],
    ) -> Tuple[str, int, int]:
        """
        Walk length defaults to 1. Reflections default to 0.25 * d for Billiard Shake and Bake, and 0 otherwise.
        """
        wl = int(walk_len) if walk_len is not None else 1

        s = str(sampler).strip().lower()
        if s in ("sb", "shake_and_bake"):
            sampler_key = "sb"
            use_bsb = False
        elif s in ("bsb", "billiard_shake_and_bake"):
            sampler_key = "bsb"
            use_bsb = True
        else:
            raise ValueError(
                'sampler must be one of {"sb","shake_and_bake","bsb","billiard_shake_and_bake"}'
            )

        if nreflections is not None:
            nref = int(nreflections)
        else:
            nref = int(0.25 * d) if use_bsb else 0

        return sampler_key, wl, nref


    @staticmethod
    def _first_k_exceeding_minESS(
        P: HPolytope,S_all: np.ndarray,ess_target: int,
    ) -> int:
        """
        Smallest k such that minESS(prefix k) >= ess_target.
        """
        S_all = np.asarray(S_all, dtype=np.float64, order="F")
        Ntot = int(S_all.shape[1])
        lo, hi = 1, Ntot

        while lo < hi:
            mid = (lo + hi) // 2
            if P.boundary_diag(S_all[:, :mid])["minESS"] >= ess_target:
                hi = mid
            else:
                lo = mid + 1
        return lo


    @staticmethod
    def boundary_sample_n(
        A,b,n: int = 1000,burn_in: int = 0,sampler: str = "sb",walk_len: Optional[int] = None,nreflections: Optional[int] = None,
    ) -> Tuple[np.ndarray, Dict]:
        """
        One boundary-sampling call for a predefined number of samples.
        """
        A = np.ascontiguousarray(A, dtype=np.float64)
        b = np.ascontiguousarray(b, dtype=np.float64)
        P = HPolytope(A, b)
        d = int(A.shape[1])

        sampler_key, wl, nref = PolytopeSampler._parse_boundary_sampler_params(d=d,sampler=sampler,walk_len=walk_len,nreflections=nreflections,)

        S = P.boundary_sample(sampler=sampler_key,number_of_points=int(n),number_of_points_to_burn=int(burn_in),walk_len=int(wl),nreflections=int(nref),)
        S = np.asarray(S, dtype=np.float64, order="F")

        diag = P.boundary_diag(S)

        info = {
            "minESS": float(diag["minESS"]),
            "maxPSRF": float(diag["maxPSRF"]),
            "N": int(diag["N"]),
            "calls": 1,
            "sampler": sampler_key,
            "walk_len": int(wl),
            "nreflections": int(nref),
        }

        return np.asarray(S), info


    @staticmethod
    def boundary_sample_ess(
        A,b,ess_target: int = 1000,chunk_n: int = 5000, burn_in_first: int = 0,sampler: str = "sb",walk_len: Optional[int] = None,nreflections: Optional[int] = None,max_calls: int = 1000,
    ) -> Tuple[np.ndarray, Dict]:
        """
        Iteratively samples the boundary of an H-polytope = in chunks until the 
        minimum Effective Sample Size (minESS) across all dimensions reaches a specified target.

        This function executes a specified random walk on the polytope boundary and accumulates 
        samples in discrete batches of size `chunk_n`. After each batch is appended to the 
        cumulative chain, MCMC diagnostics are evaluated. If the overall minESS meets or 
        exceeds `ess_target`, the function identifies the exact step index (`k_star`) where 
        the target was first achieved. The cumulative sample matrix is then strictly truncated 
        to this length to avoid returning unnecessary over-sampled points.

        Args:
            ess_target (int): The target minimum Effective Sample Size to achieve.
            chunk_n (int): Number of points to sample in each chunk (batch size).
            burn_in_first (int): Number of initial samples to discard as burn-in (applied to the first chunk only).
            max_calls (int): Maximum number of chunk iterations allowed to prevent infinite loops.

        Returns:
            S_star (np.ndarray): The truncated sample matrix of shape (d, k_star) that exactly meets the minESS target.
            info (dict): A dictionary containing the final MCMC diagnostics (minESS, maxPSRF), sampling parameters, total iterations (calls), and the truncation index (k_star).
        """
        ess_target = int(ess_target)
        A = np.ascontiguousarray(A, dtype=np.float64)
        b = np.ascontiguousarray(b, dtype=np.float64)

        P = HPolytope(A, b)
        d = int(A.shape[1])

        sampler_key, wl, nref = PolytopeSampler._parse_boundary_sampler_params(d=d,sampler=sampler,walk_len=walk_len,nreflections=nreflections,)
        S_all = np.zeros((d, 0), dtype=np.float64, order="F")
        calls = 0

        while calls < int(max_calls):
            calls += 1
            burn = int(burn_in_first) if S_all.shape[1] == 0 else 0

            S_chunk = P.boundary_sample(sampler=sampler_key,number_of_points=int(chunk_n),number_of_points_to_burn=int(burn),walk_len=int(wl),nreflections=int(nref))
            S_chunk = np.asarray(S_chunk, dtype=np.float64, order="F")
            S_all = np.asfortranarray(np.concatenate([S_all, S_chunk], axis=1))

            diag_all = P.boundary_diag(S_all)
            if diag_all["minESS"] < ess_target:
                continue

            k_star = PolytopeSampler._first_k_exceeding_minESS(P, S_all, ess_target)
            S_star = S_all[:, :k_star]
            diag_star = P.boundary_diag(S_star)

            info = {
                "minESS": float(diag_star["minESS"]),
                "maxPSRF": float(diag_star["maxPSRF"]),
                "N": int(diag_star["N"]),
                "calls": int(calls),
                "chunk_n": int(chunk_n),
                "sampler": sampler_key,
                "walk_len": int(wl),
                "nreflections": int(nref),
                "ess_target": int(ess_target),
                "k_star": int(k_star),
                "total_N": int(S_all.shape[1]),
            }

            return np.asarray(S_star), info

        last_minESS = (
            P.boundary_diag(S_all)["minESS"]
            if S_all.shape[1]
            else "N/A"
        )

        raise RuntimeError(
            f"ESS target not reached after max_calls={max_calls} "
            f"(last minESS={last_minESS})."
        )


    def boundary_diagnostics(self, S):
        """
        Diagnostics under current instance polytope.
        """
        self.get_polytope()

        S = np.asarray(S, dtype=np.float64, order="F")

        if S.shape[0] != self._A.shape[1]:
            raise ValueError(
                "S must have shape (d, N), where d = number of variables."
            )

        P = HPolytope(self._A, self._b)
        return P.boundary_diag(S)
    
    @staticmethod
    def _facet_coverage_count_from_sr(coverage_mat):
        """
        Counts the number of covered facets from a scaling ratio matrix.
        """
        cov = np.asarray(coverage_mat, dtype=float)
        if cov.ndim != 2:
            return 0
        finite_row = np.any(np.isfinite(cov), axis=1)
        return int(np.sum(finite_row))

    @staticmethod
    def boundary_sample_coverage(
        A,b,target_pcts=(10, 20, 30, 40, 50, 60, 70, 80, 90, 100),sampler="bsb",chunk_n=5000,burn_in_first=0,walk_len=1,nreflections=None,max_calls=200,sr_tol=1e-10,sr_min_ratio=0.01,
    ):
        """
        Samples an H-polytope boundary in chunks until specified percentages of facets are adequately visited. 
        A facet counts as covered only when it has enough samples to produce a finite scaling ratio.

        Returns:
            rows (list[dict]): Diagnostics and metrics saved at each target coverage milestone.
            final (dict): Total execution time and final sample counts.
        """
        A = np.asarray(A, dtype=np.float64, order="C")
        b = np.asarray(b, dtype=np.float64, order="C")
        P = HPolytope(A, b)

        m = int(A.shape[0])
        d = int(A.shape[1])

        targets = []
        for tp in target_pcts:
            tc = int(math.ceil((float(tp) / 100.0) * m)) if m > 0 else 0
            targets.append((int(tp), int(tc)))

        rows = []
        next_idx = 0
        S_acc = None
        t0 = time.perf_counter()
        calls = 0
        while calls < max_calls and next_idx < len(targets):
            calls += 1

            S_chunk = P.boundary_sample(
                sampler=sampler,
                number_of_points=int(chunk_n),
                number_of_points_to_burn=(int(burn_in_first) if calls == 1 else 0),
                walk_len=int(walk_len),
                nreflections=(0 if nreflections is None else int(nreflections)),
            )
            S_chunk = np.asarray(S_chunk, dtype=np.float64, order="F")

            if S_acc is None:
                S_acc = S_chunk
            else:
                S_acc = np.concatenate([S_acc, S_chunk], axis=1)

            # ESS/PSRF
            diag = P.boundary_diag(S_acc)
            minESS = float(diag.get("minESS", float("nan")))
            maxPSRF = float(diag.get("maxPSRF", float("nan")))
            N_samples = int(diag.get("N", S_acc.shape[1]))

            # SR + zero facets
            scale, coverage, max_dev, avg_dev, zc, zpct = P.boundary_scaling_ratio(
                S_acc, tol=float(sr_tol), min_ratio=float(sr_min_ratio)
            )

            covered_facets = PolytopeSampler._facet_coverage_count_from_sr(coverage)
            covered_pct = (100.0 * covered_facets / m) if m > 0 else 0.0

            max_dev_arr = np.asarray(max_dev, dtype=float) if max_dev is not None else None
            avg_dev_arr = np.asarray(avg_dev, dtype=float) if avg_dev is not None else None

            sr_max_dev = float(np.nanmax(max_dev_arr)) if max_dev_arr is not None else float("nan")
            sr_avg_dev_global = float(np.nanmean(avg_dev_arr)) if avg_dev_arr is not None else float("nan")

            elapsed = time.perf_counter() - t0

            while next_idx < len(targets) and covered_facets >= targets[next_idx][1]:
                tp, tc = targets[next_idx]
                rows.append({
                    "target_pct": int(tp),
                    "target_count": int(tc),
                    "calls": int(calls),
                    "chunk_n": int(chunk_n),
                    "max_calls": int(max_calls),
                    "dim": int(d),
                    "m": int(m),

                    "N_samples": int(N_samples),
                    "covered_facets": int(covered_facets),
                    "covered_pct": float(covered_pct),

                    "minESS": float(minESS),
                    "maxPSRF": float(maxPSRF),
                    "elapsed_sec": float(elapsed),

                    "sr_max_dev": float(sr_max_dev),
                    "sr_avg_dev_global": float(sr_avg_dev_global),
                    "zero_facets": int(zc),
                    "zero_facets_pct": float(zpct),
                })
                next_idx += 1

        # final state 
        final = {
            "calls": int(calls),
            "dim": int(d),
            "m": int(m),
            "N_samples": int(S_acc.shape[1]) if S_acc is not None else 0,
            "elapsed_sec": float(time.perf_counter() - t0),
        }
        return rows, final

    @staticmethod
    def round_polytope(
        A, b, method = "john_position", solver = None
    ):
        P = HPolytope(A, b)
        A, b, Tr, Tr_shift, round_value = P.rounding(method, solver)

        return A, b, Tr, Tr_shift

    @staticmethod
    def sample_from_fva_output(
        min_fluxes,
        max_fluxes,
        objective_function,
        max_objective,
        S,
        opt_percentage=100,
        ess=1000,
        psrf=False,
        parallel_mmcs=False,
        num_threads=1,
        solver = None
    ):
        """A static function to sample steady states when the output of FVA is given.

        Keyword arguments:
        min_fluxes -- minimum values of the fluxes, i.e., a n-dimensional vector
        max_fluxes -- maximum values for the fluxes, i.e., a n-dimensional vector
        objective_function -- the objective function
        max_objective -- the maximum value of the objective function
        S -- stoichiometric matrix
        opt_percentage -- consider solutions that give you at least a certain
                      percentage of the optimal solution (default is to consider
                      optimal solutions only)
        ess -- the target effective sample size
        psrf -- a boolean flag to request PSRF smaller than 1.1 for all marginal fluxes
        parallel_mmcs -- a boolean flag to request the parallel mmcs
        num_threads -- the number of threads to use for parallel mmcs
        """

        A, b, Aeq, beq = get_matrices_of_low_dim_polytope(
            S, min_fluxes, max_fluxes, opt_percentage, tol
        )

        A = np.vstack((A, -objective_function))
        b = np.append(
            b,
            -(opt_percentage / 100)
            * self._parameters["tol"]
            * math.floor(max_objective / self._parameters["tol"]),
        )

        A, b, N, N_shift = get_matrices_of_full_dim_polytope(A, b, Aeq, beq)

        P = HPolytope(A, b)

        A, b, Tr, Tr_shift, samples = P.mmcs(
            ess, psrf, parallel_mmcs, num_threads, solver
        )

        steady_states = map_samples_to_steady_states(samples, N, N_shift)

        return steady_states

    @property
    def A(self):
        return self._A

    @property
    def b(self):
        return self._b

    @property
    def T(self):
        return self._T

    @property
    def T_shift(self):
        return self._T_shift

    @property
    def N(self):
        return self._N

    @property
    def N_shift(self):
        return self._N_shift

    @property
    def metabolic_network(self):
        return self._metabolic_network

    def facet_redundancy_removal(self, value):
        self._parameters["remove_redundant_facets"] = value

    def set_solver(self, solver):
        self._parameters["solver"] = solver

    def set_distribution(self, value):

        self._parameters["distribution"] = value

    def set_nullspace_method(self, value):

        self._parameters["nullspace_method"] = value

    def set_tol(self, value):

        self._parameters["tol"] = value

    def set_opt_percentage(self, value):

        self._parameters["opt_percentage"] = value