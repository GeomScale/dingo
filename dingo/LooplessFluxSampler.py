# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2024

# Licensed under GNU LGPL.3, see LICENCE file

# =============================================================================
# LooplessFluxSampler: Non-convex sampling for thermodynamically feasible fluxes
# =============================================================================
#
# Flux sampling in metabolic networks can produce thermodynamically infeasible
# solutions: internal cycles that continuously consume/generate energy without
# net effect. These "Type III pathways" violate the second law of thermodynamics.
#
# This module implements:
#   1. Detection of internal cycles via nullspace analysis
#   2. Rejection-based loopless sampling (sample + filter)
#   3. A penalty-based approach that biases sampling away from loopy solutions
#
# References:
#   [1] Chalkis et al. - dingo: Python package for metabolic flux sampling
#   [2] De Martino - Scales and multimodal flux distributions via thermodynamics
#   [3] Saa et al. - LooplessFluxSampler
# =============================================================================

import numpy as np
from scipy.linalg import null_space
from dingo.MetabolicNetwork import MetabolicNetwork
from dingo.PolytopeSampler import PolytopeSampler


class LooplessFluxSampler:
    """
    A sampler that generates thermodynamically feasible (loopless) steady states
    from a metabolic network by detecting and rejecting internal flux loops.

    Internal loops are Type III pathways: cycles of internal reactions that carry
    flux without being driven by any exchange reaction. They violate the second
    law of thermodynamics and are biologically meaningless.

    Strategy (rejection-based):
        1. Identify internal reactions (non-exchange, non-biomass).
        2. Compute the nullspace of the internal stoichiometric sub-matrix.
           Vectors in this nullspace represent potential internal cycles.
        3. Sample from the full flux polytope using dingo's standard samplers.
        4. For each sample, project the internal flux vector onto the cycle
           nullspace. If the projection norm exceeds a tolerance, the sample
           contains an active loop and is rejected.

    Parameters
    ----------
    metabolic_network : MetabolicNetwork
        A dingo MetabolicNetwork object.
    loop_tolerance : float, optional
        Maximum allowed norm of the internal-cycle projection (default 1e-5).
    """

    def __init__(self, metabolic_network, loop_tolerance=1e-5):
        if not isinstance(metabolic_network, MetabolicNetwork):
            raise TypeError("Expected a MetabolicNetwork object.")

        self._network = metabolic_network
        self._loop_tolerance = loop_tolerance

        # Identify internal vs exchange reactions
        self._internal_indices = []
        self._exchange_indices = []
        self._classify_reactions()

        # Compute the internal cycle nullspace
        self._cycle_nullspace = None
        self._compute_cycle_nullspace()

    # ------------------------------------------------------------------
    # Reaction classification
    # ------------------------------------------------------------------
    def _classify_reactions(self):
        """
        Classify reactions as internal or exchange.

        Exchange reactions are identified as those involving only one
        metabolite (single non-zero entry in the stoichiometric column),
        or reactions listed in the model's exchanges property.
        """
        S = self._network.S
        reactions = self._network.reactions
        exchanges = set(self._network.exchanges) if self._network.exchanges else set()
        n_reactions = S.shape[1]

        for j in range(n_reactions):
            col = S[:, j]
            n_nonzero = np.count_nonzero(col)

            # A reaction is considered exchange if:
            # 1) It has a single nonzero stoichiometric coefficient, or
            # 2) It is listed in the model's exchange reactions
            if n_nonzero <= 1 or reactions[j] in exchanges:
                self._exchange_indices.append(j)
            else:
                self._internal_indices.append(j)

    # ------------------------------------------------------------------
    # Internal cycle detection via nullspace
    # ------------------------------------------------------------------
    def _compute_cycle_nullspace(self):
        """
        Compute the nullspace of the internal stoichiometric sub-matrix.

        For internal reactions only, we extract S_int (the sub-matrix of S
        with columns corresponding to internal reactions). The nullspace of
        S_int gives directions in which internal reactions can carry flux
        while still satisfying S_int * v_int = 0 – these are internal cycles.
        """
        S = self._network.S
        int_idx = self._internal_indices

        if len(int_idx) == 0:
            self._cycle_nullspace = np.empty((0, 0))
            return

        # Extract the internal sub-matrix
        S_int = S[:, int_idx]

        # Compute the nullspace of S_int
        # Each column of N_int is a potential internal cycle direction
        N_int = null_space(S_int)

        self._cycle_nullspace = N_int

    def _has_active_loop(self, flux_vector):
        """
        Check whether a flux vector contains an active internal loop.

        Projects the internal flux sub-vector onto the cycle nullspace.
        If the projection has significant magnitude, the sample carries
        an internal cycle.

        Parameters
        ----------
        flux_vector : ndarray
            A complete flux vector (dimension = number of reactions).

        Returns
        -------
        bool
            True if an active internal loop is detected.
        float
            The norm of the loop projection (for diagnostics).
        """
        if self._cycle_nullspace.size == 0:
            return False, 0.0

        # Extract internal fluxes
        v_int = flux_vector[self._internal_indices]

        # Project onto the cycle nullspace: proj = N * N^T * v_int
        N = self._cycle_nullspace
        projection = N @ (N.T @ v_int)
        loop_norm = np.linalg.norm(projection)

        return loop_norm > self._loop_tolerance, loop_norm

    # ------------------------------------------------------------------
    # Loopless sampling (rejection-based)
    # ------------------------------------------------------------------
    def sample_loopless(
        self,
        n_samples=500,
        max_attempts_factor=10,
        method="billiard_walk",
        burn_in=100,
        thinning=2,
        opt_percentage=None,
    ):
        """
        Generate loopless steady-state flux samples using rejection sampling.

        Samples are drawn from the flux polytope and those containing
        internal loops are discarded.

        Parameters
        ----------
        n_samples : int
            Desired number of loopless samples.
        max_attempts_factor : int
            Maximum total samples to draw = n_samples * max_attempts_factor.
        method : str
            MCMC sampling method (default: 'billiard_walk').
        burn_in : int
            Number of burn-in samples (default: 100).
        thinning : int
            Thinning factor for the MCMC chain (default: 2).
        opt_percentage : int or None
            If not None, set the model's opt_percentage before sampling.

        Returns
        -------
        loopless_samples : ndarray
            Steady states without internal loops (shape: n_reactions x n_accepted).
        rejection_stats : dict
            Statistics about the rejection process.
        """
        if opt_percentage is not None:
            self._network.set_opt_percentage(opt_percentage)

        sampler = PolytopeSampler(self._network)

        max_total = n_samples * max_attempts_factor

        accepted = []
        total_sampled = 0
        total_rejected = 0
        loop_norms = []

        # Draw all samples at once (more efficient than batching)
        n_draw = min(max_total, max(n_samples * 5, 1000))
        steady_states = sampler.generate_steady_states_no_multiphase(
            method=method, n=n_draw, burn_in=burn_in, thinning=thinning
        )

        total_sampled = n_draw

        # Filter samples for loops
        for j in range(steady_states.shape[1]):
            flux = steady_states[:, j]
            has_loop, norm = self._has_active_loop(flux)
            loop_norms.append(norm)

            if not has_loop:
                accepted.append(flux)
                if len(accepted) >= n_samples:
                    break
            else:
                total_rejected += 1

        # Stack accepted samples
        if len(accepted) > 0:
            loopless_samples = np.column_stack(accepted[:n_samples])
        else:
            loopless_samples = np.empty((self._network.num_of_reactions(), 0))

        rejection_stats = {
            "total_sampled": total_sampled,
            "total_accepted": len(accepted),
            "total_rejected": total_rejected,
            "acceptance_rate": len(accepted) / max(total_sampled, 1),
            "mean_loop_norm": float(np.mean(loop_norms)) if loop_norms else 0.0,
            "max_loop_norm": float(np.max(loop_norms)) if loop_norms else 0.0,
            "loop_tolerance": self._loop_tolerance,
        }

        return loopless_samples, rejection_stats

    # ------------------------------------------------------------------
    # Penalty-weighted sampling
    # ------------------------------------------------------------------
    def sample_penalized(
        self,
        n_samples=500,
        penalty_weight=10.0,
        method="billiard_walk",
        burn_in=100,
        thinning=2,
        opt_percentage=None,
    ):
        """
        Generate flux samples where loopy samples are down-weighted
        using importance sampling with a loop-penalty function.

        Instead of hard rejection, each sample receives a weight:
            w_i = exp(-penalty_weight * ||projection_i||^2)

        This produces a weighted sample set that favours loopless solutions.

        Parameters
        ----------
        n_samples : int
            Number of samples to draw (all are kept, but weighted).
        penalty_weight : float
            Strength of the penalty (higher = more bias against loops).
        method : str
            MCMC sampling method.
        burn_in : int
            Burn-in samples.
        thinning : int
            Thinning factor.
        opt_percentage : int or None
            If not None, set the model's opt_percentage before sampling.

        Returns
        -------
        samples : ndarray
            All steady state samples (shape: n_reactions x n_samples).
        weights : ndarray
            Importance weights for each sample (shape: n_samples,).
        penalty_stats : dict
            Statistics about the penalty distribution.
        """
        if opt_percentage is not None:
            self._network.set_opt_percentage(opt_percentage)

        sampler = PolytopeSampler(self._network)
        steady_states = sampler.generate_steady_states_no_multiphase(
            method=method, n=n_samples, burn_in=burn_in, thinning=thinning
        )

        weights = np.zeros(steady_states.shape[1])
        loop_norms = np.zeros(steady_states.shape[1])

        for j in range(steady_states.shape[1]):
            flux = steady_states[:, j]
            _, norm = self._has_active_loop(flux)
            loop_norms[j] = norm

        # Compute weights in log-space to avoid underflow
        # w_i = exp(-penalty_weight * norm_i^2)
        log_weights = -penalty_weight * loop_norms ** 2
        # Shift so the maximum log-weight is 0 (prevents all-zero underflow)
        log_weights -= np.max(log_weights)
        weights = np.exp(log_weights)

        # Normalise weights
        weight_sum = np.sum(weights)
        if weight_sum > 0:
            weights /= weight_sum
        else:
            # Uniform fallback if all weights are zero
            weights[:] = 1.0 / len(weights)

        penalty_stats = {
            "mean_loop_norm": float(np.mean(loop_norms)),
            "max_loop_norm": float(np.max(loop_norms)),
            "n_effectively_loopless": int(np.sum(loop_norms < self._loop_tolerance)),
            "effective_sample_size": float(1.0 / np.sum(weights ** 2)),
            "penalty_weight": penalty_weight,
        }

        return steady_states, weights, penalty_stats

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------
    def analyze_loops(self, steady_states):
        """
        Analyze the loop content of a set of steady-state samples.

        Parameters
        ----------
        steady_states : ndarray
            Steady states (shape: n_reactions x n_samples).

        Returns
        -------
        dict
            Loop analysis results including per-sample norms, fraction
            with loops, and the most common loop directions.
        """
        n_samples = steady_states.shape[1]
        loop_norms = np.zeros(n_samples)
        has_loops = np.zeros(n_samples, dtype=bool)

        for j in range(n_samples):
            has_loop, norm = self._has_active_loop(steady_states[:, j])
            loop_norms[j] = norm
            has_loops[j] = has_loop

        return {
            "loop_norms": loop_norms,
            "fraction_with_loops": float(np.mean(has_loops)),
            "mean_loop_norm": float(np.mean(loop_norms)),
            "median_loop_norm": float(np.median(loop_norms)),
            "max_loop_norm": float(np.max(loop_norms)),
            "n_internal_reactions": len(self._internal_indices),
            "n_exchange_reactions": len(self._exchange_indices),
            "cycle_nullspace_dim": self._cycle_nullspace.shape[1]
            if self._cycle_nullspace.size > 0 else 0,
        }

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------
    @property
    def internal_indices(self):
        return self._internal_indices

    @property
    def exchange_indices(self):
        return self._exchange_indices

    @property
    def cycle_nullspace(self):
        return self._cycle_nullspace

    @property
    def loop_tolerance(self):
        return self._loop_tolerance

    @loop_tolerance.setter
    def loop_tolerance(self, value):
        self._loop_tolerance = value

    @property
    def network(self):
        return self._network
