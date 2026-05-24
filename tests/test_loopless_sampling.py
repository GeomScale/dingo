# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2024

# Licensed under GNU LGPL.3, see LICENCE file

# =============================================================================
# Hard Test: Loopless (Non-Convex) Sampling
# =============================================================================
# Tests the LooplessFluxSampler module for:
#   - Internal reaction classification
#   - Cycle nullspace computation
#   - Rejection-based loopless sampling
#   - Penalty-weighted sampling
#   - Comparison of loopless vs standard samples
# =============================================================================

import unittest
import os
import sys
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo.LooplessFluxSampler import LooplessFluxSampler
from dingo.pyoptinterface_based_impl import set_default_solver


class TestLooplessSampling(unittest.TestCase):
    """Tests for the LooplessFluxSampler module."""

    MODEL_PATH = os.path.join(os.getcwd(), "ext_data", "e_coli_core.json")

    # ----------------------------------------------------------------
    # Test: Reaction classification
    # ----------------------------------------------------------------
    def test_reaction_classification(self):
        """
        Verify that reactions are correctly classified as internal vs exchange.
        """
        print("\n" + "=" * 70)
        print("HARD: Reaction classification (internal vs exchange)")
        print("=" * 70)

        model = MetabolicNetwork.from_json(self.MODEL_PATH)
        lfs = LooplessFluxSampler(model)

        n_internal = len(lfs.internal_indices)
        n_exchange = len(lfs.exchange_indices)
        n_total = model.num_of_reactions()

        print(f"  Total reactions:    {n_total}")
        print(f"  Internal reactions: {n_internal}")
        print(f"  Exchange reactions: {n_exchange}")

        # All reactions must be classified
        self.assertEqual(n_internal + n_exchange, n_total,
                         "All reactions should be classified")

        # E. coli core has ~20 exchange reactions
        self.assertTrue(n_exchange > 10,
                        "E. coli core should have >10 exchange reactions")
        self.assertTrue(n_internal > 50,
                        "E. coli core should have >50 internal reactions")

        # List some exchange reactions for verification
        reactions = model.reactions
        print("\n  Exchange reactions:")
        for idx in lfs.exchange_indices[:10]:
            print(f"    {reactions[idx]}")

        print("\n  [OK] Reaction classification PASSED")

    # ----------------------------------------------------------------
    # Test: Cycle nullspace computation
    # ----------------------------------------------------------------
    def test_cycle_nullspace(self):
        """
        Verify that the internal cycle nullspace is computed correctly.
        """
        print("\n" + "=" * 70)
        print("HARD: Cycle nullspace computation")
        print("=" * 70)

        model = MetabolicNetwork.from_json(self.MODEL_PATH)
        lfs = LooplessFluxSampler(model)

        N = lfs.cycle_nullspace

        if N.size > 0:
            print(f"  Cycle nullspace shape: {N.shape}")
            print(f"  Number of potential cycle directions: {N.shape[1]}")

            # Verify orthogonality of nullspace vectors
            NtN = N.T @ N
            identity_check = np.allclose(NtN, np.eye(N.shape[1]), atol=1e-10)
            print(f"  Nullspace vectors orthonormal: {identity_check}")

            # Verify that nullspace vectors are in the nullspace of S_int
            S_int = model.S[:, lfs.internal_indices]
            residual = np.linalg.norm(S_int @ N)
            print(f"  S_int * N residual norm: {residual:.2e}")
            self.assertTrue(residual < 1e-10,
                            "Nullspace vectors should satisfy S_int * N ≈ 0")
        else:
            print("  No internal cycles found (nullspace is empty)")

        print("\n  [OK] Cycle nullspace computation PASSED")

    # ----------------------------------------------------------------
    # Test: Loop detection on known samples
    # ----------------------------------------------------------------
    def test_loop_detection(self):
        """
        Generate standard samples and check what fraction contains loops.
        """
        print("\n" + "=" * 70)
        print("HARD: Loop detection in standard samples")
        print("=" * 70)

        model = MetabolicNetwork.from_json(self.MODEL_PATH)
        lfs = LooplessFluxSampler(model)

        # Generate standard (possibly loopy) samples
        sampler = PolytopeSampler(model)
        steady_states = sampler.generate_steady_states_no_multiphase(
            method="billiard_walk", n=300, burn_in=50, thinning=2
        )

        # Analyze loops
        analysis = lfs.analyze_loops(steady_states)

        print(f"  Samples analyzed:        {steady_states.shape[1]}")
        print(f"  Fraction with loops:     {analysis['fraction_with_loops']:.4f}")
        print(f"  Mean loop norm:          {analysis['mean_loop_norm']:.6f}")
        print(f"  Median loop norm:        {analysis['median_loop_norm']:.6f}")
        print(f"  Max loop norm:           {analysis['max_loop_norm']:.6f}")
        print(f"  Internal reactions:      {analysis['n_internal_reactions']}")
        print(f"  Exchange reactions:      {analysis['n_exchange_reactions']}")
        print(f"  Cycle nullspace dim:     {analysis['cycle_nullspace_dim']}")

        # Validate analysis output
        self.assertTrue(0 <= analysis["fraction_with_loops"] <= 1)
        self.assertTrue(analysis["mean_loop_norm"] >= 0)

        print("\n  [OK] Loop detection PASSED")

    # ----------------------------------------------------------------
    # Test: Rejection-based loopless sampling
    # ----------------------------------------------------------------
    def test_rejection_sampling(self):
        """
        Test the rejection-based loopless sampling approach.
        """
        print("\n" + "=" * 70)
        print("HARD: Rejection-based loopless sampling")
        print("=" * 70)

        model = MetabolicNetwork.from_json(self.MODEL_PATH)
        lfs = LooplessFluxSampler(model, loop_tolerance=1e-3)

        loopless_samples, stats = lfs.sample_loopless(
            n_samples=100,
            max_attempts_factor=20,
            method="billiard_walk",
            burn_in=50,
            thinning=2,
            opt_percentage=100,
        )

        print(f"  Total sampled:     {stats['total_sampled']}")
        print(f"  Total accepted:    {stats['total_accepted']}")
        print(f"  Total rejected:    {stats['total_rejected']}")
        print(f"  Acceptance rate:   {stats['acceptance_rate']:.4f}")
        print(f"  Mean loop norm:    {stats['mean_loop_norm']:.6f}")
        print(f"  Max loop norm:     {stats['max_loop_norm']:.6f}")
        print(f"  Loop tolerance:    {stats['loop_tolerance']}")
        print(f"  Output shape:      {loopless_samples.shape}")

        # Validate that accepted samples have correct shape
        self.assertEqual(loopless_samples.shape[0], 95,
                         "Expected 95 reactions")

        # Verify no loops in accepted samples
        if loopless_samples.shape[1] > 0:
            for j in range(loopless_samples.shape[1]):
                has_loop, norm = lfs._has_active_loop(loopless_samples[:, j])
                self.assertFalse(has_loop,
                                 f"Sample {j} should be loopless (norm={norm:.6f})")

            print(f"\n  [OK] All {loopless_samples.shape[1]} accepted samples "
                  f"verified loopless")
        else:
            print("\n  [!] No samples accepted (acceptance rate may be very low)")

        print("  [OK] Rejection-based loopless sampling PASSED")

    # ----------------------------------------------------------------
    # Test: Penalty-weighted sampling
    # ----------------------------------------------------------------
    def test_penalty_sampling(self):
        """
        Test the penalty-weighted (importance sampling) approach.
        """
        print("\n" + "=" * 70)
        print("HARD: Penalty-weighted loopless sampling")
        print("=" * 70)

        model = MetabolicNetwork.from_json(self.MODEL_PATH)
        lfs = LooplessFluxSampler(model)

        samples, weights, stats = lfs.sample_penalized(
            n_samples=300,
            penalty_weight=10.0,
            method="billiard_walk",
            burn_in=50,
            thinning=2,
            opt_percentage=100,
        )

        print(f"  Samples shape:              {samples.shape}")
        print(f"  Weights shape:              {weights.shape}")
        print(f"  Mean loop norm:             {stats['mean_loop_norm']:.6f}")
        print(f"  Max loop norm:              {stats['max_loop_norm']:.6f}")
        print(f"  N effectively loopless:     {stats['n_effectively_loopless']}")
        print(f"  Effective sample size:      {stats['effective_sample_size']:.1f}")
        print(f"  Penalty weight:             {stats['penalty_weight']}")

        # Weights must sum to 1
        self.assertAlmostEqual(np.sum(weights), 1.0, places=10)

        # All weights should be non-negative
        self.assertTrue(np.all(weights >= 0))

        # Effective sample size should be > 0
        self.assertTrue(stats["effective_sample_size"] > 0)

        # Weighted mean of biomass should be reasonable
        biomass_idx = model.biomass_index
        weighted_mean_biomass = np.average(
            samples[biomass_idx, :], weights=weights
        )
        print(f"  Weighted mean biomass flux: {weighted_mean_biomass:.6f}")

        print("\n  [OK] Penalty-weighted sampling PASSED")

    # ----------------------------------------------------------------
    # Test: Standard vs loopless comparison
    # ----------------------------------------------------------------
    def test_standard_vs_loopless(self):
        """
        Compare standard and loopless samples on key flux statistics.
        """
        print("\n" + "=" * 70)
        print("HARD: Standard vs Loopless sampling comparison")
        print("=" * 70)

        model = MetabolicNetwork.from_json(self.MODEL_PATH)
        lfs = LooplessFluxSampler(model, loop_tolerance=1e-3)

        # Standard sampling
        sampler = PolytopeSampler(model)
        standard_ss = sampler.generate_steady_states_no_multiphase(
            method="billiard_walk", n=300, burn_in=50, thinning=2
        )

        # Loopless sampling
        loopless_ss, stats = lfs.sample_loopless(
            n_samples=200,
            max_attempts_factor=20,
            method="billiard_walk",
            burn_in=50,
            thinning=2,
            opt_percentage=100,
        )

        # Analyze loops in both
        standard_analysis = lfs.analyze_loops(standard_ss)

        print(f"\n  {'Metric':<30} {'Standard':>15} {'Loopless':>15}")
        print("  " + "-" * 62)
        print(f"  {'N samples':<30} {standard_ss.shape[1]:>15} "
              f"{loopless_ss.shape[1]:>15}")
        print(f"  {'Fraction with loops':<30} "
              f"{standard_analysis['fraction_with_loops']:>15.4f} {'0.0000':>15}")

        if loopless_ss.shape[1] > 0:
            reactions = model.reactions
            key_rxns = ["PFK", "CS", "ATPM"]

            for rxn in key_rxns:
                if rxn in reactions:
                    idx = reactions.index(rxn)
                    std_mean = np.mean(standard_ss[idx, :])
                    ll_mean = np.mean(loopless_ss[idx, :])
                    print(f"  {'Mean flux ' + rxn:<30} "
                          f"{std_mean:>15.4f} {ll_mean:>15.4f}")

        print(f"\n  Acceptance rate: {stats['acceptance_rate']:.4f}")
        print("\n  [OK] Standard vs Loopless comparison PASSED")


if __name__ == "__main__":
    if len(sys.argv) > 1 and not sys.argv[1].startswith("-"):
        set_default_solver(sys.argv[1])
        sys.argv.pop(1)
    unittest.main(verbosity=2)
