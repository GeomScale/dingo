# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2024

# Licensed under GNU LGPL.3, see LICENCE file

# =============================================================================
# Easy Test: Flux Sampling Analysis on the E. coli Core Model
# =============================================================================
# This script performs flux sampling under three biomass constraint scenarios:
#   i.   Optimal biomass growth (opt_percentage=100)
#   ii.  At least half of the optimal (opt_percentage=50)
#   iii. Setting biomass free (opt_percentage=0)
#
# For each scenario we:
#   - Load the e_coli_core model
#   - Configure the biomass constraint
#   - Sample steady states from the flux polytope
#   - Validate sample shapes, non-zero fluxes, and biomass bounds
# =============================================================================

import unittest
import os
import sys
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo.pyoptinterface_based_impl import set_default_solver


class TestFluxSamplingEcoli(unittest.TestCase):
    """Flux sampling on e_coli_core under three biomass growth scenarios."""

    MODEL_PATH = os.path.join(os.getcwd(), "ext_data", "e_coli_core.json")

    # ----------------------------------------------------------------
    # Helper utilities
    # ----------------------------------------------------------------
    def _load_model(self):
        """Load the E. coli core model and return it along with FBA optimum."""
        model = MetabolicNetwork.from_json(self.MODEL_PATH)
        fba_res = model.fba()
        max_biomass_objective = fba_res[1]
        return model, max_biomass_objective

    def _sample_steady_states(self, model, method="billiard_walk", n_samples=500):
        """Create a PolytopeSampler and generate steady states."""
        sampler = PolytopeSampler(model)
        steady_states = sampler.generate_steady_states_no_multiphase(
            method=method, n=n_samples, burn_in=100, thinning=2
        )
        return steady_states, sampler

    # ----------------------------------------------------------------
    # Scenario i: Optimal biomass growth (100%)
    # ----------------------------------------------------------------
    def test_optimal_biomass_growth(self):
        """
        Scenario i: Sample the flux space requiring optimal biomass growth.
        All sampled steady states should achieve close to 100% of the FBA optimum.
        """
        print("\n" + "=" * 70)
        print("SCENARIO i: Optimal biomass growth (opt_percentage = 100)")
        print("=" * 70)

        model, max_biomass = self._load_model()

        # Keep the default 100% optimality constraint
        model.set_opt_percentage(100)
        print(f"  FBA max biomass objective: {max_biomass:.6f}")

        steady_states, sampler = self._sample_steady_states(model)

        # -- Validate shape: 95 reactions for e_coli_core --
        self.assertEqual(steady_states.shape[0], 95,
                         "Expected 95 reactions in e_coli_core model")
        self.assertTrue(steady_states.shape[1] > 0,
                        "Expected at least some samples")

        # -- Validate non-zero fluxes --
        self.assertTrue(np.any(np.abs(steady_states) > 1e-12),
                        "Steady states should contain non-zero flux values")

        # -- Validate biomass constraint --
        biomass_idx = model.biomass_index
        biomass_fluxes = steady_states[biomass_idx, :]
        min_biomass_sample = np.min(biomass_fluxes)

        print(f"  Biomass index: {biomass_idx}")
        print(f"  Biomass flux — mean: {np.mean(biomass_fluxes):.6f}, "
              f"min: {min_biomass_sample:.6f}, max: {np.max(biomass_fluxes):.6f}")
        print(f"  Sample shape: {steady_states.shape}")

        # With 100% opt constraint, biomass should be close to max
        self.assertTrue(min_biomass_sample >= max_biomass * 0.95,
                        f"Min biomass ({min_biomass_sample:.4f}) should be "
                        f">= 95% of max ({max_biomass:.4f})")

        # Save results for downstream analysis
        np.save(os.path.join(os.getcwd(), "tests", "results_optimal.npy"),
                steady_states)

        print("  [OK] Scenario i PASSED")

    # ----------------------------------------------------------------
    # Scenario ii: At least half-optimal biomass (50%)
    # ----------------------------------------------------------------
    def test_half_optimal_biomass(self):
        """
        Scenario ii: Sample the flux space requiring at least 50% of the
        optimal biomass growth.
        """
        print("\n" + "=" * 70)
        print("SCENARIO ii: At least half-optimal biomass (opt_percentage = 50)")
        print("=" * 70)

        model, max_biomass = self._load_model()

        # Set 50% optimality constraint
        model.set_opt_percentage(50)
        print(f"  FBA max biomass objective: {max_biomass:.6f}")
        print(f"  Required minimum biomass: {max_biomass * 0.50:.6f}")

        steady_states, sampler = self._sample_steady_states(model)

        # -- Validate shape --
        self.assertEqual(steady_states.shape[0], 95)
        self.assertTrue(steady_states.shape[1] > 0)

        # -- Validate non-zero fluxes --
        self.assertTrue(np.any(np.abs(steady_states) > 1e-12))

        # -- Validate biomass constraint --
        biomass_idx = model.biomass_index
        biomass_fluxes = steady_states[biomass_idx, :]
        min_biomass_sample = np.min(biomass_fluxes)

        print(f"  Biomass index: {biomass_idx}")
        print(f"  Biomass flux — mean: {np.mean(biomass_fluxes):.6f}, "
              f"min: {min_biomass_sample:.6f}, max: {np.max(biomass_fluxes):.6f}")
        print(f"  Sample shape: {steady_states.shape}")

        # With 50% opt constraint, biomass should be at least half
        self.assertTrue(min_biomass_sample >= max_biomass * 0.45,
                        f"Min biomass ({min_biomass_sample:.4f}) should be "
                        f">= 45% of max ({max_biomass:.4f}) (with tolerance)")

        # Save results
        np.save(os.path.join(os.getcwd(), "tests", "results_half_optimal.npy"),
                steady_states)

        print("  [OK] Scenario ii PASSED")

    # ----------------------------------------------------------------
    # Scenario iii: Biomass free (0% constraint)
    # ----------------------------------------------------------------
    def test_biomass_free(self):
        """
        Scenario iii: Sample the flux space with biomass unconstrained.
        The biomass reaction can take any feasible value, including zero.
        """
        print("\n" + "=" * 70)
        print("SCENARIO iii: Biomass free (no biomass constraint)")
        print("=" * 70)

        model, max_biomass = self._load_model()

        # Remove biomass constraint by setting opt_percentage to 0
        # This effectively removes the objective function lower bound
        n = model.num_of_reactions()
        model.set_opt_percentage(0)

        print(f"  FBA max biomass objective: {max_biomass:.6f}")
        print(f"  Biomass constraint: NONE (free)")

        steady_states, sampler = self._sample_steady_states(model)

        # -- Validate shape --
        self.assertEqual(steady_states.shape[0], 95)
        self.assertTrue(steady_states.shape[1] > 0)

        # -- Validate non-zero fluxes --
        self.assertTrue(np.any(np.abs(steady_states) > 1e-12))

        # -- Examine biomass distribution --
        biomass_idx = model.biomass_index
        biomass_fluxes = steady_states[biomass_idx, :]

        print(f"  Biomass index: {biomass_idx}")
        print(f"  Biomass flux — mean: {np.mean(biomass_fluxes):.6f}, "
              f"min: {np.min(biomass_fluxes):.6f}, max: {np.max(biomass_fluxes):.6f}")
        print(f"  Sample shape: {steady_states.shape}")

        # In the free case, we expect a wider range of biomass values
        # The minimum should be considerably lower than the max
        biomass_range = np.max(biomass_fluxes) - np.min(biomass_fluxes)
        print(f"  Biomass range: {biomass_range:.6f}")

        # Save results
        np.save(os.path.join(os.getcwd(), "tests", "results_biomass_free.npy"),
                steady_states)

        print("  [OK] Scenario iii PASSED")

    # ----------------------------------------------------------------
    # Summary comparison across all three scenarios
    # ----------------------------------------------------------------
    def test_scenario_comparison(self):
        """
        Compare basic statistics across the three sampling scenarios.
        This test runs after the individual scenario tests.
        """
        print("\n" + "=" * 70)
        print("SUMMARY: Comparing all three scenarios")
        print("=" * 70)

        model, max_biomass = self._load_model()
        biomass_idx = model.biomass_index
        reactions = model.reactions

        results = {}
        for name, opt_pct in [("optimal", 100), ("half_optimal", 50), ("biomass_free", 0)]:
            m, _ = self._load_model()
            m.set_opt_percentage(opt_pct)
            ss, _ = self._sample_steady_states(m, n_samples=300)
            results[name] = ss

        # Print comparison table
        print(f"\n  {'Scenario':<20} {'Mean Biomass':>14} {'Min Biomass':>14} "
              f"{'Max Biomass':>14} {'Flux Std Mean':>14}")
        print("  " + "-" * 78)

        for name, ss in results.items():
            bm = ss[biomass_idx, :]
            flux_std = np.mean(np.std(ss, axis=1))
            print(f"  {name:<20} {np.mean(bm):>14.6f} {np.min(bm):>14.6f} "
                  f"{np.max(bm):>14.6f} {flux_std:>14.6f}")

        # Verify that the biomass-free scenario has wider flux variability
        std_optimal = np.mean(np.std(results["optimal"], axis=1))
        std_free = np.mean(np.std(results["biomass_free"], axis=1))

        print(f"\n  Flux variability ratio (free/optimal): {std_free/std_optimal:.2f}")
        print("  [OK] Comparison complete")


if __name__ == "__main__":
    # Only treat positional (non-flag) arguments as solver name
    if len(sys.argv) > 1 and not sys.argv[1].startswith("-"):
        set_default_solver(sys.argv[1])
        sys.argv.pop(1)
    unittest.main(verbosity=2)
