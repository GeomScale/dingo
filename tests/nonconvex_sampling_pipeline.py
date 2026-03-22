# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2024

# Licensed under GNU LGPL.3, see LICENCE file

# =============================================================================
# Non-Convex Sampling Pipeline
# =============================================================================
# End-to-end pipeline demonstrating:
#   1. Flux sampling under three biomass scenarios
#   2. Loop analysis and loopless sampling
#   3. Reaction clustering & biological interpretation
#   4. Comparison of standard vs loopless sample properties
# =============================================================================

import unittest
import os
import sys
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo.LooplessFluxSampler import LooplessFluxSampler
from dingo.utils import correlated_reactions, cluster_corr_reactions
from dingo.pyoptinterface_based_impl import set_default_solver


class TestNonConvexPipeline(unittest.TestCase):
    """End-to-end non-convex sampling pipeline on E. coli core."""

    MODEL_PATH = os.path.join(os.getcwd(), "ext_data", "e_coli_core.json")

    # ----------------------------------------------------------------
    # Full pipeline test
    # ----------------------------------------------------------------
    def test_full_pipeline(self):
        """
        Complete pipeline: sample → detect loops → filter → cluster → interpret.
        """
        print("\n" + "=" * 70)
        print("PIPELINE: Non-Convex Sampling End-to-End")
        print("=" * 70)

        model = MetabolicNetwork.from_json(self.MODEL_PATH)
        reactions = model.reactions
        fba_result = model.fba()
        max_biomass = fba_result[1]
        biomass_idx = model.biomass_index

        print(f"\n  Model: E. coli core")
        print(f"  Reactions: {model.num_of_reactions()}")
        print(f"  Metabolites: {model.num_of_metabolites()}")
        print(f"  FBA optimal biomass: {max_biomass:.6f}")

        # ==============================================================
        # STEP 1: Standard sampling under 3 biomass scenarios
        # ==============================================================
        print("\n" + "-" * 50)
        print("  STEP 1: Standard flux sampling (3 scenarios)")
        print("-" * 50)

        scenarios = {
            "optimal (100%)": 100,
            "half-optimal (50%)": 50,
            "biomass-free (0%)": 0,
        }

        standard_results = {}
        for name, opt_pct in scenarios.items():
            m = MetabolicNetwork.from_json(self.MODEL_PATH)
            m.set_opt_percentage(opt_pct)
            sampler = PolytopeSampler(m)
            ss = sampler.generate_steady_states_no_multiphase(
                method="billiard_walk", n=300, burn_in=50, thinning=2
            )
            standard_results[name] = ss
            bm = ss[biomass_idx, :]
            print(f"    {name}: {ss.shape[1]} samples, "
                  f"biomass mean={np.mean(bm):.4f}")

        # ==============================================================
        # STEP 2: Loop analysis on standard samples
        # ==============================================================
        print("\n" + "-" * 50)
        print("  STEP 2: Loop analysis on standard samples")
        print("-" * 50)

        lfs = LooplessFluxSampler(model, loop_tolerance=1e-3)

        print(f"    Internal reactions: {len(lfs.internal_indices)}")
        print(f"    Exchange reactions: {len(lfs.exchange_indices)}")
        print(f"    Cycle nullspace dim: "
              f"{lfs.cycle_nullspace.shape[1] if lfs.cycle_nullspace.size > 0 else 0}")

        for name, ss in standard_results.items():
            analysis = lfs.analyze_loops(ss)
            print(f"    {name}: "
                  f"loops={analysis['fraction_with_loops']*100:.1f}%, "
                  f"mean_norm={analysis['mean_loop_norm']:.6f}")

        # ==============================================================
        # STEP 3: Loopless sampling (rejection-based)
        # ==============================================================
        print("\n" + "-" * 50)
        print("  STEP 3: Loopless sampling (rejection)")
        print("-" * 50)

        loopless_ss, rej_stats = lfs.sample_loopless(
            n_samples=200,
            max_attempts_factor=20,
            method="billiard_walk",
            opt_percentage=100,
        )

        print(f"    Accepted: {rej_stats['total_accepted']}")
        print(f"    Rejected: {rej_stats['total_rejected']}")
        print(f"    Acceptance rate: {rej_stats['acceptance_rate']:.4f}")

        # ==============================================================
        # STEP 4: Penalty-weighted sampling
        # ==============================================================
        print("\n" + "-" * 50)
        print("  STEP 4: Penalty-weighted sampling")
        print("-" * 50)

        penalized_ss, weights, pen_stats = lfs.sample_penalized(
            n_samples=300,
            penalty_weight=10.0,
            method="billiard_walk",
            opt_percentage=100,
        )

        print(f"    Effective sample size: {pen_stats['effective_sample_size']:.1f}")
        print(f"    N loopless: {pen_stats['n_effectively_loopless']}")

        # ==============================================================
        # STEP 5: Reaction clustering on standard vs loopless
        # ==============================================================
        print("\n" + "-" * 50)
        print("  STEP 5: Reaction clustering comparison")
        print("-" * 50)

        for label, ss in [("standard", standard_results["optimal (100%)"]),
                          ("loopless", loopless_ss)]:
            if ss.shape[1] < 10:
                print(f"    {label}: insufficient samples, skipping")
                continue

            corr_matrix = correlated_reactions(
                ss,
                reactions=reactions,
                pearson_cutoff=0.0,
                indicator_cutoff=0,
                cells=10,
                lower_triangle=False,
            )

            _, labels, clusters = cluster_corr_reactions(
                corr_matrix, reactions, linkage="ward", t=4.0
            )

            multi_clusters = [c for c in clusters if len(c) > 1]
            print(f"    {label}: {len(clusters)} total clusters, "
                  f"{len(multi_clusters)} non-trivial")
            for i, c in enumerate(multi_clusters[:3]):
                print(f"      Cluster {i+1}: {c[:4]}{'...' if len(c)>4 else ''} "
                      f"(n={len(c)})")

        # ==============================================================
        # STEP 6: Biological interpretation summary
        # ==============================================================
        print("\n" + "-" * 50)
        print("  STEP 6: Biological interpretation")
        print("-" * 50)

        # Compare key metabolic pathway statistics
        pathways = {
            "Glycolysis": ["PFK", "PYK", "GAPD", "PGK", "ENO", "PGM"],
            "TCA Cycle": ["CS", "ACONTa", "ACONTb", "AKGDH", "SUCOAS",
                          "FUM", "MDH", "ICDHyr"],
            "Pentose Phosphate": ["G6PDH2r", "GND", "RPI", "TKT1", "TKT2",
                                  "TALA"],
            "Respiration": ["NADH16", "CYTBD", "ATPS4r"],
        }

        standard_opt = standard_results["optimal (100%)"]

        print(f"\n    {'Pathway':<22} {'Mean±Std (Standard)':>24} "
              f"{'Mean±Std (Loopless)':>24}")
        print("    " + "-" * 72)

        for pathway_name, rxn_list in pathways.items():
            std_fluxes = []
            ll_fluxes = []
            for rxn in rxn_list:
                if rxn in reactions:
                    idx = reactions.index(rxn)
                    std_fluxes.append(np.mean(np.abs(standard_opt[idx, :])))
                    if loopless_ss.shape[1] > 0:
                        ll_fluxes.append(
                            np.mean(np.abs(loopless_ss[idx, :]))
                        )

            if std_fluxes:
                std_mean = np.mean(std_fluxes)
                std_std = np.std(std_fluxes)
                if ll_fluxes:
                    ll_mean = np.mean(ll_fluxes)
                    ll_std = np.std(ll_fluxes)
                    print(f"    {pathway_name:<22} "
                          f"{std_mean:>10.4f} ± {std_std:<10.4f} "
                          f"{ll_mean:>10.4f} ± {ll_std:<10.4f}")
                else:
                    print(f"    {pathway_name:<22} "
                          f"{std_mean:>10.4f} ± {std_std:<10.4f} "
                          f"{'N/A':>24}")

        # ==============================================================
        # STEP 7: Assessment of non-convex sampling necessity
        # ==============================================================
        print("\n" + "-" * 50)
        print("  STEP 7: Assessment — Are non-convex methods needed?")
        print("-" * 50)

        opt_analysis = lfs.analyze_loops(standard_opt)
        loop_fraction = opt_analysis["fraction_with_loops"]
        mean_norm = opt_analysis["mean_loop_norm"]

        print(f"    Fraction of standard samples with loops: "
              f"{loop_fraction:.4f}")
        print(f"    Mean loop norm: {mean_norm:.6f}")

        if loop_fraction > 0.1:
            print(f"\n    CONCLUSION: {loop_fraction*100:.1f}% of standard "
                  f"samples contain internal loops.")
            print("    Non-convex sampling IS recommended for this model.")
            print("    Loopless filtering removes thermodynamically "
                  "infeasible solutions.")
        elif loop_fraction > 0.01:
            print(f"\n    CONCLUSION: {loop_fraction*100:.1f}% of samples "
                  f"have loops — moderate impact.")
            print("    Non-convex sampling is BENEFICIAL but not critical.")
        else:
            print(f"\n    CONCLUSION: Only {loop_fraction*100:.2f}% of samples "
                  f"have loops — negligible impact.")
            print("    Standard convex sampling appears sufficient for "
                  "this model.")

        print("\n" + "=" * 70)
        print("  [OK] FULL PIPELINE COMPLETED SUCCESSFULLY")
        print("=" * 70)


if __name__ == "__main__":
    if len(sys.argv) > 1 and not sys.argv[1].startswith("-"):
        set_default_solver(sys.argv[1])
        sys.argv.pop(1)
    unittest.main(verbosity=2)
