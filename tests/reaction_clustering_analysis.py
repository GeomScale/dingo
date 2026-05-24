# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2024

# Licensed under GNU LGPL.3, see LICENCE file

# =============================================================================
# Medium Test: Reaction Clustering Analysis
# =============================================================================
# This script exploits dingo's functionalities to discover reaction clusters
# that differentiate the three biomass constraint scenarios:
#   i.   opt_percentage = 100 (optimal biomass)
#   ii.  opt_percentage = 50  (half-optimal)
#   iii. opt_percentage = 0   (biomass free)
#
# Methods used:
#   - Pearson correlation with copula-based filtering
#   - Hierarchical clustering
#   - Cross-scenario comparison of cluster composition
# =============================================================================

import unittest
import os
import sys
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo.utils import correlated_reactions, cluster_corr_reactions
from dingo.pyoptinterface_based_impl import set_default_solver


class TestReactionClustering(unittest.TestCase):
    """Reaction clustering analysis across three biomass scenarios."""

    MODEL_PATH = os.path.join(os.getcwd(), "ext_data", "e_coli_core.json")

    # ----------------------------------------------------------------
    # Helpers
    # ----------------------------------------------------------------
    def _sample_scenario(self, opt_percentage, n_samples=500):
        """Sample steady states for a given opt_percentage."""
        model = MetabolicNetwork.from_json(self.MODEL_PATH)
        model.set_opt_percentage(opt_percentage)
        sampler = PolytopeSampler(model)
        steady_states = sampler.generate_steady_states_no_multiphase(
            method="billiard_walk", n=n_samples, burn_in=100, thinning=2
        )
        return model, steady_states

    # ----------------------------------------------------------------
    # Test: Correlation matrices for each scenario
    # ----------------------------------------------------------------
    def test_correlation_matrices(self):
        """
        Compute and validate correlation matrices for all three scenarios.
        """
        print("\n" + "=" * 70)
        print("MEDIUM: Computing correlation matrices for each scenario")
        print("=" * 70)

        scenarios = {
            "optimal (100%)": 100,
            "half_optimal (50%)": 50,
            "biomass_free (0%)": 0,
        }

        for name, opt_pct in scenarios.items():
            print(f"\n  --- {name} ---")
            model, steady_states = self._sample_scenario(opt_pct)
            reactions = model.reactions

            # Compute correlation matrix without copula filtering
            corr_matrix = correlated_reactions(
                steady_states,
                reactions=reactions,
                pearson_cutoff=0.0,
                indicator_cutoff=0,
                cells=10,
                cop_coeff=0.3,
                lower_triangle=False,
            )

            # Validate the correlation matrix
            self.assertEqual(corr_matrix.shape[0], len(reactions))
            self.assertEqual(corr_matrix.shape[1], len(reactions))
            self.assertAlmostEqual(np.trace(corr_matrix), len(reactions), places=5)

            # Count highly correlated pairs (|r| > 0.9)
            upper_tri = np.triu(corr_matrix, k=1)
            n_high_corr = np.sum(np.abs(upper_tri) > 0.9)

            print(f"    Matrix shape: {corr_matrix.shape}")
            print(f"    Highly correlated pairs (|r|>0.9): {n_high_corr}")
            print(f"    Mean abs correlation: {np.mean(np.abs(upper_tri)):.4f}")

        print("\n  [OK] Correlation matrix computation PASSED for all scenarios")

    # ----------------------------------------------------------------
    # Test: Copula-based filtered correlation
    # ----------------------------------------------------------------
    def test_copula_filtered_correlation(self):
        """
        Compute correlation matrices with copula indicator filtering for
        each scenario and compare the number of truly correlated pairs.
        """
        print("\n" + "=" * 70)
        print("MEDIUM: Copula-filtered correlation analysis")
        print("=" * 70)

        scenarios = {
            "optimal (100%)": 100,
            "half_optimal (50%)": 50,
            "biomass_free (0%)": 0,
        }

        corr_results = {}

        for name, opt_pct in scenarios.items():
            print(f"\n  --- {name} ---")
            model, steady_states = self._sample_scenario(opt_pct)
            reactions = model.reactions

            # Compute correlation matrix WITH copula indicator filtering
            corr_matrix, indicator_dict = correlated_reactions(
                steady_states,
                reactions=reactions,
                pearson_cutoff=0.90,
                indicator_cutoff=5,
                cells=10,
                cop_coeff=0.3,
                lower_triangle=False,
                verbose=False,
            )

            corr_results[name] = {
                "corr_matrix": corr_matrix,
                "indicator_dict": indicator_dict,
                "reactions": reactions,
            }

            # Count positive and negative correlations
            n_positive = sum(
                1 for v in indicator_dict.values()
                if v["classification"] == "positive"
            )
            n_negative = sum(
                1 for v in indicator_dict.values()
                if v["classification"] == "negative"
            )
            n_none = sum(
                1 for v in indicator_dict.values()
                if v["classification"] == "no correlation"
            )

            print(f"    Positively correlated pairs: {n_positive}")
            print(f"    Negatively correlated pairs: {n_negative}")
            print(f"    No correlation (filtered):   {n_none}")

            self.assertEqual(corr_matrix.shape[0], len(reactions))

        print("\n  [OK] Copula-filtered correlation analysis PASSED")

    # ----------------------------------------------------------------
    # Test: Hierarchical clustering and differentiation
    # ----------------------------------------------------------------
    def test_hierarchical_clustering(self):
        """
        Perform hierarchical clustering on each scenario and compare
        which reaction clusters appear, disappear, or change across scenarios.
        """
        print("\n" + "=" * 70)
        print("MEDIUM: Hierarchical clustering & cluster differentiation")
        print("=" * 70)

        scenarios = {
            "optimal": 100,
            "half_optimal": 50,
            "biomass_free": 0,
        }

        all_clusters = {}

        for name, opt_pct in scenarios.items():
            print(f"\n  --- {name} (opt={opt_pct}%) ---")
            model, steady_states = self._sample_scenario(opt_pct, n_samples=500)
            reactions = model.reactions

            # Compute correlation matrix (unfiltered for clustering)
            corr_matrix = correlated_reactions(
                steady_states,
                reactions=reactions,
                pearson_cutoff=0.0,
                indicator_cutoff=0,
                cells=10,
                lower_triangle=False,
            )

            # Perform hierarchical clustering
            dissimilarity, labels, clusters = cluster_corr_reactions(
                corr_matrix,
                reactions,
                linkage="ward",
                t=4.0,
                correction=True,
            )

            all_clusters[name] = clusters

            print(f"    Number of clusters: {len(clusters)}")
            for i, cluster in enumerate(clusters):
                if len(cluster) > 1:
                    print(f"      Cluster {i+1} ({len(cluster)} reactions): "
                          f"{cluster[:5]}{'...' if len(cluster) > 5 else ''}")

            # Validate clustering output
            self.assertTrue(len(clusters) > 0, "Should find at least one cluster")
            total_reactions_in_clusters = sum(len(c) for c in clusters)
            self.assertEqual(total_reactions_in_clusters, len(reactions),
                             "All reactions should be assigned to clusters")

        # ---------------------------------------------------------------
        # Compare clusters across scenarios
        # ---------------------------------------------------------------
        print("\n  --- Cross-scenario cluster comparison ---")

        # Convert clusters to sets for comparison
        for name, clusters in all_clusters.items():
            cluster_sets = [frozenset(c) for c in clusters]
            all_clusters[name] = cluster_sets

        # Find clusters unique to each scenario
        for name in all_clusters:
            other_names = [n for n in all_clusters if n != name]
            unique_clusters = []
            for c in all_clusters[name]:
                if len(c) > 1:
                    is_unique = True
                    for other in other_names:
                        if c in all_clusters[other]:
                            is_unique = False
                            break
                    if is_unique:
                        unique_clusters.append(c)

            if unique_clusters:
                print(f"\n    Clusters unique to {name}:")
                for c in unique_clusters[:3]:
                    print(f"      {list(c)[:5]}{'...' if len(c) > 5 else ''} "
                          f"(size={len(c)})")
            else:
                print(f"\n    No strictly unique clusters for {name} (clusters may overlap)")

        print("\n  [OK] Hierarchical clustering & differentiation PASSED")

    # ----------------------------------------------------------------
    # Test: Key reaction flux distributions across scenarios
    # ----------------------------------------------------------------
    def test_flux_distribution_comparison(self):
        """
        Compare the flux distributions of key reactions across scenarios
        to identify biologically meaningful differences.
        """
        print("\n" + "=" * 70)
        print("MEDIUM: Flux distribution comparison for key reactions")
        print("=" * 70)

        # Key reactions in E. coli core to examine
        key_reaction_names = [
            "PFK",       # Phosphofructokinase (glycolysis)
            "PYK",       # Pyruvate kinase (glycolysis)
            "CS",        # Citrate synthase (TCA cycle)
            "AKGDH",     # Alpha-ketoglutarate dehydrogenase (TCA)
            "PPC",       # Phosphoenolpyruvate carboxylase
            "ATPM",      # ATP maintenance requirement
        ]

        scenarios = {"optimal": 100, "half_optimal": 50, "biomass_free": 0}
        scenario_results = {}

        for name, opt_pct in scenarios.items():
            model, steady_states = self._sample_scenario(opt_pct, n_samples=500)
            scenario_results[name] = {
                "model": model,
                "steady_states": steady_states,
            }

        # Compare distributions
        print(f"\n  {'Reaction':<12} {'Scenario':<18} {'Mean':>10} "
              f"{'Std':>10} {'Min':>10} {'Max':>10}")
        print("  " + "-" * 72)

        reactions = scenario_results["optimal"]["model"].reactions

        for rxn_name in key_reaction_names:
            if rxn_name in reactions:
                rxn_idx = reactions.index(rxn_name)

                for s_name, data in scenario_results.items():
                    fluxes = data["steady_states"][rxn_idx, :]
                    print(f"  {rxn_name:<12} {s_name:<18} {np.mean(fluxes):>10.4f} "
                          f"{np.std(fluxes):>10.4f} {np.min(fluxes):>10.4f} "
                          f"{np.max(fluxes):>10.4f}")

                print()  # blank line between reactions

        print("  [OK] Flux distribution comparison PASSED")


if __name__ == "__main__":
    if len(sys.argv) > 1 and not sys.argv[1].startswith("-"):
        set_default_solver(sys.argv[1])
        sys.argv.pop(1)
    unittest.main(verbosity=2)
