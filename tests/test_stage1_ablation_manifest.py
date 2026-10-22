"""Standard-library validation of the Stage 1 design manifest; runs no fits."""

import itertools
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "config" / "stage1_ablation.json"


class Stage1AblationManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = MANIFEST.read_text(encoding="utf-8")
        cls.cfg = json.loads(cls.raw)

    def test_scope_and_readiness_are_fail_closed(self):
        self.assertEqual(self.cfg["status"], "DESIGN_AND_SMOKE_ONLY_NOT_RUN")
        self.assertEqual(
            self.cfg["scientific_readiness"],
            {
                "calibration": "NOT_RUN",
                "simulation_study": "NOT_RUN",
                "figures": "NOT_READY",
                "manuscript": "NOT_READY",
                "software_promotion": "NOT_READY",
            },
        )
        prohibited = set(self.cfg["scope"]["prohibited"])
        self.assertIn("scientific_fits", prohibited)
        self.assertIn("calibration", prohibited)
        self.assertFalse(self.cfg["scope"]["inherited_source_bundle_included"])

    def test_planned_internal_helper_contract_is_exact(self):
        helpers = {item["name"]: item for item in self.cfg["internal_helpers"]}
        self.assertEqual(
            helpers[".ablation_prepare"]["planned_signature"],
            ".ablation_prepare(blocks, initial_ranks, seed)",
        )
        self.assertEqual(
            helpers[".ablation_geometry"]["planned_signature"],
            ".ablation_geometry(cache, joint_rank, representation='raw'/'orthonormal', aggregation='robust'/'classical', seed)",
        )
        self.assertEqual(
            helpers[".ablation_finish"]["planned_signature"],
            ".ablation_finish(cache, geometry, individual='original'/'constrained', seed)",
        )
        self.assertTrue(all(x["visibility"] == "internal_only" for x in helpers.values()))

    def test_arms_are_the_complete_nonconfounded_factorial(self):
        arms = self.cfg["arms"]
        observed = {(a["R"], a["S"], a["I"]) for a in arms}
        expected = set(itertools.product((0, 1), repeat=3))
        self.assertEqual(len(arms), 8)
        self.assertEqual(observed, expected)
        self.assertEqual(len({a["id"] for a in arms}), 8)
        baselines = [a for a in arms if a["baseline"]]
        self.assertEqual(baselines, [{"id": "R0S0I0", "R": 0, "S": 0, "I": 0, "baseline": True}])
        for arm in arms:
            self.assertEqual(arm["id"], f"R{arm['R']}S{arm['S']}I{arm['I']}")

    def test_representation_covers_block_and_final_joint_scores(self):
        rep = self.cfg["factor_levels"]["R"]
        self.assertEqual(rep["applies_to"], ["block_scores", "final_joint_scores"])
        levels = {x["label"]: x for x in rep["levels"]}
        self.assertIn("no hidden QR", levels["raw"]["block_score_handling"])
        self.assertIn("raw amplitudes", levels["raw"]["final_joint_score_handling"])
        self.assertIn("rank-revealing", levels["orthonormal"]["block_score_handling"])
        self.assertIn("rank-revealing", levels["orthonormal"]["final_joint_score_handling"])

    def test_stack_factor_reuses_cached_robust_block_fits(self):
        stack = self.cfg["factor_levels"]["S"]
        self.assertIn("cached robust block fits", stack["shared_input"])
        self.assertFalse(stack["classical_block_refit"])
        self.assertEqual(stack["projector_sum_valid_only_when"], "R1")
        self.assertFalse(self.cfg["scope"]["classical_block_fit_in_primary_factorial"])

    def test_individual_pair_freezes_selection_and_logs_rank_drop(self):
        indiv = self.cfg["factor_levels"]["I"]
        self.assertEqual(indiv["selected_robust_fit_count_per_dataset_block"], 1)
        self.assertEqual(indiv["rank_decision"], "freeze_before_constraint")
        self.assertFalse(indiv["post_constraint_robust_refit"])
        self.assertEqual(
            indiv["required_rank_outputs"],
            [
                "selected_rank_pre_constraint",
                "achieved_rank_post_constraint",
                "rank_drop",
            ],
        )
        pipe = self.cfg["pipeline"]
        self.assertTrue(pipe["geometry_immutable_during_finish"])
        self.assertTrue(pipe["i_pair_shares_preselected_reconstruction"])

    def test_final_canonicalization_precedes_common_l2_and_has_pre_post_outputs(self):
        pipe = self.cfg["pipeline"]
        stages = pipe["ordered_stages"]
        self.assertLess(
            stages.index("canonicalize_final_joint_scores"),
            stages.index("apply_common_l2_identifiability_filter"),
        )
        self.assertTrue(pipe["canonicalization_before_l2_filter"])
        self.assertTrue(pipe["common_l2_filter_across_arms"])
        outputs = set(pipe["required_outputs"])
        required = {
            "canonical_joint_scores_pre_l2",
            "joint_rank_pre_l2",
            "joint_gram_error_pre_l2",
            "joint_scores_post_l2",
            "joint_rank_post_l2",
            "joint_gram_error_post_l2",
            "l2_retained_direction_ids",
            "l2_dropped_direction_ids",
        }
        self.assertTrue(required.issubset(outputs))

    def test_auto_rank_and_wedin_are_separately_blocked(self):
        rank = self.cfg["rank_selection"]
        self.assertEqual(
            rank["fixed_rank_layer"]["status"], "AUTHORIZED_FOR_DESIGN_AND_SMOKE"
        )
        auto = rank["automatic_rank_layer"]
        self.assertEqual(auto["status"], "BLOCKED_PENDING_MATCHED_NULL")
        self.assertTrue(auto["independent_permutations_across_blocks"])
        self.assertTrue(auto["underlying_permutations_shared_across_arms"])
        self.assertTrue(auto["isolated_stage_seeds"])
        self.assertEqual(
            auto["no_effective_finite_threshold_status"],
            "BLOCKED_NO_EFFECTIVE_FINITE_THRESHOLD",
        )
        self.assertFalse(auto["negative_infinity_is_success"])
        wedin = rank["wedin_alternative"]
        self.assertEqual(wedin["status"], "BLOCKED_PENDING_DERIVATION_REVIEW")
        self.assertTrue(wedin["separate_from_permutation_candidate"])
        self.assertFalse(wedin["may_be_combined_as_fallback_max_or_min"])
        self.assertFalse(wedin["unverified_citations_asserted_here"])

    def test_scenarios_are_exact_centered_and_certified(self):
        scenarios = self.cfg["scenario_manifest"]
        generator = scenarios["generator"]
        self.assertEqual(generator["n"], 60)
        self.assertEqual(generator["block_features"], [30, 25, 20])
        self.assertIn("H_n", generator["center_operator"])
        self.assertFalse(generator["recenter_after_contamination"])
        self.assertIn("never silently resample", generator["failed_draw_policy"])
        self.assertTrue(scenarios["truth_certificate"]["required_post_transform"])
        signals = {x["id"]: x for x in scenarios["signals"]}
        self.assertEqual(
            set(signals),
            {"joint_rank0", "joint_rank1", "joint_rank2", "joint_rank2_weak_last"},
        )
        self.assertEqual([signals[x]["joint_rank"] for x in signals], [0, 1, 2, 2])
        self.assertEqual(signals["joint_rank2"]["joint_amplitudes"], [2.0, 1.5])
        self.assertEqual(
            signals["joint_rank2_weak_last"]["joint_amplitudes"], [2.0, 0.35]
        )
        self.assertTrue(all(x["individual_ranks"] == [1, 1, 1] for x in signals.values()))

    def test_contamination_definitions_and_pairing_are_explicit(self):
        conditions = {
            x["id"]: x
            for x in self.cfg["scenario_manifest"]["contamination_conditions"]
        }
        cell = conditions["cell_5pct_plus10_robust_scale"]
        self.assertEqual(cell["fraction_per_block"], 0.05)
        self.assertEqual(cell["cell_counts_by_block"], [90, 75, 60])
        self.assertEqual(cell["addition_in_robust_scale_units"], 10.0)
        self.assertEqual(cell["addition_sign"], 1)
        self.assertTrue(cell["mask_fixed_across_all_arms"])
        row = conditions["row_10pct_plus10_robust_scale"]
        self.assertEqual(row["row_fraction"], 0.10)
        self.assertEqual(row["row_count"], 6)
        self.assertEqual(row["status"], "OPTIONAL_DEVELOPMENT_ONLY_NOT_RUN")
        self.assertIsNone(row["confirmation_pool"])
        self.assertTrue(
            self.cfg["scenario_manifest"]["core_cross_product"][
                "all_eight_arms_per_dataset"
            ]
        )

    def test_seed_lists_are_finite_explicit_disjoint_and_untouched(self):
        seeds = self.cfg["seed_manifest"]
        self.assertEqual(seeds["smoke"]["seeds"], [1729])
        self.assertEqual(seeds["development"]["seeds"], list(range(1001, 1021)))
        confirmation = seeds["confirmation"]
        self.assertEqual(confirmation["null"], list(range(20001, 20201)))
        self.assertEqual(confirmation["rank_one"], list(range(21001, 21101)))
        self.assertEqual(confirmation["rank_two"], list(range(22001, 22101)))
        self.assertEqual(confirmation["status"], "UNTOUCHED_LOCKED_NOT_RUN")
        self.assertFalse(confirmation["available_for_development"])
        pools = [
            seeds["smoke"]["seeds"],
            seeds["development"]["seeds"],
            confirmation["null"],
            confirmation["rank_one"],
            confirmation["rank_two"],
        ]
        flattened = [seed for pool in pools for seed in pool]
        self.assertEqual(len(flattened), 421)
        self.assertEqual(len(flattened), len(set(flattened)))
        self.assertTrue(all(type(seed) is int for seed in flattened))

    def test_historical_endpoint_and_future_gates_are_unchanged(self):
        primary = self.cfg["endpoints"]["historical_primary_rank2_projector_loss"]
        self.assertEqual(primary["formula"], "norm(P_est - P_true, 'F') / sqrt(4)")
        self.assertEqual(primary["truth_rank"], 2)
        self.assertEqual(primary["candidate_degradation_margin"], 0.02)
        self.assertEqual(primary["margin_scope"], "this exact endpoint and normalization only")
        gates = self.cfg["future_gates"]
        self.assertEqual(
            gates["geometry"],
            {
                "block_gram_error_max": 1e-10,
                "final_gram_error_max": 1e-10,
                "constrained_individual_fitted_joint_leakage_max": 1e-8,
            },
        )
        self.assertEqual(gates["rank_null"], {"datasets": 200, "type_I_max": 0.10})
        self.assertEqual(
            gates["rank_one"], {"datasets": 100, "exact_recovery_min": 0.80}
        )
        self.assertEqual(
            gates["rank_two"],
            {
                "datasets": 100,
                "exact_recovery_min": 0.80,
                "over_selection_max": 0.10,
            },
        )
        self.assertEqual(
            gates["rank_gate_status"], "BLOCKED_WHILE_AUTOMATIC_RANK_IS_BLOCKED"
        )
        self.assertIn("not evidence of superiority", gates["interpretation"])

    def test_paired_analysis_multiplicity_and_failures_are_prespecified(self):
        analysis = self.cfg["paired_analysis"]
        self.assertEqual(analysis["contrasts"], ["R", "S", "I", "RS", "RI", "SI", "RSI"])
        self.assertEqual(analysis["uncertainty"]["unit"], "dataset_with_all_eight_arms_kept_together")
        self.assertEqual(analysis["uncertainty"]["resamples"], 10000)
        self.assertEqual(analysis["multiplicity"]["method"], "Holm")
        failures = analysis["failure_handling"]
        self.assertTrue(failures["all_attempts_in_denominator"])
        self.assertTrue(failures["report_complete_octet_count"])
        self.assertTrue(failures["report_missing_arm_patterns"])
        self.assertTrue(failures["analyze_failure_indicator"])
        self.assertFalse(failures["silent_complete_case_deletion"])

    def test_smoke_target_is_only_planned_and_no_sources_are_bundled(self):
        smoke = self.cfg["planned_orchestrator"]
        self.assertEqual(smoke["target"], "scripts/run_stage1_smoke.R")
        self.assertEqual(smoke["status"], "PLANNED_NOT_IMPLEMENTED_NOT_RUN")
        self.assertEqual(smoke["required_arm_count"], 8)
        boundary = self.cfg["provenance_boundary"]
        self.assertFalse(boundary["historical_manifest_reconstructed"])
        self.assertEqual(boundary["inherited_source_files"], [])
        self.assertEqual(boundary["unverified_external_citations"], [])
        self.assertIn("Wedin derivation", boundary["derivation_cross_reference"])


if __name__ == "__main__":
    unittest.main()
