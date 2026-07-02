#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 14:48:35 2026

@author: shivangayathri
"""

"""
Main runner for the ECG-HI exploratory feature discovery pipeline.

Version 1 runs:
Stage 1 — Data loading and audit
Stage 2 — Full ECG feature-set definition
Stage 3 — Remove non-informative features

Run from project root:

    python pipeline/run_pipeline.py
"""

import sys
from pathlib import Path

# Ensure pipeline directory is on the Python path
PIPELINE_DIR = Path(__file__).resolve().parent
if str(PIPELINE_DIR) not in sys.path:
    sys.path.insert(0, str(PIPELINE_DIR))

import config
from utils import (
    create_output_dirs,
    setup_logging,
    validate_input_files,
    write_run_manifest,
)

from stages import (
    stage01_data_audit,
    stage02_feature_definition,
    stage03_remove_noninformative_features,
    stage04_clinical_metadata_merge,
    stage05_segment_qc_flags,
    stage06_patient_condition_medians,
    stage07_primary_paired_statistics,
    stage08_fdr_correction,
    stage09_redundancy_filtering,
    stage10_preliminary_feature_prioritization,
    stage11_hrv_valid_subgroup_analysis,
    stage12_supportive_mixed_effects_modelling,
    stage13_sensitivity_analysis,
    stage14_leave_one_patient_out_analysis,
    stage15_patient_agreement_disagreement_analysis,
    stage16_final_candidate_scorecard,
    stage17_final_workbook_export,
)


def main() -> None:
    create_output_dirs()
    setup_logging()

    validate_input_files()
    write_run_manifest(
        extra={
            "pipeline_version": "v4_stages_01_to_10",
            "stages_run": [
                "stage01_data_audit",
                "stage02_feature_definition",
                "stage03_remove_noninformative_features",
                "stage04_clinical_metadata_merge",
                "stage05_segment_qc_flags",
                "stage06_patient_condition_medians",
                "stage07_primary_paired_statistics",
                "stage08_fdr_correction",
                "stage09_redundancy_filtering",
                "stage10_preliminary_feature_prioritization",
                "stage11_hrv_valid_subgroup_analysis",
                "stage12_supportive_mixed_effects_modelling",
                "stage13_sensitivity_analysis",
                "stage14_leave_one_patient_out_analysis",
                "stage15_patient_agreement_disagreement_analysis",
                "stage16_final_candidate_scorecard",
                "stage17_final_workbook_export",
            ],
        }
    )

    stage01_data_audit()
    stage02_feature_definition()
    stage03_remove_noninformative_features()
    stage04_clinical_metadata_merge()
    stage05_segment_qc_flags()
    stage06_patient_condition_medians()
    stage07_primary_paired_statistics()
    stage08_fdr_correction()
    stage09_redundancy_filtering()
    stage10_preliminary_feature_prioritization()
    stage11_hrv_valid_subgroup_analysis()
    stage12_supportive_mixed_effects_modelling()
    stage13_sensitivity_analysis()
    stage14_leave_one_patient_out_analysis()
    stage15_patient_agreement_disagreement_analysis()
    stage16_final_candidate_scorecard()
    stage17_final_workbook_export()

    print("\nPipeline completed successfully: Stages 1–17.")
    print(f"Outputs saved in: {config.OUTPUTS_DIR}")


if __name__ == "__main__":
    main()