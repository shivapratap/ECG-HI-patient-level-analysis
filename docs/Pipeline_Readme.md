# ECG-HI Exploratory Feature Discovery Pipeline

## Purpose

This pipeline identifies ECG-derived features that change between **Pre-HI** and **HI** periods in a 20-patient hemodynamic instability cohort.

The pipeline is designed as a **feature discovery workflow**. It does not assume the final ECG feature list in advance. Instead, it starts from the full ECG feature dataset, removes unusable features, performs patient-level paired analysis, tests robustness, examines clinical explanations for patient-level disagreement, and finally creates a candidate feature scorecard and review workbook.

The central scientific question is:

> Which ECG-derived features show consistent and interpretable PreHI-to-HI changes across patients?

The primary analysis unit is the **patient**, not the ECG segment. Segment-level ECG features are first summarized within each patient and condition before statistical testing.

---

## High-level pipeline idea

The workflow follows this logic:

1. Start with segment-level ECG features.
2. Clean and audit the dataset.
3. Define the full ECG feature set.
4. Remove non-informative features.
5. Merge patient-level clinical context.
6. Create conservative signal-quality flags.
7. Convert segment-level data into patient-level PreHI and HI medians.
8. Run paired patient-level statistics.
9. Correct for multiple testing.
10. Identify redundant feature clusters.
11. Create a preliminary feature shortlist.
12. Run supportive subgroup, mixed-effects, sensitivity, and influence analyses.
13. Study patient-level agreement and disagreement.
14. Build the final exploratory candidate scorecard.
15. Export key outputs into a final Excel workbook.

This is not a single-test script. It is a multi-stage discovery and evidence-integration pipeline.

---

## Repository files

The core pipeline is organized into four main Python files.

### `run_pipeline.py`

Main execution script.

It performs setup, validates input files, writes the run manifest, and executes Stages 1–17 in order.

Run from the project root using:

```bash
python pipeline/run_pipeline.py
```

### `config.py`

Central configuration file.

This file stores:

- project paths,
- input file locations,
- output folder locations,
- expected column names,
- condition-label mappings,
- feature inclusion thresholds,
- QC thresholds,
- FDR settings,
- sensitivity scenarios,
- leave-one-patient-out thresholds,
- final scorecard rules.

Keeping these settings in one place makes the pipeline easier to audit and reproduce.

### `utils.py`

Shared helper functions.

Examples include:

- patient ID standardization,
- condition-label standardization,
- feature-family inference,
- safe numeric conversion,
- missingness handling,
- Wilcoxon signed-rank testing,
- bootstrap confidence intervals,
- matched-pairs rank-biserial effect size,
- FDR correction,
- robust modified z-score calculation,
- CSV saving.

### `stages.py`

Main analysis file.

This file contains all 17 pipeline stages. Each stage reads outputs from previous stages, performs one logical part of the analysis, and writes structured CSV outputs.

---

## Expected project structure

The pipeline expects a project structure similar to this:

```text
project_root/
├── data/
│   └── processed/
│       └── 20PatientsFinalData.xlsx
├── metadata/
│   ├── patient_clinical_flags.csv
│   ├── patient_clinical_summary.csv
│   └── ECG_Clinical_Summary_Report.docx
├── docs/
│   ├── data_dictionary/
│   │   ├── 20PatientsFinalData_data_dictionary.csv
│   │   └── 20PatientsFinalData_data_dictionary.md
│   └── feature_documentation/
│       └── ecg_lead_ii_feature_list.md
├── pipeline/
│   ├── run_pipeline.py
│   ├── config.py
│   ├── utils.py
│   └── stages.py
└── outputs/
```

The `outputs/` folder is created automatically if it does not already exist.

---

## Required input files

The pipeline expects the following required files:

| File | Purpose |
|---|---|
| `data/processed/20PatientsFinalData.xlsx` | Main segment-level ECG feature dataset |
| `metadata/patient_clinical_flags.csv` | Patient-level clinical flags used for subgroup and disagreement analysis |
| `metadata/patient_clinical_summary.csv` | Patient-level clinical summary reference |
| `docs/data_dictionary/20PatientsFinalData_data_dictionary.csv` | Data dictionary for the main feature dataset |
| `docs/feature_documentation/ecg_lead_ii_feature_list.md` | Feature documentation |

Optional files include:

| File | Purpose |
|---|---|
| `metadata/ECG_Clinical_Summary_Report.docx` | Human-readable clinical summary report |
| `docs/data_dictionary/20PatientsFinalData_data_dictionary.md` | Markdown version of the data dictionary |

---

## Main outputs

The pipeline writes results into subfolders under `outputs/`.

Important final outputs include:

| Output | Description |
|---|---|
| `outputs/candidate_scorecard/final_candidate_scorecard.csv` | Full final evidence scorecard for candidate features |
| `outputs/candidate_scorecard/final_candidate_feature_categories.csv` | Summary of final exploratory feature categories |
| `outputs/candidate_scorecard/final_candidate_feature_summary_for_manuscript.csv` | Compact manuscript-friendly summary |
| `outputs/final_workbook/HI_ECG_candidate_feature_analysis_workbook.xlsx` | Excel workbook containing key validated outputs |
| `outputs/primary_statistics/primary_paired_feature_results_all_features_fdr.csv` | Primary paired statistics with FDR correction |
| `outputs/sensitivity/candidate_feature_robustness_matrix.csv` | Candidate robustness across sensitivity scenarios |
| `outputs/leave_one_patient_out/leave_one_patient_out_summary.csv` | Influence analysis summary |
| `outputs/agreement_disagreement/feature_agreement_summary.csv` | Feature-level patient agreement summary |
| `outputs/agreement_disagreement/patient_disagreement_summary.csv` | Patient-level disagreement summary |
| `outputs/clinical_metadata/cohort_description_table.csv` | Clinical cohort summary |

---

# Pipeline stages

## Stage 1 — Data loading and audit

### Objective

Load the main ECG feature dataset and check whether the structure is valid.

### Input

- `20PatientsFinalData.xlsx`

### What happens

- Reads the expected Excel sheet.
- Checks that required columns are present.
- Standardizes condition labels into `PreHI` and `HI`.
- Checks whether `Class` agrees with `Condition`.
- Counts segments per patient and condition.
- Checks missing values.
- Checks duplicate rows.
- Saves the cleaned segment-level dataset.

### Outputs

- `cleaned_segment_feature_dataset.csv`
- `data_integrity_summary.csv`
- `segment_counts_by_patient_condition.csv`
- `missing_values_by_column.csv`
- `class_condition_crosscheck.csv`
- `patient_condition_presence.csv`

---

## Stage 2 — Full ECG feature-set definition

### Objective

Identify all numeric ECG-derived feature columns.

### Input

- `cleaned_segment_feature_dataset.csv`

### What happens

- Removes metadata columns such as patient ID, segment number, class, condition, and episode label.
- Checks which remaining columns are numeric or convertible to numeric.
- Assigns each feature to a broad feature family.

Feature families include examples such as:

- HRV,
- ECG morphology and signal energy,
- entropy and complexity,
- entropy-profile aggregate features,
- fractal-dimension summary features,
- fractal and nonlinear dynamics,
- spectral features,
- time-domain statistics.

### Outputs

- `full_ecg_feature_list.csv`
- `feature_family_map.csv`
- `feature_family_summary.csv`
- `feature_definition_audit.csv`

---

## Stage 3 — Remove non-informative features

### Objective

Remove features that cannot support paired patient-level PreHI versus HI analysis.

### Inputs

- `cleaned_segment_feature_dataset.csv`
- `full_ecg_feature_list.csv`

### What happens

A feature is excluded if:

- it is entirely missing,
- it has only one unique non-missing value,
- it has too few paired patients,
- or it has problematic missingness.

The minimum number of paired patients is controlled in `config.py`.

### Outputs

- `informative_feature_list.csv`
- `excluded_noninformative_features.csv`
- `feature_missingness_summary.csv`
- `informative_feature_summary.csv`

---

## Stage 4 — Clinical metadata merge

### Objective

Attach patient-level clinical context to the segment-level ECG dataset.

### Inputs

- `cleaned_segment_feature_dataset.csv`
- `patient_clinical_flags.csv`
- `patient_clinical_summary.csv`

### What happens

- Cleans column names.
- Standardizes patient IDs.
- Converts clinical flags into binary 0/1 values.
- Creates derived clinical flags, including:
  - `any_af`,
  - `paced_or_bbb`,
  - `hrv_invalid_flag`,
  - `confound_count`,
  - `high_confound`,
  - `moderate_or_high_confound`,
  - `known_artifact_patient_flag`.
- Adds a stable row-level identifier named `segment_row_id`.
- Merges clinical flags onto the ECG segment dataset.

### Outputs

- `merged_segment_clinical_dataset.csv`
- `clinical_flags_standardized.csv`
- `clinical_flags_column_audit.csv`
- `clinical_merge_audit.csv`
- `patient_clinical_flag_summary.csv`
- `cohort_description_table.csv`
- `patient_clinical_summary_merged_reference.csv`

---

## Stage 5 — Segment-level QC flags

### Objective

Create conservative signal-quality flags for ECG segments.

### Inputs

- `merged_segment_clinical_dataset.csv`
- `informative_feature_list.csv`

### What happens

The stage creates flags for:

- missing informative features,
- near-zero signal energy,
- low signal variability,
- zero QRS standard deviation,
- known artifact segments,
- extreme values,
- extreme kurtosis,
- known patient-level QC cautions.

Important:

> These QC flags are not used to automatically remove segments from the primary analysis. They are mainly used for review and sensitivity analysis.

### Outputs

- `segment_qc_flags.csv`
- `segment_qc_summary_by_patient_condition.csv`
- `qc_exclusion_counts.csv`
- `qc_flagged_segments_for_review.csv`

---

## Stage 6 — Patient-condition medians and paired deltas

### Objective

Convert segment-level ECG values into patient-level PreHI versus HI comparisons.

### Inputs

- `merged_segment_clinical_dataset.csv`
- `segment_qc_flags.csv`
- `informative_feature_list.csv`

### What happens

- Merges QC flags back into the segment-level data using `segment_row_id`.
- Computes each feature's median value for every patient and condition.
- Creates full-cohort medians using all available non-missing feature values.
- Creates QC-clean medians after removing clear QC-exclusion segments.
- Calculates patient-level feature deltas:

```text
delta = HI median − PreHI median
```

- Labels each patient-level change as:
  - increase,
  - decrease,
  - no change,
  - missing.

### Outputs

- `row_id_merge_audit.csv`
- `patient_condition_medians_all_features_full.csv`
- `patient_condition_medians_all_features_qc_clean.csv`
- `patient_feature_delta_table_all_features.csv`
- `patient_level_summary_all_features.csv`

---

## Stage 7 — Primary paired statistical analysis

### Objective

Run the main patient-level paired analysis for every informative feature.

### Input

- `patient_feature_delta_table_all_features.csv`

### What happens

For each feature, the pipeline calculates:

- number of paired patients,
- median PreHI value,
- median HI value,
- median HI minus PreHI delta,
- median percent change,
- dominant direction,
- direction consistency,
- Wilcoxon signed-rank p-value,
- matched-pairs rank-biserial effect size,
- bootstrap confidence interval for the median delta.

This is the primary statistical analysis of the pipeline.

### Output

- `primary_paired_feature_results_all_features.csv`

---

## Stage 8 — FDR correction

### Objective

Correct for multiple testing.

### Input

- `primary_paired_feature_results_all_features.csv`

### What happens

- Applies global FDR correction across all informative features.
- Applies feature-family-level FDR correction.
- Adds support labels based on FDR, raw p-value, effect size, and direction consistency.

### Outputs

- `primary_paired_feature_results_all_features_fdr.csv`
- `fdr_summary_by_family.csv`

---

## Stage 9 — Redundancy filtering

### Objective

Identify features that are highly correlated and avoid selecting many redundant features.

### Inputs

- `merged_segment_clinical_dataset.csv`
- `informative_feature_list.csv`
- `primary_paired_feature_results_all_features_fdr.csv`

### What happens

- Centers feature values within each patient.
- Computes within-patient Spearman correlations.
- Identifies highly correlated feature pairs.
- Builds redundancy clusters.
- Selects one representative feature from each cluster using a transparent score.

### Outputs

- `within_subject_spearman_correlation_matrix.csv`
- `high_correlation_pairs.csv`
- `redundancy_clusters.csv`
- `cluster_representatives.csv`

---

## Stage 10 — Preliminary feature prioritization

### Objective

Create a shortlist of promising features for deeper robustness testing.

### Inputs

- `primary_paired_feature_results_all_features_fdr.csv`
- `informative_feature_list.csv`
- `redundancy_clusters.csv`
- `cluster_representatives.csv`

### What happens

Features are scored using:

- raw p-value,
- effect size,
- direction consistency,
- bootstrap support,
- FDR support,
- redundancy representative status,
- missingness.

This creates a preliminary shortlist. This is not the final candidate list.

### Outputs

- `preliminary_feature_priority_table.csv`
- `preliminary_candidate_shortlist.csv`

---

## Stage 11 — HRV-valid subgroup analysis

### Objective

Evaluate HRV features in patients where HRV interpretation is more reliable.

### Inputs

- `patient_feature_delta_table_all_features.csv`
- `patient_clinical_flag_summary.csv`
- `primary_paired_feature_results_all_features_fdr.csv`
- `feature_family_map.csv`

### What happens

- Identifies HRV features.
- Defines an HRV-valid patient subgroup.
- Excludes patients with HRV-confounding conditions such as pacing, continuous AF, or major HRV outlier status.
- Recomputes paired statistics for HRV features in this subgroup.
- Compares subgroup direction with full-cohort direction.

Important:

> This is a descriptive/supportive subgroup analysis. It is not an independent validation cohort.

### Outputs

- `hrv_valid_patient_inclusion_table.csv`
- `hrv_valid_feature_results.csv`

---

## Stage 12 — Supportive mixed-effects modelling

### Objective

Run a supportive segment-level analysis while accounting for repeated segments within each patient.

### Inputs

- `merged_segment_clinical_dataset.csv`
- `informative_feature_list.csv`
- `primary_paired_feature_results_all_features_fdr.csv`

### What happens

For each feature, the pipeline fits a mixed-effects model:

```text
feature_z ~ condition_binary + (1 | Patient_id)
```

Here:

- `feature_z` is the z-scored feature value,
- `condition_binary` is 0 for PreHI and 1 for HI,
- `Patient_id` is included as a random intercept.

This checks whether segment-level results agree with the primary patient-level direction while accounting for clustering within patients.

Important:

> This analysis is supportive only. It does not replace the patient-level paired analysis.

### Output

- `lmm_supportive_results_all_features.csv`

---

## Stage 13 — Sensitivity analysis

### Objective

Test whether shortlisted features remain stable under clinically and signal-quality motivated scenarios.

### Inputs

- `merged_segment_clinical_dataset.csv`
- `segment_qc_flags.csv`
- `preliminary_candidate_shortlist.csv`
- `primary_paired_feature_results_all_features_fdr.csv`

### What happens

The pipeline recomputes paired statistics for shortlisted features under scenarios such as:

- full cohort baseline,
- QC-clean segments only,
- excluding high-confound patients,
- excluding paced or BBB patients,
- excluding AF patients,
- excluding major HRV outlier patients,
- excluding mechanical support patients,
- excluding atypical HI patients.

It checks whether the feature direction, delta sign, direction consistency, p-value support, and bootstrap support remain stable.

### Outputs

- `sensitivity_scenario_definitions.csv`
- `sensitivity_patient_counts.csv`
- `sensitivity_analysis_by_scenario_recomputed.csv`
- `sensitivity_vs_full_cohort_summary.csv`
- `candidate_feature_robustness_matrix.csv`

---

## Stage 14 — Leave-one-patient-out influence analysis

### Objective

Check whether a feature's apparent signal is driven by one influential patient.

### Inputs

- `patient_feature_delta_table_all_features.csv`
- `preliminary_candidate_shortlist.csv`
- `primary_paired_feature_results_all_features_fdr.csv`

### What happens

For each shortlisted feature:

- calculates the baseline result using all patients,
- removes one patient at a time,
- recomputes the result,
- checks whether the dominant direction changes,
- checks whether the delta sign changes,
- checks whether direction consistency drops substantially,
- checks whether the p-value weakens strongly.

### Outputs

- `leave_one_patient_out_results.csv`
- `leave_one_patient_out_summary.csv`

---

## Stage 15 — Patient-level agreement/disagreement analysis

### Objective

Understand which patients agree or disagree with each feature's cohort-level direction.

### Inputs

- `patient_feature_delta_table_all_features.csv`
- `preliminary_candidate_shortlist.csv`
- `primary_paired_feature_results_all_features_fdr.csv`
- `patient_clinical_flag_summary.csv`
- `candidate_feature_robustness_matrix.csv`
- `leave_one_patient_out_summary.csv`

### What happens

For each shortlisted feature and each patient, the pipeline determines whether the patient's observed PreHI-to-HI change:

- agrees with the cohort direction,
- disagrees with the cohort direction,
- shows neutral/no change,
- or is missing.

It then adds clinical context, such as:

- pacing or ICD,
- BBB,
- AF,
- mechanical support,
- atypical HI,
- major HRV outlier,
- signal/QC caution.

This helps explain why some patients may not follow the dominant cohort-level feature direction.

### Outputs

- `patient_feature_agreement_table.csv`
- `feature_agreement_summary.csv`
- `patient_disagreement_summary.csv`
- `clinical_flag_disagreement_summary.csv`
- `feature_patient_disagreement_matrix.csv`

Important:

> Clinical-flag disagreement summaries are descriptive only. Because the cohort is small, these should not be interpreted as formal association tests.

---

## Stage 16 — Final candidate scorecard

### Objective

Combine all evidence into a final exploratory candidate feature scorecard.

### Inputs

- `preliminary_candidate_shortlist.csv`
- `primary_paired_feature_results_all_features_fdr.csv`
- `cluster_representatives.csv`
- `hrv_valid_feature_results.csv`
- `lmm_supportive_results_all_features.csv`
- `candidate_feature_robustness_matrix.csv`
- `leave_one_patient_out_summary.csv`
- `feature_agreement_summary.csv`
- `clinical_flag_disagreement_summary.csv`

### What happens

The final scorecard combines evidence from:

- primary patient-level paired statistics,
- global and feature-family FDR correction,
- redundancy filtering,
- HRV-valid subgroup analysis,
- supportive mixed-effects modelling,
- sensitivity analysis,
- leave-one-patient-out influence analysis,
- patient-level agreement/disagreement analysis,
- clinical context for disagreement.

Each shortlisted feature receives:

- an evidence score,
- a final candidate category,
- a human-readable interpretation note.

Final categories include:

| Category | Meaning |
|---|---|
| `primary_exploratory_candidate` | Stronger exploratory support within this cohort |
| `secondary_exploratory_candidate` | Moderate exploratory support |
| `secondary_rhythm_sensitive_candidate` | Candidate feature whose interpretation is affected by rhythm/HRV-valid subgroup results |
| `exploratory_or_downgraded` | Weaker, mixed, or less interpretable evidence |

### Outputs

- `final_candidate_scorecard.csv`
- `final_candidate_feature_categories.csv`
- `final_candidate_feature_summary_for_manuscript.csv`

---

## Stage 17 — Final workbook export

### Objective

Package important outputs into one Excel workbook for review and manuscript preparation.

### Inputs

Validated CSV outputs from previous stages.

### What happens

Creates an Excel workbook with sheets including:

- README,
- workbook index,
- final scorecard,
- final categories,
- manuscript summary,
- preliminary shortlist,
- primary statistics,
- redundancy clusters,
- cluster representatives,
- HRV subgroup results,
- supportive LMM results,
- sensitivity robustness,
- leave-one-patient-out summary,
- feature agreement,
- patient disagreement,
- clinical disagreement,
- QC counts,
- cohort description.

Important:

> Stage 17 does not create new statistical results. It only packages already generated outputs.

### Outputs

- `HI_ECG_candidate_feature_analysis_workbook.xlsx`
- `final_workbook_export_manifest.csv`

---

# Scientific interpretation

The pipeline is designed to answer several connected questions:

1. Which ECG features are usable in this dataset?
2. Which features change from PreHI to HI at the patient level?
3. Are those changes consistent across patients?
4. Are the changes statistically supported?
5. Are candidate features redundant with other features?
6. Do findings survive signal-quality and clinical sensitivity checks?
7. Are findings driven by one patient?
8. Which patients disagree with cohort-level trends?
9. Can disagreement be explained by clinical context such as AF, BBB, pacing, mechanical support, atypical HI, or artifact concerns?
10. Which features are most reasonable to report as exploratory ECG-HI candidate features?

---

# Important interpretation cautions

This is a 20-patient exploratory cohort. The results should be interpreted cautiously.

Appropriate wording:

> These ECG-derived features are exploratory candidate markers associated with PreHI-to-HI change in this cohort.

Avoid wording such as:

> These are validated diagnostic biomarkers.

The final scorecard categories are evidence summaries from this cohort. They are not external validation labels.

---

# Recommended files for manuscript preparation

The most useful outputs for manuscript tables and Results writing are:

| File | Use |
|---|---|
| `final_candidate_feature_summary_for_manuscript.csv` | Main compact feature summary table |
| `final_candidate_scorecard.csv` | Full evidence table |
| `primary_paired_feature_results_all_features_fdr.csv` | Primary statistical results |
| `candidate_feature_robustness_matrix.csv` | Sensitivity robustness table |
| `leave_one_patient_out_summary.csv` | Patient influence analysis |
| `feature_agreement_summary.csv` | Feature-level agreement/disagreement |
| `patient_disagreement_summary.csv` | Patient-level disagreement patterns |
| `clinical_flag_disagreement_summary.csv` | Clinical context for disagreement |
| `cohort_description_table.csv` | Clinical cohort summary |
| `HI_ECG_candidate_feature_analysis_workbook.xlsx` | Integrated review workbook |

---

# Notes for future development

Possible future improvements include:

1. Add manuscript-ready tables and figures as a separate notebook or pipeline module.
2. Add automated plots for final candidate features.
3. Add clearer clinical interpretation summaries for each candidate feature.
4. Add external validation when a larger independent cohort becomes available.
5. Improve feature-family mapping using the full data dictionary instead of only rule-based naming.
6. Add automated tests for input schema, expected output files, and key stage dependencies.

---

# Summary

This pipeline is a structured ECG feature discovery workflow for hemodynamic instability analysis.

It begins with raw segment-level ECG feature data and ends with a clinically interpretable exploratory candidate feature scorecard. The workflow emphasizes patient-level paired analysis, robustness testing, redundancy reduction, clinical context, and cautious interpretation.
