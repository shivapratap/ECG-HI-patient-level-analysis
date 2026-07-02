"""
Configuration file for the ECG-HI exploratory feature discovery pipeline.

All file paths, thresholds, metadata column definitions, and analysis constants
should live here so that the pipeline is reproducible and auditable.
"""

from pathlib import Path


# ---------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_ROOT / "data"
METADATA_DIR = PROJECT_ROOT / "metadata"
DOCS_DIR = PROJECT_ROOT / "docs"
OUTPUTS_DIR = PROJECT_ROOT / "outputs"


# ---------------------------------------------------------------------
# Input files
# ---------------------------------------------------------------------

DATA_FILE = DATA_DIR / "processed" / "20PatientsFinalData.xlsx"

CLINICAL_FLAGS_FILE = METADATA_DIR / "patient_clinical_flags.csv"
CLINICAL_SUMMARY_FILE = METADATA_DIR / "patient_clinical_summary.csv"
CLINICAL_REPORT_FILE = METADATA_DIR / "ECG_Clinical_Summary_Report.docx"

DATA_DICTIONARY_FILE = DOCS_DIR / "data_dictionary" / "20PatientsFinalData_data_dictionary.csv"
DATA_DICTIONARY_MD_FILE = DOCS_DIR / "data_dictionary" / "20PatientsFinalData_data_dictionary.md"

FEATURE_DOCUMENTATION_FILE = DOCS_DIR / "feature_documentation" / "ecg_lead_ii_feature_list.md"


# ---------------------------------------------------------------------
# Output folders
# ---------------------------------------------------------------------

DATA_AUDIT_DIR = OUTPUTS_DIR / "data_audit"
CLEANED_DATA_DIR = OUTPUTS_DIR / "cleaned_data"
FEATURE_DEFINITION_DIR = OUTPUTS_DIR / "feature_definition"
CLINICAL_METADATA_DIR = OUTPUTS_DIR / "clinical_metadata"
SEGMENT_QC_DIR = OUTPUTS_DIR / "segment_qc"
PATIENT_LEVEL_DIR = OUTPUTS_DIR / "patient_level"
PRIMARY_STATS_DIR = OUTPUTS_DIR / "primary_statistics"
HRV_SUBGROUP_DIR = OUTPUTS_DIR / "hrv_subgroup"
MIXED_EFFECTS_DIR = OUTPUTS_DIR / "mixed_effects"
REDUNDANCY_DIR = OUTPUTS_DIR / "redundancy"
SENSITIVITY_DIR = OUTPUTS_DIR / "sensitivity"
LOO_DIR = OUTPUTS_DIR / "leave_one_patient_out"
AGREEMENT_DIR = OUTPUTS_DIR / "agreement_disagreement"
CANDIDATE_SCORECARD_DIR = OUTPUTS_DIR / "candidate_scorecard"
MANUSCRIPT_TABLES_DIR = OUTPUTS_DIR / "manuscript_tables"
FIGURES_DIR = OUTPUTS_DIR / "figures"
FINAL_WORKBOOK_DIR = OUTPUTS_DIR / "final_workbook"

OUTPUT_DIRS = [
    DATA_AUDIT_DIR,
    CLEANED_DATA_DIR,
    FEATURE_DEFINITION_DIR,
    CLINICAL_METADATA_DIR,
    SEGMENT_QC_DIR,
    PATIENT_LEVEL_DIR,
    PRIMARY_STATS_DIR,
    HRV_SUBGROUP_DIR,
    MIXED_EFFECTS_DIR,
    REDUNDANCY_DIR,
    SENSITIVITY_DIR,
    LOO_DIR,
    AGREEMENT_DIR,
    CANDIDATE_SCORECARD_DIR,
    MANUSCRIPT_TABLES_DIR,
    FIGURES_DIR,
    FINAL_WORKBOOK_DIR,
]


# ---------------------------------------------------------------------
# Dataset structure
# ---------------------------------------------------------------------

EXPECTED_DATA_SHEET_NAME = "20PatientsFinalData"

PATIENT_ID_COL = "Patient_id"
CLASS_COL = "Class"
CONDITION_COL = "Condition"
CONDITION_STD_COL = "Condition_std"
SEGMENT_COL = "Segment_Number"
EPISODE_COL = "Episode_Label"

METADATA_COLUMNS = [
    "Sl_No",
    "Segment_Number",
    "Patient_id",
    "Class",
    "Condition",
    "Episode_Label",
    "Condition_std",
]

EXPECTED_CLASS_TO_CONDITION = {
    0: "PreHI",
    1: "HI",
}

CONDITION_LABEL_MAP = {
    "prehi": "PreHI",
    "pre-hi": "PreHI",
    "pre_hi": "PreHI",
    "pre hi": "PreHI",
    "beforehi": "PreHI",
    "before hi": "PreHI",
    "baseline": "PreHI",
    "0": "PreHI",
    0: "PreHI",

    "hi": "HI",
    "duringhi": "HI",
    "during hi": "HI",
    "during-hi": "HI",
    "during_hi": "HI",
    "hemodynamic instability": "HI",
    "haemodynamic instability": "HI",
    "1": "HI",
    1: "HI",
}

ROW_ID_COL = "segment_row_id"


# ---------------------------------------------------------------------
# Feature inclusion/exclusion thresholds
# ---------------------------------------------------------------------

MIN_PAIRED_PATIENTS_PER_FEATURE = 10
MODERATE_MISSINGNESS_THRESHOLD = 0.10
HIGH_MISSINGNESS_THRESHOLD = 0.30


# ---------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------

RANDOM_SEED = 20260702

# ---------------------------------------------------------------------
# Clinical metadata columns
# ---------------------------------------------------------------------

CLINICAL_PATIENT_ID_COL = "Patient_ID"

CLINICAL_FLAG_COLUMNS = [
    "paced_or_icd",
    "continuous_AF",
    "intermittent_or_new_AF",
    "any_AF",
    "atypical_HI",
    "major_HRV_outlier",
    "BBB",
    "antiarrhythmic_exposure",
    "mechanical_support",
    "clean_sinus_candidate",
    "exclude_from_HRV_valid_subgroup",
    "notes",
]

# Available columns used to derive confound_count.
# These are based on the actual patient_clinical_flags.csv schema.
CONFOUND_COUNT_SOURCE_COLUMNS = [
    "paced_or_icd",
    "BBB",
    "continuous_AF",
    "intermittent_or_new_AF",
    "atypical_HI",
    "major_HRV_outlier",
    "antiarrhythmic_exposure",
    "mechanical_support",
]


# ---------------------------------------------------------------------
# Segment-level QC thresholds
# ---------------------------------------------------------------------

# These are intentionally conservative. They create review/sensitivity flags,
# not automatic primary-analysis exclusions.
NEAR_ZERO_SIGNAL_ENERGY_ABS_THRESHOLD = 1e-8
LOW_VARIABILITY_STD_ABS_THRESHOLD = 1e-8
LOW_VARIABILITY_VAR_ABS_THRESHOLD = 1e-12

# Extreme-value review threshold based on robust modified z-score.
# This is for descriptive review only.
EXTREME_VALUE_MODIFIED_Z_THRESHOLD = 6.0

# Kurtosis gets a separate review flag because kurtosis itself may be a
# candidate feature and should not be automatically excluded.
EXTREME_KURTOSIS_MODIFIED_Z_THRESHOLD = 6.0

# Known segment-level artifact flags from prior waveform/artifact review.
# These will only match if Segment_Number and Condition_std align.
KNOWN_ARTIFACT_SEGMENTS = [
    {
        "Patient_id": "96879",
        "Condition_std": "HI",
        "Segment_Number": 11,
        "reason": "Prior waveform review: localized high-amplitude/spiky artifact concern",
    },
    {
        "Patient_id": "96879",
        "Condition_std": "HI",
        "Segment_Number": 29,
        "reason": "Prior waveform review: localized high-amplitude/spiky artifact concern",
    },
]

# Known patient-level caution flags. These do not remove segments.
KNOWN_PATIENT_QC_CAUTIONS = {
    "70447": "Prior review: near-zero signal-energy segment concern",
    "27245": "Prior review: early PreHI low-signal-energy concern",
    "30851": "Prior review: QRS_std = 0 segment concern",
    "906": "Prior review: HRV/R-peak outlier concern",
    "96879": "Prior review: HI segments 11 and 29 artifact concern",
}


# ---------------------------------------------------------------------
# Patient-level aggregation and statistical analysis
# ---------------------------------------------------------------------

NEAR_ZERO_PREHI_ABS_THRESHOLD = 1e-8

BOOTSTRAP_N_ITERATIONS = 5000
BOOTSTRAP_CONFIDENCE_LEVEL = 0.95

# Wilcoxon signed-rank settings
WILCOXON_ZERO_METHOD = "wilcox"
WILCOXON_ALTERNATIVE = "two-sided"

# FDR correction
FDR_METHOD = "fdr_bh"
FDR_ALPHA = 0.05

# ---------------------------------------------------------------------
# Redundancy filtering and preliminary prioritization
# ---------------------------------------------------------------------

# Correlation thresholds for redundancy review
REDUNDANCY_CORR_THRESHOLD_MODERATE = 0.85
REDUNDANCY_CORR_THRESHOLD_HIGH = 0.90

# Preliminary prioritization thresholds
PRELIM_RAW_P_THRESHOLD = 0.10
PRELIM_STRONG_RAW_P_THRESHOLD = 0.05
PRELIM_DIRECTION_CONSISTENCY_THRESHOLD = 0.65
PRELIM_STRONG_DIRECTION_CONSISTENCY_THRESHOLD = 0.70
PRELIM_EFFECT_SIZE_THRESHOLD = 0.40
PRELIM_STRONG_EFFECT_SIZE_THRESHOLD = 0.50

# Shortlist rules
MAX_PRELIMINARY_SHORTLIST_FEATURES = 15
MIN_PRELIMINARY_SHORTLIST_FEATURES = 5

# ---------------------------------------------------------------------
# Supportive mixed-effects modelling
# ---------------------------------------------------------------------

LMM_MIN_PATIENTS = 10
LMM_MIN_SEGMENTS = 50

# Features are z-scored before LMM fitting.
# Very small standard deviation means the model should be skipped.
LMM_MIN_FEATURE_STD = 1e-12

# ---------------------------------------------------------------------
# Sensitivity analysis
# ---------------------------------------------------------------------

SENSITIVITY_MIN_PAIRED_PATIENTS = 8

# Sensitivity robustness labels are descriptive, not validation labels.
SENSITIVITY_DIRECTION_CONSISTENCY_REASONABLE = 0.60
SENSITIVITY_DIRECTION_CONSISTENCY_STRONG = 0.70

# Scenarios are implemented in Stage 13 using available clinical/QC flags.
SENSITIVITY_SCENARIOS = [
    "full_cohort_baseline",
    "qc_clean_segments",
    "exclude_high_confound",
    "exclude_paced_or_bbb",
    "exclude_any_af",
    "exclude_major_hrv_outlier",
    "exclude_mechanical_support",
    "exclude_atypical_hi",
]

# ---------------------------------------------------------------------
# Leave-one-patient-out influence analysis
# ---------------------------------------------------------------------

LOO_MIN_PAIRED_PATIENTS = 8
LOO_DIRECTION_CONSISTENCY_DROP_THRESHOLD = 0.15
LOO_P_VALUE_RATIO_INFLUENCE_THRESHOLD = 5.0

# ---------------------------------------------------------------------
# Patient-level agreement/disagreement analysis
# ---------------------------------------------------------------------

AGREEMENT_MIN_FEATURES_FOR_PATIENT_SUMMARY = 1

AGREEMENT_CLINICAL_FLAGS_TO_SUMMARIZE = [
    "paced_or_icd",
    "BBB",
    "paced_or_bbb",
    "continuous_AF",
    "intermittent_or_new_AF",
    "any_af",
    "major_HRV_outlier",
    "mechanical_support",
    "atypical_HI",
    "antiarrhythmic_exposure",
    "high_confound",
    "moderate_or_high_confound",
    "known_patient_qc_caution_flag",
    "known_artifact_patient_flag",
    "hrv_invalid_flag",
    "exclude_from_HRV_valid_subgroup",
]


# ---------------------------------------------------------------------
# Final candidate scorecard
# ---------------------------------------------------------------------

# These thresholds are descriptive and cohort-specific.
# They should not be interpreted as external validation thresholds.
FINAL_SCORECARD_PRIMARY_EXPLORATORY_SCORE_THRESHOLD = 8
FINAL_SCORECARD_SECONDARY_EXPLORATORY_SCORE_THRESHOLD = 5

FINAL_SCORECARD_MIN_DIRECTION_CONSISTENCY = 0.60
FINAL_SCORECARD_STRONG_DIRECTION_CONSISTENCY = 0.70

FINAL_SCORECARD_MIN_AGREEMENT_PROPORTION = 0.60
FINAL_SCORECARD_STRONG_AGREEMENT_PROPORTION = 0.70

# A feature must have at least one form of primary patient-level support
# to be classified as a primary exploratory candidate.
FINAL_SCORECARD_REQUIRE_PRIMARY_RAW_P_OR_BOOTSTRAP_FOR_TOP_TIER = True

# ---------------------------------------------------------------------
# Final workbook export
# ---------------------------------------------------------------------

FINAL_WORKBOOK_FILE = FINAL_WORKBOOK_DIR / "HI_ECG_candidate_feature_analysis_workbook.xlsx"