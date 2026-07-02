"""
Pipeline stages for the ECG-HI exploratory feature discovery pipeline.

Version 1 implements:
Stage 1 — Data loading and audit
Stage 2 — Full ECG feature-set definition
Stage 3 — Remove non-informative features
"""

import logging

import numpy as np
import pandas as pd

import config
import statsmodels.formula.api as smf

from pathlib import Path


from utils import (
    standardize_condition_label,
    standardize_patient_id_series,
    safe_to_numeric,
    count_paired_patients,
    infer_feature_family,
    save_dataframe,
    to_binary_flag,
    robust_modified_zscore,
    clean_column_names,
    bootstrap_median_ci,
    wilcoxon_signed_rank_pvalue,
    matched_pairs_rank_biserial_effect_size,
    apply_fdr_correction,
    direction_from_delta,
    dominant_direction_from_counts,
    get_connected_components_from_edges,
    min_max_normalize,
    safe_negative_log10_pvalue,
    zscore_series,
)

# ---------------------------------------------------------------------
# Stage 1 — Data loading and audit
# ---------------------------------------------------------------------

def stage01_data_audit() -> pd.DataFrame:
    """
    Load main ECG feature dataset and perform structural audit.

    Outputs:
    - data_integrity_summary.csv
    - segment_counts_by_patient_condition.csv
    - missing_values_by_column.csv
    - class_condition_crosscheck.csv
    - duplicate_rows.csv, if any
    - cleaned_segment_feature_dataset.csv

    Important:
    Temporary audit columns are not saved into the cleaned dataset.
    """

    logging.info("Stage 1 started: Data loading and audit.")

    df = pd.read_excel(
        config.DATA_FILE,
        sheet_name=config.EXPECTED_DATA_SHEET_NAME,
    )

    if config.PATIENT_ID_COL not in df.columns:
        raise ValueError(f"Required column missing: {config.PATIENT_ID_COL}")

    if config.CONDITION_COL not in df.columns:
        raise ValueError(f"Required column missing: {config.CONDITION_COL}")

    # Standardize Patient_id
    df[config.PATIENT_ID_COL] = df[config.PATIENT_ID_COL].astype(str).str.strip()

    # Standardize condition
    df[config.CONDITION_STD_COL] = df[config.CONDITION_COL].apply(
        standardize_condition_label
    )

    unmapped_conditions = (
        df.loc[df[config.CONDITION_STD_COL].isna(), config.CONDITION_COL]
        .dropna()
        .unique()
        .tolist()
    )

    if unmapped_conditions:
        logging.warning(f"Unmapped condition labels found: {unmapped_conditions}")

    # Class-condition crosscheck
    if config.CLASS_COL in df.columns:
        expected_class_condition = df[config.CLASS_COL].map(
            config.EXPECTED_CLASS_TO_CONDITION
        )

        class_condition_match = expected_class_condition == df[config.CONDITION_STD_COL]

        class_condition_crosscheck_extended = (
            pd.DataFrame(
                {
                    config.CLASS_COL: df[config.CLASS_COL],
                    config.CONDITION_COL: df[config.CONDITION_COL],
                    config.CONDITION_STD_COL: df[config.CONDITION_STD_COL],
                    "Class_condition_expected": expected_class_condition,
                    "Class_condition_match": class_condition_match,
                }
            )
            .drop_duplicates()
            .sort_values([config.CLASS_COL, config.CONDITION_STD_COL])
        )

        save_dataframe(
            class_condition_crosscheck_extended,
            config.DATA_AUDIT_DIR / "class_condition_crosscheck.csv",
        )

        class_condition_mismatch_count = int((~class_condition_match).sum())

    else:
        class_condition_mismatch_count = np.nan
        logging.warning("Class column not found; class-condition check skipped.")

    # Duplicate rows
    duplicate_mask = df.duplicated()
    duplicate_count = int(duplicate_mask.sum())

    if duplicate_count > 0:
        save_dataframe(
            df.loc[duplicate_mask].copy(),
            config.DATA_AUDIT_DIR / "duplicate_rows.csv",
        )

    # Segment counts
    group_cols = [
        config.PATIENT_ID_COL,
        config.CONDITION_STD_COL,
    ]

    if config.CLASS_COL in df.columns:
        group_cols.append(config.CLASS_COL)

    if config.EPISODE_COL in df.columns:
        group_cols.append(config.EPISODE_COL)

    segment_counts = (
        df.groupby(group_cols, dropna=False)
        .size()
        .reset_index(name="n_segments")
        .sort_values([config.PATIENT_ID_COL, config.CONDITION_STD_COL])
    )

    save_dataframe(
        segment_counts,
        config.DATA_AUDIT_DIR / "segment_counts_by_patient_condition.csv",
    )

    # Patients with condition availability
    patient_condition_presence = (
        df.groupby([config.PATIENT_ID_COL, config.CONDITION_STD_COL])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )

    for col in ["PreHI", "HI"]:
        if col not in patient_condition_presence.columns:
            patient_condition_presence[col] = 0

    patient_condition_presence["has_PreHI"] = patient_condition_presence["PreHI"] > 0
    patient_condition_presence["has_HI"] = patient_condition_presence["HI"] > 0
    patient_condition_presence["has_both_conditions"] = (
        patient_condition_presence["has_PreHI"]
        & patient_condition_presence["has_HI"]
    )

    save_dataframe(
        patient_condition_presence,
        config.DATA_AUDIT_DIR / "patient_condition_presence.csv",
    )

    # Missing values by column
    missing_values = pd.DataFrame(
        {
            "column": df.columns,
            "missing_count": df.isna().sum().values,
            "missing_percent": (df.isna().sum().values / len(df)) * 100,
            "dtype": [str(df[col].dtype) for col in df.columns],
        }
    ).sort_values("missing_percent", ascending=False)

    save_dataframe(
        missing_values,
        config.DATA_AUDIT_DIR / "missing_values_by_column.csv",
    )

    # Candidate feature count by metadata exclusion
    metadata_present = [c for c in config.METADATA_COLUMNS if c in df.columns]
    non_metadata_cols = [c for c in df.columns if c not in metadata_present]
    numeric_non_metadata_cols = [
        c for c in non_metadata_cols
        if pd.api.types.is_numeric_dtype(df[c])
    ]

    # Data integrity summary
    data_integrity_summary = pd.DataFrame(
        [
            {
                "metric": "n_rows",
                "value": len(df),
            },
            {
                "metric": "n_columns",
                "value": df.shape[1],
            },
            {
                "metric": "n_unique_patients",
                "value": df[config.PATIENT_ID_COL].nunique(),
            },
            {
                "metric": "metadata_columns_present",
                "value": "; ".join(metadata_present),
            },
            {
                "metric": "n_non_metadata_columns",
                "value": len(non_metadata_cols),
            },
            {
                "metric": "n_numeric_non_metadata_columns",
                "value": len(numeric_non_metadata_cols),
            },
            {
                "metric": "condition_std_values",
                "value": "; ".join(
                    sorted(
                        [
                            str(x)
                            for x in df[config.CONDITION_STD_COL].dropna().unique()
                        ]
                    )
                ),
            },
            {
                "metric": "n_unmapped_condition_rows",
                "value": int(df[config.CONDITION_STD_COL].isna().sum()),
            },
            {
                "metric": "n_duplicate_rows",
                "value": duplicate_count,
            },
            {
                "metric": "n_patients_with_both_conditions",
                "value": int(patient_condition_presence["has_both_conditions"].sum()),
            },
            {
                "metric": "class_condition_mismatch_count",
                "value": class_condition_mismatch_count,
            },
        ]
    )

    save_dataframe(
        data_integrity_summary,
        config.DATA_AUDIT_DIR / "data_integrity_summary.csv",
    )

    # Save cleaned dataset without temporary audit/helper columns
    save_dataframe(
        df,
        config.CLEANED_DATA_DIR / "cleaned_segment_feature_dataset.csv",
    )

    logging.info("Stage 1 completed.")
    return df


# ---------------------------------------------------------------------
# Stage 2 — Full ECG feature-set definition
# ---------------------------------------------------------------------

def stage02_feature_definition() -> pd.DataFrame:
    """
    Define the full ECG feature set from non-metadata numeric columns.

    Outputs:
    - full_ecg_feature_list.csv
    - feature_family_map.csv
    - feature_family_summary.csv
    """

    logging.info("Stage 2 started: Full ECG feature-set definition.")

    df = pd.read_csv(config.CLEANED_DATA_DIR / "cleaned_segment_feature_dataset.csv")

    metadata_present = [c for c in config.METADATA_COLUMNS if c in df.columns]

    candidate_columns = [c for c in df.columns if c not in metadata_present]

    feature_rows = []

    for col in candidate_columns:
        numeric_version = safe_to_numeric(df[col])
        is_numeric = numeric_version.notna().any()

        feature_rows.append(
            {
                "feature": col,
                "is_numeric_or_convertible": bool(is_numeric),
                "dtype_original": str(df[col].dtype),
                "feature_family": infer_feature_family(col),
                "included_in_full_ecg_feature_set": bool(is_numeric),
                "reason": "numeric ECG-derived feature"
                if is_numeric
                else "excluded: non-numeric/non-convertible",
            }
        )

    full_feature_df = pd.DataFrame(feature_rows)

    full_ecg_feature_list = full_feature_df[
        full_feature_df["included_in_full_ecg_feature_set"] == True
    ].copy()

    save_dataframe(
        full_ecg_feature_list,
        config.FEATURE_DEFINITION_DIR / "full_ecg_feature_list.csv",
    )

    feature_family_map = full_ecg_feature_list[
        ["feature", "feature_family"]
    ].copy()

    save_dataframe(
        feature_family_map,
        config.FEATURE_DEFINITION_DIR / "feature_family_map.csv",
    )

    feature_family_summary = (
        feature_family_map.groupby("feature_family")
        .size()
        .reset_index(name="n_features")
        .sort_values("n_features", ascending=False)
    )

    save_dataframe(
        feature_family_summary,
        config.FEATURE_DEFINITION_DIR / "feature_family_summary.csv",
    )

    # Save all candidate column audit too
    save_dataframe(
        full_feature_df,
        config.FEATURE_DEFINITION_DIR / "feature_definition_audit.csv",
    )

    logging.info(
        f"Stage 2 completed. Full ECG feature count: {len(full_ecg_feature_list)}"
    )

    return full_ecg_feature_list


# ---------------------------------------------------------------------
# Stage 3 — Remove non-informative features
# ---------------------------------------------------------------------

def stage03_remove_noninformative_features() -> pd.DataFrame:
    """
    Remove ECG features that are not suitable for paired patient-level analysis.

    Exclusion criteria:
    - non-numeric after coercion
    - entirely missing
    - constant or only one unique non-missing value
    - fewer than MIN_PAIRED_PATIENTS_PER_FEATURE paired patients
    """

    logging.info("Stage 3 started: Remove non-informative features.")

    df = pd.read_csv(config.CLEANED_DATA_DIR / "cleaned_segment_feature_dataset.csv")

    full_features = pd.read_csv(
        config.FEATURE_DEFINITION_DIR / "full_ecg_feature_list.csv"
    )

    feature_names = full_features["feature"].tolist()

    informative_rows = []
    excluded_rows = []
    missingness_rows = []

    for feature in feature_names:
        numeric_series = safe_to_numeric(df[feature])

        missing_count = int(numeric_series.isna().sum())
        missing_percent = float(missing_count / len(df))
        non_missing_count = int(numeric_series.notna().sum())
        unique_non_missing = int(numeric_series.dropna().nunique())

        paired_patient_count = count_paired_patients(
            df.assign(**{feature: numeric_series}),
            feature,
        )

        exclusion_reasons = []

        if non_missing_count == 0:
            exclusion_reasons.append("entirely_missing")

        if unique_non_missing <= 1:
            exclusion_reasons.append("constant_or_single_unique_value")

        if paired_patient_count < config.MIN_PAIRED_PATIENTS_PER_FEATURE:
            exclusion_reasons.append(
                f"insufficient_paired_patients_less_than_{config.MIN_PAIRED_PATIENTS_PER_FEATURE}"
            )

        if missing_percent >= config.HIGH_MISSINGNESS_THRESHOLD:
            missingness_flag = "high_missingness"
        elif missing_percent >= config.MODERATE_MISSINGNESS_THRESHOLD:
            missingness_flag = "moderate_missingness"
        else:
            missingness_flag = "low_or_no_missingness"

        family = infer_feature_family(feature)

        missingness_rows.append(
            {
                "feature": feature,
                "feature_family": family,
                "missing_count": missing_count,
                "missing_percent": missing_percent,
                "non_missing_count": non_missing_count,
                "unique_non_missing_values": unique_non_missing,
                "paired_patient_count": paired_patient_count,
                "missingness_flag": missingness_flag,
            }
        )

        base_row = {
            "feature": feature,
            "feature_family": family,
            "missing_count": missing_count,
            "missing_percent": missing_percent,
            "non_missing_count": non_missing_count,
            "unique_non_missing_values": unique_non_missing,
            "paired_patient_count": paired_patient_count,
            "missingness_flag": missingness_flag,
        }

        if exclusion_reasons:
            excluded_rows.append(
                {
                    **base_row,
                    "included": False,
                    "exclusion_reason": "; ".join(exclusion_reasons),
                }
            )
        else:
            informative_rows.append(
                {
                    **base_row,
                    "included": True,
                    "exclusion_reason": "",
                }
            )

    informative_df = pd.DataFrame(informative_rows).sort_values(
        ["feature_family", "feature"]
    )

    excluded_df = pd.DataFrame(excluded_rows)

    if not excluded_df.empty:
        excluded_df = excluded_df.sort_values(["exclusion_reason", "feature"])

    missingness_df = pd.DataFrame(missingness_rows).sort_values(
        "missing_percent",
        ascending=False,
    )

    save_dataframe(
        informative_df,
        config.FEATURE_DEFINITION_DIR / "informative_feature_list.csv",
    )

    save_dataframe(
        excluded_df,
        config.FEATURE_DEFINITION_DIR / "excluded_noninformative_features.csv",
    )

    save_dataframe(
        missingness_df,
        config.FEATURE_DEFINITION_DIR / "feature_missingness_summary.csv",
    )

    summary = pd.DataFrame(
        [
            {
                "metric": "n_full_ecg_features",
                "value": len(feature_names),
            },
            {
                "metric": "n_informative_features",
                "value": len(informative_df),
            },
            {
                "metric": "n_excluded_features",
                "value": len(excluded_df),
            },
            {
                "metric": "min_paired_patients_required",
                "value": config.MIN_PAIRED_PATIENTS_PER_FEATURE,
            },
        ]
    )

    save_dataframe(
        summary,
        config.FEATURE_DEFINITION_DIR / "informative_feature_summary.csv",
    )

    logging.info(
        f"Stage 3 completed. Informative features retained: {len(informative_df)}"
    )

    return informative_df


# ---------------------------------------------------------------------
# Stage 4 — Clinical metadata merge
# ---------------------------------------------------------------------

def stage04_clinical_metadata_merge() -> pd.DataFrame:
    """
    Merge patient-level clinical flags onto the cleaned segment-level ECG dataset.

    This version matches the actual patient_clinical_flags.csv schema.

    Actual available columns:
    - Patient_ID
    - paced_or_icd
    - continuous_AF
    - intermittent_or_new_AF
    - any_AF
    - atypical_HI
    - major_HRV_outlier
    - BBB
    - antiarrhythmic_exposure
    - mechanical_support
    - clean_sinus_candidate
    - exclude_from_HRV_valid_subgroup
    - notes

    Outputs:
    - clinical_flags_column_audit.csv
    - clinical_flags_standardized.csv
    - clinical_merge_audit.csv
    - merged_segment_clinical_dataset.csv
    - patient_clinical_flag_summary.csv
    - cohort_description_table.csv
    - patient_clinical_summary_merged_reference.csv
    """

    logging.info("Stage 4 started: Clinical metadata merge.")

    segment_df = clean_column_names(
        pd.read_csv(config.CLEANED_DATA_DIR / "cleaned_segment_feature_dataset.csv")
    )

    clinical_flags = clean_column_names(pd.read_csv(config.CLINICAL_FLAGS_FILE))
    clinical_summary = clean_column_names(pd.read_csv(config.CLINICAL_SUMMARY_FILE))

    # -----------------------------------------------------------------
    # Clinical column audit
    # -----------------------------------------------------------------

    clinical_column_audit = pd.DataFrame(
        {
            "column": clinical_flags.columns.tolist(),
            "present": True,
            "dtype": [str(clinical_flags[c].dtype) for c in clinical_flags.columns],
            "non_missing_count": [
                int(clinical_flags[c].notna().sum()) for c in clinical_flags.columns
            ],
        }
    )

    save_dataframe(
        clinical_column_audit,
        config.CLINICAL_METADATA_DIR / "clinical_flags_column_audit.csv",
    )

    # -----------------------------------------------------------------
    # Standardize patient IDs
    # -----------------------------------------------------------------

    segment_df[config.PATIENT_ID_COL] = standardize_patient_id_series(
        segment_df[config.PATIENT_ID_COL]
    )
    
    # Create stable row-level identifier for later QC merging.
    # This avoids ambiguity when Segment_Number repeats within patient/condition.
    segment_df[config.ROW_ID_COL] = np.arange(len(segment_df))

    if config.CLINICAL_PATIENT_ID_COL in clinical_flags.columns:
        clinical_flags = clinical_flags.rename(
            columns={config.CLINICAL_PATIENT_ID_COL: config.PATIENT_ID_COL}
        )
    elif config.PATIENT_ID_COL not in clinical_flags.columns:
        raise ValueError(
            f"Could not find patient ID column in clinical flags. "
            f"Expected {config.CLINICAL_PATIENT_ID_COL} or {config.PATIENT_ID_COL}."
        )

    clinical_flags[config.PATIENT_ID_COL] = standardize_patient_id_series(
        clinical_flags[config.PATIENT_ID_COL]
    )

    if config.CLINICAL_PATIENT_ID_COL in clinical_summary.columns:
        clinical_summary = clinical_summary.rename(
            columns={config.CLINICAL_PATIENT_ID_COL: config.PATIENT_ID_COL}
        )
    elif config.PATIENT_ID_COL not in clinical_summary.columns:
        raise ValueError(
            f"Could not find patient ID column in clinical summary. "
            f"Expected {config.CLINICAL_PATIENT_ID_COL} or {config.PATIENT_ID_COL}."
        )

    clinical_summary[config.PATIENT_ID_COL] = standardize_patient_id_series(
        clinical_summary[config.PATIENT_ID_COL]
    )

    # -----------------------------------------------------------------
    # Duplicate checks
    # -----------------------------------------------------------------

    clinical_flag_duplicates = int(
        clinical_flags[config.PATIENT_ID_COL].duplicated().sum()
    )

    clinical_summary_duplicates = int(
        clinical_summary[config.PATIENT_ID_COL].duplicated().sum()
    )

    if clinical_flag_duplicates > 0:
        logging.warning(
            f"Clinical flags contains duplicated Patient_id rows: "
            f"{clinical_flag_duplicates}"
        )

    if clinical_summary_duplicates > 0:
        logging.warning(
            f"Clinical summary contains duplicated Patient_id rows: "
            f"{clinical_summary_duplicates}"
        )

    # -----------------------------------------------------------------
    # Convert actual available clinical flag columns to binary
    # -----------------------------------------------------------------

    actual_binary_cols = [
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
    ]

    for col in actual_binary_cols:
        if col in clinical_flags.columns:
            clinical_flags[col] = to_binary_flag(clinical_flags[col])
        else:
            logging.warning(f"Expected clinical flag column missing: {col}")
            clinical_flags[col] = 0

    # -----------------------------------------------------------------
    # Derived analysis-facing flags
    # -----------------------------------------------------------------

    clinical_flags["any_af"] = (
        (clinical_flags["any_AF"] == 1)
        | (clinical_flags["continuous_AF"] == 1)
        | (clinical_flags["intermittent_or_new_AF"] == 1)
    ).astype(int)

    clinical_flags["paced_or_bbb"] = (
        (clinical_flags["paced_or_icd"] == 1)
        | (clinical_flags["BBB"] == 1)
    ).astype(int)

    # Use existing mechanical_support directly.
    # Do NOT overwrite it from unavailable columns.
    clinical_flags["mechanical_support"] = clinical_flags[
        "mechanical_support"
    ].astype(int)

    # HRV invalid flag.
    # Prefer the curated exclude_from_HRV_valid_subgroup column,
    # but also enforce obvious HRV-invalid conditions.
    clinical_flags["hrv_invalid_flag"] = (
        (clinical_flags["exclude_from_HRV_valid_subgroup"] == 1)
        | (clinical_flags["paced_or_icd"] == 1)
        | (clinical_flags["continuous_AF"] == 1)
        | (clinical_flags["major_HRV_outlier"] == 1)
    ).astype(int)

    # Confound count from actual available columns.
    present_confound_cols = [
        c for c in config.CONFOUND_COUNT_SOURCE_COLUMNS if c in clinical_flags.columns
    ]

    clinical_flags["confound_count"] = clinical_flags[present_confound_cols].sum(
        axis=1
    ).astype(int)

    clinical_flags["high_confound"] = (
        clinical_flags["confound_count"] >= 3
    ).astype(int)

    clinical_flags["moderate_or_high_confound"] = (
        clinical_flags["confound_count"] >= 2
    ).astype(int)

    # Signal/artifact flags for now are based on major HRV outlier and known patient cautions.
    # More detailed artifact flags will come from segment-level QC in Stage 5.
    clinical_flags["signal_review_flag"] = (
        clinical_flags["major_HRV_outlier"] == 1
    ).astype(int)

    clinical_flags["known_patient_qc_caution"] = clinical_flags[
        config.PATIENT_ID_COL
    ].map(config.KNOWN_PATIENT_QC_CAUTIONS)

    clinical_flags["known_patient_qc_caution_flag"] = (
        clinical_flags["known_patient_qc_caution"].notna().astype(int)
    )

    clinical_flags["known_artifact_patient_flag"] = (
        (clinical_flags["signal_review_flag"] == 1)
        | (clinical_flags["known_patient_qc_caution_flag"] == 1)
    ).astype(int)

    # -----------------------------------------------------------------
    # Save standardized clinical flags
    # -----------------------------------------------------------------

    save_dataframe(
        clinical_flags,
        config.CLINICAL_METADATA_DIR / "clinical_flags_standardized.csv",
    )

    # -----------------------------------------------------------------
    # Merge onto segment-level data
    # -----------------------------------------------------------------

    merged = segment_df.merge(
        clinical_flags,
        on=config.PATIENT_ID_COL,
        how="left",
        validate="many_to_one",
        indicator=True,
    )

    merge_audit = pd.DataFrame(
        [
            {
                "metric": "n_segment_rows_before_merge",
                "value": len(segment_df),
            },
            {
                "metric": "n_segment_rows_after_merge",
                "value": len(merged),
            },
            {
                "metric": "n_unique_segment_patients",
                "value": segment_df[config.PATIENT_ID_COL].nunique(),
            },
            {
                "metric": "n_unique_clinical_flag_patients",
                "value": clinical_flags[config.PATIENT_ID_COL].nunique(),
            },
            {
                "metric": "n_unique_merged_patients",
                "value": merged[config.PATIENT_ID_COL].nunique(),
            },
            {
                "metric": "n_rows_matched_both",
                "value": int((merged["_merge"] == "both").sum()),
            },
            {
                "metric": "n_rows_left_only_no_clinical_match",
                "value": int((merged["_merge"] == "left_only").sum()),
            },
            {
                "metric": "n_clinical_flag_duplicate_patient_rows",
                "value": clinical_flag_duplicates,
            },
            {
                "metric": "n_clinical_summary_duplicate_patient_rows",
                "value": clinical_summary_duplicates,
            },
        ]
    )

    save_dataframe(
        merge_audit,
        config.CLINICAL_METADATA_DIR / "clinical_merge_audit.csv",
    )

    unmatched_patients = sorted(
        merged.loc[merged["_merge"] == "left_only", config.PATIENT_ID_COL]
        .dropna()
        .unique()
        .tolist()
    )

    if unmatched_patients:
        logging.warning(f"Patients without clinical flag match: {unmatched_patients}")

    merged = merged.drop(columns=["_merge"])

    save_dataframe(
        merged,
        config.CLINICAL_METADATA_DIR / "merged_segment_clinical_dataset.csv",
    )

    # -----------------------------------------------------------------
    # Patient-level clinical flag summary
    # -----------------------------------------------------------------

    summary_cols = [
        config.PATIENT_ID_COL,
        "paced_or_icd",
        "BBB",
        "continuous_AF",
        "intermittent_or_new_AF",
        "any_AF",
        "any_af",
        "paced_or_bbb",
        "major_HRV_outlier",
        "antiarrhythmic_exposure",
        "mechanical_support",
        "atypical_HI",
        "clean_sinus_candidate",
        "exclude_from_HRV_valid_subgroup",
        "hrv_invalid_flag",
        "signal_review_flag",
        "known_artifact_patient_flag",
        "known_patient_qc_caution_flag",
        "confound_count",
        "high_confound",
        "moderate_or_high_confound",
        "notes",
        "known_patient_qc_caution",
    ]

    summary_cols_present = [c for c in summary_cols if c in clinical_flags.columns]

    patient_clinical_flag_summary = clinical_flags[summary_cols_present].copy()

    save_dataframe(
        patient_clinical_flag_summary,
        config.CLINICAL_METADATA_DIR / "patient_clinical_flag_summary.csv",
    )

    # -----------------------------------------------------------------
    # Cohort description table
    # -----------------------------------------------------------------

    cohort_rows = []

    cohort_rows.append(
        {
            "characteristic": "Total patients",
            "value": clinical_flags[config.PATIENT_ID_COL].nunique(),
        }
    )

    cohort_rows.append(
        {
            "characteristic": "Total ECG segments",
            "value": len(segment_df),
        }
    )

    binary_summary_cols = [
        "paced_or_icd",
        "BBB",
        "continuous_AF",
        "intermittent_or_new_AF",
        "any_af",
        "paced_or_bbb",
        "major_HRV_outlier",
        "antiarrhythmic_exposure",
        "mechanical_support",
        "atypical_HI",
        "clean_sinus_candidate",
        "exclude_from_HRV_valid_subgroup",
        "hrv_invalid_flag",
        "signal_review_flag",
        "known_artifact_patient_flag",
        "high_confound",
        "moderate_or_high_confound",
    ]

    for col in binary_summary_cols:
        if col in clinical_flags.columns:
            n_flagged = int(clinical_flags[col].sum())
            n_total = int(clinical_flags[col].notna().sum())
            percent = 100 * n_flagged / n_total if n_total > 0 else np.nan

            cohort_rows.append(
                {
                    "characteristic": col,
                    "value": f"{n_flagged}/{n_total} ({percent:.1f}%)",
                }
            )

    cohort_rows.append(
        {
            "characteristic": "confound_count_median_IQR",
            "value": (
                f"{clinical_flags['confound_count'].median():.1f} "
                f"[{clinical_flags['confound_count'].quantile(0.25):.1f}, "
                f"{clinical_flags['confound_count'].quantile(0.75):.1f}]"
            ),
        }
    )

    cohort_description_table = pd.DataFrame(cohort_rows)

    save_dataframe(
        cohort_description_table,
        config.CLINICAL_METADATA_DIR / "cohort_description_table.csv",
    )

    # Save clinical summary as reference with standardized Patient_id
    save_dataframe(
        clinical_summary,
        config.CLINICAL_METADATA_DIR / "patient_clinical_summary_merged_reference.csv",
    )

    logging.info("Stage 4 completed: Clinical metadata merge.")

    return merged


# ---------------------------------------------------------------------
# Stage 5 — Segment-level QC flags
# ---------------------------------------------------------------------

def stage05_segment_qc_flags() -> pd.DataFrame:
    """
    Create conservative segment-level QC flags.

    These flags are used for sensitivity analysis and interpretation.
    They do not automatically remove segments from the primary full-cohort analysis.

    Outputs:
    - segment_qc_flags.csv
    - segment_qc_summary_by_patient_condition.csv
    - qc_exclusion_counts.csv
    - qc_flagged_segments_for_review.csv
    """

    logging.info("Stage 5 started: Segment-level QC flags.")

    merged = pd.read_csv(
        config.CLINICAL_METADATA_DIR / "merged_segment_clinical_dataset.csv"
    )

    informative_features = pd.read_csv(
        config.FEATURE_DEFINITION_DIR / "informative_feature_list.csv"
    )

    feature_names = informative_features["feature"].tolist()

    # Basic identifying columns
    id_cols = [
        config.ROW_ID_COL,
        config.PATIENT_ID_COL,
        config.SEGMENT_COL,
        config.CONDITION_STD_COL,
    ]
    
    optional_id_cols = [
        config.CLASS_COL,
        config.CONDITION_COL,
        config.EPISODE_COL,
    ]

    id_cols = id_cols + [c for c in optional_id_cols if c in merged.columns and c not in id_cols]

    qc = merged[id_cols].copy()

    qc[config.PATIENT_ID_COL] = standardize_patient_id_series(
        qc[config.PATIENT_ID_COL]
    )

    # -----------------------------------------------------------------
    # Missingness QC
    # -----------------------------------------------------------------

    available_features = [f for f in feature_names if f in merged.columns]

    feature_missing_matrix = merged[available_features].isna()

    qc["n_missing_informative_features"] = feature_missing_matrix.sum(axis=1)
    qc["missing_any_informative_feature_flag"] = (
        qc["n_missing_informative_features"] > 0
    ).astype(int)

    # -----------------------------------------------------------------
    # Near-zero signal and low-variability QC
    # -----------------------------------------------------------------

    if "signal_energy" in merged.columns:
        signal_energy = pd.to_numeric(merged["signal_energy"], errors="coerce")
        qc["near_zero_signal_flag"] = (
            signal_energy.abs() <= config.NEAR_ZERO_SIGNAL_ENERGY_ABS_THRESHOLD
        ).fillna(False).astype(int)
    else:
        qc["near_zero_signal_flag"] = 0

    if "standardDeviation" in merged.columns:
        std_values = pd.to_numeric(merged["standardDeviation"], errors="coerce")
        low_std_flag = (
            std_values.abs() <= config.LOW_VARIABILITY_STD_ABS_THRESHOLD
        ).fillna(False).astype(int)
    else:
        low_std_flag = pd.Series(np.zeros(len(merged), dtype=int), index=merged.index)

    if "variance" in merged.columns:
        var_values = pd.to_numeric(merged["variance"], errors="coerce")
        low_var_flag = (
            var_values.abs() <= config.LOW_VARIABILITY_VAR_ABS_THRESHOLD
        ).fillna(False).astype(int)
    else:
        low_var_flag = pd.Series(np.zeros(len(merged), dtype=int), index=merged.index)

    qc["low_variability_signal_flag"] = (
        (low_std_flag == 1) | (low_var_flag == 1)
    ).astype(int)

    # -----------------------------------------------------------------
    # QRS_std zero flag
    # -----------------------------------------------------------------

    if "QRS_std" in merged.columns:
        qrs_std = pd.to_numeric(merged["QRS_std"], errors="coerce")
        qc["qrs_std_zero_flag"] = (
            qrs_std.abs() <= config.LOW_VARIABILITY_STD_ABS_THRESHOLD
        ).fillna(False).astype(int)
    else:
        qc["qrs_std_zero_flag"] = 0

    # -----------------------------------------------------------------
    # Known artifact segment flag
    # -----------------------------------------------------------------

    qc["known_artifact_segment_flag"] = 0
    qc["known_artifact_segment_reason"] = ""

    for item in config.KNOWN_ARTIFACT_SEGMENTS:
        patient_id = str(item["Patient_id"])
        condition = item["Condition_std"]
        segment_number = item["Segment_Number"]
        reason = item["reason"]

        match_mask = (
            (qc[config.PATIENT_ID_COL].astype(str) == patient_id)
            & (qc[config.CONDITION_STD_COL].astype(str) == condition)
            & (pd.to_numeric(qc[config.SEGMENT_COL], errors="coerce") == segment_number)
        )

        qc.loc[match_mask, "known_artifact_segment_flag"] = 1
        qc.loc[match_mask, "known_artifact_segment_reason"] = reason

    # -----------------------------------------------------------------
    # Extreme-value review flag
    # -----------------------------------------------------------------

    review_feature_flags = []

    for feature in available_features:
        if feature == "kurtosis":
            continue

        numeric_values = pd.to_numeric(merged[feature], errors="coerce")
        modified_z = robust_modified_zscore(numeric_values)

        flag_col = f"extreme_review__{feature}"
        qc[flag_col] = (
            modified_z.abs() >= config.EXTREME_VALUE_MODIFIED_Z_THRESHOLD
        ).fillna(False).astype(int)

        review_feature_flags.append(flag_col)

    if review_feature_flags:
        qc["n_extreme_value_review_features"] = qc[review_feature_flags].sum(axis=1)
        qc["extreme_value_review_flag"] = (
            qc["n_extreme_value_review_features"] > 0
        ).astype(int)
    else:
        qc["n_extreme_value_review_features"] = 0
        qc["extreme_value_review_flag"] = 0

    # -----------------------------------------------------------------
    # Extreme kurtosis review flag
    # -----------------------------------------------------------------

    if "kurtosis" in merged.columns:
        kurtosis_values = pd.to_numeric(merged["kurtosis"], errors="coerce")
        kurtosis_z = robust_modified_zscore(kurtosis_values)
        qc["extreme_kurtosis_review_flag"] = (
            kurtosis_z.abs() >= config.EXTREME_KURTOSIS_MODIFIED_Z_THRESHOLD
        ).fillna(False).astype(int)
    else:
        qc["extreme_kurtosis_review_flag"] = 0

    # -----------------------------------------------------------------
    # Patient-level caution flag copied to segment-level
    # -----------------------------------------------------------------

    if "known_patient_qc_caution_flag" in merged.columns:
        qc["known_patient_qc_caution_flag"] = pd.to_numeric(
            merged["known_patient_qc_caution_flag"], errors="coerce"
        ).fillna(0).astype(int)
    else:
        qc["known_patient_qc_caution_flag"] = 0

    if "known_patient_qc_caution" in merged.columns:
        qc["known_patient_qc_caution"] = merged["known_patient_qc_caution"].fillna("")
    else:
        qc["known_patient_qc_caution"] = ""

    if "signal_review_flag" in merged.columns:
        qc["signal_review_flag"] = pd.to_numeric(
            merged["signal_review_flag"], errors="coerce"
        ).fillna(0).astype(int)
    else:
        qc["signal_review_flag"] = 0

    # -----------------------------------------------------------------
    # Clear segment flag for secondary QC-clean analysis
    # -----------------------------------------------------------------

    # Conservative exclusion flag for QC-clean sensitivity analysis.
    # Extreme review flags are NOT automatic exclusions.
    qc["clear_segment_qc_exclusion_flag"] = (
        (qc["near_zero_signal_flag"] == 1)
        | (qc["low_variability_signal_flag"] == 1)
        | (qc["qrs_std_zero_flag"] == 1)
        | (qc["known_artifact_segment_flag"] == 1)
    ).astype(int)

    qc["qc_clear_segment_flag"] = (
        qc["clear_segment_qc_exclusion_flag"] == 0
    ).astype(int)

    # -----------------------------------------------------------------
    # Save QC files
    # -----------------------------------------------------------------

    save_dataframe(
        qc,
        config.SEGMENT_QC_DIR / "segment_qc_flags.csv",
    )

    summary_cols = [
        "missing_any_informative_feature_flag",
        "near_zero_signal_flag",
        "low_variability_signal_flag",
        "qrs_std_zero_flag",
        "known_artifact_segment_flag",
        "extreme_value_review_flag",
        "extreme_kurtosis_review_flag",
        "known_patient_qc_caution_flag",
        "signal_review_flag",
        "clear_segment_qc_exclusion_flag",
        "qc_clear_segment_flag",
    ]

    qc_summary_by_patient_condition = (
        qc.groupby([config.PATIENT_ID_COL, config.CONDITION_STD_COL])[summary_cols]
        .agg(["sum", "mean"])
        .reset_index()
    )

    qc_summary_by_patient_condition.columns = [
        "_".join([str(x) for x in col if str(x) != ""])
        for col in qc_summary_by_patient_condition.columns
    ]

    save_dataframe(
        qc_summary_by_patient_condition,
        config.SEGMENT_QC_DIR / "segment_qc_summary_by_patient_condition.csv",
    )

    qc_exclusion_counts = pd.DataFrame(
        [
            {
                "qc_flag": col,
                "n_segments_flagged": int(qc[col].sum()),
                "percent_segments_flagged": 100 * float(qc[col].sum()) / len(qc),
            }
            for col in summary_cols
        ]
    ).sort_values("n_segments_flagged", ascending=False)

    save_dataframe(
        qc_exclusion_counts,
        config.SEGMENT_QC_DIR / "qc_exclusion_counts.csv",
    )

    flagged_for_review = qc[
        (qc["clear_segment_qc_exclusion_flag"] == 1)
        | (qc["extreme_value_review_flag"] == 1)
        | (qc["extreme_kurtosis_review_flag"] == 1)
        | (qc["known_patient_qc_caution_flag"] == 1)
        | (qc["signal_review_flag"] == 1)
    ].copy()

    save_dataframe(
        flagged_for_review,
        config.SEGMENT_QC_DIR / "qc_flagged_segments_for_review.csv",
    )

    logging.info("Stage 5 completed: Segment-level QC flags.")

    return qc

# ---------------------------------------------------------------------
# Stage 6 — Patient-condition medians and paired deltas
# ---------------------------------------------------------------------

def stage06_patient_condition_medians() -> pd.DataFrame:
    """
    Aggregate segment-level ECG features into patient-condition medians.

    Primary full-cohort analysis:
    - Uses all non-missing feature values.
    - Does not exclude QC-flagged segments.

    Secondary QC-clean medians:
    - Excludes segments where clear_segment_qc_exclusion_flag == 1.
    - Extreme-value review flags are not automatic exclusions.

    Outputs:
    - row_id_merge_audit.csv
    - patient_condition_medians_all_features_full.csv
    - patient_condition_medians_all_features_qc_clean.csv
    - patient_feature_delta_table_all_features.csv
    - patient_level_summary_all_features.csv
    """

    logging.info("Stage 6 started: Patient-condition medians and paired deltas.")

    merged = pd.read_csv(
        config.CLINICAL_METADATA_DIR / "merged_segment_clinical_dataset.csv"
    )

    qc = pd.read_csv(
        config.SEGMENT_QC_DIR / "segment_qc_flags.csv"
    )

    informative_features = pd.read_csv(
        config.FEATURE_DEFINITION_DIR / "informative_feature_list.csv"
    )

    feature_names = informative_features["feature"].tolist()

    # Standardize IDs
    merged[config.PATIENT_ID_COL] = standardize_patient_id_series(
        merged[config.PATIENT_ID_COL]
    )

    qc[config.PATIENT_ID_COL] = standardize_patient_id_series(
        qc[config.PATIENT_ID_COL]
    )

    # -----------------------------------------------------------------
    # Row-level QC merge audit
    # -----------------------------------------------------------------
    # We merge QC flags using segment_row_id, not Patient_id + Segment_Number
    # + Condition_std, because Segment_Number may repeat within a patient or
    # condition. This audit confirms the row ID is unique in both files.

    if config.ROW_ID_COL not in merged.columns:
        raise ValueError(
            f"{config.ROW_ID_COL} is missing from merged_segment_clinical_dataset.csv. "
            "Please rerun Stage 4 after adding ROW_ID_COL creation."
        )

    if config.ROW_ID_COL not in qc.columns:
        raise ValueError(
            f"{config.ROW_ID_COL} is missing from segment_qc_flags.csv. "
            "Please rerun Stage 5 after adding ROW_ID_COL to id_cols."
        )

    row_id_audit = pd.DataFrame(
        [
            {
                "dataset": "merged_segment_clinical_dataset",
                "n_rows": len(merged),
                "n_unique_segment_row_id": merged[config.ROW_ID_COL].nunique(),
                "n_duplicate_segment_row_id": int(
                    merged[config.ROW_ID_COL].duplicated().sum()
                ),
                "n_missing_segment_row_id": int(
                    merged[config.ROW_ID_COL].isna().sum()
                ),
            },
            {
                "dataset": "segment_qc_flags",
                "n_rows": len(qc),
                "n_unique_segment_row_id": qc[config.ROW_ID_COL].nunique(),
                "n_duplicate_segment_row_id": int(
                    qc[config.ROW_ID_COL].duplicated().sum()
                ),
                "n_missing_segment_row_id": int(
                    qc[config.ROW_ID_COL].isna().sum()
                ),
            },
        ]
    )

    save_dataframe(
        row_id_audit,
        config.PATIENT_LEVEL_DIR / "row_id_merge_audit.csv",
    )

    # Fail loudly if the supposedly stable row ID is not actually stable.
    if merged[config.ROW_ID_COL].duplicated().any():
        raise ValueError(
            f"Duplicate {config.ROW_ID_COL} values found in "
            "merged_segment_clinical_dataset.csv. Cannot safely merge QC flags."
        )

    if qc[config.ROW_ID_COL].duplicated().any():
        raise ValueError(
            f"Duplicate {config.ROW_ID_COL} values found in segment_qc_flags.csv. "
            "Cannot safely merge QC flags."
        )

    if merged[config.ROW_ID_COL].isna().any():
        raise ValueError(
            f"Missing {config.ROW_ID_COL} values found in "
            "merged_segment_clinical_dataset.csv."
        )

    if qc[config.ROW_ID_COL].isna().any():
        raise ValueError(
            f"Missing {config.ROW_ID_COL} values found in segment_qc_flags.csv."
        )

    # -----------------------------------------------------------------
    # Merge QC flags back onto segment data using stable row-level ID
    # -----------------------------------------------------------------

    merge_keys = [config.ROW_ID_COL]

    qc_cols_to_merge = merge_keys + [
        "clear_segment_qc_exclusion_flag",
        "qc_clear_segment_flag",
        "missing_any_informative_feature_flag",
        "near_zero_signal_flag",
        "low_variability_signal_flag",
        "qrs_std_zero_flag",
        "known_artifact_segment_flag",
        "extreme_value_review_flag",
        "extreme_kurtosis_review_flag",
        "known_patient_qc_caution_flag",
        "signal_review_flag",
    ]

    qc_cols_to_merge = [c for c in qc_cols_to_merge if c in qc.columns]

    merged_qc = merged.merge(
        qc[qc_cols_to_merge],
        on=merge_keys,
        how="left",
        validate="one_to_one",
    )

    if "clear_segment_qc_exclusion_flag" not in merged_qc.columns:
        merged_qc["clear_segment_qc_exclusion_flag"] = 0

    if "qc_clear_segment_flag" not in merged_qc.columns:
        merged_qc["qc_clear_segment_flag"] = 1

    # Ensure features are numeric
    for feature in feature_names:
        merged_qc[feature] = pd.to_numeric(merged_qc[feature], errors="coerce")

    # -----------------------------------------------------------------
    # Full patient-condition medians
    # -----------------------------------------------------------------

    full_medians = (
        merged_qc
        .groupby([config.PATIENT_ID_COL, config.CONDITION_STD_COL])[feature_names]
        .median(numeric_only=True)
        .reset_index()
    )

    segment_counts_full = (
        merged_qc
        .groupby([config.PATIENT_ID_COL, config.CONDITION_STD_COL])
        .size()
        .reset_index(name="n_segments_total")
    )

    qc_flag_counts = (
        merged_qc
        .groupby([config.PATIENT_ID_COL, config.CONDITION_STD_COL])[
            "clear_segment_qc_exclusion_flag"
        ]
        .sum()
        .reset_index(name="n_clear_qc_exclusion_segments")
    )

    full_medians = full_medians.merge(
        segment_counts_full,
        on=[config.PATIENT_ID_COL, config.CONDITION_STD_COL],
        how="left",
    )

    full_medians = full_medians.merge(
        qc_flag_counts,
        on=[config.PATIENT_ID_COL, config.CONDITION_STD_COL],
        how="left",
    )

    save_dataframe(
        full_medians,
        config.PATIENT_LEVEL_DIR / "patient_condition_medians_all_features_full.csv",
    )

    # -----------------------------------------------------------------
    # QC-clean patient-condition medians
    # -----------------------------------------------------------------

    qc_clean_df = merged_qc[
        merged_qc["clear_segment_qc_exclusion_flag"] == 0
    ].copy()

    qc_clean_medians = (
        qc_clean_df
        .groupby([config.PATIENT_ID_COL, config.CONDITION_STD_COL])[feature_names]
        .median(numeric_only=True)
        .reset_index()
    )

    segment_counts_qc_clean = (
        qc_clean_df
        .groupby([config.PATIENT_ID_COL, config.CONDITION_STD_COL])
        .size()
        .reset_index(name="n_segments_qc_clean")
    )

    qc_clean_medians = qc_clean_medians.merge(
        segment_counts_qc_clean,
        on=[config.PATIENT_ID_COL, config.CONDITION_STD_COL],
        how="left",
    )

    save_dataframe(
        qc_clean_medians,
        config.PATIENT_LEVEL_DIR / "patient_condition_medians_all_features_qc_clean.csv",
    )

    # -----------------------------------------------------------------
    # Paired delta table, one row per patient-feature
    # -----------------------------------------------------------------

    delta_rows = []

    clinical_summary_cols = [
        "paced_or_icd",
        "BBB",
        "continuous_AF",
        "intermittent_or_new_AF",
        "any_AF",
        "any_af",
        "paced_or_bbb",
        "major_HRV_outlier",
        "antiarrhythmic_exposure",
        "mechanical_support",
        "atypical_HI",
        "clean_sinus_candidate",
        "exclude_from_HRV_valid_subgroup",
        "hrv_invalid_flag",
        "signal_review_flag",
        "known_artifact_patient_flag",
        "known_patient_qc_caution_flag",
        "confound_count",
        "high_confound",
        "moderate_or_high_confound",
    ]

    clinical_summary_cols = [
        c for c in clinical_summary_cols if c in merged_qc.columns
    ]

    patient_clinical = (
        merged_qc[[config.PATIENT_ID_COL] + clinical_summary_cols]
        .drop_duplicates(subset=[config.PATIENT_ID_COL])
        .copy()
    )

    for feature in feature_names:
        tmp = full_medians[
            [
                config.PATIENT_ID_COL,
                config.CONDITION_STD_COL,
                feature,
                "n_segments_total",
                "n_clear_qc_exclusion_segments",
            ]
        ].copy()

        wide_value = tmp.pivot(
            index=config.PATIENT_ID_COL,
            columns=config.CONDITION_STD_COL,
            values=feature,
        )

        wide_n = tmp.pivot(
            index=config.PATIENT_ID_COL,
            columns=config.CONDITION_STD_COL,
            values="n_segments_total",
        )

        wide_qc = tmp.pivot(
            index=config.PATIENT_ID_COL,
            columns=config.CONDITION_STD_COL,
            values="n_clear_qc_exclusion_segments",
        )

        for patient_id in wide_value.index:
            pre_value = (
                wide_value.loc[patient_id, "PreHI"]
                if "PreHI" in wide_value.columns
                else np.nan
            )

            hi_value = (
                wide_value.loc[patient_id, "HI"]
                if "HI" in wide_value.columns
                else np.nan
            )

            pre_n = (
                wide_n.loc[patient_id, "PreHI"]
                if "PreHI" in wide_n.columns
                else np.nan
            )

            hi_n = (
                wide_n.loc[patient_id, "HI"]
                if "HI" in wide_n.columns
                else np.nan
            )

            pre_qc_n = (
                wide_qc.loc[patient_id, "PreHI"]
                if "PreHI" in wide_qc.columns
                else np.nan
            )

            hi_qc_n = (
                wide_qc.loc[patient_id, "HI"]
                if "HI" in wide_qc.columns
                else np.nan
            )

            if pd.notna(pre_value) and pd.notna(hi_value):
                delta = hi_value - pre_value

                if abs(pre_value) <= config.NEAR_ZERO_PREHI_ABS_THRESHOLD:
                    percent_change = np.nan
                    percent_change_unstable_flag = 1
                else:
                    percent_change = 100 * delta / abs(pre_value)
                    percent_change_unstable_flag = 0

                direction = direction_from_delta(delta)

            else:
                delta = np.nan
                percent_change = np.nan
                percent_change_unstable_flag = 1
                direction = "missing"

            delta_rows.append(
                {
                    config.PATIENT_ID_COL: patient_id,
                    "feature": feature,
                    "feature_family": infer_feature_family(feature),
                    "PreHI_median": pre_value,
                    "HI_median": hi_value,
                    "delta_HI_minus_PreHI": delta,
                    "percent_change": percent_change,
                    "percent_change_unstable_flag": percent_change_unstable_flag,
                    "observed_direction": direction,
                    "n_segments_PreHI": pre_n,
                    "n_segments_HI": hi_n,
                    "n_clear_qc_exclusion_segments_PreHI": pre_qc_n,
                    "n_clear_qc_exclusion_segments_HI": hi_qc_n,
                    "paired_observation_available": int(
                        pd.notna(pre_value) and pd.notna(hi_value)
                    ),
                }
            )

    delta_df = pd.DataFrame(delta_rows)

    delta_df = delta_df.merge(
        patient_clinical,
        on=config.PATIENT_ID_COL,
        how="left",
    )

    save_dataframe(
        delta_df,
        config.PATIENT_LEVEL_DIR / "patient_feature_delta_table_all_features.csv",
    )

    # -----------------------------------------------------------------
    # Feature-level patient summary
    # -----------------------------------------------------------------

    summary_rows = []

    for feature in feature_names:
        sub = delta_df[delta_df["feature"] == feature].copy()
        paired = sub[sub["paired_observation_available"] == 1].copy()

        n_paired = len(paired)
        n_increase = int((paired["observed_direction"] == "increase").sum())
        n_decrease = int((paired["observed_direction"] == "decrease").sum())
        n_no_change = int((paired["observed_direction"] == "no_change").sum())

        dominant_direction = dominant_direction_from_counts(
            n_increase,
            n_decrease,
            n_no_change,
        )

        max_direction_count = (
            max(n_increase, n_decrease, n_no_change) if n_paired > 0 else 0
        )

        direction_consistency = (
            max_direction_count / n_paired if n_paired > 0 else np.nan
        )

        summary_rows.append(
            {
                "feature": feature,
                "feature_family": infer_feature_family(feature),
                "n_paired_patients": n_paired,
                "n_increase": n_increase,
                "n_decrease": n_decrease,
                "n_no_change": n_no_change,
                "dominant_direction": dominant_direction,
                "direction_consistency": direction_consistency,
                "median_PreHI": paired["PreHI_median"].median(),
                "median_HI": paired["HI_median"].median(),
                "median_delta_HI_minus_PreHI": paired[
                    "delta_HI_minus_PreHI"
                ].median(),
                "median_percent_change": paired["percent_change"].median(),
                "n_percent_change_unstable": int(
                    paired["percent_change_unstable_flag"].sum()
                ),
            }
        )

    patient_level_summary = pd.DataFrame(summary_rows).sort_values(
        ["direction_consistency", "n_paired_patients"],
        ascending=[False, False],
    )

    save_dataframe(
        patient_level_summary,
        config.PATIENT_LEVEL_DIR / "patient_level_summary_all_features.csv",
    )

    logging.info("Stage 6 completed: Patient-condition medians and paired deltas.")

    return delta_df

# ---------------------------------------------------------------------
# Stage 7 — Primary paired statistical analysis
# ---------------------------------------------------------------------

def stage07_primary_paired_statistics() -> pd.DataFrame:
    """
    Run primary patient-level paired statistical analysis for all informative features.

    Primary statistical unit:
    - Patient

    For each feature:
    - Wilcoxon signed-rank p-value
    - matched-pairs rank-biserial effect size
    - bootstrap 95% CI for median paired delta
    - direction consistency

    Output:
    - primary_paired_feature_results_all_features.csv
    """

    logging.info("Stage 7 started: Primary paired statistical analysis.")

    delta_df = pd.read_csv(
        config.PATIENT_LEVEL_DIR / "patient_feature_delta_table_all_features.csv"
    )

    informative_features = pd.read_csv(
        config.FEATURE_DEFINITION_DIR / "informative_feature_list.csv"
    )

    feature_names = informative_features["feature"].tolist()

    results = []

    for feature in feature_names:
        sub = delta_df[
            (delta_df["feature"] == feature)
            & (delta_df["paired_observation_available"] == 1)
        ].copy()

        n_paired = len(sub)

        n_increase = int((sub["observed_direction"] == "increase").sum())
        n_decrease = int((sub["observed_direction"] == "decrease").sum())
        n_no_change = int((sub["observed_direction"] == "no_change").sum())

        dominant_direction = dominant_direction_from_counts(
            n_increase,
            n_decrease,
            n_no_change,
        )

        max_direction_count = max(n_increase, n_decrease, n_no_change) if n_paired > 0 else 0
        direction_consistency = max_direction_count / n_paired if n_paired > 0 else np.nan

        if n_paired >= 2:
            p_value = wilcoxon_signed_rank_pvalue(
                sub["PreHI_median"],
                sub["HI_median"],
            )

            effect_size = matched_pairs_rank_biserial_effect_size(
                sub["PreHI_median"],
                sub["HI_median"],
            )

            ci_low, ci_high = bootstrap_median_ci(
                sub["delta_HI_minus_PreHI"],
                n_iterations=config.BOOTSTRAP_N_ITERATIONS,
                confidence_level=config.BOOTSTRAP_CONFIDENCE_LEVEL,
                random_seed=config.RANDOM_SEED,
            )

        else:
            p_value = np.nan
            effect_size = np.nan
            ci_low = np.nan
            ci_high = np.nan

        bootstrap_ci_excludes_zero = (
            pd.notna(ci_low)
            and pd.notna(ci_high)
            and ((ci_low > 0 and ci_high > 0) or (ci_low < 0 and ci_high < 0))
        )

        results.append(
            {
                "feature": feature,
                "feature_family": infer_feature_family(feature),
                "n_paired_patients": n_paired,
                "median_PreHI": sub["PreHI_median"].median(),
                "median_HI": sub["HI_median"].median(),
                "median_delta_HI_minus_PreHI": sub[
                    "delta_HI_minus_PreHI"
                ].median(),
                "median_percent_change": sub["percent_change"].median(),
                "n_increase": n_increase,
                "n_decrease": n_decrease,
                "n_no_change": n_no_change,
                "dominant_direction": dominant_direction,
                "direction_consistency": direction_consistency,
                "wilcoxon_p": p_value,
                "rank_biserial_effect_size": effect_size,
                "bootstrap_median_delta_ci_low": ci_low,
                "bootstrap_median_delta_ci_high": ci_high,
                "bootstrap_ci_excludes_zero": int(bootstrap_ci_excludes_zero),
                "n_percent_change_unstable": int(
                    sub["percent_change_unstable_flag"].sum()
                ),
            }
        )

    results_df = pd.DataFrame(results)

    results_df["abs_rank_biserial_effect_size"] = results_df[
        "rank_biserial_effect_size"
    ].abs()

    results_df = results_df.sort_values(
        [
            "wilcoxon_p",
            "direction_consistency",
            "abs_rank_biserial_effect_size",
        ],
        ascending=[True, False, False],
    )

    save_dataframe(
        results_df,
        config.PRIMARY_STATS_DIR / "primary_paired_feature_results_all_features.csv",
    )

    logging.info("Stage 7 completed: Primary paired statistical analysis.")

    return results_df


# ---------------------------------------------------------------------
# Stage 8 — Multiple-testing correction
# ---------------------------------------------------------------------

def stage08_fdr_correction() -> pd.DataFrame:
    """
    Apply global and feature-family-level FDR correction.

    Outputs:
    - primary_paired_feature_results_all_features_fdr.csv
    - fdr_summary_by_family.csv
    """

    logging.info("Stage 8 started: FDR correction.")

    results_df = pd.read_csv(
        config.PRIMARY_STATS_DIR / "primary_paired_feature_results_all_features.csv"
    )

    # Global FDR across all informative features
    results_df["wilcoxon_p_fdr_global"] = apply_fdr_correction(
        results_df["wilcoxon_p"],
        method=config.FDR_METHOD,
    )

    # Family-level FDR
    results_df["wilcoxon_p_fdr_family"] = np.nan

    for family, idx in results_df.groupby("feature_family").groups.items():
        family_p = results_df.loc[idx, "wilcoxon_p"]
        results_df.loc[idx, "wilcoxon_p_fdr_family"] = apply_fdr_correction(
            family_p,
            method=config.FDR_METHOD,
        )

    results_df["global_fdr_significant"] = (
        results_df["wilcoxon_p_fdr_global"] < config.FDR_ALPHA
    ).astype(int)

    results_df["family_fdr_significant"] = (
        results_df["wilcoxon_p_fdr_family"] < config.FDR_ALPHA
    ).astype(int)

    results_df["raw_p_less_0_05"] = (
        results_df["wilcoxon_p"] < 0.05
    ).astype(int)

    results_df["candidate_support_tier_preliminary"] = "exploratory"

    results_df.loc[
        (results_df["global_fdr_significant"] == 1)
        & (results_df["direction_consistency"] >= 0.70),
        "candidate_support_tier_preliminary",
    ] = "strong_global_fdr"

    results_df.loc[
        (results_df["candidate_support_tier_preliminary"] == "exploratory")
        & (results_df["family_fdr_significant"] == 1)
        & (results_df["direction_consistency"] >= 0.70),
        "candidate_support_tier_preliminary",
    ] = "family_fdr_supported"

    results_df.loc[
        (results_df["candidate_support_tier_preliminary"] == "exploratory")
        & (results_df["raw_p_less_0_05"] == 1)
        & (results_df["direction_consistency"] >= 0.70)
        & (results_df["abs_rank_biserial_effect_size"] >= 0.50),
        "candidate_support_tier_preliminary",
    ] = "raw_p_effect_direction_supported"

    results_df = results_df.sort_values(
        [
            "wilcoxon_p_fdr_global",
            "wilcoxon_p",
            "direction_consistency",
            "abs_rank_biserial_effect_size",
        ],
        ascending=[True, True, False, False],
    )

    save_dataframe(
        results_df,
        config.PRIMARY_STATS_DIR / "primary_paired_feature_results_all_features_fdr.csv",
    )

    fdr_summary_rows = []

    for family, sub in results_df.groupby("feature_family"):
        fdr_summary_rows.append(
            {
                "feature_family": family,
                "n_features": len(sub),
                "n_raw_p_less_0_05": int(sub["raw_p_less_0_05"].sum()),
                "n_global_fdr_significant": int(sub["global_fdr_significant"].sum()),
                "n_family_fdr_significant": int(sub["family_fdr_significant"].sum()),
                "min_raw_p": sub["wilcoxon_p"].min(),
                "min_global_fdr_p": sub["wilcoxon_p_fdr_global"].min(),
                "min_family_fdr_p": sub["wilcoxon_p_fdr_family"].min(),
            }
        )

    fdr_summary_by_family = pd.DataFrame(fdr_summary_rows).sort_values(
        "min_raw_p",
        ascending=True,
    )

    save_dataframe(
        fdr_summary_by_family,
        config.PRIMARY_STATS_DIR / "fdr_summary_by_family.csv",
    )

    logging.info("Stage 8 completed: FDR correction.")

    return results_df

# ---------------------------------------------------------------------
# Stage 9 — Redundancy filtering
# ---------------------------------------------------------------------

def stage09_redundancy_filtering() -> pd.DataFrame:
    """
    Perform within-patient centered Spearman correlation analysis.

    Purpose:
    - Identify highly correlated ECG features.
    - Build redundancy clusters.
    - Select transparent cluster representatives using primary paired results.

    Outputs:
    - within_subject_spearman_correlation_matrix.csv
    - high_correlation_pairs.csv
    - redundancy_clusters.csv
    - cluster_representatives.csv
    """

    logging.info("Stage 9 started: Redundancy filtering.")

    merged = pd.read_csv(
        config.CLINICAL_METADATA_DIR / "merged_segment_clinical_dataset.csv"
    )

    informative_features = pd.read_csv(
        config.FEATURE_DEFINITION_DIR / "informative_feature_list.csv"
    )

    primary_results = pd.read_csv(
        config.PRIMARY_STATS_DIR / "primary_paired_feature_results_all_features_fdr.csv"
    )

    feature_names = informative_features["feature"].tolist()

    merged[config.PATIENT_ID_COL] = standardize_patient_id_series(
        merged[config.PATIENT_ID_COL]
    )

    # Ensure all informative features are numeric
    available_features = [f for f in feature_names if f in merged.columns]

    for feature in available_features:
        merged[feature] = pd.to_numeric(merged[feature], errors="coerce")

    # -----------------------------------------------------------------
    # Within-patient centering
    # -----------------------------------------------------------------
    # We center within patient to reduce between-patient baseline differences.
    # Correlation is then estimated on within-patient deviations.

    centered = merged[[config.PATIENT_ID_COL] + available_features].copy()

    for feature in available_features:
        patient_median = centered.groupby(config.PATIENT_ID_COL)[feature].transform(
            "median"
        )
        centered[feature] = centered[feature] - patient_median

    # Drop features that became entirely missing after centering
    valid_centered_features = [
        f for f in available_features if centered[f].notna().sum() > 1
    ]

    corr_matrix = centered[valid_centered_features].corr(method="spearman")

    save_dataframe(
        corr_matrix.reset_index().rename(columns={"index": "feature"}),
        config.REDUNDANCY_DIR / "within_subject_spearman_correlation_matrix.csv",
    )

    # -----------------------------------------------------------------
    # High-correlation pairs
    # -----------------------------------------------------------------

    pair_rows = []

    for i, feature_a in enumerate(valid_centered_features):
        for feature_b in valid_centered_features[i + 1:]:
            rho = corr_matrix.loc[feature_a, feature_b]

            if pd.isna(rho):
                continue

            abs_rho = abs(rho)

            if abs_rho >= config.REDUNDANCY_CORR_THRESHOLD_MODERATE:
                if abs_rho >= config.REDUNDANCY_CORR_THRESHOLD_HIGH:
                    redundancy_level = "high_abs_rho_ge_0_90"
                else:
                    redundancy_level = "moderate_abs_rho_ge_0_85"

                pair_rows.append(
                    {
                        "feature_a": feature_a,
                        "feature_b": feature_b,
                        "spearman_rho_within_patient_centered": rho,
                        "abs_rho": abs_rho,
                        "redundancy_level": redundancy_level,
                    }
                )

    high_pairs = pd.DataFrame(pair_rows)

    if not high_pairs.empty:
        high_pairs = high_pairs.sort_values("abs_rho", ascending=False)

    save_dataframe(
        high_pairs,
        config.REDUNDANCY_DIR / "high_correlation_pairs.csv",
    )

    # -----------------------------------------------------------------
    # Build redundancy clusters from abs rho >= 0.85 pairs
    # -----------------------------------------------------------------

    if not high_pairs.empty:
        edges = list(zip(high_pairs["feature_a"], high_pairs["feature_b"]))
        components = get_connected_components_from_edges(edges)
    else:
        components = []

    clustered_features = set()
    cluster_rows = []

    for cluster_idx, component in enumerate(components, start=1):
        cluster_id = f"cluster_{cluster_idx:03d}"
        clustered_features.update(component)

        for feature in component:
            cluster_rows.append(
                {
                    "cluster_id": cluster_id,
                    "feature": feature,
                    "cluster_size": len(component),
                    "cluster_type": "correlated_cluster_abs_rho_ge_0_85",
                }
            )

    # Features not in any redundancy cluster are singleton clusters
    singleton_idx = len(components) + 1

    for feature in valid_centered_features:
        if feature not in clustered_features:
            cluster_rows.append(
                {
                    "cluster_id": f"singleton_{singleton_idx:03d}",
                    "feature": feature,
                    "cluster_size": 1,
                    "cluster_type": "singleton",
                }
            )
            singleton_idx += 1

    redundancy_clusters = pd.DataFrame(cluster_rows)

    save_dataframe(
        redundancy_clusters,
        config.REDUNDANCY_DIR / "redundancy_clusters.csv",
    )

    # -----------------------------------------------------------------
    # Select cluster representatives
    # -----------------------------------------------------------------

    rep_input = redundancy_clusters.merge(
        primary_results,
        on="feature",
        how="left",
    )

    # Representative score is deliberately transparent and simple.
    # Higher score means stronger candidate to represent the redundancy cluster.

    rep_input["score_raw_p"] = min_max_normalize(
        safe_negative_log10_pvalue(rep_input["wilcoxon_p"]),
        higher_is_better=True,
    )

    rep_input["score_effect_size"] = min_max_normalize(
        rep_input["abs_rank_biserial_effect_size"],
        higher_is_better=True,
    )

    rep_input["score_direction_consistency"] = min_max_normalize(
        rep_input["direction_consistency"],
        higher_is_better=True,
    )

    rep_input["score_bootstrap"] = pd.to_numeric(
        rep_input.get("bootstrap_ci_excludes_zero", 0),
        errors="coerce",
    ).fillna(0)

    rep_input["score_fdr"] = 0

    if "family_fdr_significant" in rep_input.columns:
        rep_input["score_fdr"] += rep_input["family_fdr_significant"].fillna(0)

    if "global_fdr_significant" in rep_input.columns:
        rep_input["score_fdr"] += 2 * rep_input["global_fdr_significant"].fillna(0)

    rep_input["representative_priority_score"] = (
        2.0 * rep_input["score_raw_p"]
        + 2.0 * rep_input["score_effect_size"]
        + 2.0 * rep_input["score_direction_consistency"]
        + 1.0 * rep_input["score_bootstrap"]
        + 1.0 * rep_input["score_fdr"]
    )

    representatives = []

    for cluster_id, sub in rep_input.groupby("cluster_id"):
        sub_sorted = sub.sort_values(
            [
                "representative_priority_score",
                "direction_consistency",
                "abs_rank_biserial_effect_size",
                "wilcoxon_p",
            ],
            ascending=[False, False, False, True],
        )

        chosen = sub_sorted.iloc[0].copy()

        representatives.append(
            {
                "cluster_id": cluster_id,
                "representative_feature": chosen["feature"],
                "cluster_size": int(chosen["cluster_size"]),
                "cluster_type": chosen["cluster_type"],
                "representative_feature_family": chosen.get(
                    "feature_family", infer_feature_family(chosen["feature"])
                ),
                "representative_priority_score": chosen[
                    "representative_priority_score"
                ],
                "wilcoxon_p": chosen.get("wilcoxon_p", np.nan),
                "wilcoxon_p_fdr_global": chosen.get(
                    "wilcoxon_p_fdr_global", np.nan
                ),
                "wilcoxon_p_fdr_family": chosen.get(
                    "wilcoxon_p_fdr_family", np.nan
                ),
                "direction_consistency": chosen.get(
                    "direction_consistency", np.nan
                ),
                "rank_biserial_effect_size": chosen.get(
                    "rank_biserial_effect_size", np.nan
                ),
                "abs_rank_biserial_effect_size": chosen.get(
                    "abs_rank_biserial_effect_size", np.nan
                ),
                "dominant_direction": chosen.get("dominant_direction", ""),
                "bootstrap_ci_excludes_zero": chosen.get(
                    "bootstrap_ci_excludes_zero", np.nan
                ),
                "cluster_members": "; ".join(
                    sorted(sub["feature"].astype(str).tolist())
                ),
            }
        )

    cluster_representatives = pd.DataFrame(representatives).sort_values(
        ["representative_priority_score", "cluster_size"],
        ascending=[False, False],
    )

    save_dataframe(
        cluster_representatives,
        config.REDUNDANCY_DIR / "cluster_representatives.csv",
    )

    logging.info("Stage 9 completed: Redundancy filtering.")

    return cluster_representatives


# ---------------------------------------------------------------------
# Stage 10 — Preliminary feature prioritization
# ---------------------------------------------------------------------

def stage10_preliminary_feature_prioritization() -> pd.DataFrame:
    """
    Create a preliminary feature-priority table and shortlist.

    This is not the final candidate list.

    Purpose:
    - Rank features for robustness testing.
    - Integrate primary paired evidence, FDR, direction consistency,
      bootstrap support, missingness, and redundancy representative status.

    Outputs:
    - preliminary_feature_priority_table.csv
    - preliminary_candidate_shortlist.csv
    """

    logging.info("Stage 10 started: Preliminary feature prioritization.")

    primary_results = pd.read_csv(
        config.PRIMARY_STATS_DIR / "primary_paired_feature_results_all_features_fdr.csv"
    )

    informative_features = pd.read_csv(
        config.FEATURE_DEFINITION_DIR / "informative_feature_list.csv"
    )

    redundancy_clusters = pd.read_csv(
        config.REDUNDANCY_DIR / "redundancy_clusters.csv"
    )

    cluster_representatives = pd.read_csv(
        config.REDUNDANCY_DIR / "cluster_representatives.csv"
    )

    priority = primary_results.merge(
        informative_features[
            [
                "feature",
                "missing_count",
                "missing_percent",
                "paired_patient_count",
                "missingness_flag",
            ]
        ],
        on="feature",
        how="left",
    )

    priority = priority.merge(
        redundancy_clusters[
            ["cluster_id", "feature", "cluster_size", "cluster_type"]
        ],
        on="feature",
        how="left",
    )

    representative_features = set(
        cluster_representatives["representative_feature"].dropna().astype(str)
    )

    priority["is_cluster_representative"] = priority["feature"].isin(
        representative_features
    ).astype(int)

    # Singleton features should also be considered representatives of themselves
    priority.loc[
        priority["cluster_type"] == "singleton",
        "is_cluster_representative",
    ] = 1

    # -----------------------------------------------------------------
    # Transparent preliminary evidence flags
    # -----------------------------------------------------------------

    priority["evidence_raw_p_0_05"] = (
        priority["wilcoxon_p"] < config.PRELIM_STRONG_RAW_P_THRESHOLD
    ).astype(int)

    priority["evidence_raw_p_0_10"] = (
        priority["wilcoxon_p"] < config.PRELIM_RAW_P_THRESHOLD
    ).astype(int)

    priority["evidence_direction_consistency"] = (
        priority["direction_consistency"]
        >= config.PRELIM_DIRECTION_CONSISTENCY_THRESHOLD
    ).astype(int)

    priority["evidence_strong_direction_consistency"] = (
        priority["direction_consistency"]
        >= config.PRELIM_STRONG_DIRECTION_CONSISTENCY_THRESHOLD
    ).astype(int)

    priority["evidence_effect_size"] = (
        priority["abs_rank_biserial_effect_size"]
        >= config.PRELIM_EFFECT_SIZE_THRESHOLD
    ).astype(int)

    priority["evidence_strong_effect_size"] = (
        priority["abs_rank_biserial_effect_size"]
        >= config.PRELIM_STRONG_EFFECT_SIZE_THRESHOLD
    ).astype(int)

    priority["evidence_bootstrap_ci_excludes_zero"] = pd.to_numeric(
        priority["bootstrap_ci_excludes_zero"],
        errors="coerce",
    ).fillna(0).astype(int)

    priority["evidence_family_fdr"] = pd.to_numeric(
        priority["family_fdr_significant"],
        errors="coerce",
    ).fillna(0).astype(int)

    priority["evidence_global_fdr"] = pd.to_numeric(
        priority["global_fdr_significant"],
        errors="coerce",
    ).fillna(0).astype(int)

    priority["evidence_redundancy_representative"] = priority[
        "is_cluster_representative"
    ].astype(int)

    # -----------------------------------------------------------------
    # Numeric preliminary priority score
    # -----------------------------------------------------------------

    priority["score_raw_p"] = min_max_normalize(
        safe_negative_log10_pvalue(priority["wilcoxon_p"]),
        higher_is_better=True,
    )

    priority["score_effect_size"] = min_max_normalize(
        priority["abs_rank_biserial_effect_size"],
        higher_is_better=True,
    )

    priority["score_direction_consistency"] = min_max_normalize(
        priority["direction_consistency"],
        higher_is_better=True,
    )

    priority["score_bootstrap"] = priority[
        "evidence_bootstrap_ci_excludes_zero"
    ].astype(float)

    priority["score_fdr"] = (
        2.0 * priority["evidence_global_fdr"]
        + 1.0 * priority["evidence_family_fdr"]
    )

    priority["score_redundancy"] = priority[
        "evidence_redundancy_representative"
    ].astype(float)

    priority["score_missingness"] = min_max_normalize(
        priority["missing_percent"],
        higher_is_better=False,
    )

    priority["preliminary_priority_score"] = (
        2.0 * priority["score_raw_p"]
        + 2.0 * priority["score_effect_size"]
        + 2.0 * priority["score_direction_consistency"]
        + 1.0 * priority["score_bootstrap"]
        + 1.0 * priority["score_fdr"]
        + 1.0 * priority["score_redundancy"]
        + 0.5 * priority["score_missingness"]
    )

    # -----------------------------------------------------------------
    # Preliminary tier labels
    # -----------------------------------------------------------------

    priority["preliminary_priority_tier"] = "not_prioritized"

    strict_like_mask = (
        (priority["evidence_raw_p_0_05"] == 1)
        & (priority["evidence_strong_direction_consistency"] == 1)
        & (priority["evidence_strong_effect_size"] == 1)
        & (priority["evidence_redundancy_representative"] == 1)
    )

    moderate_like_mask = (
        (priority["preliminary_priority_tier"] == "not_prioritized")
        & (priority["evidence_raw_p_0_10"] == 1)
        & (priority["evidence_direction_consistency"] == 1)
        & (priority["evidence_effect_size"] == 1)
        & (priority["evidence_redundancy_representative"] == 1)
    )

    exploratory_like_mask = (
        (priority["preliminary_priority_tier"] == "not_prioritized")
        & (
            (
                (priority["evidence_raw_p_0_10"] == 1)
                & (priority["evidence_effect_size"] == 1)
            )
            | (
                (priority["evidence_direction_consistency"] == 1)
                & (priority["evidence_bootstrap_ci_excludes_zero"] == 1)
            )
        )
    )

    priority.loc[strict_like_mask, "preliminary_priority_tier"] = (
        "high_priority_for_robustness_testing"
    )

    priority.loc[moderate_like_mask, "preliminary_priority_tier"] = (
        "moderate_priority_for_robustness_testing"
    )

    priority.loc[exploratory_like_mask, "preliminary_priority_tier"] = (
        "exploratory_priority_for_review"
    )

    # -----------------------------------------------------------------
    # Shortlist construction
    # -----------------------------------------------------------------
    # Shortlist is not final. It determines which features go to
    # sensitivity and leave-one-patient-out analysis.

    priority = priority.sort_values(
        [
            "preliminary_priority_score",
            "wilcoxon_p",
            "direction_consistency",
            "abs_rank_biserial_effect_size",
        ],
        ascending=[False, True, False, False],
    )

    priority["preliminary_rank"] = np.arange(1, len(priority) + 1)

    shortlist_candidate_mask = (
        priority["preliminary_priority_tier"].isin(
            [
                "high_priority_for_robustness_testing",
                "moderate_priority_for_robustness_testing",
                "exploratory_priority_for_review",
            ]
        )
        & (priority["is_cluster_representative"] == 1)
    )
    
    shortlist = priority[shortlist_candidate_mask].copy()

    # If shortlist is too small, fill with top-ranked representative features.
    if len(shortlist) < config.MIN_PRELIMINARY_SHORTLIST_FEATURES:
        needed = config.MIN_PRELIMINARY_SHORTLIST_FEATURES - len(shortlist)
    
        filler = priority[
            (~priority["feature"].isin(shortlist["feature"]))
            & (priority["is_cluster_representative"] == 1)
        ].sort_values(
            [
                "preliminary_priority_score",
                "wilcoxon_p",
                "direction_consistency",
                "abs_rank_biserial_effect_size",
            ],
            ascending=[False, True, False, False],
        ).head(needed)
    
        shortlist = pd.concat([shortlist, filler], ignore_index=True)

    # If shortlist is too large, keep the top N.
    shortlist = shortlist.sort_values(
        [
            "preliminary_priority_score",
            "wilcoxon_p",
            "direction_consistency",
            "abs_rank_biserial_effect_size",
        ],
        ascending=[False, True, False, False],
    ).head(config.MAX_PRELIMINARY_SHORTLIST_FEATURES)

    shortlist["included_in_preliminary_shortlist"] = 1

    priority["included_in_preliminary_shortlist"] = priority["feature"].isin(
        shortlist["feature"]
    ).astype(int)

    # Save outputs
    save_dataframe(
        priority,
        config.CANDIDATE_SCORECARD_DIR / "preliminary_feature_priority_table.csv",
    )

    save_dataframe(
        shortlist,
        config.CANDIDATE_SCORECARD_DIR / "preliminary_candidate_shortlist.csv",
    )

    logging.info(
        f"Stage 10 completed: Preliminary prioritization. "
        f"Shortlist size = {len(shortlist)}"
    )

    return shortlist

# ---------------------------------------------------------------------
# Stage 11 — HRV-valid subgroup analysis
# ---------------------------------------------------------------------

def stage11_hrv_valid_subgroup_analysis() -> pd.DataFrame:
    """
    Analyze HRV features in the HRV-valid patient subgroup.

    Purpose:
    - HRV features are especially sensitive to pacing, continuous AF,
      and R-peak detection outliers.
    - This subgroup analysis is for physiological interpretation only.
    - It is not treated as an independent validation cohort.

    HRV-invalid patients are identified using:
    - exclude_from_HRV_valid_subgroup, if available
    - hrv_invalid_flag, if available
    - paced_or_icd
    - continuous_AF
    - major_HRV_outlier

    Outputs:
    - hrv_valid_patient_inclusion_table.csv
    - hrv_valid_feature_results.csv
    """

    logging.info("Stage 11 started: HRV-valid subgroup analysis.")

    delta_df = pd.read_csv(
        config.PATIENT_LEVEL_DIR / "patient_feature_delta_table_all_features.csv"
    )

    patient_flags = pd.read_csv(
        config.CLINICAL_METADATA_DIR / "patient_clinical_flag_summary.csv"
    )

    primary_results = pd.read_csv(
        config.PRIMARY_STATS_DIR / "primary_paired_feature_results_all_features_fdr.csv"
    )

    feature_family_map = pd.read_csv(
        config.FEATURE_DEFINITION_DIR / "feature_family_map.csv"
    )

    # Standardize patient IDs
    delta_df[config.PATIENT_ID_COL] = standardize_patient_id_series(
        delta_df[config.PATIENT_ID_COL]
    )

    patient_flags[config.PATIENT_ID_COL] = standardize_patient_id_series(
        patient_flags[config.PATIENT_ID_COL]
    )

    # -----------------------------------------------------------------
    # Define HRV features
    # -----------------------------------------------------------------

    hrv_features = (
        feature_family_map.loc[
            feature_family_map["feature_family"] == "HRV",
            "feature",
        ]
        .dropna()
        .astype(str)
        .tolist()
    )

    hrv_features = [f for f in hrv_features if f in delta_df["feature"].unique()]

    if len(hrv_features) == 0:
        logging.warning("No HRV features found for HRV-valid subgroup analysis.")

    # -----------------------------------------------------------------
    # Build HRV inclusion/exclusion table
    # -----------------------------------------------------------------

    flag_cols = [
        "exclude_from_HRV_valid_subgroup",
        "hrv_invalid_flag",
        "paced_or_icd",
        "continuous_AF",
        "major_HRV_outlier",
        "intermittent_or_new_AF",
        "any_af",
        "BBB",
        "mechanical_support",
        "clean_sinus_candidate",
        "confound_count",
        "notes",
    ]

    flag_cols_present = [c for c in flag_cols if c in patient_flags.columns]

    hrv_inclusion = patient_flags[
        [config.PATIENT_ID_COL] + flag_cols_present
    ].copy()

    # Ensure binary columns are numeric if present
    binary_cols = [
        "exclude_from_HRV_valid_subgroup",
        "hrv_invalid_flag",
        "paced_or_icd",
        "continuous_AF",
        "major_HRV_outlier",
        "intermittent_or_new_AF",
        "any_af",
        "BBB",
        "mechanical_support",
        "clean_sinus_candidate",
    ]

    for col in binary_cols:
        if col in hrv_inclusion.columns:
            hrv_inclusion[col] = pd.to_numeric(
                hrv_inclusion[col],
                errors="coerce",
            ).fillna(0).astype(int)

    # Conservative HRV-invalid rule.
    # Any one of these makes HRV interpretation unreliable.
    invalid_sources = []

    for col in [
        "exclude_from_HRV_valid_subgroup",
        "hrv_invalid_flag",
        "paced_or_icd",
        "continuous_AF",
        "major_HRV_outlier",
    ]:
        if col in hrv_inclusion.columns:
            invalid_sources.append(col)

    if invalid_sources:
        hrv_inclusion["hrv_invalid_computed"] = (
            hrv_inclusion[invalid_sources].sum(axis=1) > 0
        ).astype(int)
    else:
        hrv_inclusion["hrv_invalid_computed"] = 0

    hrv_inclusion["hrv_valid_subgroup"] = (
        hrv_inclusion["hrv_invalid_computed"] == 0
    ).astype(int)

    def _hrv_exclusion_reason(row):
        reasons = []

        if row.get("exclude_from_HRV_valid_subgroup", 0) == 1:
            reasons.append("curated_exclude_from_HRV_valid_subgroup")

        if row.get("paced_or_icd", 0) == 1:
            reasons.append("paced_or_icd")

        if row.get("continuous_AF", 0) == 1:
            reasons.append("continuous_AF")

        if row.get("major_HRV_outlier", 0) == 1:
            reasons.append("major_HRV_outlier")

        if row.get("hrv_invalid_flag", 0) == 1:
            reasons.append("hrv_invalid_flag")

        if not reasons:
            return ""

        # Remove duplicates while preserving order
        unique_reasons = []
        for reason in reasons:
            if reason not in unique_reasons:
                unique_reasons.append(reason)

        return "; ".join(unique_reasons)

    hrv_inclusion["hrv_exclusion_reason"] = hrv_inclusion.apply(
        _hrv_exclusion_reason,
        axis=1,
    )

    hrv_inclusion = hrv_inclusion.sort_values(
        [config.PATIENT_ID_COL]
    ).reset_index(drop=True)

    save_dataframe(
        hrv_inclusion,
        config.HRV_SUBGROUP_DIR / "hrv_valid_patient_inclusion_table.csv",
    )

    # -----------------------------------------------------------------
    # Run paired statistics for HRV features in HRV-valid subgroup
    # -----------------------------------------------------------------

    hrv_valid_patient_ids = set(
        hrv_inclusion.loc[
            hrv_inclusion["hrv_valid_subgroup"] == 1,
            config.PATIENT_ID_COL,
        ].astype(str)
    )

    results = []

    for feature in hrv_features:
        sub = delta_df[
            (delta_df["feature"] == feature)
            & (delta_df[config.PATIENT_ID_COL].astype(str).isin(hrv_valid_patient_ids))
            & (delta_df["paired_observation_available"] == 1)
        ].copy()

        n_paired = len(sub)

        n_increase = int((sub["observed_direction"] == "increase").sum())
        n_decrease = int((sub["observed_direction"] == "decrease").sum())
        n_no_change = int((sub["observed_direction"] == "no_change").sum())

        dominant_direction = dominant_direction_from_counts(
            n_increase,
            n_decrease,
            n_no_change,
        )

        max_direction_count = (
            max(n_increase, n_decrease, n_no_change) if n_paired > 0 else 0
        )

        direction_consistency = (
            max_direction_count / n_paired if n_paired > 0 else np.nan
        )

        if n_paired >= 2:
            p_value = wilcoxon_signed_rank_pvalue(
                sub["PreHI_median"],
                sub["HI_median"],
            )

            effect_size = matched_pairs_rank_biserial_effect_size(
                sub["PreHI_median"],
                sub["HI_median"],
            )

            ci_low, ci_high = bootstrap_median_ci(
                sub["delta_HI_minus_PreHI"],
                n_iterations=config.BOOTSTRAP_N_ITERATIONS,
                confidence_level=config.BOOTSTRAP_CONFIDENCE_LEVEL,
                random_seed=config.RANDOM_SEED,
            )

        else:
            p_value = np.nan
            effect_size = np.nan
            ci_low = np.nan
            ci_high = np.nan

        bootstrap_ci_excludes_zero = (
            pd.notna(ci_low)
            and pd.notna(ci_high)
            and ((ci_low > 0 and ci_high > 0) or (ci_low < 0 and ci_high < 0))
        )

        # Pull full-cohort result for comparison
        full_match = primary_results[primary_results["feature"] == feature]

        if not full_match.empty:
            full_row = full_match.iloc[0]
            full_n = full_row.get("n_paired_patients", np.nan)
            full_p = full_row.get("wilcoxon_p", np.nan)
            full_direction = full_row.get("dominant_direction", "")
            full_direction_consistency = full_row.get(
                "direction_consistency", np.nan
            )
            full_effect = full_row.get("rank_biserial_effect_size", np.nan)
        else:
            full_n = np.nan
            full_p = np.nan
            full_direction = ""
            full_direction_consistency = np.nan
            full_effect = np.nan

        results.append(
            {
                "feature": feature,
                "feature_family": "HRV",
                "analysis_group": "HRV_valid_subgroup",
                "n_hrv_valid_patients_total": len(hrv_valid_patient_ids),
                "n_paired_patients": n_paired,
                "median_PreHI": sub["PreHI_median"].median(),
                "median_HI": sub["HI_median"].median(),
                "median_delta_HI_minus_PreHI": sub[
                    "delta_HI_minus_PreHI"
                ].median(),
                "median_percent_change": sub["percent_change"].median(),
                "n_increase": n_increase,
                "n_decrease": n_decrease,
                "n_no_change": n_no_change,
                "dominant_direction": dominant_direction,
                "direction_consistency": direction_consistency,
                "wilcoxon_p": p_value,
                "rank_biserial_effect_size": effect_size,
                "bootstrap_median_delta_ci_low": ci_low,
                "bootstrap_median_delta_ci_high": ci_high,
                "bootstrap_ci_excludes_zero": int(bootstrap_ci_excludes_zero),
                "full_cohort_n_paired_patients": full_n,
                "full_cohort_wilcoxon_p": full_p,
                "full_cohort_dominant_direction": full_direction,
                "full_cohort_direction_consistency": full_direction_consistency,
                "full_cohort_rank_biserial_effect_size": full_effect,
                "direction_agrees_with_full_cohort": int(
                    dominant_direction == full_direction
                    if dominant_direction not in ["tie_or_mixed", "missing"]
                    and full_direction not in ["tie_or_mixed", "missing", ""]
                    else 0
                ),
                "interpretation_note": (
                    "HRV-valid subgroup analysis is descriptive/supportive only; "
                    "it is not an independent validation cohort."
                ),
            }
        )

    hrv_results = pd.DataFrame(results)

    if not hrv_results.empty:
        hrv_results = hrv_results.sort_values(
            [
                "wilcoxon_p",
                "direction_consistency",
                "rank_biserial_effect_size",
            ],
            ascending=[True, False, False],
        )

    save_dataframe(
        hrv_results,
        config.HRV_SUBGROUP_DIR / "hrv_valid_feature_results.csv",
    )

    # -----------------------------------------------------------------
    # Summary logging
    # -----------------------------------------------------------------

    n_valid = int(hrv_inclusion["hrv_valid_subgroup"].sum())
    n_total = int(len(hrv_inclusion))
    logging.info(
        f"Stage 11 completed: HRV-valid subgroup analysis. "
        f"HRV-valid patients: {n_valid}/{n_total}."
    )

    return hrv_results

# ---------------------------------------------------------------------
# Stage 12 — Supportive mixed-effects modelling
# ---------------------------------------------------------------------

def stage12_supportive_mixed_effects_modelling() -> pd.DataFrame:
    """
    Fit supportive segment-level mixed-effects models for all informative features.

    Model:
        feature_z ~ condition_binary + (1 | Patient_id)

    Purpose:
    - This analysis uses segment-level observations while accounting for
      patient-level clustering using a random intercept.
    - It is supportive only.
    - It does not replace the primary patient-level paired analysis.

    Outputs:
    - lmm_supportive_results_all_features.csv
    """

    logging.info("Stage 12 started: Supportive mixed-effects modelling.")

    merged = pd.read_csv(
        config.CLINICAL_METADATA_DIR / "merged_segment_clinical_dataset.csv"
    )

    informative_features = pd.read_csv(
        config.FEATURE_DEFINITION_DIR / "informative_feature_list.csv"
    )

    primary_results = pd.read_csv(
        config.PRIMARY_STATS_DIR / "primary_paired_feature_results_all_features_fdr.csv"
    )

    feature_names = informative_features["feature"].tolist()

    merged[config.PATIENT_ID_COL] = standardize_patient_id_series(
        merged[config.PATIENT_ID_COL]
    )

    # -----------------------------------------------------------------
    # Define condition_binary
    # -----------------------------------------------------------------

    if config.CONDITION_STD_COL not in merged.columns:
        raise ValueError(
            f"{config.CONDITION_STD_COL} is missing from merged dataset."
        )

    condition_map = {
        "PreHI": 0,
        "HI": 1,
    }

    merged["condition_binary"] = merged[config.CONDITION_STD_COL].map(condition_map)

    if merged["condition_binary"].isna().any():
        bad_conditions = (
            merged.loc[merged["condition_binary"].isna(), config.CONDITION_STD_COL]
            .dropna()
            .unique()
            .tolist()
        )
        raise ValueError(
            f"Unmapped Condition_std values for LMM: {bad_conditions}"
        )

    merged["condition_binary"] = merged["condition_binary"].astype(int)

    # -----------------------------------------------------------------
    # Fit LMM feature by feature
    # -----------------------------------------------------------------

    results = []

    for feature in feature_names:
        if feature not in merged.columns:
            results.append(
                {
                    "feature": feature,
                    "feature_family": infer_feature_family(feature),
                    "lmm_status": "skipped_feature_missing",
                    "n_segments": 0,
                    "n_patients": 0,
                    "beta_condition_HI": np.nan,
                    "p_condition_HI": np.nan,
                    "converged": 0,
                    "lmm_direction": "",
                    "primary_dominant_direction": "",
                    "lmm_agrees_with_primary_direction": np.nan,
                    "notes": "Feature not present in merged segment dataset.",
                }
            )
            continue

        model_df = merged[
            [
                config.PATIENT_ID_COL,
                "condition_binary",
                feature,
            ]
        ].copy()

        model_df[feature] = pd.to_numeric(model_df[feature], errors="coerce")
        model_df = model_df.dropna(subset=[feature, "condition_binary", config.PATIENT_ID_COL])

        n_segments = len(model_df)
        n_patients = model_df[config.PATIENT_ID_COL].nunique()

        if n_segments < config.LMM_MIN_SEGMENTS or n_patients < config.LMM_MIN_PATIENTS:
            results.append(
                {
                    "feature": feature,
                    "feature_family": infer_feature_family(feature),
                    "lmm_status": "skipped_insufficient_data",
                    "n_segments": n_segments,
                    "n_patients": n_patients,
                    "beta_condition_HI": np.nan,
                    "p_condition_HI": np.nan,
                    "converged": 0,
                    "lmm_direction": "",
                    "primary_dominant_direction": "",
                    "lmm_agrees_with_primary_direction": np.nan,
                    "notes": (
                        f"Skipped because n_segments={n_segments}, "
                        f"n_patients={n_patients}."
                    ),
                }
            )
            continue

        feature_std = model_df[feature].std(skipna=True)

        if pd.isna(feature_std) or feature_std <= config.LMM_MIN_FEATURE_STD:
            results.append(
                {
                    "feature": feature,
                    "feature_family": infer_feature_family(feature),
                    "lmm_status": "skipped_near_constant_feature",
                    "n_segments": n_segments,
                    "n_patients": n_patients,
                    "beta_condition_HI": np.nan,
                    "p_condition_HI": np.nan,
                    "converged": 0,
                    "lmm_direction": "",
                    "primary_dominant_direction": "",
                    "lmm_agrees_with_primary_direction": np.nan,
                    "notes": (
                        f"Skipped because feature standard deviation is {feature_std}."
                    ),
                }
            )
            continue

        model_df["feature_z"] = zscore_series(model_df[feature])
        model_df = model_df.dropna(subset=["feature_z"])

        try:
            model = smf.mixedlm(
                "feature_z ~ condition_binary",
                data=model_df,
                groups=model_df[config.PATIENT_ID_COL],
            )

            fit = model.fit(
                reml=False,
                method="lbfgs",
                maxiter=200,
                disp=False,
            )

            beta = fit.params.get("condition_binary", np.nan)
            p_value = fit.pvalues.get("condition_binary", np.nan)

            converged = int(getattr(fit, "converged", False))

            if pd.isna(beta):
                lmm_direction = ""
            elif beta > 0:
                lmm_direction = "increase"
            elif beta < 0:
                lmm_direction = "decrease"
            else:
                lmm_direction = "no_change"

            lmm_status = "fit_success" if converged == 1 else "fit_not_converged"
            notes = ""

        except Exception as e:
            beta = np.nan
            p_value = np.nan
            converged = 0
            lmm_direction = ""
            lmm_status = "fit_failed"
            notes = str(e)

        # Pull primary paired direction for comparison
        primary_match = primary_results[primary_results["feature"] == feature]

        if not primary_match.empty:
            primary_row = primary_match.iloc[0]
            primary_direction = primary_row.get("dominant_direction", "")
            primary_p = primary_row.get("wilcoxon_p", np.nan)
            primary_effect = primary_row.get("rank_biserial_effect_size", np.nan)
            primary_direction_consistency = primary_row.get(
                "direction_consistency", np.nan
            )
        else:
            primary_direction = ""
            primary_p = np.nan
            primary_effect = np.nan
            primary_direction_consistency = np.nan

        if lmm_direction in ["increase", "decrease", "no_change"] and primary_direction in [
            "increase",
            "decrease",
            "no_change",
        ]:
            agrees = int(lmm_direction == primary_direction)
        else:
            agrees = np.nan

        results.append(
            {
                "feature": feature,
                "feature_family": infer_feature_family(feature),
                "lmm_status": lmm_status,
                "n_segments": n_segments,
                "n_patients": n_patients,
                "beta_condition_HI": beta,
                "p_condition_HI": p_value,
                "converged": converged,
                "lmm_direction": lmm_direction,
                "primary_dominant_direction": primary_direction,
                "lmm_agrees_with_primary_direction": agrees,
                "primary_wilcoxon_p": primary_p,
                "primary_rank_biserial_effect_size": primary_effect,
                "primary_direction_consistency": primary_direction_consistency,
                "notes": notes,
            }
        )

    lmm_results = pd.DataFrame(results)

    # FDR correction across successfully fit LMM p-values
    if "p_condition_HI" in lmm_results.columns:
        lmm_results["p_condition_HI_fdr"] = apply_fdr_correction(
            lmm_results["p_condition_HI"],
            method=config.FDR_METHOD,
        )
    else:
        lmm_results["p_condition_HI_fdr"] = np.nan

    lmm_results["lmm_raw_p_less_0_05"] = (
        lmm_results["p_condition_HI"] < 0.05
    ).fillna(False).astype(int)

    lmm_results["lmm_fdr_significant"] = (
        lmm_results["p_condition_HI_fdr"] < config.FDR_ALPHA
    ).fillna(False).astype(int)

    lmm_results = lmm_results.sort_values(
        [
            "p_condition_HI",
            "lmm_agrees_with_primary_direction",
            "primary_wilcoxon_p",
        ],
        ascending=[True, False, True],
    )

    save_dataframe(
        lmm_results,
        config.MIXED_EFFECTS_DIR / "lmm_supportive_results_all_features.csv",
    )

    logging.info("Stage 12 completed: Supportive mixed-effects modelling.")

    return lmm_results

# ---------------------------------------------------------------------
# Stage 13 — Sensitivity analysis
# ---------------------------------------------------------------------

def stage13_sensitivity_analysis() -> pd.DataFrame:
    """
    Recompute patient-level paired results for shortlisted features under
    clinically and signal-quality motivated sensitivity scenarios.

    This stage uses only the preliminary candidate shortlist from Stage 10.

    Scenarios:
    - full_cohort_baseline
    - qc_clean_segments
    - exclude_high_confound
    - exclude_paced_or_bbb
    - exclude_any_af
    - exclude_major_hrv_outlier
    - exclude_mechanical_support
    - exclude_atypical_hi

    Outputs:
    - sensitivity_scenario_definitions.csv
    - sensitivity_patient_counts.csv
    - sensitivity_analysis_by_scenario_recomputed.csv
    - sensitivity_vs_full_cohort_summary.csv
    - candidate_feature_robustness_matrix.csv
    """

    logging.info("Stage 13 started: Sensitivity analysis.")

    merged = pd.read_csv(
        config.CLINICAL_METADATA_DIR / "merged_segment_clinical_dataset.csv"
    )

    qc = pd.read_csv(
        config.SEGMENT_QC_DIR / "segment_qc_flags.csv"
    )

    shortlist = pd.read_csv(
        config.CANDIDATE_SCORECARD_DIR / "preliminary_candidate_shortlist.csv"
    )

    full_primary = pd.read_csv(
        config.PRIMARY_STATS_DIR / "primary_paired_feature_results_all_features_fdr.csv"
    )

    shortlisted_features = shortlist["feature"].dropna().astype(str).tolist()

    merged[config.PATIENT_ID_COL] = standardize_patient_id_series(
        merged[config.PATIENT_ID_COL]
    )

    qc[config.PATIENT_ID_COL] = standardize_patient_id_series(
        qc[config.PATIENT_ID_COL]
    )

    # -----------------------------------------------------------------
    # Merge QC flags using stable row ID
    # -----------------------------------------------------------------

    if config.ROW_ID_COL not in merged.columns:
        raise ValueError(
            f"{config.ROW_ID_COL} missing from merged segment dataset. "
            "Please rerun Stage 4."
        )

    if config.ROW_ID_COL not in qc.columns:
        raise ValueError(
            f"{config.ROW_ID_COL} missing from QC flags. "
            "Please rerun Stage 5."
        )

    qc_cols_to_merge = [
        config.ROW_ID_COL,
        "clear_segment_qc_exclusion_flag",
        "qc_clear_segment_flag",
        "near_zero_signal_flag",
        "low_variability_signal_flag",
        "qrs_std_zero_flag",
        "known_artifact_segment_flag",
        "extreme_value_review_flag",
        "extreme_kurtosis_review_flag",
    ]

    qc_cols_to_merge = [c for c in qc_cols_to_merge if c in qc.columns]

    merged_qc = merged.merge(
        qc[qc_cols_to_merge],
        on=config.ROW_ID_COL,
        how="left",
        validate="one_to_one",
    )

    if "clear_segment_qc_exclusion_flag" not in merged_qc.columns:
        merged_qc["clear_segment_qc_exclusion_flag"] = 0

    # Ensure clinical flags exist
    required_flag_defaults = {
        "high_confound": 0,
        "paced_or_bbb": 0,
        "any_af": 0,
        "major_HRV_outlier": 0,
        "mechanical_support": 0,
        "atypical_HI": 0,
    }

    for col, default in required_flag_defaults.items():
        if col not in merged_qc.columns:
            logging.warning(
                f"{col} not found in merged dataset. Defaulting to {default}."
            )
            merged_qc[col] = default

        merged_qc[col] = pd.to_numeric(
            merged_qc[col],
            errors="coerce",
        ).fillna(default).astype(int)

    # Ensure features numeric
    for feature in shortlisted_features:
        merged_qc[feature] = pd.to_numeric(merged_qc[feature], errors="coerce")

    # -----------------------------------------------------------------
    # Scenario definitions
    # -----------------------------------------------------------------

    scenario_definitions = pd.DataFrame(
        [
            {
                "scenario": "full_cohort_baseline",
                "segment_filter": "none",
                "patient_filter": "none",
                "purpose": "Reference analysis using all available non-missing feature values.",
            },
            {
                "scenario": "qc_clean_segments",
                "segment_filter": "clear_segment_qc_exclusion_flag == 0",
                "patient_filter": "none",
                "purpose": "Tests whether findings survive removal of clear segment-level QC issues.",
            },
            {
                "scenario": "exclude_high_confound",
                "segment_filter": "none",
                "patient_filter": "high_confound == 0",
                "purpose": "Tests whether findings are driven by heavily confounded patients.",
            },
            {
                "scenario": "exclude_paced_or_bbb",
                "segment_filter": "none",
                "patient_filter": "paced_or_bbb == 0",
                "purpose": "Tests sensitivity to pacing/BBB, especially morphology/QRS features.",
            },
            {
                "scenario": "exclude_any_af",
                "segment_filter": "none",
                "patient_filter": "any_af == 0",
                "purpose": "Tests sensitivity to AF/rhythm irregularity.",
            },
            {
                "scenario": "exclude_major_hrv_outlier",
                "segment_filter": "none",
                "patient_filter": "major_HRV_outlier == 0",
                "purpose": "Tests sensitivity to major HRV/R-peak outlier patients.",
            },
            {
                "scenario": "exclude_mechanical_support",
                "segment_filter": "none",
                "patient_filter": "mechanical_support == 0",
                "purpose": "Tests sensitivity to mechanical support/device-related influence.",
            },
            {
                "scenario": "exclude_atypical_hi",
                "segment_filter": "none",
                "patient_filter": "atypical_HI == 0",
                "purpose": "Tests sensitivity to atypical HI mechanisms.",
            },
        ]
    )

    save_dataframe(
        scenario_definitions,
        config.SENSITIVITY_DIR / "sensitivity_scenario_definitions.csv",
    )

    # -----------------------------------------------------------------
    # Helper function for one scenario and one feature
    # -----------------------------------------------------------------

    def _apply_scenario_filter(df: pd.DataFrame, scenario: str) -> pd.DataFrame:
        out = df.copy()

        if scenario == "full_cohort_baseline":
            return out

        if scenario == "qc_clean_segments":
            return out[out["clear_segment_qc_exclusion_flag"] == 0].copy()

        if scenario == "exclude_high_confound":
            return out[out["high_confound"] == 0].copy()

        if scenario == "exclude_paced_or_bbb":
            return out[out["paced_or_bbb"] == 0].copy()

        if scenario == "exclude_any_af":
            return out[out["any_af"] == 0].copy()

        if scenario == "exclude_major_hrv_outlier":
            return out[out["major_HRV_outlier"] == 0].copy()

        if scenario == "exclude_mechanical_support":
            return out[out["mechanical_support"] == 0].copy()

        if scenario == "exclude_atypical_hi":
            return out[out["atypical_HI"] == 0].copy()

        raise ValueError(f"Unknown sensitivity scenario: {scenario}")

    def _compute_feature_result(
        df: pd.DataFrame,
        feature: str,
        scenario: str,
    ) -> dict:
        scenario_df = _apply_scenario_filter(df, scenario)

        # Segment and patient counts after filtering
        n_segments = len(scenario_df)
        n_patients_total = scenario_df[config.PATIENT_ID_COL].nunique()

        # Patient-condition medians
        medians = (
            scenario_df
            .groupby([config.PATIENT_ID_COL, config.CONDITION_STD_COL])[feature]
            .median()
            .reset_index()
        )

        wide = medians.pivot(
            index=config.PATIENT_ID_COL,
            columns=config.CONDITION_STD_COL,
            values=feature,
        )

        if "PreHI" not in wide.columns:
            wide["PreHI"] = np.nan

        if "HI" not in wide.columns:
            wide["HI"] = np.nan

        paired = wide[["PreHI", "HI"]].dropna().copy()
        paired["delta_HI_minus_PreHI"] = paired["HI"] - paired["PreHI"]
        paired["observed_direction"] = paired["delta_HI_minus_PreHI"].apply(
            direction_from_delta
        )

        n_paired = len(paired)

        n_increase = int((paired["observed_direction"] == "increase").sum())
        n_decrease = int((paired["observed_direction"] == "decrease").sum())
        n_no_change = int((paired["observed_direction"] == "no_change").sum())

        dominant_direction = dominant_direction_from_counts(
            n_increase,
            n_decrease,
            n_no_change,
        )

        max_direction_count = (
            max(n_increase, n_decrease, n_no_change) if n_paired > 0 else 0
        )

        direction_consistency = (
            max_direction_count / n_paired if n_paired > 0 else np.nan
        )

        if n_paired >= 2:
            p_value = wilcoxon_signed_rank_pvalue(
                paired["PreHI"],
                paired["HI"],
            )

            effect_size = matched_pairs_rank_biserial_effect_size(
                paired["PreHI"],
                paired["HI"],
            )

            ci_low, ci_high = bootstrap_median_ci(
                paired["delta_HI_minus_PreHI"],
                n_iterations=config.BOOTSTRAP_N_ITERATIONS,
                confidence_level=config.BOOTSTRAP_CONFIDENCE_LEVEL,
                random_seed=config.RANDOM_SEED,
            )

        else:
            p_value = np.nan
            effect_size = np.nan
            ci_low = np.nan
            ci_high = np.nan

        bootstrap_ci_excludes_zero = (
            pd.notna(ci_low)
            and pd.notna(ci_high)
            and ((ci_low > 0 and ci_high > 0) or (ci_low < 0 and ci_high < 0))
        )

        if n_paired < config.SENSITIVITY_MIN_PAIRED_PATIENTS:
            scenario_status = "limited_sample_size"
        else:
            scenario_status = "adequate_sample_size"

        return {
            "scenario": scenario,
            "feature": feature,
            "feature_family": infer_feature_family(feature),
            "n_segments": n_segments,
            "n_patients_total_after_filter": n_patients_total,
            "n_paired_patients": n_paired,
            "median_PreHI": paired["PreHI"].median(),
            "median_HI": paired["HI"].median(),
            "median_delta_HI_minus_PreHI": paired[
                "delta_HI_minus_PreHI"
            ].median(),
            "n_increase": n_increase,
            "n_decrease": n_decrease,
            "n_no_change": n_no_change,
            "dominant_direction": dominant_direction,
            "direction_consistency": direction_consistency,
            "wilcoxon_p": p_value,
            "rank_biserial_effect_size": effect_size,
            "bootstrap_median_delta_ci_low": ci_low,
            "bootstrap_median_delta_ci_high": ci_high,
            "bootstrap_ci_excludes_zero": int(bootstrap_ci_excludes_zero),
            "scenario_status": scenario_status,
        }

    # -----------------------------------------------------------------
    # Recompute all scenario-feature results
    # -----------------------------------------------------------------

    results = []

    for scenario in config.SENSITIVITY_SCENARIOS:
        for feature in shortlisted_features:
            results.append(
                _compute_feature_result(
                    merged_qc,
                    feature,
                    scenario,
                )
            )

    sensitivity_results = pd.DataFrame(results)

    save_dataframe(
        sensitivity_results,
        config.SENSITIVITY_DIR / "sensitivity_analysis_by_scenario_recomputed.csv",
    )

    # -----------------------------------------------------------------
    # Patient counts by scenario
    # -----------------------------------------------------------------

    patient_count_rows = []

    for scenario in config.SENSITIVITY_SCENARIOS:
        scenario_df = _apply_scenario_filter(merged_qc, scenario)

        patient_condition_presence = (
            scenario_df
            .groupby([config.PATIENT_ID_COL, config.CONDITION_STD_COL])
            .size()
            .unstack(fill_value=0)
        )

        for col in ["PreHI", "HI"]:
            if col not in patient_condition_presence.columns:
                patient_condition_presence[col] = 0

        n_patients_total = scenario_df[config.PATIENT_ID_COL].nunique()
        n_patients_with_prehi = int((patient_condition_presence["PreHI"] > 0).sum())
        n_patients_with_hi = int((patient_condition_presence["HI"] > 0).sum())
        n_patients_with_both = int(
            (
                (patient_condition_presence["PreHI"] > 0)
                & (patient_condition_presence["HI"] > 0)
            ).sum()
        )

        patient_count_rows.append(
            {
                "scenario": scenario,
                "n_segments": len(scenario_df),
                "n_patients_total": n_patients_total,
                "n_patients_with_PreHI": n_patients_with_prehi,
                "n_patients_with_HI": n_patients_with_hi,
                "n_patients_with_both_conditions": n_patients_with_both,
            }
        )

    sensitivity_patient_counts = pd.DataFrame(patient_count_rows)

    save_dataframe(
        sensitivity_patient_counts,
        config.SENSITIVITY_DIR / "sensitivity_patient_counts.csv",
    )

    # -----------------------------------------------------------------
    # Compare each scenario against full cohort baseline
    # -----------------------------------------------------------------

    baseline = sensitivity_results[
        sensitivity_results["scenario"] == "full_cohort_baseline"
    ][
        [
            "feature",
            "dominant_direction",
            "median_delta_HI_minus_PreHI",
            "direction_consistency",
            "wilcoxon_p",
            "rank_biserial_effect_size",
        ]
    ].copy()

    baseline = baseline.rename(
        columns={
            "dominant_direction": "baseline_dominant_direction",
            "median_delta_HI_minus_PreHI": "baseline_median_delta_HI_minus_PreHI",
            "direction_consistency": "baseline_direction_consistency",
            "wilcoxon_p": "baseline_wilcoxon_p",
            "rank_biserial_effect_size": "baseline_rank_biserial_effect_size",
        }
    )

    comparison = sensitivity_results.merge(
        baseline,
        on="feature",
        how="left",
    )

    comparison["direction_matches_baseline"] = (
        comparison["dominant_direction"]
        == comparison["baseline_dominant_direction"]
    ).astype(int)

    comparison.loc[
        comparison["scenario"] == "full_cohort_baseline",
        "direction_matches_baseline",
    ] = 1

    comparison["median_delta_sign_matches_baseline"] = (
        np.sign(comparison["median_delta_HI_minus_PreHI"])
        == np.sign(comparison["baseline_median_delta_HI_minus_PreHI"])
    ).astype(int)

    comparison.loc[
        comparison["scenario"] == "full_cohort_baseline",
        "median_delta_sign_matches_baseline",
    ] = 1

    comparison["p_value_ratio_vs_baseline"] = (
        comparison["wilcoxon_p"] / comparison["baseline_wilcoxon_p"]
    )

    save_dataframe(
        comparison,
        config.SENSITIVITY_DIR / "sensitivity_vs_full_cohort_summary.csv",
    )

    # -----------------------------------------------------------------
    # Robustness matrix by candidate feature
    # -----------------------------------------------------------------

    robustness_rows = []

    for feature, sub in comparison.groupby("feature"):
        non_baseline = sub[sub["scenario"] != "full_cohort_baseline"].copy()

        n_scenarios = len(non_baseline)
        n_adequate = int(
            (non_baseline["scenario_status"] == "adequate_sample_size").sum()
        )

        n_direction_retained = int(
            (
                (non_baseline["direction_matches_baseline"] == 1)
                & (non_baseline["scenario_status"] == "adequate_sample_size")
            ).sum()
        )

        n_delta_sign_retained = int(
            (
                (non_baseline["median_delta_sign_matches_baseline"] == 1)
                & (non_baseline["scenario_status"] == "adequate_sample_size")
            ).sum()
        )

        n_raw_p_less_0_05 = int(
            (
                (non_baseline["wilcoxon_p"] < 0.05)
                & (non_baseline["scenario_status"] == "adequate_sample_size")
            ).sum()
        )

        n_bootstrap_ci_excludes_zero = int(
            (
                (non_baseline["bootstrap_ci_excludes_zero"] == 1)
                & (non_baseline["scenario_status"] == "adequate_sample_size")
            ).sum()
        )

        median_direction_consistency = non_baseline.loc[
            non_baseline["scenario_status"] == "adequate_sample_size",
            "direction_consistency",
        ].median()

        if n_adequate > 0:
            prop_direction_retained = n_direction_retained / n_adequate
            prop_delta_sign_retained = n_delta_sign_retained / n_adequate
        else:
            prop_direction_retained = np.nan
            prop_delta_sign_retained = np.nan

        baseline_row = sub[sub["scenario"] == "full_cohort_baseline"].iloc[0]

        # Descriptive robustness label.
        if (
            n_adequate >= 4
            and prop_direction_retained >= 0.80
            and prop_delta_sign_retained >= 0.80
            and median_direction_consistency
            >= config.SENSITIVITY_DIRECTION_CONSISTENCY_REASONABLE
            and (
                n_raw_p_less_0_05 >= 2
                or n_bootstrap_ci_excludes_zero >= 2
            )
        ):
            robustness_label = "robust"

        elif (
            n_adequate >= 3
            and prop_direction_retained >= 0.60
            and prop_delta_sign_retained >= 0.60
        ):
            robustness_label = "moderate"

        else:
            robustness_label = "weak_or_exploratory"

        robustness_rows.append(
            {
                "feature": feature,
                "feature_family": infer_feature_family(feature),
                "baseline_dominant_direction": baseline_row[
                    "baseline_dominant_direction"
                ],
                "baseline_median_delta_HI_minus_PreHI": baseline_row[
                    "baseline_median_delta_HI_minus_PreHI"
                ],
                "baseline_wilcoxon_p": baseline_row["baseline_wilcoxon_p"],
                "baseline_rank_biserial_effect_size": baseline_row[
                    "baseline_rank_biserial_effect_size"
                ],
                "n_sensitivity_scenarios_non_baseline": n_scenarios,
                "n_scenarios_with_adequate_sample_size": n_adequate,
                "n_direction_retained": n_direction_retained,
                "prop_direction_retained": prop_direction_retained,
                "n_delta_sign_retained": n_delta_sign_retained,
                "prop_delta_sign_retained": prop_delta_sign_retained,
                "median_direction_consistency_across_scenarios": median_direction_consistency,
                "n_scenarios_raw_p_less_0_05": n_raw_p_less_0_05,
                "n_scenarios_bootstrap_ci_excludes_zero": n_bootstrap_ci_excludes_zero,
                "sensitivity_robustness_label": robustness_label,
            }
        )

    robustness_matrix = pd.DataFrame(robustness_rows).sort_values(
        [
            "sensitivity_robustness_label",
            "prop_direction_retained",
            "n_scenarios_raw_p_less_0_05",
            "n_scenarios_bootstrap_ci_excludes_zero",
        ],
        ascending=[True, False, False, False],
    )

    # Make label sort more intuitive
    label_order = {
        "robust": 1,
        "moderate": 2,
        "weak_or_exploratory": 3,
    }

    robustness_matrix["robustness_sort_order"] = robustness_matrix[
        "sensitivity_robustness_label"
    ].map(label_order)

    robustness_matrix = robustness_matrix.sort_values(
        [
            "robustness_sort_order",
            "prop_direction_retained",
            "n_scenarios_raw_p_less_0_05",
            "n_scenarios_bootstrap_ci_excludes_zero",
        ],
        ascending=[True, False, False, False],
    ).drop(columns=["robustness_sort_order"])

    save_dataframe(
        robustness_matrix,
        config.SENSITIVITY_DIR / "candidate_feature_robustness_matrix.csv",
    )

    logging.info("Stage 13 completed: Sensitivity analysis.")

    return sensitivity_results

# ---------------------------------------------------------------------
# Stage 14 — Leave-one-patient-out influence analysis
# ---------------------------------------------------------------------

def stage14_leave_one_patient_out_analysis() -> pd.DataFrame:
    """
    Perform leave-one-patient-out influence analysis for preliminary shortlisted features.

    Purpose:
    - Determine whether each candidate feature depends heavily on one patient.
    - Recompute patient-level paired statistics after removing one patient at a time.
    - Identify whether dominant direction changes, statistical support weakens,
      or direction consistency drops substantially.

    Outputs:
    - leave_one_patient_out_results.csv
    - leave_one_patient_out_summary.csv
    """

    logging.info("Stage 14 started: Leave-one-patient-out influence analysis.")

    delta_df = pd.read_csv(
        config.PATIENT_LEVEL_DIR / "patient_feature_delta_table_all_features.csv"
    )

    shortlist = pd.read_csv(
        config.CANDIDATE_SCORECARD_DIR / "preliminary_candidate_shortlist.csv"
    )

    primary_results = pd.read_csv(
        config.PRIMARY_STATS_DIR / "primary_paired_feature_results_all_features_fdr.csv"
    )

    shortlisted_features = shortlist["feature"].dropna().astype(str).tolist()

    delta_df[config.PATIENT_ID_COL] = standardize_patient_id_series(
        delta_df[config.PATIENT_ID_COL]
    )

    patient_ids = sorted(
        delta_df[config.PATIENT_ID_COL].dropna().astype(str).unique().tolist()
    )

    # -----------------------------------------------------------------
    # Helper function to compute one feature result after excluding patients
    # -----------------------------------------------------------------

    def _compute_feature_result_from_delta(
        feature: str,
        excluded_patient_id: str | None = None,
    ) -> dict:
        sub = delta_df[
            (delta_df["feature"] == feature)
            & (delta_df["paired_observation_available"] == 1)
        ].copy()

        if excluded_patient_id is not None:
            sub = sub[sub[config.PATIENT_ID_COL].astype(str) != str(excluded_patient_id)]

        n_paired = len(sub)

        n_increase = int((sub["observed_direction"] == "increase").sum())
        n_decrease = int((sub["observed_direction"] == "decrease").sum())
        n_no_change = int((sub["observed_direction"] == "no_change").sum())

        dominant_direction = dominant_direction_from_counts(
            n_increase,
            n_decrease,
            n_no_change,
        )

        max_direction_count = (
            max(n_increase, n_decrease, n_no_change) if n_paired > 0 else 0
        )

        direction_consistency = (
            max_direction_count / n_paired if n_paired > 0 else np.nan
        )

        if n_paired >= 2:
            p_value = wilcoxon_signed_rank_pvalue(
                sub["PreHI_median"],
                sub["HI_median"],
            )

            effect_size = matched_pairs_rank_biserial_effect_size(
                sub["PreHI_median"],
                sub["HI_median"],
            )

            ci_low, ci_high = bootstrap_median_ci(
                sub["delta_HI_minus_PreHI"],
                n_iterations=config.BOOTSTRAP_N_ITERATIONS,
                confidence_level=config.BOOTSTRAP_CONFIDENCE_LEVEL,
                random_seed=config.RANDOM_SEED,
            )

        else:
            p_value = np.nan
            effect_size = np.nan
            ci_low = np.nan
            ci_high = np.nan

        bootstrap_ci_excludes_zero = (
            pd.notna(ci_low)
            and pd.notna(ci_high)
            and ((ci_low > 0 and ci_high > 0) or (ci_low < 0 and ci_high < 0))
        )

        if n_paired < config.LOO_MIN_PAIRED_PATIENTS:
            loo_status = "limited_sample_size"
        else:
            loo_status = "adequate_sample_size"

        return {
            "feature": feature,
            "feature_family": infer_feature_family(feature),
            "excluded_patient_id": excluded_patient_id if excluded_patient_id is not None else "",
            "n_paired_patients": n_paired,
            "median_PreHI": sub["PreHI_median"].median(),
            "median_HI": sub["HI_median"].median(),
            "median_delta_HI_minus_PreHI": sub["delta_HI_minus_PreHI"].median(),
            "n_increase": n_increase,
            "n_decrease": n_decrease,
            "n_no_change": n_no_change,
            "dominant_direction": dominant_direction,
            "direction_consistency": direction_consistency,
            "wilcoxon_p": p_value,
            "rank_biserial_effect_size": effect_size,
            "bootstrap_median_delta_ci_low": ci_low,
            "bootstrap_median_delta_ci_high": ci_high,
            "bootstrap_ci_excludes_zero": int(bootstrap_ci_excludes_zero),
            "loo_status": loo_status,
        }

    # -----------------------------------------------------------------
    # Compute baseline and leave-one-out results
    # -----------------------------------------------------------------

    loo_rows = []

    for feature in shortlisted_features:
        baseline = _compute_feature_result_from_delta(
            feature=feature,
            excluded_patient_id=None,
        )

        baseline["analysis_type"] = "baseline_all_patients"
        loo_rows.append(baseline)

        for patient_id in patient_ids:
            result = _compute_feature_result_from_delta(
                feature=feature,
                excluded_patient_id=patient_id,
            )

            result["analysis_type"] = "leave_one_patient_out"
            loo_rows.append(result)

    loo_results = pd.DataFrame(loo_rows)

    # -----------------------------------------------------------------
    # Compare LOO rows against baseline rows
    # -----------------------------------------------------------------

    baseline_for_merge = loo_results[
        loo_results["analysis_type"] == "baseline_all_patients"
    ][
        [
            "feature",
            "dominant_direction",
            "median_delta_HI_minus_PreHI",
            "direction_consistency",
            "wilcoxon_p",
            "rank_biserial_effect_size",
            "bootstrap_ci_excludes_zero",
        ]
    ].copy()

    baseline_for_merge = baseline_for_merge.rename(
        columns={
            "dominant_direction": "baseline_dominant_direction",
            "median_delta_HI_minus_PreHI": "baseline_median_delta_HI_minus_PreHI",
            "direction_consistency": "baseline_direction_consistency",
            "wilcoxon_p": "baseline_wilcoxon_p",
            "rank_biserial_effect_size": "baseline_rank_biserial_effect_size",
            "bootstrap_ci_excludes_zero": "baseline_bootstrap_ci_excludes_zero",
        }
    )

    loo_results = loo_results.merge(
        baseline_for_merge,
        on="feature",
        how="left",
    )

    loo_results["direction_matches_baseline"] = (
        loo_results["dominant_direction"]
        == loo_results["baseline_dominant_direction"]
    ).astype(int)

    loo_results.loc[
        loo_results["analysis_type"] == "baseline_all_patients",
        "direction_matches_baseline",
    ] = 1

    loo_results["median_delta_sign_matches_baseline"] = (
        np.sign(loo_results["median_delta_HI_minus_PreHI"])
        == np.sign(loo_results["baseline_median_delta_HI_minus_PreHI"])
    ).astype(int)

    loo_results.loc[
        loo_results["analysis_type"] == "baseline_all_patients",
        "median_delta_sign_matches_baseline",
    ] = 1

    loo_results["direction_consistency_drop"] = (
        loo_results["baseline_direction_consistency"]
        - loo_results["direction_consistency"]
    )

    loo_results["p_value_ratio_vs_baseline"] = (
        loo_results["wilcoxon_p"] / loo_results["baseline_wilcoxon_p"]
    )

    loo_results["potentially_influential_patient_flag"] = 0

    loo_mask = loo_results["analysis_type"] == "leave_one_patient_out"

    loo_results.loc[
        loo_mask
        & (
            (loo_results["direction_matches_baseline"] == 0)
            | (loo_results["median_delta_sign_matches_baseline"] == 0)
            | (
                loo_results["direction_consistency_drop"]
                >= config.LOO_DIRECTION_CONSISTENCY_DROP_THRESHOLD
            )
            | (
                loo_results["p_value_ratio_vs_baseline"]
                >= config.LOO_P_VALUE_RATIO_INFLUENCE_THRESHOLD
            )
        ),
        "potentially_influential_patient_flag",
    ] = 1

    save_dataframe(
        loo_results,
        config.LOO_DIR / "leave_one_patient_out_results.csv",
    )

    # -----------------------------------------------------------------
    # Feature-level LOO summary
    # -----------------------------------------------------------------

    summary_rows = []

    for feature, sub in loo_results.groupby("feature"):
        baseline_row = sub[sub["analysis_type"] == "baseline_all_patients"].iloc[0]
        loo_sub = sub[sub["analysis_type"] == "leave_one_patient_out"].copy()

        adequate_loo = loo_sub[loo_sub["loo_status"] == "adequate_sample_size"].copy()

        n_loo_tests = len(loo_sub)
        n_adequate_loo = len(adequate_loo)

        n_direction_changes = int(
            (
                (adequate_loo["direction_matches_baseline"] == 0)
            ).sum()
        )

        n_delta_sign_changes = int(
            (
                (adequate_loo["median_delta_sign_matches_baseline"] == 0)
            ).sum()
        )

        n_large_direction_consistency_drop = int(
            (
                adequate_loo["direction_consistency_drop"]
                >= config.LOO_DIRECTION_CONSISTENCY_DROP_THRESHOLD
            ).sum()
        )

        n_large_p_value_ratio = int(
            (
                adequate_loo["p_value_ratio_vs_baseline"]
                >= config.LOO_P_VALUE_RATIO_INFLUENCE_THRESHOLD
            ).sum()
        )

        n_potentially_influential_patients = int(
            adequate_loo["potentially_influential_patient_flag"].sum()
        )

        if n_adequate_loo > 0:
            max_direction_consistency_drop = adequate_loo[
                "direction_consistency_drop"
            ].max()

            max_p_value_ratio = adequate_loo[
                "p_value_ratio_vs_baseline"
            ].replace([np.inf, -np.inf], np.nan).max()

            min_loo_direction_consistency = adequate_loo[
                "direction_consistency"
            ].min()

            max_loo_p_value = adequate_loo["wilcoxon_p"].max()

        else:
            max_direction_consistency_drop = np.nan
            max_p_value_ratio = np.nan
            min_loo_direction_consistency = np.nan
            max_loo_p_value = np.nan

        influential_patient_ids = (
            adequate_loo.loc[
                adequate_loo["potentially_influential_patient_flag"] == 1,
                "excluded_patient_id",
            ]
            .dropna()
            .astype(str)
            .tolist()
        )

        if (
            n_direction_changes == 0
            and n_delta_sign_changes == 0
            and n_potentially_influential_patients == 0
        ):
            loo_stability_label = "stable"

        elif (
            n_direction_changes == 0
            and n_delta_sign_changes == 0
            and n_potentially_influential_patients <= 2
        ):
            loo_stability_label = "moderately_stable"

        else:
            loo_stability_label = "influential_patient_sensitive"

        summary_rows.append(
            {
                "feature": feature,
                "feature_family": infer_feature_family(feature),
                "baseline_dominant_direction": baseline_row[
                    "baseline_dominant_direction"
                ],
                "baseline_median_delta_HI_minus_PreHI": baseline_row[
                    "baseline_median_delta_HI_minus_PreHI"
                ],
                "baseline_direction_consistency": baseline_row[
                    "baseline_direction_consistency"
                ],
                "baseline_wilcoxon_p": baseline_row["baseline_wilcoxon_p"],
                "baseline_rank_biserial_effect_size": baseline_row[
                    "baseline_rank_biserial_effect_size"
                ],
                "baseline_bootstrap_ci_excludes_zero": baseline_row[
                    "baseline_bootstrap_ci_excludes_zero"
                ],
                "n_loo_tests": n_loo_tests,
                "n_adequate_loo_tests": n_adequate_loo,
                "n_direction_changes": n_direction_changes,
                "n_delta_sign_changes": n_delta_sign_changes,
                "n_large_direction_consistency_drop": n_large_direction_consistency_drop,
                "n_large_p_value_ratio": n_large_p_value_ratio,
                "n_potentially_influential_patients": n_potentially_influential_patients,
                "max_direction_consistency_drop": max_direction_consistency_drop,
                "min_loo_direction_consistency": min_loo_direction_consistency,
                "max_p_value_ratio_vs_baseline": max_p_value_ratio,
                "max_loo_wilcoxon_p": max_loo_p_value,
                "potentially_influential_patient_ids": "; ".join(
                    influential_patient_ids
                ),
                "loo_stability_label": loo_stability_label,
            }
        )

    loo_summary = pd.DataFrame(summary_rows)

    label_order = {
        "stable": 1,
        "moderately_stable": 2,
        "influential_patient_sensitive": 3,
    }

    loo_summary["loo_stability_sort_order"] = loo_summary[
        "loo_stability_label"
    ].map(label_order)

    loo_summary = loo_summary.sort_values(
        [
            "loo_stability_sort_order",
            "n_potentially_influential_patients",
            "max_direction_consistency_drop",
            "baseline_wilcoxon_p",
        ],
        ascending=[True, True, True, True],
    ).drop(columns=["loo_stability_sort_order"])

    save_dataframe(
        loo_summary,
        config.LOO_DIR / "leave_one_patient_out_summary.csv",
    )

    logging.info("Stage 14 completed: Leave-one-patient-out influence analysis.")

    return loo_results

# ---------------------------------------------------------------------
# Stage 15 — Patient-level agreement/disagreement analysis
# ---------------------------------------------------------------------

def stage15_patient_agreement_disagreement_analysis() -> pd.DataFrame:
    """
    Analyze patient-level agreement and disagreement for shortlisted features.

    Purpose:
    - For each shortlisted feature, determine whether each patient's observed
      PreHI -> HI change agrees with the feature's cohort-level dominant direction.
    - Summarize which patients repeatedly disagree across features.
    - Examine whether disagreement is enriched among clinically confounded
      patients or signal-quality caution patients.

    Outputs:
    - patient_feature_agreement_table.csv
    - feature_agreement_summary.csv
    - patient_disagreement_summary.csv
    - clinical_flag_disagreement_summary.csv
    - feature_patient_disagreement_matrix.csv
    """

    logging.info("Stage 15 started: Patient-level agreement/disagreement analysis.")

    delta_df = pd.read_csv(
        config.PATIENT_LEVEL_DIR / "patient_feature_delta_table_all_features.csv"
    )

    shortlist = pd.read_csv(
        config.CANDIDATE_SCORECARD_DIR / "preliminary_candidate_shortlist.csv"
    )

    primary_results = pd.read_csv(
        config.PRIMARY_STATS_DIR / "primary_paired_feature_results_all_features_fdr.csv"
    )

    patient_flags = pd.read_csv(
        config.CLINICAL_METADATA_DIR / "patient_clinical_flag_summary.csv"
    )

    sensitivity_summary = pd.read_csv(
        config.SENSITIVITY_DIR / "candidate_feature_robustness_matrix.csv"
    )

    loo_summary = pd.read_csv(
        config.LOO_DIR / "leave_one_patient_out_summary.csv"
    )

    # -----------------------------------------------------------------
    # Standardize patient IDs
    # -----------------------------------------------------------------

    delta_df[config.PATIENT_ID_COL] = standardize_patient_id_series(
        delta_df[config.PATIENT_ID_COL]
    )

    patient_flags[config.PATIENT_ID_COL] = standardize_patient_id_series(
        patient_flags[config.PATIENT_ID_COL]
    )

    shortlisted_features = shortlist["feature"].dropna().astype(str).tolist()

    # -----------------------------------------------------------------
    # Build cohort-level direction reference
    # -----------------------------------------------------------------

    direction_reference = primary_results[
        primary_results["feature"].isin(shortlisted_features)
    ][
        [
            "feature",
            "feature_family",
            "dominant_direction",
            "direction_consistency",
            "wilcoxon_p",
            "rank_biserial_effect_size",
            "abs_rank_biserial_effect_size",
            "bootstrap_ci_excludes_zero",
            "candidate_support_tier_preliminary",
        ]
    ].copy()

    direction_reference = direction_reference.rename(
        columns={
            "dominant_direction": "cohort_dominant_direction",
            "direction_consistency": "cohort_direction_consistency",
            "wilcoxon_p": "cohort_wilcoxon_p",
            "rank_biserial_effect_size": "cohort_rank_biserial_effect_size",
            "abs_rank_biserial_effect_size": "cohort_abs_rank_biserial_effect_size",
            "bootstrap_ci_excludes_zero": "cohort_bootstrap_ci_excludes_zero",
        }
    )

    # Add sensitivity and LOO summaries for context
    sensitivity_cols = [
        "feature",
        "sensitivity_robustness_label",
        "prop_direction_retained",
        "prop_delta_sign_retained",
        "n_scenarios_raw_p_less_0_05",
        "n_scenarios_bootstrap_ci_excludes_zero",
    ]

    sensitivity_cols = [
        c for c in sensitivity_cols if c in sensitivity_summary.columns
    ]

    direction_reference = direction_reference.merge(
        sensitivity_summary[sensitivity_cols],
        on="feature",
        how="left",
    )

    loo_cols = [
        "feature",
        "loo_stability_label",
        "n_direction_changes",
        "n_delta_sign_changes",
        "n_potentially_influential_patients",
        "potentially_influential_patient_ids",
    ]

    loo_cols = [c for c in loo_cols if c in loo_summary.columns]

    direction_reference = direction_reference.merge(
        loo_summary[loo_cols],
        on="feature",
        how="left",
    )

    # -----------------------------------------------------------------
    # Patient-feature agreement table
    # -----------------------------------------------------------------

    agreement = delta_df[
        delta_df["feature"].isin(shortlisted_features)
    ].copy()

    agreement = agreement.merge(
        direction_reference,
        on="feature",
        how="left",
    )

    # Agreement is defined against the primary cohort dominant direction.
    # If the cohort dominant direction is no_change, agreement means patient
    # also shows no_change. This mainly affects fd_minimum.

    def _agreement_status(row):
        observed = row.get("observed_direction", "")
        cohort = row.get("cohort_dominant_direction", "")

        if row.get("paired_observation_available", 0) != 1:
            return "missing_pair"

        if cohort in ["", "missing", "tie_or_mixed"] or pd.isna(cohort):
            return "cohort_direction_unavailable"

        if observed == cohort:
            return "agree"

        if observed == "no_change":
            return "neutral_no_change"

        return "disagree"

    agreement["agreement_status"] = agreement.apply(
        _agreement_status,
        axis=1,
    )

    agreement["agrees_with_cohort_direction"] = (
        agreement["agreement_status"] == "agree"
    ).astype(int)

    agreement["disagrees_with_cohort_direction"] = (
        agreement["agreement_status"] == "disagree"
    ).astype(int)

    agreement["neutral_no_change_against_direction"] = (
        agreement["agreement_status"] == "neutral_no_change"
    ).astype(int)

    # Add clinical flags to every patient-feature row
    agreement = agreement.merge(
        patient_flags,
        on=config.PATIENT_ID_COL,
        how="left",
        suffixes=("", "_clinical"),
    )

    # Ensure configured clinical flags exist and are numeric
    for flag in config.AGREEMENT_CLINICAL_FLAGS_TO_SUMMARIZE:
        if flag not in agreement.columns:
            agreement[flag] = 0

        agreement[flag] = pd.to_numeric(
            agreement[flag],
            errors="coerce",
        ).fillna(0).astype(int)

    # Helpful human-readable patient-level interpretation note
    def _patient_feature_note(row):
        feature = row["feature"]
        status = row["agreement_status"]
        observed = row.get("observed_direction", "")
        cohort = row.get("cohort_dominant_direction", "")

        clinical_reasons = []

        if row.get("paced_or_icd", 0) == 1:
            clinical_reasons.append("paced/ICD")

        if row.get("BBB", 0) == 1:
            clinical_reasons.append("BBB")

        if row.get("any_af", 0) == 1:
            clinical_reasons.append("AF")

        if row.get("mechanical_support", 0) == 1:
            clinical_reasons.append("mechanical support")

        if row.get("atypical_HI", 0) == 1:
            clinical_reasons.append("atypical HI")

        if row.get("major_HRV_outlier", 0) == 1:
            clinical_reasons.append("major HRV outlier")

        if row.get("known_patient_qc_caution_flag", 0) == 1:
            clinical_reasons.append("signal/QC caution")

        if status == "agree":
            base = (
                f"Patient change agrees with cohort-level {feature} "
                f"direction ({cohort})."
            )
        elif status == "disagree":
            base = (
                f"Patient change disagrees with cohort-level {feature} "
                f"direction; observed {observed}, cohort {cohort}."
            )
        elif status == "neutral_no_change":
            base = (
                f"Patient shows no change for {feature} while cohort direction "
                f"is {cohort}."
            )
        else:
            base = f"Agreement status for {feature}: {status}."

        if clinical_reasons:
            return base + " Clinical context: " + "; ".join(clinical_reasons) + "."

        return base

    agreement["patient_feature_interpretation_note"] = agreement.apply(
        _patient_feature_note,
        axis=1,
    )

    # Sort for readability
    agreement = agreement.sort_values(
        [
            "feature",
            "agreement_status",
            config.PATIENT_ID_COL,
        ]
    ).reset_index(drop=True)

    save_dataframe(
        agreement,
        config.AGREEMENT_DIR / "patient_feature_agreement_table.csv",
    )

    # -----------------------------------------------------------------
    # Feature-level agreement summary
    # -----------------------------------------------------------------

    feature_rows = []

    for feature, sub in agreement.groupby("feature"):
        n_total = len(sub)
        n_available = int((sub["paired_observation_available"] == 1).sum())
        n_agree = int((sub["agreement_status"] == "agree").sum())
        n_disagree = int((sub["agreement_status"] == "disagree").sum())
        n_neutral = int((sub["agreement_status"] == "neutral_no_change").sum())
        n_missing = int((sub["agreement_status"] == "missing_pair").sum())

        if n_available > 0:
            prop_agree = n_agree / n_available
            prop_disagree = n_disagree / n_available
            prop_neutral = n_neutral / n_available
        else:
            prop_agree = np.nan
            prop_disagree = np.nan
            prop_neutral = np.nan

        direction_row = direction_reference[
            direction_reference["feature"] == feature
        ]

        if not direction_row.empty:
            direction_row = direction_row.iloc[0]
            cohort_direction = direction_row.get("cohort_dominant_direction", "")
            cohort_p = direction_row.get("cohort_wilcoxon_p", np.nan)
            cohort_consistency = direction_row.get(
                "cohort_direction_consistency", np.nan
            )
            sensitivity_label = direction_row.get(
                "sensitivity_robustness_label", ""
            )
            loo_label = direction_row.get("loo_stability_label", "")
        else:
            cohort_direction = ""
            cohort_p = np.nan
            cohort_consistency = np.nan
            sensitivity_label = ""
            loo_label = ""

        disagreeing_patients = (
            sub.loc[
                sub["agreement_status"] == "disagree",
                config.PATIENT_ID_COL,
            ]
            .dropna()
            .astype(str)
            .tolist()
        )

        neutral_patients = (
            sub.loc[
                sub["agreement_status"] == "neutral_no_change",
                config.PATIENT_ID_COL,
            ]
            .dropna()
            .astype(str)
            .tolist()
        )

        feature_rows.append(
            {
                "feature": feature,
                "feature_family": infer_feature_family(feature),
                "cohort_dominant_direction": cohort_direction,
                "cohort_wilcoxon_p": cohort_p,
                "cohort_direction_consistency": cohort_consistency,
                "n_patient_pairs_available": n_available,
                "n_agree": n_agree,
                "n_disagree": n_disagree,
                "n_neutral_no_change": n_neutral,
                "n_missing_pairs": n_missing,
                "prop_agree": prop_agree,
                "prop_disagree": prop_disagree,
                "prop_neutral_no_change": prop_neutral,
                "disagreeing_patient_ids": "; ".join(disagreeing_patients),
                "neutral_no_change_patient_ids": "; ".join(neutral_patients),
                "sensitivity_robustness_label": sensitivity_label,
                "loo_stability_label": loo_label,
            }
        )

    feature_agreement_summary = pd.DataFrame(feature_rows).sort_values(
        [
            "prop_agree",
            "cohort_wilcoxon_p",
        ],
        ascending=[False, True],
    )

    save_dataframe(
        feature_agreement_summary,
        config.AGREEMENT_DIR / "feature_agreement_summary.csv",
    )

    # -----------------------------------------------------------------
    # Patient-level disagreement summary
    # -----------------------------------------------------------------

    patient_rows = []

    for patient_id, sub in agreement.groupby(config.PATIENT_ID_COL):
        n_features_available = int((sub["paired_observation_available"] == 1).sum())
        n_agree = int((sub["agreement_status"] == "agree").sum())
        n_disagree = int((sub["agreement_status"] == "disagree").sum())
        n_neutral = int((sub["agreement_status"] == "neutral_no_change").sum())

        if n_features_available > 0:
            prop_agree = n_agree / n_features_available
            prop_disagree = n_disagree / n_features_available
            prop_neutral = n_neutral / n_features_available
        else:
            prop_agree = np.nan
            prop_disagree = np.nan
            prop_neutral = np.nan

        patient_flag_row = patient_flags[
            patient_flags[config.PATIENT_ID_COL].astype(str) == str(patient_id)
        ]

        if not patient_flag_row.empty:
            patient_flag_row = patient_flag_row.iloc[0]
        else:
            patient_flag_row = pd.Series(dtype="object")

        disagreeing_features = (
            sub.loc[sub["agreement_status"] == "disagree", "feature"]
            .dropna()
            .astype(str)
            .tolist()
        )

        neutral_features = (
            sub.loc[sub["agreement_status"] == "neutral_no_change", "feature"]
            .dropna()
            .astype(str)
            .tolist()
        )

        clinical_context = []

        for flag in config.AGREEMENT_CLINICAL_FLAGS_TO_SUMMARIZE:
            if flag in patient_flag_row.index:
                try:
                    if int(patient_flag_row.get(flag, 0)) == 1:
                        clinical_context.append(flag)
                except Exception:
                    pass

        patient_rows.append(
            {
                config.PATIENT_ID_COL: patient_id,
                "n_features_available": n_features_available,
                "n_agree": n_agree,
                "n_disagree": n_disagree,
                "n_neutral_no_change": n_neutral,
                "prop_agree": prop_agree,
                "prop_disagree": prop_disagree,
                "prop_neutral_no_change": prop_neutral,
                "disagreeing_features": "; ".join(disagreeing_features),
                "neutral_no_change_features": "; ".join(neutral_features),
                "clinical_flags_present": "; ".join(clinical_context),
            }
        )

    patient_disagreement_summary = pd.DataFrame(patient_rows).sort_values(
        [
            "n_disagree",
            "prop_disagree",
            config.PATIENT_ID_COL,
        ],
        ascending=[False, False, True],
    )

    save_dataframe(
        patient_disagreement_summary,
        config.AGREEMENT_DIR / "patient_disagreement_summary.csv",
    )

    # -----------------------------------------------------------------
    # Clinical flag disagreement summary
    # -----------------------------------------------------------------

    clinical_rows = []

    for feature in shortlisted_features:
        feature_sub = agreement[
            (agreement["feature"] == feature)
            & (agreement["paired_observation_available"] == 1)
        ].copy()

        for flag in config.AGREEMENT_CLINICAL_FLAGS_TO_SUMMARIZE:
            if flag not in feature_sub.columns:
                continue

            flag_positive = feature_sub[feature_sub[flag] == 1]
            flag_negative = feature_sub[feature_sub[flag] == 0]

            n_flag_positive = len(flag_positive)
            n_flag_negative = len(flag_negative)

            n_disagree_flag_positive = int(
                (flag_positive["agreement_status"] == "disagree").sum()
            )

            n_disagree_flag_negative = int(
                (flag_negative["agreement_status"] == "disagree").sum()
            )

            prop_disagree_flag_positive = (
                n_disagree_flag_positive / n_flag_positive
                if n_flag_positive > 0
                else np.nan
            )

            prop_disagree_flag_negative = (
                n_disagree_flag_negative / n_flag_negative
                if n_flag_negative > 0
                else np.nan
            )

            clinical_rows.append(
                {
                    "feature": feature,
                    "feature_family": infer_feature_family(feature),
                    "clinical_flag": flag,
                    "n_flag_positive": n_flag_positive,
                    "n_flag_negative": n_flag_negative,
                    "n_disagree_flag_positive": n_disagree_flag_positive,
                    "n_disagree_flag_negative": n_disagree_flag_negative,
                    "prop_disagree_flag_positive": prop_disagree_flag_positive,
                    "prop_disagree_flag_negative": prop_disagree_flag_negative,
                    "difference_in_disagreement_proportion_positive_minus_negative": (
                        prop_disagree_flag_positive - prop_disagree_flag_negative
                        if pd.notna(prop_disagree_flag_positive)
                        and pd.notna(prop_disagree_flag_negative)
                        else np.nan
                    ),
                    "note": (
                        "Descriptive only; small subgroup counts should not be "
                        "interpreted as formal association testing."
                    ),
                }
            )

    clinical_flag_disagreement_summary = pd.DataFrame(clinical_rows).sort_values(
        [
            "feature",
            "difference_in_disagreement_proportion_positive_minus_negative",
        ],
        ascending=[True, False],
    )

    save_dataframe(
        clinical_flag_disagreement_summary,
        config.AGREEMENT_DIR / "clinical_flag_disagreement_summary.csv",
    )

    # -----------------------------------------------------------------
    # Feature-patient disagreement matrix
    # -----------------------------------------------------------------

    matrix_input = agreement[
        [
            config.PATIENT_ID_COL,
            "feature",
            "agreement_status",
        ]
    ].copy()

    matrix_input["agreement_code"] = matrix_input["agreement_status"].map(
        {
            "agree": 1,
            "neutral_no_change": 0,
            "disagree": -1,
            "missing_pair": np.nan,
            "cohort_direction_unavailable": np.nan,
        }
    )

    feature_patient_matrix = matrix_input.pivot(
        index=config.PATIENT_ID_COL,
        columns="feature",
        values="agreement_code",
    ).reset_index()

    save_dataframe(
        feature_patient_matrix,
        config.AGREEMENT_DIR / "feature_patient_disagreement_matrix.csv",
    )

    logging.info("Stage 15 completed: Patient-level agreement/disagreement analysis.")

    return agreement

# ---------------------------------------------------------------------
# Stage 16 — Final candidate scorecard
# ---------------------------------------------------------------------

def stage16_final_candidate_scorecard() -> pd.DataFrame:
    """
    Build the final candidate feature scorecard.

    Purpose:
    - Combine evidence from primary paired analysis, FDR correction,
      redundancy filtering, HRV-valid subgroup analysis, supportive LMM,
      sensitivity analysis, leave-one-patient-out analysis, and
      patient-level agreement/disagreement analysis.
    - Classify shortlisted features into transparent exploratory categories.

    Important:
    - This is not external validation.
    - Final categories are evidence summaries for this 20-patient cohort.
    - HRV features are treated cautiously if weakened in the HRV-valid subgroup.

    Outputs:
    - final_candidate_scorecard.csv
    - final_candidate_feature_categories.csv
    - final_candidate_feature_summary_for_manuscript.csv
    """

    logging.info("Stage 16 started: Final candidate scorecard.")

    # -----------------------------------------------------------------
    # Load inputs
    # -----------------------------------------------------------------

    shortlist = pd.read_csv(
        config.CANDIDATE_SCORECARD_DIR / "preliminary_candidate_shortlist.csv"
    )

    primary = pd.read_csv(
        config.PRIMARY_STATS_DIR / "primary_paired_feature_results_all_features_fdr.csv"
    )

    cluster_reps = pd.read_csv(
        config.REDUNDANCY_DIR / "cluster_representatives.csv"
    )

    hrv_results = pd.read_csv(
        config.HRV_SUBGROUP_DIR / "hrv_valid_feature_results.csv"
    )

    lmm_results = pd.read_csv(
        config.MIXED_EFFECTS_DIR / "lmm_supportive_results_all_features.csv"
    )

    sensitivity = pd.read_csv(
        config.SENSITIVITY_DIR / "candidate_feature_robustness_matrix.csv"
    )

    loo_summary = pd.read_csv(
        config.LOO_DIR / "leave_one_patient_out_summary.csv"
    )

    feature_agreement = pd.read_csv(
        config.AGREEMENT_DIR / "feature_agreement_summary.csv"
    )

    clinical_disagreement = pd.read_csv(
        config.AGREEMENT_DIR / "clinical_flag_disagreement_summary.csv"
    )

    shortlisted_features = shortlist["feature"].dropna().astype(str).tolist()

    # -----------------------------------------------------------------
    # Start scorecard from shortlist
    # -----------------------------------------------------------------

    scorecard = shortlist.copy()

    keep_shortlist_cols = [
        "feature",
        "feature_family",
        "preliminary_rank",
        "preliminary_priority_score",
        "preliminary_priority_tier",
        "included_in_preliminary_shortlist",
        "cluster_id",
        "cluster_size",
        "cluster_type",
        "is_cluster_representative",
    ]

    keep_shortlist_cols = [c for c in keep_shortlist_cols if c in scorecard.columns]

    scorecard = scorecard[keep_shortlist_cols].copy()

    # -----------------------------------------------------------------
    # Merge primary paired evidence
    # -----------------------------------------------------------------

    primary_cols = [
        "feature",
        "n_paired_patients",
        "median_PreHI",
        "median_HI",
        "median_delta_HI_minus_PreHI",
        "median_percent_change",
        "n_increase",
        "n_decrease",
        "n_no_change",
        "dominant_direction",
        "direction_consistency",
        "wilcoxon_p",
        "wilcoxon_p_fdr_global",
        "wilcoxon_p_fdr_family",
        "global_fdr_significant",
        "family_fdr_significant",
        "rank_biserial_effect_size",
        "abs_rank_biserial_effect_size",
        "bootstrap_median_delta_ci_low",
        "bootstrap_median_delta_ci_high",
        "bootstrap_ci_excludes_zero",
        "candidate_support_tier_preliminary",
    ]

    primary_cols = [c for c in primary_cols if c in primary.columns]

    scorecard = scorecard.merge(
        primary[primary_cols],
        on="feature",
        how="left",
    )

    # -----------------------------------------------------------------
    # Merge redundancy representative information
    # -----------------------------------------------------------------

    cluster_reps_cols = [
        "cluster_id",
        "representative_feature",
        "cluster_members",
    ]

    cluster_reps_cols = [c for c in cluster_reps_cols if c in cluster_reps.columns]

    if "cluster_id" in scorecard.columns and "cluster_id" in cluster_reps.columns:
        scorecard = scorecard.merge(
            cluster_reps[cluster_reps_cols],
            on="cluster_id",
            how="left",
        )

    # -----------------------------------------------------------------
    # Merge HRV-valid subgroup results
    # -----------------------------------------------------------------

    hrv_cols = [
        "feature",
        "n_hrv_valid_patients_total",
        "n_paired_patients",
        "dominant_direction",
        "direction_consistency",
        "wilcoxon_p",
        "rank_biserial_effect_size",
        "bootstrap_ci_excludes_zero",
        "direction_agrees_with_full_cohort",
    ]

    hrv_cols = [c for c in hrv_cols if c in hrv_results.columns]

    hrv_rename = {
        "n_paired_patients": "hrv_valid_n_paired_patients",
        "dominant_direction": "hrv_valid_dominant_direction",
        "direction_consistency": "hrv_valid_direction_consistency",
        "wilcoxon_p": "hrv_valid_wilcoxon_p",
        "rank_biserial_effect_size": "hrv_valid_rank_biserial_effect_size",
        "bootstrap_ci_excludes_zero": "hrv_valid_bootstrap_ci_excludes_zero",
        "direction_agrees_with_full_cohort": "hrv_valid_direction_agrees_with_full_cohort",
    }

    if not hrv_results.empty:
        scorecard = scorecard.merge(
            hrv_results[hrv_cols].rename(columns=hrv_rename),
            on="feature",
            how="left",
        )

    # -----------------------------------------------------------------
    # Merge LMM evidence
    # -----------------------------------------------------------------

    lmm_cols = [
        "feature",
        "lmm_status",
        "n_segments",
        "n_patients",
        "beta_condition_HI",
        "p_condition_HI",
        "p_condition_HI_fdr",
        "lmm_raw_p_less_0_05",
        "lmm_fdr_significant",
        "lmm_direction",
        "lmm_agrees_with_primary_direction",
        "converged",
    ]

    lmm_cols = [c for c in lmm_cols if c in lmm_results.columns]

    lmm_rename = {
        "n_segments": "lmm_n_segments",
        "n_patients": "lmm_n_patients",
    }

    scorecard = scorecard.merge(
        lmm_results[lmm_cols].rename(columns=lmm_rename),
        on="feature",
        how="left",
    )

    # -----------------------------------------------------------------
    # Merge sensitivity evidence
    # -----------------------------------------------------------------

    sensitivity_cols = [
        "feature",
        "sensitivity_robustness_label",
        "n_sensitivity_scenarios_non_baseline",
        "n_scenarios_with_adequate_sample_size",
        "n_direction_retained",
        "prop_direction_retained",
        "n_delta_sign_retained",
        "prop_delta_sign_retained",
        "median_direction_consistency_across_scenarios",
        "n_scenarios_raw_p_less_0_05",
        "n_scenarios_bootstrap_ci_excludes_zero",
    ]

    sensitivity_cols = [c for c in sensitivity_cols if c in sensitivity.columns]

    scorecard = scorecard.merge(
        sensitivity[sensitivity_cols],
        on="feature",
        how="left",
    )

    # -----------------------------------------------------------------
    # Merge leave-one-patient-out evidence
    # -----------------------------------------------------------------

    loo_cols = [
        "feature",
        "loo_stability_label",
        "n_loo_tests",
        "n_adequate_loo_tests",
        "n_direction_changes",
        "n_delta_sign_changes",
        "n_potentially_influential_patients",
        "max_direction_consistency_drop",
        "min_loo_direction_consistency",
        "max_p_value_ratio_vs_baseline",
        "max_loo_wilcoxon_p",
        "potentially_influential_patient_ids",
    ]

    loo_cols = [c for c in loo_cols if c in loo_summary.columns]

    scorecard = scorecard.merge(
        loo_summary[loo_cols],
        on="feature",
        how="left",
    )

    # -----------------------------------------------------------------
    # Merge feature agreement evidence
    # -----------------------------------------------------------------

    agreement_cols = [
        "feature",
        "n_patient_pairs_available",
        "n_agree",
        "n_disagree",
        "n_neutral_no_change",
        "prop_agree",
        "prop_disagree",
        "prop_neutral_no_change",
        "disagreeing_patient_ids",
        "neutral_no_change_patient_ids",
    ]

    agreement_cols = [c for c in agreement_cols if c in feature_agreement.columns]

    scorecard = scorecard.merge(
        feature_agreement[agreement_cols],
        on="feature",
        how="left",
    )

    # -----------------------------------------------------------------
    # Clinical disagreement context
    # -----------------------------------------------------------------

    # For each feature, summarize the clinical flags where disagreement appears
    # most enriched. This is descriptive only.

    clinical_context_rows = []

    if not clinical_disagreement.empty:
        for feature in shortlisted_features:
            sub = clinical_disagreement[
                clinical_disagreement["feature"] == feature
            ].copy()

            if sub.empty:
                clinical_context_rows.append(
                    {
                        "feature": feature,
                        "top_clinical_disagreement_context": "",
                    }
                )
                continue

            sub = sub[
                sub["difference_in_disagreement_proportion_positive_minus_negative"]
                .notna()
            ].copy()

            sub = sub[
                sub["n_flag_positive"] > 0
            ].copy()

            sub = sub.sort_values(
                [
                    "difference_in_disagreement_proportion_positive_minus_negative",
                    "n_disagree_flag_positive",
                ],
                ascending=[False, False],
            )

            top_items = []

            for _, row in sub.head(3).iterrows():
                flag = row["clinical_flag"]
                n_pos = int(row["n_flag_positive"])
                n_dis_pos = int(row["n_disagree_flag_positive"])
                prop_pos = row["prop_disagree_flag_positive"]

                top_items.append(
                    f"{flag}: {n_dis_pos}/{n_pos} disagree"
                )

            clinical_context_rows.append(
                {
                    "feature": feature,
                    "top_clinical_disagreement_context": "; ".join(top_items),
                }
            )

    clinical_context = pd.DataFrame(clinical_context_rows)

    scorecard = scorecard.merge(
        clinical_context,
        on="feature",
        how="left",
    )

    # -----------------------------------------------------------------
    # Evidence scoring
    # -----------------------------------------------------------------

    # Convert relevant columns safely
    numeric_like_cols = [
        "global_fdr_significant",
        "family_fdr_significant",
        "bootstrap_ci_excludes_zero",
        "direction_consistency",
        "abs_rank_biserial_effect_size",
        "wilcoxon_p",
        "lmm_agrees_with_primary_direction",
        "lmm_raw_p_less_0_05",
        "lmm_fdr_significant",
        "hrv_valid_direction_agrees_with_full_cohort",
        "hrv_valid_wilcoxon_p",
        "prop_direction_retained",
        "prop_delta_sign_retained",
        "n_scenarios_raw_p_less_0_05",
        "n_scenarios_bootstrap_ci_excludes_zero",
        "n_direction_changes",
        "n_delta_sign_changes",
        "n_potentially_influential_patients",
        "prop_agree",
        "prop_disagree",
    ]

    for col in numeric_like_cols:
        if col in scorecard.columns:
            scorecard[col] = pd.to_numeric(scorecard[col], errors="coerce")

    scorecard["final_evidence_score"] = 0

    # Primary patient-level evidence
    scorecard["score_primary_raw_p"] = (
        scorecard["wilcoxon_p"] < 0.05
    ).fillna(False).astype(int) * 2

    scorecard["score_primary_near_p"] = (
        (scorecard["wilcoxon_p"] >= 0.05)
        & (scorecard["wilcoxon_p"] < 0.10)
    ).fillna(False).astype(int) * 1

    scorecard["score_global_fdr"] = (
        scorecard["global_fdr_significant"] == 1
    ).fillna(False).astype(int) * 3

    scorecard["score_family_fdr"] = (
        scorecard["family_fdr_significant"] == 1
    ).fillna(False).astype(int) * 2

    scorecard["score_effect_size"] = (
        scorecard["abs_rank_biserial_effect_size"] >= 0.50
    ).fillna(False).astype(int) * 1

    scorecard["score_direction_consistency"] = (
        scorecard["direction_consistency"]
        >= config.FINAL_SCORECARD_STRONG_DIRECTION_CONSISTENCY
    ).fillna(False).astype(int) * 1

    scorecard["score_bootstrap_primary"] = (
        scorecard["bootstrap_ci_excludes_zero"] == 1
    ).fillna(False).astype(int) * 1

    # Redundancy
    if "is_cluster_representative" in scorecard.columns:
        scorecard["score_redundancy_representative"] = (
            scorecard["is_cluster_representative"] == 1
        ).fillna(False).astype(int) * 1
    else:
        scorecard["score_redundancy_representative"] = 0

    # LMM supportive evidence
    scorecard["score_lmm_direction"] = (
        scorecard["lmm_agrees_with_primary_direction"] == 1
    ).fillna(False).astype(int) * 1

    scorecard["score_lmm_fdr"] = (
        scorecard["lmm_fdr_significant"] == 1
    ).fillna(False).astype(int) * 1

    # Sensitivity
    scorecard["score_sensitivity_direction"] = (
        (scorecard["prop_direction_retained"] >= 0.80)
        & (scorecard["prop_delta_sign_retained"] >= 0.80)
    ).fillna(False).astype(int) * 1

    scorecard["score_sensitivity_statistical"] = (
        (scorecard["n_scenarios_raw_p_less_0_05"] >= 2)
        | (scorecard["n_scenarios_bootstrap_ci_excludes_zero"] >= 2)
    ).fillna(False).astype(int) * 1

    # Leave-one-out
    scorecard["score_loo_stable"] = (
        (scorecard["loo_stability_label"] == "stable")
        & (scorecard["n_direction_changes"] == 0)
        & (scorecard["n_delta_sign_changes"] == 0)
        & (scorecard["n_potentially_influential_patients"] == 0)
    ).fillna(False).astype(int) * 1

    # Patient agreement
    scorecard["score_patient_agreement"] = (
        scorecard["prop_agree"]
        >= config.FINAL_SCORECARD_MIN_AGREEMENT_PROPORTION
    ).fillna(False).astype(int) * 1

    # HRV-valid subgroup
    # Only applies to HRV features. Non-HRV features receive 0, not a penalty.
    scorecard["score_hrv_valid_support"] = 0

    hrv_mask = scorecard["feature_family"] == "HRV"

    scorecard.loc[
        hrv_mask
        & (scorecard["hrv_valid_direction_agrees_with_full_cohort"] == 1)
        & (scorecard["hrv_valid_wilcoxon_p"] < 0.10),
        "score_hrv_valid_support",
    ] = 1

    # Penalties / caution flags
    scorecard["penalty_hrv_valid_weakened"] = 0

    scorecard.loc[
        hrv_mask
        & (
            (scorecard["hrv_valid_direction_agrees_with_full_cohort"] != 1)
            | (scorecard["hrv_valid_wilcoxon_p"] >= 0.10)
            | (scorecard["hrv_valid_wilcoxon_p"].isna())
        ),
        "penalty_hrv_valid_weakened",
    ] = -1

    scorecard["penalty_low_agreement"] = (
        scorecard["prop_agree"] < config.FINAL_SCORECARD_MIN_AGREEMENT_PROPORTION
    ).fillna(False).astype(int) * -1

    scorecard["penalty_no_change_dominant_direction"] = (
        scorecard["dominant_direction"] == "no_change"
    ).fillna(False).astype(int) * -2

    score_components = [
        "score_primary_raw_p",
        "score_primary_near_p",
        "score_global_fdr",
        "score_family_fdr",
        "score_effect_size",
        "score_direction_consistency",
        "score_bootstrap_primary",
        "score_redundancy_representative",
        "score_lmm_direction",
        "score_lmm_fdr",
        "score_sensitivity_direction",
        "score_sensitivity_statistical",
        "score_loo_stable",
        "score_patient_agreement",
        "score_hrv_valid_support",
        "penalty_hrv_valid_weakened",
        "penalty_low_agreement",
        "penalty_no_change_dominant_direction",
    ]

    scorecard["final_evidence_score"] = scorecard[score_components].sum(axis=1)

    # -----------------------------------------------------------------
    # Final categories
    # -----------------------------------------------------------------
    def _final_category(row):
        """
        Assign cautious final exploratory categories.

        Important:
        These are cohort-level exploratory evidence categories, not validation labels.
        """

        feature_family = row.get("feature_family", "")
        score = row.get("final_evidence_score", np.nan)
        dominant_direction = row.get("dominant_direction", "")
        prop_agree = row.get("prop_agree", np.nan)

        hrv_weakened = row.get("penalty_hrv_valid_weakened", 0) < 0
        no_change = dominant_direction == "no_change"

        primary_raw_p_supported = (
            pd.notna(row.get("wilcoxon_p", np.nan))
            and row.get("wilcoxon_p", np.nan) < 0.05
        )

        primary_bootstrap_supported = (
            row.get("bootstrap_ci_excludes_zero", 0) == 1
        )

        primary_patient_level_supported = (
            primary_raw_p_supported or primary_bootstrap_supported
        )

        if feature_family == "HRV" and hrv_weakened:
            return "secondary_rhythm_sensitive_candidate"

        if (
            pd.notna(score)
            and score >= config.FINAL_SCORECARD_PRIMARY_EXPLORATORY_SCORE_THRESHOLD
            and not no_change
            and pd.notna(prop_agree)
            and prop_agree >= config.FINAL_SCORECARD_MIN_AGREEMENT_PROPORTION
            and primary_patient_level_supported
        ):
            return "primary_exploratory_candidate"

        if (
            pd.notna(score)
            and score >= config.FINAL_SCORECARD_SECONDARY_EXPLORATORY_SCORE_THRESHOLD
            and not no_change
        ):
            return "secondary_exploratory_candidate"

        return "exploratory_or_downgraded"

    scorecard["final_candidate_category"] = scorecard.apply(
        _final_category,
        axis=1,
    )

    # -----------------------------------------------------------------
    # Interpretation notes
    # -----------------------------------------------------------------



    # -----------------------------------------------------------------
    # Interpretation notes
    # -----------------------------------------------------------------

    def _interpretation_note(row):
        notes = []

        feature = row["feature"]
        direction = row.get("dominant_direction", "")
        category = row.get("final_candidate_category", "")
        family = row.get("feature_family", "")

        notes.append(
            f"{feature} showed a cohort-level {direction} direction and was classified as {category} in this exploratory cohort."
        )

        if row.get("wilcoxon_p", np.nan) < 0.05:
            notes.append("Primary paired test had raw p < 0.05.")

        if row.get("global_fdr_significant", 0) == 1:
            notes.append("Global FDR support was present.")
        elif row.get("family_fdr_significant", 0) == 1:
            notes.append("Family-level FDR support was present.")
        else:
            notes.append(
                "No global or family-level FDR support was present; interpretation remains exploratory."
            )

        if row.get("bootstrap_ci_excludes_zero", 0) == 1:
            notes.append("Primary bootstrap median-delta CI excluded zero.")

        if row.get("lmm_agrees_with_primary_direction", 0) == 1:
            notes.append("Supportive LMM direction agreed with the primary direction.")

        if row.get("sensitivity_robustness_label", ""):
            notes.append(
                f"Sensitivity label: {row.get('sensitivity_robustness_label')}."
            )

        if row.get("loo_stability_label", ""):
            notes.append(
                f"Leave-one-patient-out label: {row.get('loo_stability_label')}."
            )

        if family == "HRV":
            if row.get("penalty_hrv_valid_weakened", 0) < 0:
                notes.append(
                    "HRV-valid subgroup weakened statistical support; interpret as rhythm-sensitive."
                )
            else:
                notes.append(
                    "HRV-valid subgroup did not weaken direction/statistical support."
                )

        clinical_context = row.get("top_clinical_disagreement_context", "")
        if isinstance(clinical_context, str) and clinical_context.strip():
            notes.append(
                f"Disagreement context: {clinical_context}."
            )

        if direction == "no_change":
            notes.append(
                "Dominant no-change direction limits physiological interpretability."
            )

        return " ".join(notes)

    scorecard["final_interpretation_note"] = scorecard.apply(
        _interpretation_note,
        axis=1,
    )

    # -----------------------------------------------------------------
    # Sort and save final scorecard
    # -----------------------------------------------------------------

    category_order = {
        "primary_exploratory_candidate": 1,
        "secondary_exploratory_candidate": 2,
        "secondary_rhythm_sensitive_candidate": 3,
        "exploratory_or_downgraded": 4,
    }

    scorecard["category_sort_order"] = scorecard["final_candidate_category"].map(
        category_order
    )

    scorecard = scorecard.sort_values(
        [
            "category_sort_order",
            "final_evidence_score",
            "wilcoxon_p",
            "prop_agree",
        ],
        ascending=[True, False, True, False],
    ).drop(columns=["category_sort_order"])

    save_dataframe(
        scorecard,
        config.CANDIDATE_SCORECARD_DIR / "final_candidate_scorecard.csv",
    )

    # -----------------------------------------------------------------
    # Category summary
    # -----------------------------------------------------------------

    category_summary = (
        scorecard
        .groupby("final_candidate_category")
        .agg(
            n_features=("feature", "count"),
            features=("feature", lambda x: "; ".join(x.astype(str))),
        )
        .reset_index()
    )

    category_summary["category_sort_order"] = category_summary[
        "final_candidate_category"
    ].map(category_order)

    category_summary = category_summary.sort_values(
        "category_sort_order"
    ).drop(columns=["category_sort_order"])

    save_dataframe(
        category_summary,
        config.CANDIDATE_SCORECARD_DIR / "final_candidate_feature_categories.csv",
    )

    # -----------------------------------------------------------------
    # Manuscript-friendly compact summary
    # -----------------------------------------------------------------

    manuscript_cols = [
        "feature",
        "feature_family",
        "final_candidate_category",
        "dominant_direction",
        "n_paired_patients",
        "direction_consistency",
        "wilcoxon_p",
        "wilcoxon_p_fdr_global",
        "rank_biserial_effect_size",
        "bootstrap_ci_excludes_zero",
        "sensitivity_robustness_label",
        "loo_stability_label",
        "prop_agree",
        "n_agree",
        "n_disagree",
        "top_clinical_disagreement_context",
        "final_interpretation_note",
    ]

    manuscript_cols = [c for c in manuscript_cols if c in scorecard.columns]

    manuscript_summary = scorecard[manuscript_cols].copy()

    save_dataframe(
        manuscript_summary,
        config.CANDIDATE_SCORECARD_DIR / "final_candidate_feature_summary_for_manuscript.csv",
    )

    logging.info("Stage 16 completed: Final candidate scorecard.")

    return scorecard


# ---------------------------------------------------------------------
# Stage 17 — Final workbook export
# ---------------------------------------------------------------------

def stage17_final_workbook_export() -> Path:
    """
    Export key validated pipeline outputs into one Excel workbook.

    Purpose:
    - Package final candidate scorecard outputs.
    - Include primary statistics, redundancy, HRV subgroup, LMM,
      sensitivity, leave-one-patient-out, and agreement/disagreement summaries.
    - Provide a manuscript-review workbook without creating new statistical results.

    Output:
    - HI_ECG_candidate_feature_analysis_workbook.xlsx
    """

    logging.info("Stage 17 started: Final workbook export.")

    workbook_path = config.FINAL_WORKBOOK_FILE
    workbook_path.parent.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------
    # Helper functions
    # -----------------------------------------------------------------

    def _load_csv(path: Path, required: bool = True) -> pd.DataFrame:
        if path.exists():
            return pd.read_csv(path)

        message = f"Workbook export input missing: {path}"

        if required:
            raise FileNotFoundError(message)

        logging.warning(message)
        return pd.DataFrame({"note": [message]})

    def _safe_sheet_name(name: str) -> str:
        """
        Excel sheet names must be <= 31 characters.
        """
        invalid_chars = ["[", "]", "*", "?", "/", "\\", ":"]
        cleaned = name

        for ch in invalid_chars:
            cleaned = cleaned.replace(ch, "_")

        return cleaned[:31]

    def _write_sheet(
        writer: pd.ExcelWriter,
        df: pd.DataFrame,
        sheet_name: str,
        freeze_panes: tuple[int, int] = (1, 0),
    ) -> None:
        """
        Write a dataframe to an Excel sheet with basic formatting.

        Formatting is intentionally simple to avoid dependency-heavy workbook logic.
        """
        sheet_name = _safe_sheet_name(sheet_name)

        df.to_excel(
            writer,
            sheet_name=sheet_name,
            index=False,
        )

        worksheet = writer.sheets[sheet_name]

        # Freeze header row
        try:
            worksheet.freeze_panes(*freeze_panes)
        except Exception:
            pass

        # Auto-width with caps
        try:
            for idx, col in enumerate(df.columns):
                series = df[col].astype(str)
                max_len = max(
                    [len(str(col))]
                    + series.head(500).map(len).tolist()
                )

                width = min(max(max_len + 2, 10), 45)
                worksheet.set_column(idx, idx, width)
        except Exception:
            pass

        # Header formatting
        try:
            workbook = writer.book
            header_format = workbook.add_format(
                {
                    "bold": True,
                    "text_wrap": True,
                    "valign": "top",
                    "border": 1,
                }
            )

            for col_num, value in enumerate(df.columns.values):
                worksheet.write(0, col_num, value, header_format)
        except Exception:
            pass

        # Apply autofilter
        try:
            if len(df.columns) > 0 and len(df) > 0:
                worksheet.autofilter(0, 0, len(df), len(df.columns) - 1)
        except Exception:
            pass

    # -----------------------------------------------------------------
    # Workbook index / README
    # -----------------------------------------------------------------

    workbook_index = pd.DataFrame(
        [
            {
                "sheet_name": "README",
                "description": "Workbook overview and interpretation cautions.",
            },
            {
                "sheet_name": "Final_Scorecard",
                "description": "Final candidate feature scorecard from Stage 16.",
            },
            {
                "sheet_name": "Final_Categories",
                "description": "Final candidate feature categories.",
            },
            {
                "sheet_name": "Manuscript_Summary",
                "description": "Compact manuscript-friendly feature summary.",
            },
            {
                "sheet_name": "Prelim_Shortlist",
                "description": "Stage 10 preliminary candidate shortlist.",
            },
            {
                "sheet_name": "Primary_Stats",
                "description": "Stage 8 primary paired statistics with FDR correction.",
            },
            {
                "sheet_name": "Redundancy_Clusters",
                "description": "Stage 9 redundancy cluster membership.",
            },
            {
                "sheet_name": "Cluster_Reps",
                "description": "Stage 9 selected cluster representatives.",
            },
            {
                "sheet_name": "HRV_Subgroup",
                "description": "Stage 11 HRV-valid subgroup feature results.",
            },
            {
                "sheet_name": "LMM_Supportive",
                "description": "Stage 12 supportive mixed-effects model results.",
            },
            {
                "sheet_name": "Sensitivity_Robustness",
                "description": "Stage 13 candidate robustness matrix.",
            },
            {
                "sheet_name": "LOO_Summary",
                "description": "Stage 14 leave-one-patient-out summary.",
            },
            {
                "sheet_name": "Feature_Agreement",
                "description": "Stage 15 feature-level patient agreement summary.",
            },
            {
                "sheet_name": "Patient_Disagreement",
                "description": "Stage 15 patient-level disagreement summary.",
            },
            {
                "sheet_name": "Clinical_Disagreement",
                "description": "Stage 15 clinical-flag disagreement summary.",
            },
            {
                "sheet_name": "QC_Counts",
                "description": "Stage 5 QC exclusion counts.",
            },
            {
                "sheet_name": "Cohort_Description",
                "description": "Stage 4 cohort clinical flag summary.",
            },
        ]
    )

    readme = pd.DataFrame(
        [
            {
                "item": "Workbook purpose",
                "description": (
                    "This workbook packages validated outputs from the exploratory "
                    "20-patient ECG hemodynamic instability feature analysis pipeline."
                ),
            },
            {
                "item": "Primary analysis unit",
                "description": (
                    "Patient-level paired PreHI vs HI summaries are the primary "
                    "analysis unit."
                ),
            },
            {
                "item": "Interpretation caution",
                "description": (
                    "Final categories are exploratory cohort-level evidence categories, "
                    "not externally validated diagnostic biomarkers."
                ),
            },
            {
                "item": "FDR interpretation",
                "description": (
                    "No candidate should be described as validated solely from this "
                    "workbook. Global and family-level FDR columns should be reviewed."
                ),
            },
            {
                "item": "HRV interpretation",
                "description": (
                    "HRV features should be interpreted as rhythm-sensitive because "
                    "support weakened in the HRV-valid subgroup."
                ),
            },
            {
                "item": "Clinical disagreement",
                "description": (
                    "Clinical-flag disagreement summaries are descriptive only. "
                    "Small subgroup counts should not be interpreted as formal "
                    "association testing."
                ),
            },
        ]
    )

    # -----------------------------------------------------------------
    # Load validated outputs
    # -----------------------------------------------------------------

    tables = {
        "Final_Scorecard": _load_csv(
            config.CANDIDATE_SCORECARD_DIR / "final_candidate_scorecard.csv"
        ),
        "Final_Categories": _load_csv(
            config.CANDIDATE_SCORECARD_DIR / "final_candidate_feature_categories.csv"
        ),
        "Manuscript_Summary": _load_csv(
            config.CANDIDATE_SCORECARD_DIR / "final_candidate_feature_summary_for_manuscript.csv"
        ),
        "Prelim_Shortlist": _load_csv(
            config.CANDIDATE_SCORECARD_DIR / "preliminary_candidate_shortlist.csv"
        ),
        "Primary_Stats": _load_csv(
            config.PRIMARY_STATS_DIR / "primary_paired_feature_results_all_features_fdr.csv"
        ),
        "Redundancy_Clusters": _load_csv(
            config.REDUNDANCY_DIR / "redundancy_clusters.csv"
        ),
        "Cluster_Reps": _load_csv(
            config.REDUNDANCY_DIR / "cluster_representatives.csv"
        ),
        "HRV_Subgroup": _load_csv(
            config.HRV_SUBGROUP_DIR / "hrv_valid_feature_results.csv"
        ),
        "LMM_Supportive": _load_csv(
            config.MIXED_EFFECTS_DIR / "lmm_supportive_results_all_features.csv"
        ),
        "Sensitivity_Robustness": _load_csv(
            config.SENSITIVITY_DIR / "candidate_feature_robustness_matrix.csv"
        ),
        "LOO_Summary": _load_csv(
            config.LOO_DIR / "leave_one_patient_out_summary.csv"
        ),
        "Feature_Agreement": _load_csv(
            config.AGREEMENT_DIR / "feature_agreement_summary.csv"
        ),
        "Patient_Disagreement": _load_csv(
            config.AGREEMENT_DIR / "patient_disagreement_summary.csv"
        ),
        "Clinical_Disagreement": _load_csv(
            config.AGREEMENT_DIR / "clinical_flag_disagreement_summary.csv"
        ),
        "QC_Counts": _load_csv(
            config.SEGMENT_QC_DIR / "qc_exclusion_counts.csv"
        ),
        "Cohort_Description": _load_csv(
            config.CLINICAL_METADATA_DIR / "cohort_description_table.csv"
        ),
    }

    # -----------------------------------------------------------------
    # Write workbook
    # -----------------------------------------------------------------

    with pd.ExcelWriter(workbook_path, engine="xlsxwriter") as writer:
        _write_sheet(writer, readme, "README")
        _write_sheet(writer, workbook_index, "Workbook_Index")

        for sheet_name, df in tables.items():
            _write_sheet(writer, df, sheet_name)

    # -----------------------------------------------------------------
    # Export manifest
    # -----------------------------------------------------------------

    export_manifest = pd.DataFrame(
        [
            {
                "output_file": str(workbook_path),
                "n_sheets": len(tables) + 2,
                "export_note": (
                    "Final workbook created from validated pipeline outputs. "
                    "No new statistical analysis was performed in Stage 17."
                ),
            }
        ]
    )

    save_dataframe(
        export_manifest,
        config.FINAL_WORKBOOK_DIR / "final_workbook_export_manifest.csv",
    )

    logging.info(f"Stage 17 completed: Workbook exported to {workbook_path}")

    return workbook_path