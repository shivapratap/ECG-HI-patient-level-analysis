
"""
HI ECG exploratory re-analysis pipeline.

Primary statistical unit: subject.
Primary analysis: paired within-subject changes (before_HI vs during_HI).
Supportive analysis: linear mixed-effects models with patient random intercept.

This module intentionally does NOT perform segment-wise random cross-validation,
because windows from the same subject are repeated and temporally dependent.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence
import warnings

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, rankdata
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt


META_COLUMNS = ("subject_id", "label")


@dataclass(frozen=True)
class AnalysisConfig:
    subject_col: str = "subject_id"
    label_col: str = "label"
    before_label: int = 0
    during_label: int = 1
    alpha: float = 0.05
    aggregation: str = "median"


def load_processed_data(path: str | Path, config: AnalysisConfig = AnalysisConfig()) -> pd.DataFrame:
    """Load processed feature data and standardize labels."""
    df = pd.read_csv(path)

    # Support either numeric or string labels.
    if df[config.label_col].dtype == object:
        mapping = {
            "before_HI": config.before_label,
            "during_HI": config.during_label,
            "0": config.before_label,
            "1": config.during_label,
        }
        df[config.label_col] = df[config.label_col].astype(str).map(mapping)

    df[config.subject_col] = df[config.subject_col].astype(str)
    return df


def feature_columns(df: pd.DataFrame, config: AnalysisConfig = AnalysisConfig()) -> list[str]:
    """Return numeric feature columns, excluding subject and label."""
    excluded = {config.subject_col, config.label_col}
    return [
        c for c in df.columns
        if c not in excluded and pd.api.types.is_numeric_dtype(df[c])
    ]


def validate_data(df: pd.DataFrame, config: AnalysisConfig = AnalysisConfig()) -> dict:
    """Run hard checks and return a compact quality-control report."""
    required = {config.subject_col, config.label_col}
    missing_required = required.difference(df.columns)
    if missing_required:
        raise ValueError(f"Missing required columns: {sorted(missing_required)}")

    allowed = {config.before_label, config.during_label}
    observed = set(df[config.label_col].dropna().unique())
    if not observed.issubset(allowed):
        raise ValueError(f"Unexpected labels: {sorted(observed - allowed)}")

    features = feature_columns(df, config)
    if not features:
        raise ValueError("No numeric feature columns found.")

    counts = pd.crosstab(df[config.subject_col], df[config.label_col])
    for lab in (config.before_label, config.during_label):
        if lab not in counts.columns:
            counts[lab] = 0

    subjects_missing_a_state = counts[
        (counts[config.before_label] == 0) | (counts[config.during_label] == 0)
    ].index.tolist()

    constant_features = [
        c for c in features if df[c].nunique(dropna=True) <= 1
    ]

    report = {
        "n_rows": len(df),
        "n_subjects": df[config.subject_col].nunique(),
        "n_features": len(features),
        "label_counts": df[config.label_col].value_counts().sort_index().to_dict(),
        "missing_values_total": int(df.isna().sum().sum()),
        "duplicate_rows": int(df.duplicated().sum()),
        "subjects_missing_a_state": subjects_missing_a_state,
        "constant_features": constant_features,
        "subject_by_label": counts.sort_index(),
    }

    if subjects_missing_a_state:
        raise ValueError(
            "Every subject must contribute both states for paired analysis. "
            f"Problem subjects: {subjects_missing_a_state}"
        )
    return report


def cohort_summary(df: pd.DataFrame, config: AnalysisConfig = AnalysisConfig()) -> pd.DataFrame:
    """Number of windows per subject and condition."""
    out = (
        df.groupby([config.subject_col, config.label_col])
          .size()
          .unstack(fill_value=0)
          .rename(columns={
              config.before_label: "n_before_HI",
              config.during_label: "n_during_HI",
          })
    )
    out["n_total"] = out.sum(axis=1)
    return out.reset_index()


def aggregate_subject_state(
    df: pd.DataFrame,
    config: AnalysisConfig = AnalysisConfig(),
    features: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Aggregate repeated windows to one value per subject per condition.

    Median is the default because many ECG features are skewed and contain outliers.
    """
    if features is None:
        features = feature_columns(df, config)

    agg = config.aggregation
    if agg not in {"median", "mean"}:
        raise ValueError("aggregation must be 'median' or 'mean'")

    return (
        df.groupby([config.subject_col, config.label_col], as_index=False)[list(features)]
          .agg(agg)
    )


def _matched_rank_biserial(differences: np.ndarray) -> float:
    """
    Matched-pairs rank-biserial correlation.
    Positive values indicate larger values during HI.
    """
    d = np.asarray(differences, dtype=float)
    d = d[np.isfinite(d) & (d != 0)]
    if len(d) == 0:
        return np.nan
    ranks = rankdata(np.abs(d), method="average")
    total = ranks.sum()
    return float((ranks[d > 0].sum() - ranks[d < 0].sum()) / total)

def paired_feature_analysis(
    subject_state: pd.DataFrame,
    config: AnalysisConfig = AnalysisConfig(),
    features: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Primary patient-level paired analysis.

    For each feature, this function compares one aggregated before-HI value
    and one aggregated during-HI value per subject.

    Reported outputs include:
    - median before-HI value
    - median during-HI value
    - median paired difference
    - median percentage change, where defined
    - numbers of subjects increasing, decreasing, or unchanged
    - direction consistency
    - matched-pairs rank-biserial effect size
    - Wilcoxon signed-rank test
    - Benjamini-Hochberg FDR-adjusted p-value

    Notes
    -----
    Percentage change is undefined when the before-HI value is zero.
    Wilcoxon testing is set to p = 1.0 when all paired differences are zero.
    """

    if features is None:
        features = [
            column
            for column in subject_state.columns
            if column not in {config.subject_col, config.label_col}
        ]

    rows = []

    for feature in features:
        # Convert the long-format subject/state data into paired wide format.
        wide = (
            subject_state.pivot(
                index=config.subject_col,
                columns=config.label_col,
                values=feature,
            )
            .dropna(
                subset=[
                    config.before_label,
                    config.during_label,
                ]
            )
        )

        before = wide[config.before_label].to_numpy(dtype=float)
        during = wide[config.during_label].to_numpy(dtype=float)

        # Paired within-subject change.
        diff = during - before

        # ---------------------------------------------------------
        # Safe Wilcoxon signed-rank test
        # ---------------------------------------------------------
        finite_diff = diff[np.isfinite(diff)]
        nonzero_diff = finite_diff[finite_diff != 0]

        if len(finite_diff) == 0:
            wilcoxon_statistic = np.nan
            p_value = np.nan

        elif len(nonzero_diff) == 0:
            # Every subject has exactly zero change.
            wilcoxon_statistic = 0.0
            p_value = 1.0

        else:
            try:
                wilcoxon_statistic, p_value = wilcoxon(
                    during,
                    before,
                    alternative="two-sided",
                    zero_method="wilcox",
                    method="auto",
                )
            except ValueError:
                wilcoxon_statistic = np.nan
                p_value = np.nan

        # ---------------------------------------------------------
        # Safe percentage-change calculation
        # ---------------------------------------------------------
        pct_change = np.full(diff.shape, np.nan, dtype=float)

        valid_baseline = (
            np.isfinite(before)
            & np.isfinite(diff)
            & (np.abs(before) > np.finfo(float).eps)
        )

        pct_change[valid_baseline] = (
            100.0
            * diff[valid_baseline]
            / np.abs(before[valid_baseline])
        )

        if np.any(np.isfinite(pct_change)):
            median_percent_change = float(np.nanmedian(pct_change))
        else:
            median_percent_change = np.nan

        # ---------------------------------------------------------
        # Direction counts and consistency
        # ---------------------------------------------------------
        n_increased = int(np.sum(diff > 0))
        n_decreased = int(np.sum(diff < 0))
        n_unchanged = int(np.sum(diff == 0))

        n_valid_differences = int(np.sum(np.isfinite(diff)))

        if n_valid_differences > 0:
            direction_consistency = max(
                n_increased / n_valid_differences,
                n_decreased / n_valid_differences,
            )
        else:
            direction_consistency = np.nan

        rows.append(
            {
                "feature": feature,
                "n_subjects": len(wide),
                "median_before": (
                    float(np.nanmedian(before))
                    if np.any(np.isfinite(before))
                    else np.nan
                ),
                "median_during": (
                    float(np.nanmedian(during))
                    if np.any(np.isfinite(during))
                    else np.nan
                ),
                "median_difference": (
                    float(np.nanmedian(diff))
                    if np.any(np.isfinite(diff))
                    else np.nan
                ),
                "median_percent_change": median_percent_change,
                "n_increased": n_increased,
                "n_decreased": n_decreased,
                "n_unchanged": n_unchanged,
                "n_nonzero_differences": int(np.sum(diff != 0)),
                "n_valid_percent_changes": int(
                    np.sum(np.isfinite(pct_change))
                ),
                "direction_consistency": direction_consistency,
                "rank_biserial": _matched_rank_biserial(diff),
                "wilcoxon_statistic": wilcoxon_statistic,
                "p_value": p_value,
            }
        )

    result = pd.DataFrame(rows)

    # Use p = 1.0 for FDR correction when a p-value is undefined.
    p_values_for_fdr = result["p_value"].fillna(1.0)

    result["p_fdr_bh"] = multipletests(
        p_values_for_fdr,
        alpha=config.alpha,
        method="fdr_bh",
    )[1]

    result["significant_fdr"] = (
        result["p_fdr_bh"] < config.alpha
    )

    result = result.sort_values(
        by=[
            "significant_fdr",
            "direction_consistency",
            "p_fdr_bh",
        ],
        ascending=[
            False,
            False,
            True,
        ],
    ).reset_index(drop=True)

    return result

def fit_lmm_per_feature(
    df: pd.DataFrame,
    config: AnalysisConfig = AnalysisConfig(),
    features: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Supportive segment-level mixed-effects analysis.

    Model for each feature:
        standardized_feature ~ HI_condition + (1 | subject)

    The outcome is z-standardized to make beta coefficients comparable.
    With only six subjects, treat LMM inference as supportive and report
    convergence diagnostics.
    """
    if features is None:
        features = feature_columns(df, config)

    rows = []
    work = df.copy()
    work["_condition"] = (work[config.label_col] == config.during_label).astype(int)

    for feature in features:
        tmp = work[[config.subject_col, "_condition", feature]].dropna().copy()
        sd = tmp[feature].std(ddof=0)

        if not np.isfinite(sd) or sd == 0:
            rows.append({
                "feature": feature,
                "beta_condition_std": np.nan,
                "se": np.nan,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p_value": np.nan,
                "converged": False,
                "warning": "constant_or_invalid_feature",
            })
            continue

        tmp["_z"] = (tmp[feature] - tmp[feature].mean()) / sd

        warning_text = ""
        converged = False
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                model = smf.mixedlm(
                    "_z ~ _condition",
                    data=tmp,
                    groups=tmp[config.subject_col],
                    re_formula="1",
                )
                fit = model.fit(reml=False, method="lbfgs", maxiter=500, disp=False)
                warning_text = " | ".join(str(w.message) for w in caught)
            beta = fit.params.get("_condition", np.nan)
            se = fit.bse.get("_condition", np.nan)
            ci = fit.conf_int().loc["_condition"].to_numpy()
            p = fit.pvalues.get("_condition", np.nan)
            converged = bool(getattr(fit, "converged", False))
        except Exception as exc:
            beta = se = p = np.nan
            ci = np.array([np.nan, np.nan])
            warning_text = repr(exc)

        rows.append({
            "feature": feature,
            "beta_condition_std": beta,
            "se": se,
            "ci_low": ci[0],
            "ci_high": ci[1],
            "p_value": p,
            "converged": converged,
            "warning": warning_text,
        })

    result = pd.DataFrame(rows)
    valid_p = result["p_value"].fillna(1.0)
    result["p_fdr_bh"] = multipletests(
        valid_p, alpha=config.alpha, method="fdr_bh"
    )[1]
    result["significant_fdr"] = result["p_fdr_bh"] < config.alpha
    return result.sort_values(
        ["significant_fdr", "p_fdr_bh"],
        ascending=[False, True],
    ).reset_index(drop=True)


def within_subject_centered_correlations(
    df: pd.DataFrame,
    config: AnalysisConfig = AnalysisConfig(),
    features: Sequence[str] | None = None,
    method: str = "spearman",
) -> pd.DataFrame:
    """
    Correlation after removing each subject's mean.

    This is preferable to pooled raw correlation when the goal is to examine
    within-subject redundancy. It still treats windows as repeated observations,
    so use it for descriptive redundancy mapping, not formal inference.
    """
    if features is None:
        features = feature_columns(df, config)

    centered = df[list(features)].copy()
    centered = centered - df.groupby(config.subject_col)[list(features)].transform("mean")
    return centered.corr(method=method)


def plot_paired_feature(
    subject_state: pd.DataFrame,
    feature: str,
    output_path: str | Path | None = None,
    config: AnalysisConfig = AnalysisConfig(),
):
    """Patient-linked before-vs-during plot."""
    wide = subject_state.pivot(
        index=config.subject_col,
        columns=config.label_col,
        values=feature,
    ).dropna()

    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for subject, row in wide.iterrows():
        ax.plot([0, 1], [row[config.before_label], row[config.during_label]],
                marker="o", linewidth=1)
    ax.set_xticks([0, 1], ["Before HI", "During HI"])
    ax.set_ylabel(feature)
    ax.set_title(f"Within-subject change: {feature}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig, ax


def save_tables(
    output_dir: str | Path,
    cohort: pd.DataFrame,
    paired: pd.DataFrame,
    lmm: pd.DataFrame | None = None,
) -> None:
    """Save analysis tables as CSV files."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    cohort.to_csv(out / "cohort_window_counts.csv", index=False)
    paired.to_csv(out / "paired_subject_level_feature_results.csv", index=False)
    if lmm is not None:
        lmm.to_csv(out / "supportive_lmm_results.csv", index=False)
