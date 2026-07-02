"""
Utility functions for the ECG-HI exploratory feature discovery pipeline.
"""

import hashlib
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Iterable, Dict, Any, List

import numpy as np
import pandas as pd

from scipy import stats
from statsmodels.stats.multitest import multipletests

import config


# ---------------------------------------------------------------------
# Logging and folder setup
# ---------------------------------------------------------------------

def create_output_dirs() -> None:
    """Create all configured output directories."""
    for folder in config.OUTPUT_DIRS:
        folder.mkdir(parents=True, exist_ok=True)


def setup_logging() -> None:
    """Set up console and file logging."""
    create_output_dirs()

    log_file = config.OUTPUTS_DIR / "analysis_run_log.txt"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode="w"),
            logging.StreamHandler(),
        ],
    )

    logging.info("Logging initialized.")
    logging.info(f"Project root: {config.PROJECT_ROOT}")


# ---------------------------------------------------------------------
# File validation and manifest
# ---------------------------------------------------------------------

def file_checksum(path: Path) -> str:
    """Return SHA256 checksum for a file."""
    sha256 = hashlib.sha256()

    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8192), b""):
            sha256.update(block)

    return sha256.hexdigest()


def validate_input_files() -> pd.DataFrame:
    """
    Validate required and optional input files.

    Returns
    -------
    pd.DataFrame
        File validation table.
    """

    required_files = [
        config.DATA_FILE,
        config.CLINICAL_FLAGS_FILE,
        config.CLINICAL_SUMMARY_FILE,
        config.DATA_DICTIONARY_FILE,
        config.FEATURE_DOCUMENTATION_FILE,
    ]

    optional_files = [
        config.CLINICAL_REPORT_FILE,
        config.DATA_DICTIONARY_MD_FILE,
    ]

    rows = []

    for path in required_files:
        rows.append(
            {
                "file": str(path),
                "required": True,
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else np.nan,
                "sha256": file_checksum(path) if path.exists() else None,
            }
        )

    for path in optional_files:
        rows.append(
            {
                "file": str(path),
                "required": False,
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else np.nan,
                "sha256": file_checksum(path) if path.exists() else None,
            }
        )

    validation_df = pd.DataFrame(rows)
    validation_df.to_csv(
        config.DATA_AUDIT_DIR / "input_file_validation.csv",
        index=False,
    )

    missing_required = validation_df[
        (validation_df["required"] == True) & (validation_df["exists"] == False)
    ]

    if not missing_required.empty:
        missing_files = missing_required["file"].tolist()
        raise FileNotFoundError(
            "Missing required input files:\n" + "\n".join(missing_files)
        )

    logging.info("Input file validation completed.")
    return validation_df


def write_run_manifest(extra: Dict[str, Any] | None = None) -> None:
    """Write a simple run manifest."""
    manifest = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "project_root": str(config.PROJECT_ROOT),
        "data_file": str(config.DATA_FILE),
        "clinical_flags_file": str(config.CLINICAL_FLAGS_FILE),
        "clinical_summary_file": str(config.CLINICAL_SUMMARY_FILE),
        "data_dictionary_file": str(config.DATA_DICTIONARY_FILE),
        "feature_documentation_file": str(config.FEATURE_DOCUMENTATION_FILE),
        "random_seed": config.RANDOM_SEED,
        "min_paired_patients_per_feature": config.MIN_PAIRED_PATIENTS_PER_FEATURE,
        "moderate_missingness_threshold": config.MODERATE_MISSINGNESS_THRESHOLD,
        "high_missingness_threshold": config.HIGH_MISSINGNESS_THRESHOLD,
    }

    if extra:
        manifest.update(extra)

    with open(config.OUTPUTS_DIR / "run_manifest.json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=4)

    logging.info("Run manifest written.")


# ---------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------

def standardize_condition_label(value) -> str | None:
    """
    Convert condition labels to standardized values:
    PreHI or HI.

    Returns None if label cannot be mapped.
    """

    if pd.isna(value):
        return None

    if value in config.CONDITION_LABEL_MAP:
        return config.CONDITION_LABEL_MAP[value]

    key = str(value).strip().lower()
    key = key.replace("_", " ").replace("-", " ")
    key = " ".join(key.split())

    compact_key = key.replace(" ", "")

    if key in config.CONDITION_LABEL_MAP:
        return config.CONDITION_LABEL_MAP[key]

    if compact_key in config.CONDITION_LABEL_MAP:
        return config.CONDITION_LABEL_MAP[compact_key]

    return None


def safe_to_numeric(series: pd.Series) -> pd.Series:
    """Convert a pandas Series to numeric while preserving non-convertible as NaN."""
    return pd.to_numeric(series, errors="coerce")


def count_paired_patients(
    df: pd.DataFrame,
    feature: str,
    patient_col: str = config.PATIENT_ID_COL,
    condition_col: str = config.CONDITION_STD_COL,
) -> int:
    """
    Count patients with both PreHI and HI non-missing medians for a feature.
    """

    tmp = (
        df[[patient_col, condition_col, feature]]
        .dropna(subset=[feature])
        .groupby([patient_col, condition_col])[feature]
        .median()
        .reset_index()
    )

    wide = tmp.pivot(index=patient_col, columns=condition_col, values=feature)

    if "PreHI" not in wide.columns or "HI" not in wide.columns:
        return 0

    paired = wide[["PreHI", "HI"]].dropna()
    return int(len(paired))


def infer_feature_family(feature_name: str) -> str:
    """
    Infer broad ECG feature family from feature name.

    This is intentionally transparent and rule-based.
    Later we can improve this using the data dictionary if needed.
    """

    f = feature_name.lower()

    if f.startswith("hrv_"):
        return "HRV"

    if f in {"r_amp_mean", "r_amp_std", "qrs_mean", "qrs_std", "signal_energy"}:
        return "ECG morphology and signal energy"

    if f.startswith("entropyprofiled_"):
        return "entropy-profile aggregate features"

    if f.startswith("fd_"):
        return "fractal-dimension summary features"

    entropy_terms = [
        "entropy",
        "lempelziv",
        "complexity",
        "approximateentropy",
        "sampleentropy",
        "permutationentropy",
        "fuzzyentropy",
        "distributionentropy",
        "shannonentropy",
        "renyientropy",
        "singularvaluedecompositionentropy",
    ]

    if any(term in f for term in entropy_terms):
        return "entropy and complexity"

    fractal_terms = [
        "hjorth",
        "fisher",
        "petrosian",
        "katz",
        "higuchi",
        "detrended",
    ]

    if any(term in f for term in fractal_terms):
        return "fractal and nonlinear dynamics"

    spectral_terms = [
        "spectralentropy",
        "bandpower",
    ]

    if any(term in f for term in spectral_terms):
        return "spectral features"

    time_domain_terms = [
        "maximum",
        "minimum",
        "mean",
        "median",
        "standarddeviation",
        "variance",
        "kurtosis",
        "skewness",
        "numberofzerocrossing",
        "positivetonegativesampleratio",
        "positivetonegativepeakratio",
        "meanabsolutevalue",
    ]

    if f in time_domain_terms:
        return "time-domain statistics"

    return "unknown/unmapped"


def save_dataframe(df: pd.DataFrame, path: Path) -> None:
    """Save DataFrame as CSV with folder creation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    logging.info(f"Saved: {path}")
    
def standardize_patient_id_series(series: pd.Series) -> pd.Series:
    """
    Standardize patient IDs as clean strings.

    Handles values that may be read as integers, floats, or strings.
    Example:
        491      -> "491"
        491.0    -> "491"
        "491 "   -> "491"
    """

    def _clean_one(x):
        if pd.isna(x):
            return None

        text = str(x).strip()

        if text.endswith(".0"):
            text = text[:-2]

        return text

    return series.apply(_clean_one)


def to_binary_flag(series: pd.Series) -> pd.Series:
    """
    Convert a column to conservative binary 0/1 flags.

    Recognizes 1/0, True/False, yes/no, y/n.
    Missing values become 0.
    """

    def _convert(x):
        if pd.isna(x):
            return 0

        if isinstance(x, bool):
            return int(x)

        text = str(x).strip().lower()

        if text in {"1", "1.0", "true", "yes", "y"}:
            return 1

        if text in {"0", "0.0", "false", "no", "n", "none", "nan", ""}:
            return 0

        # Non-empty text in a flag-like column is treated cautiously as present.
        return 1

    return series.apply(_convert).astype(int)


def robust_modified_zscore(series: pd.Series) -> pd.Series:
    """
    Compute robust modified z-score using median absolute deviation.

    Returns zeros if MAD is zero or undefined.
    """

    x = pd.to_numeric(series, errors="coerce")
    median = x.median(skipna=True)
    mad = (x - median).abs().median(skipna=True)

    if pd.isna(mad) or mad == 0:
        return pd.Series(np.zeros(len(series)), index=series.index)

    return 0.6745 * (x - median) / mad

def clean_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return a copy of a DataFrame with stripped column names.

    This prevents hidden errors from trailing spaces or accidental whitespace
    in CSV headers.
    """
    out = df.copy()
    out.columns = [str(c).strip() for c in out.columns]
    return out

def bootstrap_median_ci(
    values: pd.Series,
    n_iterations: int = config.BOOTSTRAP_N_ITERATIONS,
    confidence_level: float = config.BOOTSTRAP_CONFIDENCE_LEVEL,
    random_seed: int = config.RANDOM_SEED,
) -> tuple[float, float]:
    """
    Bootstrap confidence interval for the median of paired deltas.

    Parameters
    ----------
    values : pd.Series
        Numeric paired delta values.
    n_iterations : int
        Number of bootstrap resamples.
    confidence_level : float
        Confidence level, e.g. 0.95.
    random_seed : int
        Random seed for reproducibility.

    Returns
    -------
    tuple
        Lower and upper confidence interval bounds.
    """

    x = pd.to_numeric(values, errors="coerce").dropna().values

    if len(x) == 0:
        return np.nan, np.nan

    if len(x) == 1:
        return float(x[0]), float(x[0])

    rng = np.random.default_rng(random_seed)

    boot_medians = np.empty(n_iterations)

    for i in range(n_iterations):
        sample = rng.choice(x, size=len(x), replace=True)
        boot_medians[i] = np.median(sample)

    alpha = 1 - confidence_level
    lower = np.quantile(boot_medians, alpha / 2)
    upper = np.quantile(boot_medians, 1 - alpha / 2)

    return float(lower), float(upper)


def wilcoxon_signed_rank_pvalue(
    pre_values: pd.Series,
    hi_values: pd.Series,
) -> float:
    """
    Compute Wilcoxon signed-rank p-value for paired PreHI and HI values.

    Returns NaN if test cannot be computed.
    """

    pre = pd.to_numeric(pre_values, errors="coerce")
    hi = pd.to_numeric(hi_values, errors="coerce")

    paired = pd.DataFrame({"PreHI": pre, "HI": hi}).dropna()

    if len(paired) < 2:
        return np.nan

    deltas = paired["HI"] - paired["PreHI"]

    # Wilcoxon cannot run if all paired differences are zero.
    if np.allclose(deltas.values, 0):
        return 1.0

    try:
        result = stats.wilcoxon(
            paired["HI"],
            paired["PreHI"],
            zero_method=config.WILCOXON_ZERO_METHOD,
            alternative=config.WILCOXON_ALTERNATIVE,
            mode="auto",
        )
        return float(result.pvalue)

    except Exception:
        return np.nan


def matched_pairs_rank_biserial_effect_size(
    pre_values: pd.Series,
    hi_values: pd.Series,
) -> float:
    """
    Compute matched-pairs rank-biserial effect size.

    This is calculated from signed ranks of non-zero paired differences:

        r_rb = (sum_positive_ranks - sum_negative_ranks) / total_rank_sum

    Positive values indicate HI > PreHI.
    Negative values indicate HI < PreHI.
    """

    pre = pd.to_numeric(pre_values, errors="coerce")
    hi = pd.to_numeric(hi_values, errors="coerce")

    paired = pd.DataFrame({"PreHI": pre, "HI": hi}).dropna()

    if len(paired) < 2:
        return np.nan

    deltas = paired["HI"] - paired["PreHI"]
    nonzero = deltas[deltas != 0]

    if len(nonzero) == 0:
        return 0.0

    abs_ranks = stats.rankdata(np.abs(nonzero.values))

    positive_rank_sum = abs_ranks[nonzero.values > 0].sum()
    negative_rank_sum = abs_ranks[nonzero.values < 0].sum()
    total_rank_sum = abs_ranks.sum()

    if total_rank_sum == 0:
        return 0.0

    effect = (positive_rank_sum - negative_rank_sum) / total_rank_sum

    return float(effect)


def apply_fdr_correction(
    pvalues: pd.Series,
    method: str = config.FDR_METHOD,
) -> pd.Series:
    """
    Apply FDR correction to a vector of p-values.

    Missing p-values remain missing.
    """

    p = pd.to_numeric(pvalues, errors="coerce")
    adjusted = pd.Series(np.nan, index=p.index, dtype=float)

    valid_mask = p.notna()

    if valid_mask.sum() == 0:
        return adjusted

    _, p_adjusted, _, _ = multipletests(
        p.loc[valid_mask],
        alpha=config.FDR_ALPHA,
        method=method,
    )

    adjusted.loc[valid_mask] = p_adjusted

    return adjusted


def direction_from_delta(delta: float, tolerance: float = 0.0) -> str:
    """
    Convert numeric delta to direction label.
    """

    if pd.isna(delta):
        return "missing"

    if delta > tolerance:
        return "increase"

    if delta < -tolerance:
        return "decrease"

    return "no_change"


def dominant_direction_from_counts(
    n_increase: int,
    n_decrease: int,
    n_no_change: int = 0,
) -> str:
    """
    Determine dominant direction from patient-level delta counts.
    """

    if n_increase > n_decrease and n_increase > n_no_change:
        return "increase"

    if n_decrease > n_increase and n_decrease > n_no_change:
        return "decrease"

    if n_no_change > n_increase and n_no_change > n_decrease:
        return "no_change"

    return "tie_or_mixed"


def get_connected_components_from_edges(edges: list[tuple[str, str]]) -> list[list[str]]:
    """
    Build connected components from an edge list.

    This avoids requiring networkx.
    Each edge is a tuple: (feature_a, feature_b).
    """

    adjacency = {}

    for a, b in edges:
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)

    visited = set()
    components = []

    for node in adjacency:
        if node in visited:
            continue

        stack = [node]
        component = []

        while stack:
            current = stack.pop()

            if current in visited:
                continue

            visited.add(current)
            component.append(current)

            for neighbor in adjacency.get(current, []):
                if neighbor not in visited:
                    stack.append(neighbor)

        components.append(sorted(component))

    return components


def min_max_normalize(series: pd.Series, higher_is_better: bool = True) -> pd.Series:
    """
    Min-max normalize a numeric series to 0-1.

    Missing values become 0 after normalization.
    """

    x = pd.to_numeric(series, errors="coerce")

    if x.notna().sum() == 0:
        return pd.Series(0.0, index=series.index)

    min_val = x.min()
    max_val = x.max()

    if pd.isna(min_val) or pd.isna(max_val) or max_val == min_val:
        normalized = pd.Series(0.0, index=series.index)
    else:
        normalized = (x - min_val) / (max_val - min_val)

    if not higher_is_better:
        normalized = 1 - normalized

    return normalized.fillna(0.0)


def safe_negative_log10_pvalue(series: pd.Series) -> pd.Series:
    """
    Convert p-values to -log10(p), safely.

    Smaller p-values get larger scores.
    """

    p = pd.to_numeric(series, errors="coerce")
    p = p.clip(lower=1e-300)

    out = -np.log10(p)
    out = out.replace([np.inf, -np.inf], np.nan).fillna(0.0)

    return out

def zscore_series(series: pd.Series) -> pd.Series:
    """
    Z-score a numeric pandas Series.

    Returns NaN if standard deviation is zero or undefined.
    """

    x = pd.to_numeric(series, errors="coerce")
    mean = x.mean(skipna=True)
    std = x.std(skipna=True)

    if pd.isna(std) or std <= 0:
        return pd.Series(np.nan, index=series.index)

    return (x - mean) / std

