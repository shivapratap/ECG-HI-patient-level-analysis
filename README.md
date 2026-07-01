# ECG-Derived Candidate Features Associated With Hemodynamic Instability in ICU Patients

This repository contains the analysis notebook, processed ECG feature dataset, curated clinical metadata, data dictionaries, intermediate analysis outputs, and figures supporting the manuscript:

**ECG-Derived Candidate Features Associated With Hemodynamic Instability in ICU Patients: A Patient-Level Exploratory Analysis of Ultra-Short Segments**

## Overview

This project investigates ECG-derived candidate features associated with hemodynamic instability (HI) in critically ill ICU patients.

The analysis uses ultra-short 2-minute Lead II ECG segments from a curated 20-patient ICU cohort. Each ECG segment is labelled as either:

- **Pre-HI**: segment preceding hemodynamic instability
- **During-HI**: segment recorded during hemodynamic instability

The primary goal is **exploratory candidate feature discovery**, not prospective prediction or clinical validation.

The patient, rather than the individual ECG segment, is treated as the primary statistical unit. For each ECG-derived feature, patient-level median values are computed separately for Pre-HI and During-HI segments, and paired patient-level differences are analysed.

## Repository structure

```text
.
├── README.md
├── LICENSE
├── requirements.txt
├── ecg_hi_patient_level_exploratory_analysis.ipynb
│
├── data/
│   └── processed/
│       └── 20PatientsFinalData.xlsx
│
├── metadata/
│   ├── README.md
│   ├── ECG_Clinical_Summary_Report.docx
│   ├── patient_clinical_summary.csv
│   └── patient_clinical_flags.csv
│
├── docs/
│   ├── data_dictionary/
│   │   ├── 20PatientsFinalData_data_dictionary.md
│   │   └── 20PatientsFinalData_data_dictionary.csv
│   │
│   └── feature_documentation/
│       └── ecg_lead_ii_feature_list.md
│
└── analysis_outputs/
    ├── tables/
    └── figures/
```

## Main notebook

The complete analysis is provided in:

```text
ecg_hi_patient_level_exploratory_analysis.ipynb
```

The notebook performs the full analysis workflow, including:

1. data loading and integrity checks,
2. feature-family mapping,
3. clinical and signal-quality flagging,
4. cohort description,
5. primary patient-level paired analysis,
6. HRV-valid subgroup analysis,
7. supportive mixed-effects modelling,
8. redundancy analysis,
9. sensitivity analyses,
10. descriptive Episode_Label analysis,
11. candidate feature scoring and ranking,
12. generation of manuscript-facing tables and figures.

## Data

The main processed input dataset is stored in:

```text
data/processed/20PatientsFinalData.xlsx
```

This file contains the ECG-derived feature dataset used by the notebook.

Dataset summary:

| Item | Value |
|---|---:|
| Patients | 20 |
| ECG segment length | 2 minutes |
| ECG lead | Lead II |
| Rows | 1,260 |
| Columns | 62 |
| Metadata columns | 6 |
| ECG-derived feature columns | 56 |

Each row corresponds to one 2-minute ECG segment from one patient.

The six metadata columns are:

```text
Sl_No
Segment_Number
Patient_id
Class
Condition
Episode_Label
```

The remaining 56 columns are ECG-derived features.

A full data dictionary is provided in:

```text
docs/data_dictionary/20PatientsFinalData_data_dictionary.md
docs/data_dictionary/20PatientsFinalData_data_dictionary.csv
```

## ECG feature documentation

The ECG feature list is documented in:

```text
docs/feature_documentation/ecg_lead_ii_feature_list.md
```

The 56 ECG-derived features include feature families such as:

- time-domain statistical features,
- entropy and complexity features,
- fractal and nonlinear dynamic features,
- entropy-profile aggregate features,
- fractal-dimension summary features,
- spectral features,
- HRV features,
- ECG morphology and signal-energy features.

## Clinical metadata

Curated patient-level clinical metadata are stored in:

```text
metadata/
```

This folder contains:

| File | Description |
|---|---|
| `patient_clinical_summary.csv` | Curated patient-level clinical summary table |
| `patient_clinical_flags.csv` | Machine-readable patient-level binary clinical flags |
| `ECG_Clinical_Summary_Report.docx` | Human-readable clinical reference report |
| `README.md` | Description of the clinical metadata files |

The file `patient_clinical_flags.csv` contains patient-level flags used for interpretation and sensitivity analyses, including:

```text
paced_or_icd
continuous_AF
intermittent_or_new_AF
any_AF
atypical_HI
major_HRV_outlier
BBB
antiarrhythmic_exposure
mechanical_support
clean_sinus_candidate
exclude_from_HRV_valid_subgroup
```

These annotations support analyses such as:

- HRV-valid subgroup analysis,
- clean sinus subgroup analysis,
- exclusion of atypical HI patients,
- exclusion of major HRV outlier patient,
- BBB-specific interpretation of QRS morphology,
- mechanical-support sensitivity analysis,
- antiarrhythmic exposure sensitivity analysis.

### Note on clinical data

Raw MIMIC-III discharge summaries and raw clinical notes are **not included** in this repository.

Only curated, derived clinical metadata are provided to support transparency and reproducibility of the analysis. Users requiring access to the original MIMIC-III clinical records must obtain access through the official MIMIC data access process and comply with the relevant data-use agreement.

## Analysis outputs

Generated outputs are stored in:

```text
analysis_outputs/
├── tables/
└── figures/
```

The `tables/` folder contains intermediate and final CSV files generated during the analysis. These include patient-level summaries, primary paired results, HRV-valid subgroup results, mixed-effects model results, redundancy summaries, sensitivity analyses, candidate scorecards, and manuscript-facing interpretation tables.

The `figures/` folder contains analysis and manuscript figures generated by the notebook.

## Statistical analysis summary

The primary analysis is a patient-level paired comparison.

For each ECG-derived feature:

1. the median value is computed for each patient during Pre-HI segments,
2. the median value is computed for each patient during During-HI segments,
3. the paired difference is computed as:

```text
During-HI median − Pre-HI median
```

The primary patient-level analysis includes:

- Wilcoxon signed-rank test,
- matched-pairs rank-biserial effect size,
- bootstrap confidence intervals for the median paired difference,
- direction consistency across patients,
- global and feature-family-level false discovery rate correction.

Supportive analyses include:

- HRV-valid subgroup analysis,
- segment-level linear mixed-effects models,
- within-subject centered Spearman correlation analysis,
- redundancy clustering,
- leave-one-patient-out analysis,
- signal-quality filtered analysis,
- clinically motivated sensitivity analyses,
- descriptive Episode_Label analysis.

The mixed-effects model used as a supportive analysis is:

```text
feature_z ~ Class + (1 | Patient_id)
```

The mixed-effects analysis is supportive rather than primary because multiple ECG segments are available per patient, while the cohort contains 20 patients.

## Final candidate features

The final scorecard selected six ECG-derived candidate features.

### Strict candidates

| Feature | Direction during HI | Interpretation |
|---|---|---|
| `kurtosis` | Increased | Leading patient-level distributional candidate |
| `QRS_mean` | Increased | Leading morphology/conduction candidate |

### Secondary candidates

| Feature | Direction during HI | Interpretation |
|---|---|---|
| `fuzzyEntropy` | Decreased | Representative entropy/complexity candidate |
| `fd_mean` | Decreased | Fractal-dimension summary candidate |
| `entropyProfiled_maximum_sampleEntropy` | Decreased | Derived entropy-profile candidate |
| `fd_median` | Decreased | Secondary fractal-dimension summary candidate |

These features should be interpreted as **exploratory candidate markers**, not validated biomarkers or clinical predictors.

## Important interpretation cautions

This repository supports an exploratory analysis. The results should not be interpreted as establishing:

- clinical diagnostic validity,
- prospective predictive performance,
- causal physiological mechanisms,
- deployment readiness,
- generalizability beyond this cohort.

Important limitations include:

- small cohort size,
- retrospective design,
- clinical heterogeneity,
- rhythm and pacing confounding,
- antiarrhythmic medication effects,
- possible signal-quality artefacts,
- no external validation cohort,
- no prospective prediction model.

## Reproducing the analysis

### 1. Clone the repository

```bash
git clone <repository-url>
cd <repository-name>
```

### 2. Create a Python environment

Using `venv`:

```bash
python -m venv .venv
source .venv/bin/activate
```

On Windows:

```bash
python -m venv .venv
.venv\Scripts\activate
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

The main dependencies are:

```text
numpy
pandas
scipy
statsmodels
scikit-learn
matplotlib
seaborn
openpyxl
jupyter
ipykernel
```

### 4. Run the notebook

```bash
jupyter notebook ecg_hi_patient_level_exploratory_analysis.ipynb
```

Run the notebook cells sequentially from top to bottom.

The notebook expects the input dataset at:

```text
data/processed/20PatientsFinalData.xlsx
```

and writes generated outputs to:

```text
analysis_outputs/tables/
analysis_outputs/figures/
```

## Citation

If you use this repository, please cite the associated manuscript:

```text
ECG-Derived Candidate Features Associated With Hemodynamic Instability in ICU Patients:
A Patient-Level Exploratory Analysis of Ultra-Short Segments
```

A complete citation will be added after publication.

## License

See the `LICENSE` file for reuse terms.

## Contact

For questions about the analysis or repository, please contact the corresponding author of the associated manuscript.
