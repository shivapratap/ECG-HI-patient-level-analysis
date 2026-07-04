# ECG-HI Patient-Level Analysis

**ECG-Derived Candidate Features Associated With Hemodynamic Instability in ICU Patients: A Patient-Level Exploratory Analysis of Ultra-Short Segments**

This repository contains the analysis pipeline, processed ECG feature tables, clinical-context annotations, signal-quality flags, and documentation for an exploratory study of ECG-derived candidate features associated with hemodynamic instability (HI) in ICU patients.

The analysis uses a clinically curated, cardiogenic-shock–enriched MIMIC-III cohort of 20 ICU patients. The central question is whether continuously available **Lead-II ECG** contains morphology, rhythm, distributional, entropy, or nonlinear-complexity features that change within patients between **Pre-HI** and **During-HI** conditions.

This repository is intended for reproducible research and transparent review of the feature-discovery workflow. It is **not** a clinical prediction model, a real-time HI detector, or a validated ECG biomarker package.

---

## Study overview

Hemodynamic instability can be difficult to recognize early in ICU practice because clinically meaningful deterioration may occur between intermittent blood pressure measurements and may be partly masked by compensatory mechanisms. ECG is continuously monitored in many ICU patients and is therefore attractive for studying rapid physiological change. However, ECG does **not** directly measure blood pressure, cardiac output, vascular tone, or tissue perfusion.

For this reason, this project uses a cautious patient-level feature-discovery design. Segment-level ECG features are first summarized into paired patient-level medians, and statistical inference is performed at the **patient level**, not the segment level.

### Cohort and analysis summary

| Item | Value |
|---|---:|
| Data source | MIMIC-III Clinical Database and matched Waveform Database |
| Cohort | Clinically curated cardiogenic-shock–enriched ICU cohort |
| Patients | 20 |
| ECG lead | Lead II |
| ECG sampling rate | 125 Hz |
| Segment length | 2 minutes |
| Total ECG segments | 1,260 |
| Pre-HI segments | 593 |
| During-HI segments | 667 |
| Initial ECG-derived features | 56 |
| Informative ECG-derived features analyzed | 52 |
| Primary statistical unit | Patient |

---

## What is included in this repository

This repository supports the patient-level exploratory analysis reported in the associated manuscript. It includes:

- processed ECG-derived feature data for the curated 20-patient cohort;
- patient-level clinical-context and rhythm flags;
- signal-quality and artifact-review flags;
- the full 17-stage patient-level analysis pipeline;
- intermediate and final analysis outputs;
- documentation for understanding and rerunning the workflow.

The pipeline starts from a processed segment-level ECG feature workbook. It does not download MIMIC-III, extract raw waveforms, or regenerate ECG features from raw waveform files.

---

## Repository structure

The expected project layout is:

```text
ECG-HI-patient-level-analysis/
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
│   ├── feature_documentation/
│   │   └── ecg_lead_ii_feature_list.md
│   └── Pipeline_Readme.md
├── pipeline/
│   ├── config.py
│   ├── utils.py
│   ├── stages.py
│   └── run_pipeline.py
├── outputs/
├── requirements.txt
├── LICENSE
└── README.md
```

The `outputs/` directory is created automatically when the pipeline is run.

---

## Required input files

The pipeline expects the following required files:

| File | Purpose |
|---|---|
| `data/processed/20PatientsFinalData.xlsx` | Main segment-level ECG feature dataset |
| `metadata/patient_clinical_flags.csv` | Patient-level rhythm, clinical-context, and confound flags |
| `metadata/patient_clinical_summary.csv` | Patient-level clinical summary reference |
| `docs/data_dictionary/20PatientsFinalData_data_dictionary.csv` | Data dictionary for the processed feature dataset |
| `docs/feature_documentation/ecg_lead_ii_feature_list.md` | ECG feature documentation |

Optional/reference files:

| File | Purpose |
|---|---|
| `metadata/ECG_Clinical_Summary_Report.docx` | Human-readable clinical summary report |
| `docs/data_dictionary/20PatientsFinalData_data_dictionary.md` | Markdown version of the data dictionary |

---

## Installation

Clone the repository:

```bash
git clone https://github.com/shivapratap/ECG-HI-patient-level-analysis.git
cd ECG-HI-patient-level-analysis
```

Create and activate a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
```

On Windows, use:

```bash
.venv\Scripts\activate
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Python 3.10 or later is recommended.

---

## How to run the analysis

Run the full pipeline from the repository root:

```bash
python pipeline/run_pipeline.py
```

The runner executes all pipeline stages in order, from data audit through final workbook export. When the run completes successfully, outputs are written under:

```text
outputs/
```

The final workbook is written to:

```text
outputs/final_workbook/HI_ECG_candidate_feature_analysis_workbook.xlsx
```

---

## Pipeline overview

The pipeline contains 17 stages:

| Stage | Purpose |
|---:|---|
| 1 | Load and audit the ECG feature dataset |
| 2 | Define the full ECG-derived feature set |
| 3 | Remove non-informative features |
| 4 | Merge clinical-context metadata |
| 5 | Create segment-level signal-quality flags |
| 6 | Aggregate segment-level features into patient-condition medians |
| 7 | Run primary patient-level paired statistics |
| 8 | Apply global and feature-family FDR correction |
| 9 | Map redundancy using within-patient centered Spearman correlation |
| 10 | Create a preliminary candidate-feature shortlist |
| 11 | Run HRV-valid subgroup analysis |
| 12 | Run supportive mixed-effects models |
| 13 | Run sensitivity analyses |
| 14 | Run leave-one-patient-out influence analysis |
| 15 | Analyze patient-level agreement and disagreement |
| 16 | Build the final exploratory candidate scorecard |
| 17 | Export the final review workbook |

For a more detailed explanation of each stage, see:

```text
docs/Pipeline_Readme.md
```

---

## Main analysis outputs

Important outputs include:

| Output | Description |
|---|---|
| `outputs/data_audit/data_integrity_summary.csv` | Dataset structure and integrity summary |
| `outputs/data_audit/segment_counts_by_patient_condition.csv` | Segment counts by patient and condition |
| `outputs/feature_definition/informative_feature_list.csv` | Features retained for analysis |
| `outputs/clinical_metadata/cohort_description_table.csv` | Clinical/rhythm-context cohort summary |
| `outputs/patient_level/patient_feature_delta_table_all_features.csv` | Patient-level Pre-HI vs During-HI feature deltas |
| `outputs/primary_statistics/primary_paired_feature_results_all_features_fdr.csv` | Primary paired statistics with FDR correction |
| `outputs/hrv_subgroup/hrv_valid_feature_results.csv` | HRV-valid subgroup results |
| `outputs/mixed_effects/lmm_supportive_results_all_features.csv` | Supportive mixed-effects model results |
| `outputs/sensitivity/candidate_feature_robustness_matrix.csv` | Candidate robustness across sensitivity scenarios |
| `outputs/agreement_disagreement/feature_agreement_summary.csv` | Feature-level agreement summary |
| `outputs/agreement_disagreement/patient_disagreement_summary.csv` | Patient-level disagreement summary |
| `outputs/candidate_scorecard/final_candidate_scorecard.csv` | Full candidate evidence scorecard |
| `outputs/candidate_scorecard/final_candidate_feature_categories.csv` | Final candidate category summary |
| `outputs/final_workbook/HI_ECG_candidate_feature_analysis_workbook.xlsx` | Consolidated workbook of key results |

---

## Main result summary

The primary patient-level paired analysis did **not** identify any ECG-derived feature that survived global or feature-family false discovery rate correction. Therefore, all candidate features should be interpreted as exploratory and hypothesis-generating.

Within this FDR-null constraint, the final scorecard prioritized a nine-feature exploratory structure:

| Candidate tier | Features |
|---|---|
| Primary exploratory candidates | `kurtosis`, `QRS_mean`, `entropyProfiled_maximum_sampleEntropy` |
| Secondary exploratory candidates | `fuzzyEntropy`, `fd_median`, `entropyProfiled_standardDeviation_sampleEntropy` |
| Secondary rhythm-sensitive candidates | `HRV_pNN50`, `HRV_MeanNN` |
| Exploratory/downgraded | `fd_minimum` |

The HRV-derived features are treated separately because HRV interpretation depends on reliable beat-to-beat interval estimation and is sensitive to pacing, atrial fibrillation, R-peak detection errors, and rhythm irregularity.

---

## Interpretation cautions

This project is an exploratory feature-discovery analysis. The results should be described as:

> ECG-derived candidate features associated with Pre-HI to During-HI change in this small, clinically curated ICU cohort.

They should **not** be described as:

- validated biomarkers;
- a deployable HI detector;
- a blood-pressure estimation method;
- a clinical decision-support tool;
- externally validated predictive features.

The final scorecard organizes exploratory evidence. It does not override the null FDR-corrected primary analysis.

---

## Data and licensing

The processed feature files and analysis code in this repository are available under the MIT License.

This project was derived from MIMIC-III clinical and waveform resources. Users who wish to access or reproduce the original raw clinical or waveform data must obtain access through the appropriate MIMIC/PhysioNet credentialing and data-use process. Raw MIMIC-III records are not redistributed here.

Users are responsible for complying with all applicable data-use agreements and must not attempt to re-identify patients.

---

## Citation

A citation will be added here when available.

```bibtex
# Citation placeholder
```

---

## Contact

Repository maintainer:

- GitHub: [@shivapratap](https://github.com/shivapratap)
- Email: shivapg@am.amrita.edu

---

## License

This repository is released under the MIT License. See `LICENSE` for details.
