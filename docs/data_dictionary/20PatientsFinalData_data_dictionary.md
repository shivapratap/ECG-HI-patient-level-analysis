# Data Dictionary: `20PatientsFinalData.xlsx`

This document describes the main processed dataset used by the notebook:

`data/processed/20PatientsFinalData.xlsx`

Each row corresponds to one ultra-short 2-minute ECG Lead II segment from one ICU patient. Each segment is labelled as either Pre-HI or During-HI.

## Dataset overview

| Item | Value |
|---|---:|
| Rows | 1,260 |
| Columns | 62 |
| Metadata columns | 6 |
| ECG-derived feature columns | 56 |
| Patients | 20 |
| ECG lead | Lead II |
| Segment length | 2 minutes |

## Column groups

| Feature family | Number of columns |
| --- | --- |
| ECG morphology and signal-energy features | 5 |
| Entropy and complexity features | 9 |
| Entropy-profile aggregate features | 9 |
| Fractal and nonlinear dynamic features | 7 |
| Fractal-dimension summary features | 8 |
| HRV features | 4 |
| Metadata | 6 |
| Spectral features | 2 |
| Time-domain statistical features | 12 |

## Metadata columns

| column_name | column_type | description | allowed_values_or_units | role_in_analysis | notes |
| --- | --- | --- | --- | --- | --- |
| Sl_No | Integer | Row serial number in the processed dataset. | Sequential integer | Traceability | Not used as a statistical variable. |
| Segment_Number | Integer | Segment number assigned to the 2-minute ECG segment within the processed dataset context. | Positive integer | Traceability | Used for row/segment traceability; not used as the primary statistical unit. |
| Patient_id | Identifier | De-identified patient identifier. | De-identified patient ID | Grouping variable | Primary statistical unit for paired patient-level analysis. |
| Class | Categorical integer | Binary condition label for each ECG segment. | 0 = Pre-HI; 1 = During-HI | Primary condition variable | Used to define paired Pre-HI versus During-HI comparisons. |
| Condition | Categorical text | Human-readable segment condition. | PreHI; HI | Condition label | Should correspond to Class values. |
| Episode_Label | Categorical text | Detailed episode label for Pre-HI or HI episode grouping. | PreHI; HI; HI1; HI2; HI3; HI4; HI5; HI6 | Descriptive episode grouping | Used descriptively only; not used as a primary model term. |

## ECG-derived feature columns by family

### HRV features

| column_name | column_type | description | allowed_values_or_units | role_in_analysis | notes |
| --- | --- | --- | --- | --- | --- |
| HRV_MeanNN | Numeric | Mean normal-to-normal interval estimated from detected beats within the segment. | milliseconds (ms) | ECG-derived feature | HRV feature; not prioritized after HRV-valid subgroup analysis. |
| HRV_SDNN | Numeric | Standard deviation of normal-to-normal intervals within the segment. | milliseconds (ms) | ECG-derived feature | HRV feature; not prioritized after HRV-valid subgroup analysis. |
| HRV_RMSSD | Numeric | Root mean square of successive normal-to-normal interval differences within the segment. | milliseconds (ms) | ECG-derived feature | HRV feature; not prioritized after HRV-valid subgroup analysis. |
| HRV_pNN50 | Numeric | Percentage of successive normal-to-normal intervals differing by more than 50 ms. | percent (%) | ECG-derived feature | HRV feature; not prioritized after HRV-valid subgroup analysis. |

### ECG morphology and signal-energy features

| column_name | column_type | description | allowed_values_or_units | role_in_analysis | notes |
| --- | --- | --- | --- | --- | --- |
| R_amp_mean | Numeric | Mean R-wave amplitude within the ECG segment. | ECG amplitude units from processed signal | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| R_amp_std | Numeric | Standard deviation of R-wave amplitude within the ECG segment. | ECG amplitude units from processed signal | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| QRS_mean | Numeric | Mean QRS duration within the ECG segment. | seconds (s) | Candidate ECG feature | Final strict candidate; increased during HI; rhythm/conduction sensitive. |
| QRS_std | Numeric | Standard deviation of QRS duration within the ECG segment. | seconds (s) | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| signal_energy | Numeric | Energy of the ECG signal over the segment. | signal-energy units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |

### Time-domain statistical features

| column_name | column_type | description | allowed_values_or_units | role_in_analysis | notes |
| --- | --- | --- | --- | --- | --- |
| maximum | Numeric | Maximum ECG signal amplitude within the segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| minimum | Numeric | Minimum ECG signal amplitude within the segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| mean | Numeric | Mean ECG signal amplitude within the segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. This feature was later identified as non-informative/constant in this cohort and excluded from informative-feature analyses. |
| median | Numeric | Median ECG signal amplitude within the segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| standardDeviation | Numeric | Standard deviation of ECG signal amplitude within the segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| variance | Numeric | Variance of ECG signal amplitude within the segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| kurtosis | Numeric | Kurtosis of the ECG signal amplitude distribution within the segment. | Numeric value; unitless or feature-specific processed-signal units | Candidate ECG feature | Final strict candidate; increased during HI in the integrated scorecard. |
| skewness | Numeric | Skewness of the ECG signal amplitude distribution within the segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| numberOfZeroCrossing | Numeric | Number of times the ECG signal crosses zero within the segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| positiveToNegativeSampleRatio | Numeric | Ratio of positive to negative ECG signal samples within the segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| positiveToNegativePeakRatio | Numeric | Ratio of positive to negative ECG signal peaks within the segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. This feature was later identified as non-informative/constant in this cohort and excluded from informative-feature analyses. |
| meanAbsoluteValue | Numeric | Mean absolute ECG signal amplitude within the segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |

### Entropy and complexity features

| column_name | column_type | description | allowed_values_or_units | role_in_analysis | notes |
| --- | --- | --- | --- | --- | --- |
| approximateEntropy | Numeric | Approximate entropy of the ECG segment, representing signal regularity/complexity. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| sampleEntropy | Numeric | Sample entropy of the ECG segment, representing signal regularity/complexity. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| permutationEntropy | Numeric | Permutation entropy of the ECG segment, representing ordinal complexity. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| singularValueDecompositionEntropy | Numeric | Entropy derived from singular value decomposition of the ECG segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| fuzzyEntropy | Numeric | Fuzzy entropy of the ECG segment, representing signal irregularity/complexity. | Numeric value; unitless or feature-specific processed-signal units | Candidate ECG feature | Final secondary candidate; representative entropy/complexity feature; decreased during HI. |
| distributionEntropy | Numeric | Distribution entropy of the ECG segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| shannonEntropy | Numeric | Shannon entropy of the ECG segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. This feature was later identified as non-informative/constant in this cohort and excluded from informative-feature analyses. |
| renyiEntropy | Numeric | Renyi entropy of the ECG segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| lempelZivComplexity | Numeric | Lempel-Ziv complexity of the ECG segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |

### Fractal and nonlinear dynamic features

| column_name | column_type | description | allowed_values_or_units | role_in_analysis | notes |
| --- | --- | --- | --- | --- | --- |
| hjorthMobility | Numeric | Hjorth mobility parameter describing signal frequency/mobility characteristics. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| hjorthComplexity | Numeric | Hjorth complexity parameter describing signal shape complexity. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| fisherInfo | Numeric | Fisher information measure derived from the ECG segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| petrosianFd | Numeric | Petrosian fractal dimension of the ECG segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| katzFd | Numeric | Katz fractal dimension of the ECG segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| higuchiFd | Numeric | Higuchi fractal dimension of the ECG segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| detrendedFluctuation | Numeric | Detrended fluctuation analysis measure for long-range correlation/fractal behavior. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |

### Entropy-profile aggregate features

| column_name | column_type | description | allowed_values_or_units | role_in_analysis | notes |
| --- | --- | --- | --- | --- | --- |
| entropyProfiled_total_sampleEntropy | Numeric | Total sample entropy across entropy-profile windows. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| entropyProfiled_average_sampleEntropy | Numeric | Average sample entropy across entropy-profile windows. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| entropyProfiled_maximum_sampleEntropy | Numeric | Maximum sample entropy across entropy-profile windows. | Numeric value; unitless or feature-specific processed-signal units | Candidate ECG feature | Final secondary candidate; derived entropy-profile feature; decreased during HI. |
| entropyProfiled_minimum_sampleEntropy | Numeric | Minimum sample entropy across entropy-profile windows. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. This feature was later identified as non-informative/constant in this cohort and excluded from informative-feature analyses. |
| entropyProfiled_median_sampleEntropy | Numeric | Median sample entropy across entropy-profile windows. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| entropyProfiled_standardDeviation_sampleEntropy | Numeric | Standard deviation of sample entropy across entropy-profile windows. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| entropyProfiled_variance_sampleEntropy | Numeric | Variance of sample entropy across entropy-profile windows. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| entropyProfiled_kurtosis_sampleEntropy | Numeric | Kurtosis of sample entropy values across entropy-profile windows. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| entropyProfiled_skewness_sampleEntropy | Numeric | Skewness of sample entropy values across entropy-profile windows. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |

### Frequency Domain summary features

| column_name | column_type | description | allowed_values_or_units | role_in_analysis | notes |
| --- | --- | --- | --- | --- | --- |
| fd_maximum | Numeric | Maximum value of the frequency domain profile. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| fd_minimum | Numeric | Minimum value of the frequency domain profile. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| fd_mean | Numeric | Mean value of the frequency domain profile. | Numeric value; unitless or feature-specific processed-signal units | Candidate ECG feature | Final secondary candidate; frequency domain summary feature; decreased during HI. |
| fd_median | Numeric | Median value of the frequency domain profile. | Numeric value; unitless or feature-specific processed-signal units | Candidate ECG feature | Final secondary candidate; frequency domain summary feature; decreased during HI; may be redundant with fd_mean. |
| fd_standardDeviation | Numeric | Standard deviation of the frequency domain profile. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| fd_variance | Numeric | Variance of the frequency domain profile. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| fd_kurtosis | Numeric | Kurtosis of the frequency domain profile. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| fd_skewness | Numeric | Skewness of the frequency domain profile. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |

### Spectral features

| column_name | column_type | description | allowed_values_or_units | role_in_analysis | notes |
| --- | --- | --- | --- | --- | --- |
| spectralEntropy | Numeric | Spectral entropy of the ECG segment. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |
| fd_bandPower | Numeric | Band-power feature associated with the spectral/fractal representation. | Numeric value; unitless or feature-specific processed-signal units | ECG-derived feature | Extracted from ultra-short 2-minute Lead II ECG segment. |

## Interpretation notes

- `Patient_id` is used as the primary grouping variable for patient-level paired analysis.
- `Class` is the primary binary condition variable, where `0 = Pre-HI` and `1 = During-HI`.
- `Condition` is a human-readable version of the binary class label.
- `Episode_Label` is retained for descriptive analysis only and is not used as a primary model term.
- HRV features require cautious interpretation because they may be affected by atrial fibrillation, pacing, ectopy, and rhythm irregularity.
- Morphology features such as `QRS_mean` may be affected by pacing, bundle branch block, antiarrhythmic exposure, ischemia, electrolyte abnormalities, and signal quality.
- The final manuscript candidate features were selected through patient-level paired analysis, supportive mixed-effects modelling, redundancy filtering, sensitivity analyses, and a transparent candidate scorecard.

## Final candidate feature notes

| Tier | Feature | Direction during HI | Note |
|---|---|---|---|
| Strict candidate | `kurtosis` | Increased | Leading patient-level distributional candidate |
| Strict candidate | `QRS_mean` | Increased | Leading morphology/conduction candidate |
| Secondary candidate | `fuzzyEntropy` | Decreased | Representative entropy/complexity candidate |
| Secondary candidate | `fd_mean` | Decreased | Fractal-dimension summary candidate |
| Secondary candidate | `entropyProfiled_maximum_sampleEntropy` | Decreased | Derived entropy-profile candidate |
| Secondary candidate | `fd_median` | Decreased | Secondary frequency domain summary candidate |

## Reproducibility

This Markdown data dictionary is intended for human-readable documentation. The companion machine-readable file is:

`docs/data_dictionary/20PatientsFinalData_data_dictionary.csv`
