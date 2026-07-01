# ECG Lead II Feature List

This document describes the ECG-derived features extracted from Lead II ECG signals for the hemodynamic instability analysis.

A total of 56 ECG-derived features were extracted from ultra-short 2-minute ECG segments. These features were grouped into time-domain, entropy/complexity, fractal/nonlinear, spectral, HRV, morphology, and signal-energy feature families.

## 1. Time-domain statistical features

| Feature | Description |
|---|---|
| maximum | Maximum signal amplitude within the ECG segment |
| minimum | Minimum signal amplitude within the ECG segment |
| mean | Mean signal amplitude |
| median | Median signal amplitude |
| standardDeviation | Standard deviation of signal amplitude |
| variance | Variance of signal amplitude |
| kurtosis | Kurtosis of the signal amplitude distribution |
| skewness | Skewness of the signal amplitude distribution |
| numberOfZeroCrossing | Number of zero crossings in the ECG segment |
| positiveToNegativeSampleRatio | Ratio of positive to negative signal samples |
| positiveToNegativePeakRatio | Ratio of positive to negative signal peaks |
| meanAbsoluteValue | Mean absolute signal amplitude |

## 2. Entropy and complexity features

| Feature | Description |
|---|---|
| approximateEntropy | Approximate entropy of the ECG segment |
| sampleEntropy | Sample entropy of the ECG segment |
| permutationEntropy | Permutation entropy of the ECG segment |
| singularValueDecompositionEntropy | Singular value decomposition entropy |
| fuzzyEntropy | Fuzzy entropy of the ECG segment |
| distributionEntropy | Distribution entropy of the ECG segment |
| shannonEntropy | Shannon entropy |
| renyiEntropy | Rényi entropy |
| lempelZivComplexity | Lempel-Ziv complexity |

## 3. Fractal and nonlinear dynamic features

| Feature | Description |
|---|---|
| hjorthMobility | Hjorth mobility parameter |
| hjorthComplexity | Hjorth complexity parameter |
| fisherInfo | Fisher information measure |
| petrosianFd | Petrosian fractal dimension |
| katzFd | Katz fractal dimension |
| higuchiFd | Higuchi fractal dimension |
| detrendedFluctuation | Detrended fluctuation analysis measure |

## 4. Entropy-profile aggregate features

| Feature | Description |
|---|---|
| entropyProfiled_total_sampleEntropy | Total sample entropy across entropy-profile windows |
| entropyProfiled_average_sampleEntropy | Average sample entropy across entropy-profile windows |
| entropyProfiled_maximum_sampleEntropy | Maximum sample entropy across entropy-profile windows |
| entropyProfiled_minimum_sampleEntropy | Minimum sample entropy across entropy-profile windows |
| entropyProfiled_median_sampleEntropy | Median sample entropy across entropy-profile windows |
| entropyProfiled_standardDeviation_sampleEntropy | Standard deviation of sample entropy across entropy-profile windows |
| entropyProfiled_variance_sampleEntropy | Variance of sample entropy across entropy-profile windows |
| entropyProfiled_kurtosis_sampleEntropy | Kurtosis of sample entropy across entropy-profile windows |
| entropyProfiled_skewness_sampleEntropy | Skewness of sample entropy across entropy-profile windows |

## 5. Fractal-dimension summary features

| Feature | Description |
|---|---|
| fd_maximum | Maximum value of fractal-dimension profile |
| fd_minimum | Minimum value of fractal-dimension profile |
| fd_mean | Mean value of fractal-dimension profile |
| fd_median | Median value of fractal-dimension profile |
| fd_standardDeviation | Standard deviation of fractal-dimension profile |
| fd_variance | Variance of fractal-dimension profile |
| fd_kurtosis | Kurtosis of fractal-dimension profile |
| fd_skewness | Skewness of fractal-dimension profile |

## 6. Spectral features

| Feature | Description |
|---|---|
| spectralEntropy | Spectral entropy of the ECG segment |
| fd_bandPower | Band-power measure associated with the fractal-dimension/spectral representation |

## 7. HRV features

| Feature | Description |
|---|---|
| HRV_MeanNN | Mean normal-to-normal interval |
| HRV_SDNN | Standard deviation of normal-to-normal intervals |
| HRV_RMSSD | Root mean square of successive NN interval differences |
| HRV_pNN50 | Percentage of successive NN intervals differing by more than 50 ms |

## 8. ECG morphology and signal-energy features

| Feature | Description |
|---|---|
| R_amp_mean | Mean R-wave amplitude |
| R_amp_std | Standard deviation of R-wave amplitude |
| QRS_mean | Mean QRS duration |
| QRS_std | Standard deviation of QRS duration |
| signal_energy | Signal energy of the ECG segment |

## Notes

The feature extraction was performed on Lead II ECG signals using ultra-short 2-minute segments.

Some features may be sensitive to rhythm abnormalities, pacing, ectopy, filtering, signal quality, and waveform artefacts. Therefore, downstream analyses include clinical flagging, signal-quality checks, HRV-valid subgroup analysis, redundancy filtering, and sensitivity analyses.