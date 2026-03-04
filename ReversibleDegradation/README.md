# PERMANENT_DataAnalysis
Data Analysis for PERMANENT PEM fuel cell degradation.

This repository aims at sharing the R files for reproducible research and collaboration.

The folder `ReversibleDegradation` includes all R scripts to reproduce the results
obtained in the paper *"A Stochastic Recursive Nonlinear Regression Model for
Recoverable Degradation in PEM Fuel Cells"* by Parzer et al.

The corresponding dataset can be downloaded from
[https://doi.org/10.5281/zenodo.17223247](https://doi.org/10.5281/zenodo.17223247).
Large data files in the folder `data/` need to be stored locally and are ignored
by git due to their size.

---

## File Structure

### `RFNR_help_functions.R`
The mathematical foundation of the RFNR model. Contains the linear predictor,
its analytical gradient (used for inference), the MSE objective function with
its gradient, and evaluation metrics (RMSE, MAE). This file is
sourced by all other scripts.

### `RFNR_model_functions.R`
Implements the RFNR model as an R S3 class. The main function `RFNR()` optimises
the model parameters via `opm()`, which tries multiple optimisation algorithms
and selects the best solution. It then computes fitted values, residuals, and
estimates a HAC covariance matrix for inference. Standard S3 methods (`predict`,
`plot`, `summary`, `coef`, `residuals`) allow the model to be used with familiar
R syntax.

### `benchmark_model_functions.R`
Four benchmark models sharing the same S3 interface as RFNR:
- **WienerDeg**: Wiener process with linear drift on log(y)
- **GammaDeg**: Gamma process on the degradation increments
- **EKF_JD**: Extended Kalman Filter with deterministic jump inputs
- **SparseGP**: Sparse Gaussian Process via the `laGP` package

### `auto_tune_final.R`
Automatic selection of the jump detection threshold τ via the empirical quantile
of the absolute symmetric gradient of log(current):

$$\tau = Q_{p}\|\left(\left|\nabla \log y\right|\right), \quad p = 0.987$$

The two RH100 datasets (i=10, i=23) use
manually set thresholds due to their distinct sampling frequency and signal
structure.

### `main_analysis.R`
The main script orchestrating the full analysis in three phases:
- **Phase 1**: fits all four RFNR variants in parallel across all 26 datasets,
  with automatic resume (skips already completed datasets)
- **Phase 2a**: fits WienerDeg, GammaDeg, and EKF_JD in parallel
- **Phase 2b**: fits SparseGP sequentially (memory-intensive, not safe in PSOCK clusters)
- **Phase 3**: aggregates all results and produces the RMSE comparison plot


---

## Key Design Choices

All models share the same data preparation pipeline (identical train/test split,
subsampling index `ndx`, and jump indicator vector `dum`), ensuring that
performance comparisons are fair and directly interpretable. Results are saved
incrementally after each dataset, so the analysis can be interrupted and resumed
without losing progress.

The script saves a full `sessionInfo()` report at the end of the analysis,
recording the exact R version, operating system, and package versions used,
ensuring full computational reproducibility.