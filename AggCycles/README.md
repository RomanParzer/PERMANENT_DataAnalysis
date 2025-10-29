# PERMANENT_DataAnalysis
Data Analysis for PERMANENT PEM fuel cell degradation 


This repository aims at sharing the .R-files for reproducable research and collaboration,
specifically for exploration of the average current densities over cycles of specific voltage profiles.

It includes all R-scripts to reproduce all results collected in the file Cycle_protocol.tex on overleaf.
This folder includes the files 'RFNR_model_functions_v2.R' implementing the fitting and utility 
functions to fit the proposed model (with option of using the standard covariance matrix instead of the HAC; see folder ReversibleDegradation), 'Protocol_aggregatedCycles_clean.R' implementing the application
to all data sets using parallel computing and saving the results in a compressed R-file in the folder 'saved_results',
and 'Results_agg_cycles.R' reading in the saved R-files and producing all figures and Tables used in the paper.
