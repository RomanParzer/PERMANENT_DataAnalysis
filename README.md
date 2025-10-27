# PERMANENT_DataAnalysis
Data Analysis for PERMANENT PEM fuel cell degradation

This repository aims at sharing the .R-files for reproducable research and collaboration.
The folder 'ReversibleDegradation' includes all R-scripts to reproduce all results obtained in the paper 
'A Recursive Nonlinear Regression Model for Recoverable Degradation in PEM Fuel Cells' submitted at JRSSC.
This folder includes the files 'RFNR_model_functions.R' implementing the fitting and utility 
functions to fit the proposed model, 'Protocol_RFNR_reversible_all.R' implementing the application
to all data sets using parallel computing and saving the results in a compressed R-file in the folder 'saved_results',
and 'Results_RFNR_model_all.R' reading in the saved R-files and producing all figures and Tables used in the paper.

Large data files in an folder ReversibleDegradation/data/ need to be stored locally and are ignored by git due to their large size.
For the paper 'A Recursive Nonlinear Regression Model for Recoverable Degradation in PEM Fuel Cells' by Parzer et. al., the corresponding data set can be downloaded from https://doi.org/10.5281/zenodo.17223247.
