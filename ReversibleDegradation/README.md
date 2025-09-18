# PERMANENT_DataAnalysis
Data Analysis for PERMANENT PEM fuel cell degradation

This repository aims at sharing the .R-files for reproducable research and collaboration.
The folder 'ReversibleDegradation' includes all R-scripts to reproduce all results obtained in the paper 
'A Recursive Nonlinear Regression Model for Recoverable Degradation in PEM Fuel Cells' submitted at JRSSC.
This folder includes the files 'RFNR_model_functions.R' implementing the fitting and utility 
functions to fit the proposed model, 'Protocol_RFNR_reversible_all.R' implementing the application
to all data sets using parallel computing and saving the results in a compressed R-file in the folder 'saved_results',
and 'Results_RFNR_model_all.R' reading in the saved R-files and producing all figures and Tables used in the paper.

Large data files in an folder data/ (in parent repository, not in ReversibleDegradation/) need to be stored locally and are ignored by git due to their large size.

For reproducability, here is the sessionInfo() (used packages and their versions) when executing the file 'Protocol_RFNR_reversible_all.R' to obtain the results presented in the paper:

R version 4.2.1 (2022-06-23)

Platform: aarch64-apple-darwin20 (64-bit)

Running under: macOS 14.3.1

Matrix products: default

LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

locale:

[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#### attached base packages:

[1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

#### other attached packages:

 [1] foreach_1.5.2     readxl_1.4.5      lubridate_1.9.3   forcats_1.0.0     stringr_1.5.1     purrr_1.0.2      
 [7] readr_2.1.4       tibble_3.2.1      tidyverse_2.0.0   ggplotify_0.1.2   optimx_2023-10.21 kableExtra_1.3.4 
[13] ggrepel_0.9.4     knitr_1.44        ggplot2_3.5.2     tidyr_1.3.1       dplyr_1.1.4      

#### loaded via a namespace (and not attached):

 [1] httr_1.4.7          viridisLite_0.4.2   ucminf_1.2.0        highr_0.9           yulab.utils_0.2.1  
 [6] cellranger_1.1.0    numDeriv_2016.8-1.1 pillar_1.10.2       lattice_0.20-45     glue_1.8.0         
[11] quadprog_1.5-8      digest_0.6.37       rvest_1.0.4         minqa_1.2.4         colorspace_2.1-1   
[16] sandwich_3.1-1      Matrix_1.6-1.1      dfoptim_2023.1.0    htmltools_0.5.7     pkgconfig_2.0.3    
[21] scales_1.3.0        webshot_0.5.3       svglite_2.1.0       marqLevAlg_2.0.8    tzdb_0.4.0         
[26] pracma_2.3.8        lbfgsb3c_2020-3.3   timechange_0.3.0    generics_0.1.3      farver_2.1.2       
[31] ellipsis_0.3.2      pacman_0.5.1        withr_3.0.2         cli_3.6.5           magrittr_2.0.3     
[36] evaluate_1.0.0      fs_1.6.4            MASS_7.3-60         doParallel_1.0.17   xml2_1.3.3         
[41] textshaping_0.3.6   tools_4.2.1         hms_1.1.2           lifecycle_1.0.4     munsell_0.5.1      
[46] compiler_4.2.1      gridGraphics_0.5-1  systemfonts_1.0.4   rlang_1.1.6         nloptr_2.0.3       
[51] iterators_1.0.14    rstudioapi_0.17.1   rappdirs_0.3.3      subplex_1.9         labeling_0.4.3     
[56] rmarkdown_2.25      codetools_0.2-18    gtable_0.3.6        BB_2019.10-1        R6_2.6.1           
[61] zoo_1.8-12          fastmap_1.1.1       ragg_1.2.5          stringi_1.8.4       Rcpp_1.0.14        
[66] vctrs_0.6.5         tidyselect_1.2.1    xfun_0.42          
