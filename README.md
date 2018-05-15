# Evaluation procedures for forecasting with spatio-temporal data

This repository contains the research compendium of "Evaluation procedures for forecasting with spatio-temporal data", authored by Mariana Oliveira, Luis Torgo, and Vitor Santos Costa, and submitted to ECML2018.


## Prerequisites

To install necessary packages, run:

```
for(p in c("dplyr", "performanceEstimation", "starma",
"spdep", "reshape2", "assertthat", "wavethresh", "MASS",
"stringr", "doParallel", "foreach", "ranger")){
	if(!(p in installed.packages()))
		install.packages(p)
	}
```

To replicate figures, installing the following package is also necessary:

```
if(!("ggplot2" in installed.packages())) install.packages("ggplot2")
```

## Reproducing experiments

To run experiments, run the following lines from the main directory:

```
source("analysis/step1_gen_data.R")
source("analysis/step2_artificial_experiments.R")
source("analysis/step3_real_experiments.R")
```

## Contents

The repository is organized as follows:
* **inst/analysis** - Contains scripts for generating data and launching experiments
* **data** - Contains the real-world data sets used in the experiments in an appropriate format
* **man** - Contains function documentation
* **R** - Contains reusable functions
* **inst/results** - Contains files of the results presented in the paper

#### R

* analyse\_utils.R - functions that can be used to summarize produced by experiments
* calc_metrics.R - functions to calculate error metrics
* experiment_utils.R - functions designed to launch experiments
* fold_alloc.R - functions for allocating observations to cross-validation folds
* gen\_data\_utils.R - functions needed to generate lists of STARMA-generated datasets, depends on starma_utils.R
* methods.R - functions implementing different error estimation methods, depends on fold_alloc.R, utils.R
* starma_utils.R - functions needed to generate STARMA artificial data and check coefficient stationarity
* st_indicators.R - functions designed to calculate spatio-temporal indicators
* utils.R - general utility functions used inside some evaluation methods
* workflow.R - functions needed for the multiple steps in experiments, including learning, generating and evaluating predictions


#### inst/analysis

* step1_gen_data.R - data generation script
* step2_artificial_experiments.R - running experiments on artificially generated datasets
* step3_real_experiments.R - calculating spatio-temporal indicators and running experiments on real data sets


## Session Info

```
R version 3.4.0 (2017-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 14.04.5 LTS

Matrix products: default
BLAS: /usr/lib/libblas/libblas.so.3.0
LAPACK: /usr/lib/lapack/liblapack.so.3.0

locale:
 [1] LC_CTYPE=pt_PT.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=pt_PT.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=pt_PT.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=pt_PT.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=pt_PT.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] wavethresh_4.6.8            MASS_7.3-48                
 [3] stringr_1.2.0               assertthat_0.2.0           
 [5] doParallel_1.0.10           iterators_1.0.8            
 [7] foreach_1.4.3               rpart_4.1-11               
 [9] ranger_0.9.0                performanceEstimation_1.1.0
[11] dplyr_0.7.4                

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.11     codetools_0.2-15 lattice_0.20-35  grid_3.4.0      
 [5] R6_2.1.2         magrittr_1.5     stringi_1.0-1    rlang_0.1.6     
 [9] bindrcpp_0.2     Matrix_1.2-10    tools_3.4.0      glue_1.1.1      
[13] compiler_3.4.0   pkgconfig_2.0.1  bindr_0.1        tibble_1.3.3
```
