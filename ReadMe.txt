
This folder includes data, R codes, and outputs for the manuscript “A multi-modal study of microbiomes and metabolomes reveals a system-wide dysbiosis preceding HIV-1 infection”. 

The main script to run all the analyses is "code/data_analysis_workflow.R". All the required R packages were installed in RStudio. 
Correlation analyses were performed using SECOM (from the package ANCOMBC, https://github.com/FrederickHuangLin/ANCOMBC) with modifications to regularize variances and to account for covariates. The scripts can be found at code/ANCOMBC_mod_covariates.zip. Please unzip the folder code/ANCOMBC_mod_covariates.zip before running code/data_analysis_workflow.R.
 

###################################################
            Software and Platform
###################################################


sessionInfo()

# Version information about R and attached packages:
RStudio 2023.09.0+463
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS Ventura 13.6.8

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods  
[9] base     

other attached packages:
 [1] mixOmics_6.28.0                 MASS_7.3-60.2                  
 [3] microViz_0.12.4                 kableExtra_1.4.0               
 [5] readxl_1.4.3                    RColorBrewer_1.1-3             
 [7] NetCoMi_1.1.0                   SpiecEasi_1.1.3                
 [9] rstatix_0.7.2                   vegan_2.6-6                    
[11] lattice_0.22-6                  permute_0.9-7                  
[13] lubridate_1.9.3                 forcats_1.0.0                  
[15] stringr_1.5.1                   dplyr_1.1.4                    
[17] purrr_1.0.2                     readr_2.1.5                    
[19] tidyr_1.3.1                     tidyverse_2.0.0                
[21] tibble_3.2.1                    viridis_0.6.5                  
[23] viridisLite_0.4.2               circlize_0.4.16                
[25] ComplexHeatmap_2.20.0           gridExtra_2.3                  
[27] ggpubr_0.6.0                    ggrepel_0.9.5                  
[29] ggplot2_3.5.1                   phyloseq_1.48.0                
[31] TreeSummarizedExperiment_2.12.0 Biostrings_2.72.0              
[33] XVector_0.44.0                  SingleCellExperiment_1.26.0    
[35] SummarizedExperiment_1.34.0     Biobase_2.64.0                 
[37] GenomicRanges_1.56.0            GenomeInfoDb_1.40.0            
[39] IRanges_2.38.0                  S4Vectors_0.42.0               
[41] BiocGenerics_0.50.0             MatrixGenerics_1.16.0          
[43] matrixStats_1.3.0              

loaded via a namespace (and not attached):
  [1] vroom_1.6.5                 gld_2.6.6                   urlchecker_1.0.1           
  [4] rARPACK_0.11-0              nnet_7.3-19                 TH.data_1.1-2              
  [7] vctrs_0.6.5                 energy_1.7-11               digest_0.6.35              
 [10] png_0.1-8                   corpcor_1.6.10              shape_1.4.6.1              
 [13] proxy_0.4-27                Exact_3.2                   pcaPP_2.0-4                
 [16] corrplot_0.92               SPRING_1.0.4                reshape2_1.4.4             
 [19] foreach_1.5.2               httpuv_1.6.15               withr_3.0.0                
 [22] psych_2.4.3                 xfun_0.44                   ellipsis_0.3.2             
 [25] survival_3.6-4              doRNG_1.8.6                 memoise_2.0.1              
 [28] ggbeeswarm_0.7.2            profvis_0.3.8               gmp_0.7-4                  
 [31] systemfonts_1.0.6           tidytree_0.4.6              zoo_1.8-12                 
 [34] GlobalOptions_0.1.2         gtools_3.9.5                pbapply_1.7-2              
 [37] Formula_1.2-5               ellipse_0.5.0               KEGGREST_1.44.0            
 [40] promises_1.3.0              httr_1.4.7                  rhdf5filters_1.16.0        
 [43] rhdf5_2.48.0                rstudioapi_0.16.0           UCSC.utils_1.0.0           
 [46] miniUI_0.1.1.1              generics_0.1.3              base64enc_0.1-3            
 [49] zlibbioc_1.50.0             ScaledMatrix_1.12.0         doSNOW_1.0.20              
 [52] GenomeInfoDbData_1.2.12     quadprog_1.5-8              SparseArray_1.4.3          
 [55] xtable_1.8-4                ade4_1.7-22                 doParallel_1.0.17          
 [58] evaluate_0.23               S4Arrays_1.4.0              preprocessCore_1.66.0      
 [61] hms_1.1.3                   glmnet_4.1-8                pulsar_0.3.11              
 [64] irlba_2.3.5.1               colorspace_2.1-0            magrittr_2.0.3             
 [67] later_1.3.2                 DECIPHER_3.0.0              cowplot_1.1.3              
 [70] scuttle_1.14.0              class_7.3-22                Hmisc_5.1-2                
 [73] pillar_1.9.0                nlme_3.1-164                huge_1.3.5                 
 [76] iterators_1.0.14            decontam_1.24.0             compiler_4.4.1             
 [79] beachmat_2.20.0             RSpectra_0.16-1             stringi_1.8.4              
 [82] biomformat_1.32.0           DescTools_0.99.57           minqa_1.2.6                
 [85] devtools_2.4.5              plyr_1.8.9                  crayon_1.5.2               
 [88] abind_1.4-5                 scater_1.32.0               bit_4.0.5                  
 [91] mia_1.12.0                  rootSolve_1.8.2.4           mixedCCA_1.6.2             
 [94] sandwich_3.1-0              fastcluster_1.2.6           codetools_0.2-20           
 [97] multcomp_1.4-25             BiocSingular_1.20.0         e1071_1.7-14               
[100] lmom_3.0                    GetoptLong_1.0.5            mime_0.12                  
[103] multtest_2.60.0             MultiAssayExperiment_1.30.1 splines_4.4.1              
[106] Rcpp_1.0.12                 sparseMatrixStats_1.16.0    qgraph_1.9.8               
[109] cellranger_1.1.0            knitr_1.46                  blob_1.2.4                 
[112] utf8_1.2.4                  clue_0.3-65                 lme4_1.1-35.3              
[115] pbivnorm_0.6.0              fs_1.6.4                    checkmate_2.3.1            
[118] DelayedMatrixStats_1.26.0   Rdpack_2.6                  pkgbuild_1.4.4             
[121] expm_0.999-9                gsl_2.1-8                   ggsignif_0.6.4             
[124] lavaan_0.6-17               Matrix_1.7-0                tzdb_0.4.0                 
[127] svglite_2.1.3               pkgconfig_2.0.3             tools_4.4.1                
[130] cachem_1.0.8                rbibutils_2.2.16            RSQLite_2.3.6              
[133] DBI_1.2.2                   numDeriv_2016.8-1.1         impute_1.78.0              
[136] fastmap_1.1.1               rmarkdown_2.26              scales_1.3.0               
[139] usethis_2.2.3               broom_1.0.5                 xlsx_0.6.5                 
[142] carData_3.0-5               rpart_4.1.23                farver_2.1.2               
[145] snow_0.4-4                  mgcv_1.9-1                  VGAM_1.1-10                
[148] foreign_0.8-86              cli_3.6.2                   lifecycle_1.0.4            
[151] mvtnorm_1.2-4               sessioninfo_1.2.2           bluster_1.14.0             
[154] backports_1.4.1             BiocParallel_1.38.0         timechange_0.3.0           
[157] gtable_0.3.5                rjson_0.2.21                ANCOMBC_2.6.0              
[160] parallel_4.4.1              ape_5.8                     CVXR_1.0-12                
[163] jsonlite_1.8.8              bit64_4.0.5                 Rtsne_0.17                 
[166] glasso_1.11                 yulab.utils_0.1.8           BiocNeighbors_1.22.0       
[169] orca_1.1-2                  lazyeval_0.2.2              shiny_1.8.1.1              
[172] dynamicTreeCut_1.63-1       htmltools_0.5.8.1           rJava_1.0-11               
[175] GO.db_3.19.1                glue_1.7.0                  treeio_1.28.0              
[178] mnormt_2.1.1                jpeg_0.1-10                 boot_1.3-30                
[181] igraph_2.0.3                R6_2.5.1                    fdrtool_1.2.17             
[184] labeling_0.4.3              Rmpfr_0.9-5                 xlsxjars_0.6.1             
[187] cluster_2.1.6               rngtools_1.5.2              pkgload_1.3.4              
[190] Rhdf5lib_1.26.0             nloptr_2.0.3                DirichletMultinomial_1.46.0
[193] DelayedArray_0.30.1         tidyselect_1.2.1            vipor_0.4.7                
[196] WGCNA_1.72-5                htmlTable_2.4.2             microbiome_1.26.0          
[199] xml2_1.3.6                  car_3.1-2                   AnnotationDbi_1.66.0       
[202] filematrix_1.3              rsvd_1.0.5                  munsell_0.5.1              
[205] data.table_1.15.4           htmlwidgets_1.6.4           rlang_1.1.3                
[208] remotes_2.5.0               lmerTest_3.1-3              fansi_1.0.6                
[211] beeswarm_0.4.0  


########### Notes:
ANCOMBC 2.6.0 was used for differential analyses.

Permutation tests were performed at High-Performance Computing Cluster with x86_64-pc-linux-gnu (64-bit), Rocky Linux 8.8 (Green Obsidian), and R version 4.3.2 (2023-10-31).




