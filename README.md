
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Microbiome_Metabolomics_HIV

<!-- badges: start -->
<!-- badges: end -->

This repository includes data, R codes, and outputs for the following
paper:

**Fouladi et al. A taxon-specific measurement of disruption in a
multi-modal study of microbiomes and metabolomes reveals system-wide
dysbiosis preceding HIV-1 infection. 2025. Nature Communications.**

Please refer to *data availability* in the paper to obtain the metadata
for this project

The main script to run all the analyses is
*code/data_analysis_workflow.R*

**Platform and Software**

``` r
sessionInfo()
#> R version 4.4.1 (2024-06-14)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS 15.5
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: America/New_York
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] grid      stats4    stats     graphics  grDevices utils     datasets 
#> [8] methods   base     
#> 
#> other attached packages:
#>  [1] microViz_0.12.4                 kableExtra_1.4.0               
#>  [3] readxl_1.4.3                    RColorBrewer_1.1-3             
#>  [5] NetCoMi_1.1.0                   SpiecEasi_1.1.3                
#>  [7] rstatix_0.7.2                   vegan_2.6-6                    
#>  [9] lattice_0.22-6                  permute_0.9-7                  
#> [11] lubridate_1.9.3                 forcats_1.0.0                  
#> [13] stringr_1.5.1                   dplyr_1.1.4                    
#> [15] purrr_1.0.2                     readr_2.1.5                    
#> [17] tidyr_1.3.1                     tidyverse_2.0.0                
#> [19] tibble_3.2.1                    viridis_0.6.5                  
#> [21] viridisLite_0.4.2               circlize_0.4.16                
#> [23] ComplexHeatmap_2.20.0           gridExtra_2.3                  
#> [25] ggpubr_0.6.0                    ggrepel_0.9.5                  
#> [27] ggplot2_3.5.1                   phyloseq_1.48.0                
#> [29] TreeSummarizedExperiment_2.12.0 Biostrings_2.72.0              
#> [31] XVector_0.44.0                  SingleCellExperiment_1.26.0    
#> [33] SummarizedExperiment_1.34.0     Biobase_2.64.0                 
#> [35] GenomicRanges_1.56.0            GenomeInfoDb_1.40.0            
#> [37] IRanges_2.38.0                  S4Vectors_0.42.0               
#> [39] BiocGenerics_0.50.0             MatrixGenerics_1.16.0          
#> [41] matrixStats_1.3.0              
#> 
#> loaded via a namespace (and not attached):
#>   [1] splines_4.4.1           cellranger_1.1.0        preprocessCore_1.66.0  
#>   [4] rpart_4.1.23            lifecycle_1.0.4         Rdpack_2.6             
#>   [7] fastcluster_1.2.6       doParallel_1.0.17       MASS_7.3-60.2          
#>  [10] backports_1.4.1         magrittr_2.0.3          Hmisc_5.1-2            
#>  [13] rmarkdown_2.26          yaml_2.3.8              qgraph_1.9.8           
#>  [16] pbapply_1.7-2           SPRING_1.0.4            DBI_1.2.2              
#>  [19] ade4_1.7-22             abind_1.4-5             zlibbioc_1.50.0        
#>  [22] quadprog_1.5-8          yulab.utils_0.1.8       nnet_7.3-19            
#>  [25] GenomeInfoDbData_1.2.12 irlba_2.3.5.1           tidytree_0.4.6         
#>  [28] svglite_2.1.3           codetools_0.2-20        DelayedArray_0.30.1    
#>  [31] xml2_1.3.6              tidyselect_1.2.1        shape_1.4.6.1          
#>  [34] UCSC.utils_1.0.0        dynamicTreeCut_1.63-1   base64enc_0.1-3        
#>  [37] jsonlite_1.8.8          GetoptLong_1.0.5        multtest_2.60.0        
#>  [40] Formula_1.2-5           survival_3.6-4          iterators_1.0.14       
#>  [43] systemfonts_1.0.6       foreach_1.5.2           tools_4.4.1            
#>  [46] treeio_1.28.0           snow_0.4-4              filematrix_1.3         
#>  [49] Rcpp_1.0.12             glue_1.8.0              mnormt_2.1.1           
#>  [52] SparseArray_1.4.3       xfun_0.44               mgcv_1.9-1             
#>  [55] withr_3.0.0             fastmap_1.1.1           rhdf5filters_1.16.0    
#>  [58] fansi_1.0.6             digest_0.6.35           timechange_0.3.0       
#>  [61] R6_2.5.1                colorspace_2.1-0        GO.db_3.19.1           
#>  [64] gtools_3.9.5            jpeg_0.1-10             RSQLite_2.3.6          
#>  [67] utf8_1.2.4              generics_0.1.3          corpcor_1.6.10         
#>  [70] data.table_1.15.4       pulsar_0.3.11           httr_1.4.7             
#>  [73] htmlwidgets_1.6.4       S4Arrays_1.4.0          pkgconfig_2.0.3        
#>  [76] gtable_0.3.5            blob_1.2.4              impute_1.78.0          
#>  [79] pcaPP_2.0-4             htmltools_0.5.8.1       lavaan_0.6-17          
#>  [82] carData_3.0-5           biomformat_1.32.0       clue_0.3-65            
#>  [85] scales_1.3.0            orca_1.1-2              png_0.1-8              
#>  [88] doSNOW_1.0.20           corrplot_0.92           knitr_1.46             
#>  [91] rstudioapi_0.16.0       tzdb_0.4.0              reshape2_1.4.4         
#>  [94] rjson_0.2.21            checkmate_2.3.1         nlme_3.1-164           
#>  [97] cachem_1.0.8            rhdf5_2.48.0            GlobalOptions_0.1.2    
#> [100] rootSolve_1.8.2.4       parallel_4.4.1          foreign_0.8-86         
#> [103] AnnotationDbi_1.66.0    pillar_1.9.0            vctrs_0.6.5            
#> [106] VGAM_1.1-10             car_3.1-2               huge_1.3.5             
#> [109] cluster_2.1.6           htmlTable_2.4.2         pbivnorm_0.6.0         
#> [112] evaluate_0.23           mvtnorm_1.2-4           cli_3.6.5              
#> [115] compiler_4.4.1          rlang_1.1.6             crayon_1.5.2           
#> [118] ggsignif_0.6.4          fdrtool_1.2.17          plyr_1.8.9             
#> [121] fs_1.6.4                psych_2.4.3             stringi_1.8.4          
#> [124] WGCNA_1.72-5            BiocParallel_1.38.0     munsell_0.5.1          
#> [127] lazyeval_0.2.2          glmnet_4.1-8            Matrix_1.7-0           
#> [130] glasso_1.11             hms_1.1.3               bit64_4.0.5            
#> [133] Rhdf5lib_1.26.0         KEGGREST_1.44.0         rbibutils_2.2.16       
#> [136] igraph_2.0.3            broom_1.0.8             memoise_2.0.1          
#> [139] mixedCCA_1.6.2          bit_4.0.5               ape_5.8
```
