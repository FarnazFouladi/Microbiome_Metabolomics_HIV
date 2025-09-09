# Set your working directory at Microbiome_Metabolomics_HIV
# setwd("~/Microbiome_Metabolomics_HIV")

##############################################
#Prepare data for the down-stream analysis
##############################################
source("code/gut_prepare_data.R")
source("code/oral_prepare_data.R")


##############################################
#Summary table
##############################################
source("code/summary_table.R")


###################################################
#Differential Abundance Testing for Pre-HIV vs. Non-HIV
###################################################
# Microbiome 
system("Rscript code/differential_abundance.R Gut")
system("Rscript code/differential_abundance.R Oral")

# Metabolomics
system("Rscript code/metabolomics_differential_abundance.R Gut")
system("Rscript code/metabolomics_differential_abundance.R Oral")
system("Rscript code/metabolomics_differential_abundance.R Plasma")

# SCFA
source("code/SCFA_differential_analysis.R")

###################################################
#Trend analysis for Sexual Activity
###################################################
# Note: Trend analysis using ANCOM-BC2 with 1000 bootstrapping takes some time to run.
# Microbiome 
system("Rscript code/differential_abundance_SA.R Gut")
system("Rscript code/differential_abundance_SA.R Oral")

# Metabolomics
system("Rscript code/metabolomics_differential_abundance_SA.R Gut")
system("Rscript code/metabolomics_differential_abundance_SA.R Oral")
system("Rscript code/metabolomics_differential_abundance_SA.R Plasma")

# Create a heatmap to compare results from differential abundance versus trend analysis
source("code/differential_abundance_comparisons.R")

##############################################
#DIABLO
##############################################
source("code/diablo_gut.R")


##############################################
#Differential Correlations and DYSCO
##############################################
source("code/process_data.R")
# Differential correlation
source("code/differential_correlation.R")
# DYSCO
system("Rscript code/dysco_heatmap.R") 
# Subnetworks
system("Rscript code/subset_networks_gut.R")
system("Rscript code/subset_networks_oral.R")

##############################################
#Correlations trend analysis
##############################################
system("Rscript code/correlation_trend_analysis.R")
system("Rscript code/correlation_trend_heatmap.R")


##############################################
#External datasets
##############################################
source("code/external_data_Fulcher_2022_preprocess.R")
source("code/external_data_Fulcher_2022_dysco.R")

source("code/external_data_Garcia_2024_preprocs.R")
source("code/external_data_Garcia_2024_dysco.R")

source("code/external_data_Rocafort_2024_dysco.R")
source("code/external_data_Armstrong_2018.R")

source("code/dysco_heatmap_across_studies.R")

##############################################
#DYSCO validation
##############################################

source("code/differential_correlations_validation.R")

##############################################
#Supplementary Tables
##############################################

source("code/Supplementary_Tables.R")

