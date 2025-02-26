# Set your working directory
# setwd("~/HIV_Microbiome_2025")

##############################################
#Prepare data for the down-stream analysis
##############################################
source("code/gut_prepare_data.R")
source("code/oral_prepare_data.R")

##############################################
#Differential Abundance Testing
##############################################
# Microbiome 
system("Rscript code/differential_abundance.R Gut")
system("Rscript code/differential_abundance.R Oral")

# Metabolomics
system("Rscript code/metabolomics_differential_abundance.R Gut")
system("Rscript code/metabolomics_differential_abundance.R Oral")
system("Rscript code/metabolomics_differential_abundance.R Plasma")

# SCFA
source("code/SCFA_differential_analysis.R")

##############################################
#DIABLO
##############################################
source("code/diablo_gut.R")

##############################################
#Differential Correlations and Dysbiosos Score
##############################################
# Gut
source("code/dysbiosis_score_gut.R")
source("code/networks_gut.R")
source("code/subset_networks_gut.R")
# Controlled for sexual activity (number of sexual partners)
source("code/dysbiosis_score_gut_controlled_for_sa.R")

## Permutation test (n = 1000), performed at NIEHS HPC using the following scripts
# Note: p-value results can be found at:
  # Without controlling for sexual activity:
    # Outputs/Dysbiosis/Gut/Scores/pvals_p_diff.csv
  # With controlling for sexual activity:
    # Outputs/Dysbiosis_Controlled_for_SA/Gut/Scores/pvals_p_diff.csv

#Rscript dysbiosis_score_permutation_prepare_gut.R
#Rscript code/dysbiosis_score_permutation_gut.R Gut i # run by swarm, i is the permutation number, see code/R_script_permutation_gut.swarm
#Rscript dysbiosis_score_permutation_prepare_gut_controlled_for_sa.R
#Rscript code/dysbiosis_score_permutation_gut_controlled_for_sa.R Gut i # run by swarm, i is the permutation number, see code/R_script_permutation_gut_controlled_for_sa.swarm 
#Rscript code/dysbiosis_score_permutation_pvals.R Gut Dysbiosis_Permutation Dysbiosis # Calculate p-values
#Rscript code/dysbiosis_score_permutation_pvals.R Gut Dysbiosis_Permutation_Controlled_for_SA Dysbiosis_Controlled_for_SA # Calculate p-values


# Oral
source("code/dysbiosis_score_oral.R")
source("code/networks_oral.R")
source("code/subset_networks_oral.R")
source("code/dysbiosis_score_oral_controlled_for_sa.R")
## Permutation test (n = 1000), performed at NIEHS HPC
# Note: p-value results can be found at:
  # Without controlling for sexual activity:
    # Outputs/Dysbiosis/Oral/Scores/pvals_p_diff.csv
  # With controlling for sexual activity:
    # Outputs/Dysbiosis_Controlled_for_SA/Oral/Scores/pvals_p_diff.csv

#Rscript dysbiosis_score_permutation_prepare_oral.R
#Rscript code/dysbiosis_score_permutation_oral.R Oral i # run by swarm, i is the permutation number, see code/R_script_permutation_oral.swarm
#Rscript dysbiosis_score_permutation_prepare_oral_controlled_for_SA.R
#Rscript code/dysbiosis_score_permutation_oral_controlled_for_SA.R Oral i # run by swarm, i is the permutation number, see code/R_script_permutation_oral_controlled_for_sa.swarm 
#Rscript code/dysbiosis_score_permutation_pvals.R Oral Dysbiosis_Permutation Dysbiosis # Calculate p-values
#Rscript code/dysbiosis_score_permutation_pvals.R Oral Dysbiosis_Permutation_Controlled_for_SA Dysbiosis_Controlled_for_SA # Calculate p-values

# Figures for dysbiosis scores
system("Rscript dysbiosis_score_scatterplots.R Gut")
system("Rscript dysbiosis_score_scatterplots.R Oral")
system("Rscript dysbiosis_score_heatmap.R Gut")
system("Rscript dysbiosis_score_heatmap.R Oral")


##############################################
#Correlation between microbiome and cytokines
##############################################
source("code/Cytokine_Microbiome_Correlation.R")


##############################################
#Supplementary Tables
##############################################

source("code/Supplementary_Tables.R")
