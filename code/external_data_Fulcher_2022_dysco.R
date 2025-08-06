# Project PRJNA83636
rm(list = ls())
source("code/load.R")
source("code/bias_correction.R")
source("code/differential_correlations_functions.R")
source("code/create_heatmap.R")

method <- "spearman"
cov_cat <- c("geo_loc_name_country")
cov_num <- NULL
assay_name <- "counts"
tax_level <- "Species"
filt_th <- 0.001
prv_cut <- 0.5
threshold <- 0.3
pval_t <- 0.01
groups <- c("Control", "Case")
sample_type <- "Gut"


for (hiv in c("Negative", "Positive")) {
  
  output_dir <- file.path("Outputs/External_Dataset/Fulcher_2022", hiv)
  dir.create(output_dir, recursive = T, showWarnings = F)
  name <- paste0("Fulcher_2022_",hiv)
  
  # Load data
  input_dir <- "data/external_dataset/Fulcher_2022/processed_data"
  phy_obj_species <- readRDS(file.path(
    input_dir,
    paste0("LKT", "_phy.rds")
  ))
  phy_obj_go <- readRDS(file.path(
    input_dir,
    paste0("GO", "_phy.rds")
  ))

  data_list <- list(phy_obj_species, phy_obj_go)

  # Subset to case (g2) and control (g1)
  data_list_g2 <- lapply(data_list, function(x) subset_samples(x, Group == "Case" & HIVSerostatus == hiv))
  data_list_g1 <- lapply(data_list, function(x) subset_samples(x, Group == "Control"))

  # Bias correction
  g2_bias_corrected <- lapply(data_list_g2, function(x) {
    if ("Species" %in% colnames(tax_table(x))) {
      bias_correction(list(x),
        tax_level = "Species", pseudo = 0, prv_cut = prv_cut,
        lib_cut = 0, cov_cat = cov_cat, cov_num = cov_num
      )$residuals
    } else {
      bias_correction(list(x),
        tax_level = NULL, pseudo = 0, prv_cut = prv_cut,
        lib_cut = 0, cov_cat = cov_cat, cov_num = cov_num
      )$residuals
    }
  })

  g1_bias_corrected <- lapply(data_list_g1, function(x) {
    if ("Species" %in% colnames(tax_table(x))) {
      bias_correction(list(x),
        tax_level = "Species", pseudo = 0, prv_cut = prv_cut,
        lib_cut = 0, cov_cat = cov_cat, cov_num = cov_num
      )$residuals
    } else {
      bias_correction(list(x),
        tax_level = NULL, pseudo = 0, prv_cut = prv_cut,
        lib_cut = 0, cov_cat = cov_cat, cov_num = cov_num
      )$residuals
    }
  })

  # Gut species
  df.g1 <- t(g1_bias_corrected[[1]])
  df.g2 <- t(g2_bias_corrected[[1]])
  common_features <- intersect(colnames(df.g1), colnames(df.g2))
  common_features <- common_features[!grepl("s__sp.", common_features)]
  df.g1 <- df.g1[, common_features]
  df.g2 <- df.g2[, common_features]
  
  species_res <- dysco_score(
    data1 = list(df.g1), data2 = list(df.g2), 
    output_dir =output_dir ,prefix = paste0(name,"_gut_species")
  )

  # Gut Go
  go.g1 <- t(g1_bias_corrected[[2]])
  go.g2 <- t(g2_bias_corrected[[2]])
  common_features <- intersect(colnames(go.g1), colnames(go.g2))
  go.g1 <- go.g1[, common_features]
  go.g2 <- go.g2[, common_features]

  go_res <- dysco_score(
    data1 = list(df.g1, go.g1), data2 = list(df.g2, go.g2), 
    output_dir = output_dir,prefix = paste0(name,"_gut_go")
  )

}

# Heatmap
Fulcher_neg_species <- read.csv("Outputs/External_Dataset/Fulcher_2022/Negative/Fulcher_2022_Negative_gut_species_dysco.csv",check.names = F)
Fulcher_neg_go <- read.csv("Outputs/External_Dataset/Fulcher_2022/Negative/Fulcher_2022_Negative_gut_go_dysco.csv",check.names = F)
Fulcher_pos_species <- read.csv("Outputs/External_Dataset/Fulcher_2022/Positive/Fulcher_2022_Positive_gut_species_dysco.csv",check.names = F)
Fulcher_pos_go <- read.csv("Outputs/External_Dataset/Fulcher_2022/Positive/Fulcher_2022_Positive_gut_go_dysco.csv",check.names = F)

list_res <- list(Fulcher_neg_species,Fulcher_neg_go,Fulcher_pos_species,Fulcher_pos_go)
names(list_res) <- c("Pre-HIV Gut Species","Pre-HIV Gut GO",
                     "Post-HIV Gut Species","Post-HIV Gut GO")

create_heatmap (list_res, p_alpha = 0.01,bh_alpha = 0.1,
                col_names = names(list_res),top_n = 10,
                sample_type = "Gut",output_dir = "Outputs/External_Dataset/Fulcher_2022")


