rm(list = ls())
source("code/load.R")
source("code/bias_correction.R")
source("code/differential_correlations_functions.R")
source("code/create_heatmap.R")

method <- "spearman"
cov_cat <- c("Ethnicity", "HIV_Profile")
name <- "Garcia_2024"
cov_num <- NULL
assay_name <- "counts"
tax_level <- "Species"
filt_th <- 0.001
prv_cut <- 0.5
threshold <- 0.3
pval_t <- 0.01
filt_th <- 0.001
groups <- c("negative", "positive") # HIV_serostatus

sample_type <- "Gut"
output_dir <- file.path("Outputs/External_Dataset/Garcia_2024")
dir.create(output_dir, recursive = T, showWarnings = F)
# Load data
input_dir <- "data/external_dataset/Garcia_2024/processed_data"
phy_obj_species <- readRDS(file.path(
  input_dir,
  paste0("LKT", "_phy.rds")
))
phy_obj_go <- readRDS(file.path(
  input_dir,
  paste0("GO", "_phy.rds")
))

data_list <- list(phy_obj_species, phy_obj_go)

# Subset to control (g1) and HIV positive (g2)
data_list_g2 <- lapply(data_list, function(x) subset_samples(x, HIV_serostatus == groups[2]))
data_list_g1 <- lapply(data_list, function(x) subset_samples(x, HIV_serostatus == groups[1]))


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

n2 <- nrow(df.g2)
n1 <- nrow(df.g1)
species_res <- dysco_score(
  data1 = list(df.g1), data2 = list(df.g2),
  output_dir = output_dir,prefix = paste0(name,"_gut_species")
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

# Heatmap
Garcia_species <- read.csv("Outputs/External_Dataset/Garcia_2024/Garcia_2024_gut_species_dysco.csv",check.names = F)
Garcia_go <- read.csv("Outputs/External_Dataset/Garcia_2024/Garcia_2024_gut_go_dysco.csv",check.names = F)

list_res <- list(Garcia_species,Garcia_go)
names(list_res) <- c("Post-HIV Gut Species","Post-HIV Gut GO")

create_heatmap (list_res, p_alpha = 0.01,bh_alpha = 0.1,
                col_names = names(list_res),top_n = 10,
                sample_type = "Gut",output_dir = "Outputs/External_Dataset/Garcia_2024")

