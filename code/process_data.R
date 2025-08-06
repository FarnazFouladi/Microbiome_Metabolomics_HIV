source("code/load.R")
source("code/bias_correction.R")

output_dir <- file.path("Outputs", "Processed_data")
cov_cat <- c("loc", "abx_use")
cov_num <- c("age")

dir.create(output_dir, showWarnings = F, recursive = T)
prv_cut <- 0.5

# Load data
phy_species <- readRDS("data/processed_data/gut/LKT_phy.rds")
phy_go <- readRDS("data/processed_data/gut/GO_phy.rds")
phy_met_pos <- readRDS("data/processed_data/gut/metabolites_positive_phy.rds")
phy_met_neg <- readRDS("data/processed_data/gut/metabolites_negative_phy.rds")
phy_plasma_met_pos <- readRDS("data/processed_data/plasma/metabolites_positive_phy.rds")
phy_plasma_met_neg <- readRDS("data/processed_data/plasma/metabolites_negative_phy.rds")
phy_species_o <- readRDS("data/processed_data/oral/LKT_phy.rds")
phy_go_o <- readRDS("data/processed_data/oral/GO_phy.rds")
phy_met_pos_o <- readRDS("data/processed_data/oral/metabolites_positive_phy.rds")
phy_met_neg_o <- readRDS("data/processed_data/oral/metabolites_negative_phy.rds")

data_list <- list(
  phy_species, phy_go, phy_met_pos, phy_met_neg, phy_plasma_met_pos, phy_plasma_met_neg,
  phy_species_o, phy_go_o, phy_met_pos_o, phy_met_neg_o
)
names(data_list) <- c(
  "gut_species", "gut_go", "gut_mtb_p", "gut_mtb_n",
  "plasma_mtb_p", "plasma_mtb_n", "oral_species", "oral_go", "oral_mtb_p", "oral_mtb_n"
)

# Set sample names to subject ID
data_list <- lapply(data_list, function(x) {
  sample_names(x) <- sample_data(x)$subjid
  return(x)
})

# Subset to NC (Non-HIV) and SC (Pre-HIV)
data_list_sc <- lapply(data_list, function(x) subset_samples(x, status == "sc"))
data_list_nc <- lapply(data_list, function(x) subset_samples(x, status == "nc"))

# Bias correction
sc_bias_corrected <- lapply(data_list_sc, function(x) {
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

nc_bias_corrected <- lapply(data_list_nc, function(x) {
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
write_rds(sc_bias_corrected, file.path(output_dir, "processed_pre_hiv.rds"))
write_rds(nc_bias_corrected, file.path(output_dir, "processed_non_hiv.rds"))


# Subset to g1, g2, g3, g4 (4 sexual activity groups)
data_list_g1 <- lapply(data_list, function(x) subset_samples(x, group1 == "g1"))
data_list_g2 <- lapply(data_list, function(x) subset_samples(x, group1 == "g2"))
data_list_g3 <- lapply(data_list, function(x) subset_samples(x, group1 == "g3"))
data_list_g4 <- lapply(data_list, function(x) subset_samples(x, group1 == "g4"))

# Bias correction
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

g3_bias_corrected <- lapply(data_list_g3, function(x) {
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
g4_bias_corrected <- lapply(data_list_g4, function(x) {
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



write_rds(g1_bias_corrected, file.path(output_dir, "processed_g1.rds"))
write_rds(g2_bias_corrected, file.path(output_dir, "processed_g2.rds"))
write_rds(g3_bias_corrected, file.path(output_dir, "processed_g3.rds"))
write_rds(g4_bias_corrected, file.path(output_dir, "processed_g4.rds"))
