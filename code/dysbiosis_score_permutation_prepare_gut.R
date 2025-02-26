rm(list = ls())
sample_type <- "Gut"
source("code/load.R")
output_dir <- file.path(output_dir, "Dysbiosis_Permutation", sample_type)
dir.create(output_dir, recursive = T, showWarnings = F)
devtools::load_all("code/ANCOMBC_mod_covariates")
set.seed(17385)

# Load data
phy_mic <- readRDS(file.path("data", "processed_data", tolower(sample_type), "LKT_phy.rds"))
phy_go <- readRDS(file.path("data", "processed_data", tolower(sample_type), "GO_phy.rds"))
phy_met_pos <- readRDS(file.path("data/processed_data", tolower(sample_type), paste0("metabolites", "_", "positive", "_phy.rds")))
phy_met_neg <- readRDS(file.path("data/processed_data", tolower(sample_type), paste0("metabolites", "_", "negative", "_phy.rds")))
phy_plasma_met_pos <- readRDS(file.path("data/processed_data", "plasma", paste0("metabolites", "_", "positive", "_phy.rds")))
phy_plasma_met_neg <- readRDS(file.path("data/processed_data", "plasma", paste0("metabolites", "_", "negative", "_phy.rds")))

data_list <- list(phy_mic, phy_go, phy_met_pos, phy_met_neg, phy_plasma_met_pos, phy_plasma_met_neg)
data_list_v1 <- lapply(data_list, function(x) subset_samples(x, visit == "v1"))
data_list_v1 <- lapply(data_list_v1, function(x) {
  sample_names(x) <- sample_data(x)$subjid
  return(x)
})

samples <- lapply(data_list_v1, function(x) sample_names(x))
comm_samples <- Reduce(intersect,samples)

# Subset to common sample and make sure the orders are the same
data_list_comm <- lapply(data_list_v1, function(x){
  x <- subset_samples(x, subjid %in% comm_samples)
  x <- x %>% ps_reorder(comm_samples)
  return(x)
})

meta <- sample_data(data_list_comm[[1]])
# Shuffle the labels 
samples_perm <- lapply(1:1000, function(i){
 sample(meta$status,size = nrow(meta),replace = F)
})

# species names
tse_mic <- mia::makeTreeSummarizedExperimentFromPhyloseq(phy_mic)
tse_species <- mia::agglomerateByRank(tse_mic, rank = "Species")
species <- rownames(rowData(tse_species))
species <- species[!grepl("__Missing|__sp.", species)]
species_df <- data.frame(Taxa = rep(species, each = 2), group = rep(c("cor1", "cor2")))

saveRDS(samples_perm, file.path(output_dir,"perms.rds"))
saveRDS(species_df, file.path(output_dir,"species_df.rds"))
saveRDS(data_list_comm, file.path(output_dir,"data_list_comm.rds"))

