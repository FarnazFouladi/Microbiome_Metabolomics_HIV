# DYSCO scores for oral samples
rm(list = ls())
sample_type <- "Oral"
source("code/load.R")
output_dir <- file.path("Outputs", "Dysbiosis", sample_type)
dir.create(output_dir, recursive = T, showWarnings = F)
devtools::load_all("code/ANCOMBC_mod_covariates")
source("code/dysbiosis_score_functions.R")
method <- "spearman"
cov_cat <- c("abx_use", "loc")
cov_num <- c("age")
assay_name <- c("counts", "counts")
tax_level <- c(NULL, NULL)
prv_cut <- 0.5
threshold <- 0.3
pval_t <- 0.01

# Load data
phy_mic <- readRDS(file.path("data", "processed_data", tolower(sample_type), "LKT_phy.rds"))
phy_go <- readRDS(file.path("data", "processed_data", tolower(sample_type), "GO_phy.rds"))
phy_met_pos <- readRDS(file.path("data/processed_data", tolower(sample_type), paste0("metabolites", "_", "positive", "_phy.rds")))
phy_met_neg <- readRDS(file.path("data/processed_data", tolower(sample_type), paste0("metabolites", "_", "negative", "_phy.rds")))
phy_plasma_met_pos <- readRDS(file.path("data/processed_data", "plasma", paste0("metabolites", "_", "positive", "_phy.rds")))
phy_plasma_met_neg <- readRDS(file.path("data/processed_data", "plasma", paste0("metabolites", "_", "negative", "_phy.rds")))
phy_mic1 <- readRDS(file.path("data", "processed_data", "gut", "LKT_phy.rds"))

data_list <- list(phy_mic, phy_go, phy_met_pos, phy_met_neg, phy_plasma_met_pos, phy_plasma_met_neg, phy_mic1)
data_list_v1 <- lapply(data_list, function(x) subset_samples(x, visit == "v1"))
data_list_v1 <- lapply(data_list_v1, function(x) {
  sample_names(x) <- sample_data(x)$subjid
  return(x)
})

# Common samples
samples <- lapply(data_list_v1, function(x) sample_names(x))
comm_samples <- Reduce(intersect, samples)

# Subset to common sample and make sure the orders are the same
data_list_comm <- lapply(data_list_v1, function(x) {
  x <- subset_samples(x, subjid %in% comm_samples)
  x <- x %>% ps_reorder(comm_samples)
  return(x)
})


data_list_sc <- lapply(data_list_comm, function(x) subset_samples(x, status == "sc"))
data_list_nc <- lapply(data_list_comm, function(x) subset_samples(x, status == "nc"))

# species
tse_mic <- mia::makeTreeSummarizedExperimentFromPhyloseq(phy_mic)
tse_species <- mia::agglomerateByRank(tse_mic, rank = "Species")
species <- rownames(rowData(tse_species))
species <- species[!grepl("__Missing|__sp.", species)]
species_df <- data.frame(Taxa = rep(species, each = 2), group = rep(c("cor1", "cor2")))

data_list_sc_tse <- lapply(data_list_sc, function(x) mia::makeTreeSummarizedExperimentFromPhyloseq(x))
data_list_nc_tse <- lapply(data_list_nc, function(x) mia::makeTreeSummarizedExperimentFromPhyloseq(x))

# Oral species
species <- dysbiosis_score(
  data_g1 = data_list_nc_tse[c(1)],
  data_g2 = data_list_sc_tse[c(1)],
  assay_name, tax_level = "Species", prv_cut, method, cov_cat, cov_num, threshold, pval_t, name = "species", two_datasets = FALSE, output_dir = file.path(output_dir, paste0(sample_type, "_Species"))
)


# Gut species and oral species
species2 <- dysbiosis_score(
  data_g1 = data_list_nc_tse[c(1, 7)],
  data_g2 = data_list_sc_tse[c(1, 7)],
  assay_name, tax_level = c("Species", "Species"), prv_cut, method, cov_cat, cov_num, threshold, pval_t, name = "species2", two_datasets = TRUE, output_dir = file.path(output_dir, paste0(ifelse(sample_type == "Gut", "Oral", "Gut"), "_Species"))
)


# Convert to Species
data_list_sc_tse[[1]] <- mia::agglomerateByRank(data_list_sc_tse[[1]], rank = "Species")
data_list_sc_tse[[7]] <- mia::agglomerateByRank(data_list_sc_tse[[7]], rank = "Species")
data_list_nc_tse[[1]] <- mia::agglomerateByRank(data_list_nc_tse[[1]], rank = "Species")
data_list_nc_tse[[7]] <- mia::agglomerateByRank(data_list_nc_tse[[7]], rank = "Species")

# Oral species and oral GO
species_go <- dysbiosis_score(
  data_g1 = data_list_nc_tse[c(1, 2)],
  data_g2 = data_list_sc_tse[c(1, 2)],
  assay_name, tax_level = c(NULL, NULL), prv_cut, method, cov_cat, cov_num, threshold, pval_t, name = "go", two_datasets = TRUE, output_dir = file.path(output_dir, paste0(sample_type, "_GO"))
)

# Oral species and oral metabolites positive channel
species_mtb_pos <- dysbiosis_score(
  data_g1 = data_list_nc_tse[c(1, 3)],
  data_g2 = data_list_sc_tse[c(1, 3)],
  assay_name, tax_level = c(NULL, NULL), prv_cut, method, cov_cat, cov_num, threshold, pval_t, name = "mtb_pos", two_datasets = TRUE, output_dir = file.path(output_dir, paste0(sample_type, "_MTB_Pos"))
)


# Oral species and oral metabolites negative channel
species_mtb_neg <- dysbiosis_score(
  data_g1 = data_list_nc_tse[c(1, 4)],
  data_g2 = data_list_sc_tse[c(1, 4)],
  assay_name, tax_level = c(NULL, NULL), prv_cut, method, cov_cat, cov_num, threshold, pval_t, name = "mtb_neg", two_datasets = TRUE, output_dir = file.path(output_dir, paste0(sample_type, "_MTB_Neg"))
)


# Oral species and plasma metabolites positive channel
species_plasma_mtb_pos <- dysbiosis_score(
  data_g1 = data_list_nc_tse[c(1, 5)],
  data_g2 = data_list_sc_tse[c(1, 5)],
  assay_name, tax_level = c(NULL, NULL), prv_cut, method, cov_cat, cov_num, threshold, pval_t, name = "plasma_mtb_pos", two_datasets = TRUE, output_dir = file.path(output_dir, "Plasma_MTB_Pos")
)


# Oral species and plasma metabolites negative channel
species_plasma_mtb_neg <- dysbiosis_score(
  data_g1 = data_list_nc_tse[c(1, 6)],
  data_g2 = data_list_sc_tse[c(1, 6)],
  assay_name, tax_level = c(NULL, NULL), prv_cut, method, cov_cat, cov_num, threshold, pval_t, name = "plasma_mtb_neg", two_datasets = TRUE, output_dir = file.path(output_dir, "Plasma_MTB_Neg")
)

# Number of differential correlations
species_all_p <- species_df %>%
  left_join(species[[1]]) %>%
  left_join(species_go[[1]]) %>%
  left_join(species_mtb_pos[[1]]) %>%
  left_join(species_mtb_neg[[1]]) %>%
  left_join(species_plasma_mtb_pos[[1]]) %>%
  left_join(species_plasma_mtb_neg[[1]]) %>%
  left_join(species2[[1]])

# Number of total features
species_fm <- species_df %>%
  left_join(species[[2]]) %>%
  left_join(species_go[[2]]) %>%
  left_join(species_mtb_pos[[2]]) %>%
  left_join(species_mtb_neg[[2]]) %>%
  left_join(species_plasma_mtb_pos[[2]]) %>%
  left_join(species_plasma_mtb_neg[[2]]) %>%
  left_join(species2[[2]])


species_all_p[is.na(species_all_p)] <- 0
stopifnot(all.equal(species_all_p$Taxa, species_fm$Taxa))
stopifnot(all.equal(colnames(species_all_p), colnames(species_fm)))

# Get the proportion of differential correlations
species_all_p2 <- (species_all_p[, 3:ncol(species_all_p)] / species_fm[, 3:ncol(species_fm)]) * 100
species_all_p2[is.na(species_all_p2)] <- 0
species_all_p2$Taxa <- species_all_p$Taxa
species_all_p2$group <- species_all_p$group


species_all_p2 <- species_all_p2 %>% relocate(c("Taxa", "group"), .before = species)

# Average between the two modes of metabolites
species_all_p3 <- data.frame(
  Taxa = species_all_p2$Taxa,
  group = species_all_p2$group,
  species = species_all_p2$species,
  go = species_all_p2$go,
  mtb = rowMeans(species_all_p2[, c("mtb_pos", "mtb_neg")], na.rm = T),
  plasma_mtb = rowMeans(species_all_p2[, c("plasma_mtb_pos", "plasma_mtb_neg")], na.rm = T),
  species2 = species_all_p2$species2
)

species_all_p3[is.na(species_all_p3)] <- 0

species_p_wide <- species_all_p3 %>%
  mutate(group = ifelse(group == "cor1", "NC", "SC")) %>%
  pivot_wider(names_from = group, values_from = colnames(species_all_p3)[3:ncol(species_all_p3)])

species_p_wide_nc <- species_p_wide %>%
  column_to_rownames("Taxa") %>%
  dplyr::select(grep("_NC", colnames(.)))

species_p_wide_sc <- species_p_wide %>%
  column_to_rownames("Taxa") %>%
  dplyr::select(grep("_SC", colnames(.)))

# Difference in scores between SC and NC
species_p_wide_diff <- species_p_wide_sc - species_p_wide_nc

species_p_long <- species_all_p3 %>%
  mutate(group = ifelse(group == "cor1", "NC", "SC")) %>%
  filter(Taxa %in% species_p_wide$Taxa)


fig_dir <- file.path(output_dir, "Figures")
score_dir <- file.path(output_dir, "Scores")
dir.create(fig_dir, showWarnings = F, recursive = T)
dir.create(score_dir, showWarnings = F, recursive = T)

write.csv(species_p_long, file.path(score_dir, paste0("scores_long_p.csv")), row.names = F, quote = F)
write.csv(species_p_wide, file.path(score_dir, paste0("scores_wide_p.csv")), row.names = F, quote = F)
write.csv(species_p_wide_diff, file.path(score_dir, paste0("scores_wide_diff_p.csv")), row.names = T, quote = F)
