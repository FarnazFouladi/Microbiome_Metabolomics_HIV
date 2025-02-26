args <- commandArgs(TRUE)
sample_type <-args[1]
i <- as.numeric(args[2])

devtools::load_all("code/ANCOMBC_mod_covariates")
source("code/load.R")
source("code/dysbiosis_score_functions_permutation.R")
sample_type = sample_type
output_dir <- file.path("Outputs/Dysbiosis_Permutation",sample_type)
samples_perm <- readRDS(file.path(output_dir, "perms.rds"))
species_df <- readRDS(file.path(output_dir, "species_df.rds"))
data_list_comm <- readRDS(file.path(output_dir, "data_list_comm.rds"))

print(i)
method <- "spearman"
cov_cat <- c("abx_use", "loc")
cov_num <- c("age")
assay_name <- c("counts", "counts")
tax_level <- c(NULL, NULL)
prv_cut <- 0.5
threshold <- 0.3
pval_t <- 0.01

score_dir <- file.path(output_dir, "Scores")
dir.create(score_dir, showWarnings = F, recursive = T)


data_list_perm <- lapply(1:length(data_list_comm), function(x) {
  ps <- data_list_comm[[x]] %>% ps_mutate(status = samples_perm[[i]])
  return(ps)
})

data_list_sc <- lapply(data_list_perm, function(x) subset_samples(x, status == "sc"))
data_list_nc <- lapply(data_list_perm, function(x) subset_samples(x, status == "nc"))


data_list_sc_tse <- lapply(data_list_sc, function(x) mia::makeTreeSummarizedExperimentFromPhyloseq(x))
data_list_nc_tse <- lapply(data_list_nc, function(x) mia::makeTreeSummarizedExperimentFromPhyloseq(x))

# Gut taxa
species <- dysbiosis_score(
  data_g1 = data_list_nc_tse[c(1)],
  data_g2 = data_list_sc_tse[c(1)],
  assay_name, tax_level = "Species", prv_cut, method, cov_cat, cov_num, threshold, pval_t, name = "species", two_datasets = FALSE, output_dir = file.path(output_dir, paste0(sample_type, "_Species")),perm = i)


# Convert to Species
data_list_sc_tse[[1]] <- mia::agglomerateByRank(data_list_sc_tse[[1]], rank = "Species")
data_list_nc_tse[[1]] <- mia::agglomerateByRank(data_list_nc_tse[[1]], rank = "Species")

# Gut species and gut GO
species_go <- dysbiosis_score(
  data_g1 = data_list_nc_tse[c(1, 2)],
  data_g2 = data_list_sc_tse[c(1, 2)],
  assay_name, tax_level = c(NULL, NULL), prv_cut, method, cov_cat, cov_num, threshold, pval_t, name = "go", two_datasets = TRUE, output_dir = file.path(output_dir, paste0(sample_type, "_GO")),perm = i
)


# Gut species and gut metabolites positive channel
species_mtb_pos <- dysbiosis_score(
  data_g1 = data_list_nc_tse[c(1, 3)],
  data_g2 = data_list_sc_tse[c(1, 3)],
  assay_name, tax_level = c(NULL, NULL), prv_cut, method, cov_cat, cov_num, threshold, pval_t, name = "mtb_pos", two_datasets = TRUE, output_dir = file.path(output_dir, paste0(sample_type, "_MTB_Pos")),perm = i
)

# Gut species and gut metabolites negative channel
species_mtb_neg <- dysbiosis_score(
  data_g1 = data_list_nc_tse[c(1, 4)],
  data_g2 = data_list_sc_tse[c(1, 4)],
  assay_name, tax_level = c(NULL, NULL), prv_cut, method, cov_cat, cov_num, threshold, pval_t, name = "mtb_neg", two_datasets = TRUE, output_dir = file.path(output_dir, paste0(sample_type, "_MTB_Neg")),perm = i
)

# Gut species and plasma metabolites positive channel
species_plasma_mtb_pos <- dysbiosis_score(
  data_g1 = data_list_nc_tse[c(1, 5)],
  data_g2 = data_list_sc_tse[c(1, 5)],
  assay_name, tax_level = c(NULL, NULL), prv_cut, method, cov_cat, cov_num, threshold, pval_t, name = "plasma_mtb_pos", two_datasets = TRUE, output_dir = file.path(output_dir, "Plasma_MTB_Pos"),perm = i
)

# Gut species and plasma metabolites negative channel
species_plasma_mtb_neg <- dysbiosis_score(
  data_g1 = data_list_nc_tse[c(1, 6)],
  data_g2 = data_list_sc_tse[c(1, 6)],
  assay_name, tax_level = c(NULL, NULL), prv_cut, method, cov_cat, cov_num, threshold, pval_t, name = "plasma_mtb_neg", two_datasets = TRUE, output_dir = file.path(output_dir, "Plasma_MTB_Neg"),perm = i
)

# Get the number of species
species_fm <- species_df %>%
  left_join(species[[2]]) %>%
  left_join(species_go[[2]]) %>%
  left_join(species_mtb_pos[[2]]) %>%
  left_join(species_mtb_neg[[2]]) %>%
  left_join(species_plasma_mtb_pos[[2]]) %>%
  left_join(species_plasma_mtb_neg[[2]])


# get the number of differential correlations per species
species_all_p <- species_df %>%
  left_join(species[[1]]) %>%
  left_join(species_go[[1]]) %>%
  left_join(species_mtb_pos[[1]]) %>%
  left_join(species_mtb_neg[[1]]) %>%
  left_join(species_plasma_mtb_pos[[1]]) %>%
  left_join(species_plasma_mtb_neg[[1]]) 


species_all_p[is.na(species_all_p)] <- 0
stopifnot(all.equal(species_all_p$Taxa, species_fm$Taxa))
stopifnot(all.equal(colnames(species_all_p), colnames(species_fm)))

# Get the proportions of differential correlations
species_all_p2 <- (species_all_p[, 3:ncol(species_all_p)] / species_fm[, 3:ncol(species_fm)]) * 100
species_all_p2[is.na(species_all_p2)] <- 0
species_all_p2$Taxa <- species_all_p$Taxa
species_all_p2$group <- species_all_p$group

species_all_p2 <- species_all_p2 %>% relocate(c("Taxa", "group"), .before = species)

# Average scores across the two channels of metabolomics
species_all_p3 <- data.frame(
  Taxa = species_all_p2$Taxa,
  group = species_all_p2$group,
  species = species_all_p2$species,
  go = species_all_p2$go,
  mtb = rowMeans(species_all_p2[, c("mtb_pos", "mtb_neg")], na.rm = T),
  plasma_mtb = rowMeans(species_all_p2[, c("plasma_mtb_pos", "plasma_mtb_neg")], na.rm = T)
)

species_all_p3[is.na(species_all_p3)] <- 0


# Convert to wide format
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

# Convert to long format
species_p_long <- species_all_p3 %>%
  mutate(group = ifelse(group == "cor1", "NC", "SC")) %>%
  filter(Taxa %in% species_p_wide$Taxa)


write.csv(species_p_long, file.path(score_dir, paste0("scores_long_p_", i, ".csv")), row.names = F, quote = F)
write.csv(species_p_wide, file.path(score_dir, paste0("scores_wide_p_", i, ".csv")), row.names = F, quote = F)
write.csv(species_p_wide_diff, file.path(score_dir, paste0("scores_wide_diff_p_", i, ".csv")), row.names = T, quote = F)

