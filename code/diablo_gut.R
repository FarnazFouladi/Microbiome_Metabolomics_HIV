library(tidyverse)
library(readr)
library(mixOmics)
library(phyloseq)

# Parameters
microbiome_sample_type <- "Gut"
microbiome_sample_id <- "nci_metagenomic_stool_id"
nrepeat <- 10
nfold <- 10
ncore <- 6
output_dir <- file.path("Outputs/Diablo/Gut")
dir.create(output_dir, recursive = T)

# Load metadata
meta <- read_csv("data/metadata/metadata_c.csv")
meta$plasma_sample_id <- as.character(meta$plasma_sample_id)
meta$gut_sample_id <- as.character(meta$gut_sample_id)

# Species, GO, and metabolites table, bias-corrected abundances obtained from ANCOM-BC2
df_taxa <- read.csv(file.path(
  "Outputs/Differential_Abundance_Microbiome", microbiome_sample_type,
  "Species", paste0("Species", "_ANCOMBC_bias_corrected_data.txt")
), sep = "\t", check.names = F, quote = "", comment.char = "")
df_go <- read.csv(file.path(
  "Outputs/Differential_Abundance_Microbiome", microbiome_sample_type,
  "GO", paste0("GO", "_ANCOMBC_bias_corrected_data.txt")
), sep = "\t", check.names = F, quote = "", comment.char = "")
df_met_plasma_neg <- read.csv(file.path(
  "Outputs/Differential_Abundance_Metabolomics", "Plasma", "negative",
  "metabolites", paste0("negative_metabolites", "_ANCOMBC_bias_corrected_data.txt")
), sep = "\t", check.names = F, quote = "", comment.char = "")
df_met_plasma_pos <- read.csv(file.path(
  "Outputs/Differential_Abundance_Metabolomics", "Plasma", "positive",
  "metabolites", paste0("positive_metabolites", "_ANCOMBC_bias_corrected_data.txt")
), sep = "\t", check.names = F, quote = "", comment.char = "")
df_met_gut_neg <- read.csv(file.path(
  "Outputs/Differential_Abundance_Metabolomics", "gut", "negative",
  "metabolites", paste0("negative_metabolites", "_ANCOMBC_bias_corrected_data.txt")
), sep = "\t", check.names = F, quote = "", comment.char = "")
df_met_gut_pos <- read.csv(file.path(
  "Outputs/Differential_Abundance_Metabolomics", "gut", "positive",
  "metabolites", paste0("positive_metabolites", "_ANCOMBC_bias_corrected_data.txt")
), sep = "\t", check.names = F, quote = "", comment.char = "")


# Transpose and change sample IDs to the common names (subjid) across all datasets
df_taxa1 <- t(df_taxa) %>%
  as.data.frame() %>%
  rownames_to_column(microbiome_sample_id) %>%
  left_join(meta) %>%
  column_to_rownames("subjid") %>%
  dplyr::select(rownames(df_taxa)) %>%
  t() %>%
  as.data.frame()


df_go1 <- t(df_go) %>%
  as.data.frame() %>%
  rownames_to_column(microbiome_sample_id) %>%
  left_join(meta) %>%
  column_to_rownames("subjid") %>%
  dplyr::select(rownames(df_go)) %>%
  t() %>%
  as.data.frame()

df_met_plasma_neg1 <- t(df_met_plasma_neg) %>%
  as.data.frame() %>%
  rownames_to_column("plasma_sample_id") %>%
  left_join(meta) %>%
  column_to_rownames("subjid") %>%
  dplyr::select(rownames(df_met_plasma_neg)) %>%
  t() %>%
  as.data.frame()

df_met_plasma_pos1 <- t(df_met_plasma_pos) %>%
  as.data.frame() %>%
  rownames_to_column("plasma_sample_id") %>%
  left_join(meta) %>%
  column_to_rownames("subjid") %>%
  dplyr::select(rownames(df_met_plasma_pos)) %>%
  t() %>%
  as.data.frame()


df_met_gut_neg1 <- t(df_met_gut_neg) %>%
  as.data.frame() %>%
  rownames_to_column("gut_sample_id") %>%
  left_join(meta) %>%
  column_to_rownames("subjid") %>%
  dplyr::select(rownames(df_met_gut_neg)) %>%
  t() %>%
  as.data.frame()

df_met_gut_pos1 <- t(df_met_gut_pos) %>%
  as.data.frame() %>%
  rownames_to_column("gut_sample_id") %>%
  left_join(meta) %>%
  column_to_rownames("subjid") %>%
  dplyr::select(rownames(df_met_gut_pos)) %>%
  t() %>%
  as.data.frame()

# Subset to common samples across all datasets
common_samples <- Reduce(intersect, list(colnames(df_taxa1), colnames(df_met_plasma_neg1), colnames(df_met_gut_neg1)))
df_met_plasma_neg1 <- df_met_plasma_neg1[, common_samples]
df_met_plasma_pos1 <- df_met_plasma_pos1[, common_samples]
df_met_gut_neg1 <- df_met_gut_neg1[, common_samples]
df_met_gut_pos1 <- df_met_gut_pos1[, common_samples]
df_taxa1 <- df_taxa1[, common_samples]
df_go1 <- df_go1[, common_samples]

# Replace GO terms with their annotations
phy_obj <- readRDS(file.path("data/processed_data", "gut", "GO_phy.rds"))
ann_tmp <- tax_table(phy_obj) %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  dplyr::select(c("Feature", "Description"))
df_go1 <- df_go1 %>%
  rownames_to_column("Feature") %>%
  left_join(ann_tmp) %>%
  mutate(Feature = paste0(Feature, ":", Description)) %>%
  column_to_rownames("Feature") %>%
  dplyr::select(-c("Description"))

meta_sub <- meta %>%
  filter(subjid %in% common_samples) %>%
  arrange(order(match(common_samples, subjid)))
stopifnot(all.equal(meta_sub$subjid, common_samples))

# Replace metabolite names with their compound names
mtb_ann <- read.table(file.path("data/metabolomics_annotation/metabolites_annotation.txt"), sep = "\t", comment.char = "", quote = "", check.names = F, header = T)
mtb_ann <- mtb_ann[, c(1, 2)]
df_met_plasma_neg2 <- df_met_plasma_neg1 %>%
  rownames_to_column("accession") %>%
  left_join(mtb_ann) %>%
  dplyr::select(-c("accession")) %>%
  column_to_rownames("Compound_name")
df_met_plasma_pos2 <- df_met_plasma_pos1 %>%
  rownames_to_column("accession") %>%
  left_join(mtb_ann) %>%
  dplyr::select(-c("accession")) %>%
  column_to_rownames("Compound_name")
df_met_gut_neg2 <- df_met_gut_neg1 %>%
  rownames_to_column("accession") %>%
  left_join(mtb_ann) %>%
  dplyr::select(-c("accession")) %>%
  column_to_rownames("Compound_name")
df_met_gut_pos2 <- df_met_gut_pos1 %>%
  rownames_to_column("accession") %>%
  left_join(mtb_ann) %>%
  dplyr::select(-c("accession")) %>%
  column_to_rownames("Compound_name")

# Get residuals, control for age, antibiotics, and location
get_residuals <- function(df) {
  df <- as.data.frame(t(df))
  stopifnot(all.equal(rownames(df), meta_sub$subjid))
  l_tmp <- lapply(seq_len(ncol(df)), function(i) {
    df_tmp <- data.frame(y = df[, i], meta_sub)
    df_tmp$loc <- factor(df_tmp$loc)
    fit <- stats::lm(as.formula(paste0("y ~", paste(c("age", "loc", "abx_use"), collapse = "+"))), data = df_tmp)
    er <- resid(fit)
    return(er)
  })

  res_df <- bind_rows(l_tmp) %>% as.data.frame()
  rownames(res_df) <- colnames(df)
  colnames(res_df) <- rownames(df)
  return(res_df)
}
df_taxa_res <- get_residuals(df = df_taxa1)
df_go_res <- get_residuals(df = df_go1)
df_met_gut_neg_res <- get_residuals(df = df_met_gut_neg2)
df_met_gut_pos_res <- get_residuals(df = df_met_gut_pos2)
df_met_plasma_neg_res <- get_residuals(df = df_met_plasma_neg2)
df_met_plasma_pos_res <- get_residuals(df = df_met_plasma_pos2)

# Remove unclassified species
df_taxa_res <- df_taxa_res[!grepl("s__sp.", rownames(df_taxa_res)), ]

# Create a list of datasets
X <- list(
  Gut_Species = as.data.frame(t(df_taxa_res)),
  GO = as.data.frame(t(df_go_res)),
  MTB_Gut_Neg = as.data.frame(t(df_met_gut_neg_res)),
  MTB_Gut_Pos = as.data.frame(t(df_met_gut_pos_res)),
  MTB_Plasma_Neg = as.data.frame(t(df_met_plasma_neg_res)),
  MTB_Plasma_Pos = as.data.frame(t(df_met_plasma_pos_res))
)

# Outcome
Y <- meta_sub$status
summary(Y)

# design matrix
design <- matrix(0.1,
  ncol = length(X), nrow = length(X),
  dimnames = list(names(X), names(X))
)
diag(design) <- 0
design

diablo.obj <- block.plsda(X, Y, ncomp = 5, design = design)

set.seed(123)
perf.diablo.obj <- perf(diablo.obj, validation = "Mfold", folds = 10, nrepeat = 10)

# Plot of the error rates based on weighted vote, Choose the number of components based on minimum error
plot(perf.diablo.obj)
perf.diablo.obj$choice.ncomp$WeightedVote
ncomp <- perf.diablo.obj$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

set.seed(123) # for reproducibility

# Number of features to keep for each dataset
test.keepX <- list(
  Gut_Species = c(5, 10, 15),
  GO = c(5, 10, 15),
  MTB_Gut_Neg = c(3, 5, 10, 15),
  MTB_Gut_Pos = c(3, 5, 10, 15),
  MTB_Plasma_Neg = c(3, 5, 10, 15),
  MTB_Plasma_Pos = c(3, 5, 10, 15)
)

start.time <- Sys.time()
tune.diablo.obj <- tune.block.splsda(X, Y,
  ncomp = ncomp,
  test.keepX = test.keepX, design = design,
  validation = "Mfold", folds = nfold, nrepeat = nrepeat,
  # BPPARAM = BiocParallel::SnowParam(workers = ncore),
  BPPARAM = BiocParallel::MulticoreParam(workers = ncore),
  dist = "centroids.dist"
)
end.time <- Sys.time()
end.time - start.time
# Time difference of 18.97746 hours
saveRDS(tune.diablo.obj, file.path(output_dir, "tune_diablo_obj.rds"))

list.keepX <- tune.diablo.obj$choice.keepX

diablo.obj <- block.splsda(X, Y,
  ncomp = ncomp,
  keepX = list.keepX, design = design
)

saveRDS(diablo.obj, file.path(output_dir, "diablo_obj.rds"))


diablo.obj <- readRDS(file.path(output_dir, "diablo_obj.rds"))

############### Plotting

status_cols <- c("#0072b5ff", "#bc3c29ff")
block_cols <- c("#A6CEE3", "#1F78B4", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")
cor_cols <- c("green", "red")

pdf(file.path(output_dir, "circosplot_com12_0.8.pdf"), 20, 20)

c_plot <- circosPlot(diablo.obj,
  comp = 1:2,
  cutoff = 0.8,
  line = TRUE,
  size.variables = 1,
  var.adj = -1, # inner
  block.labels.adj = -2,
  color.blocks = block_cols,
  color.cor = cor_cols,
  size.labels = 1
)


dev.off()

# Extract correlation coefficients
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  diag(cormat) <- NA
  return(cormat)
}

df1 <- get_upper_tri(c_plot)
df2 <- df1 %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  pivot_longer(cols = -Feature) %>%
  filter(!is.na(value)) %>%
  mutate(d1 = ifelse(grepl("MTB_Plasma", Feature), "MTB_Plasma",
    ifelse(grepl("MTB_Gut", Feature), "MTB_Gut",
      ifelse(grepl("s__", Feature), "Gut_Species", "Gut_GO")
    )
  )) %>%
  mutate(d2 = ifelse(grepl("MTB_Plasma", name), "MTB_Plasma",
    ifelse(grepl("MTB_Gut", name), "MTB_Gut",
      ifelse(grepl("s__", name), "Gut_Species", "Gut_GO")
    )
  )) %>%
  filter(d1 != d2) %>%
  arrange(desc(abs(value))) %>%
  dplyr::rename(coefficient = value) %>%
  mutate(Feature = gsub("_Gut|_MTB_Plasma|_MTB_Gut", "", Feature)) %>%
  mutate(name = gsub("_Gut|_MTB_Plasma|_MTB_Gut", "", name)) %>%
  rename(feature1 = Feature, feature2 = name, data_mdality1 = d1, data_mdality2 = d2 )

xlsx::write.xlsx(df2, file = file.path( "Supplementary_Tables", "Supplementary_Table2.xlsx"))
