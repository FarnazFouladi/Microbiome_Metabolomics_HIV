rm(list = ls())
sample_type <- "Oral"
source("code/load.R")
source("code/network_function.R")
threshold <- 0.3
pval_t <- 0.01
res_dir <- file.path("Dysbiosis/Oral")


############################################# Oral Species
name <- "Species"
sub_dir <- file.path(output_dir, "Dysbiosis/Oral", "Oral_Species", "sub-network")
dir.create(sub_dir)
res_g1 <- readRDS(file.path(output_dir, res_dir, "Oral_Species", "species_g1.rds"))
res_g2 <- readRDS(file.path(output_dir, res_dir, "Oral_Species", "species_g2.rds"))
cor_diff <- read.csv(file.path(output_dir, res_dir, "Oral_Species", "species_differential_correlations.csv"))
sig_diff <- cor_diff[cor_diff$diff_cor > threshold & cor_diff$p < pval_t, ]
sig_diff_hold <- sig_diff %>% filter(grepl("Rothia|TM7", Taxa) | grepl("Rothia|TM7", name))
species <- unique(c(sig_diff_hold$Taxa, sig_diff_hold$name))
g1_mat <- res_g1$corr
g1_mat_cooccur <- res_g1$mat_cooccur
g2_mat <- res_g2$corr
g2_mat_cooccur <- res_g2$mat_cooccur

g1_mat_sub <- g1_mat[rownames(g1_mat) %in% species, colnames(g1_mat) %in% species]
g2_mat_sub <- g2_mat[rownames(g2_mat) %in% species, colnames(g2_mat) %in% species]
g1_mat_cooccur_sub <- g1_mat_cooccur[rownames(g1_mat_cooccur) %in% species, colnames(g1_mat_cooccur) %in% species]
g2_mat_cooccur_sub <- g2_mat_cooccur[rownames(g2_mat_cooccur) %in% species, colnames(g2_mat_cooccur) %in% species]
phylum_levels <- c(
  "p__Firmicutes", "p__Bacteroidota", "p__Actinobacteria", "p__Proteobacteria", "p__Fusobacteria",
  "p__Uroviricota", "p__Candidatus_Saccharibacteria", "p__Unclassified", "p__Missing"
)
names(phylum_levels) <- mycols[1:length(phylum_levels)]
taxonomies <- read_rds("data/processed_data/oral/taxonomies.rds")

create_networks(
  g1_mat = g1_mat_sub, g1_mat_cooccur = g1_mat_cooccur_sub,
  g2_mat = g2_mat_sub, g2_mat_cooccur = g2_mat_cooccur_sub,
  sig_diff = sig_diff_hold,
  # taxonomies = taxonomies,
  datasets = c("Oral Species", "Oral Species"),
  t = "Species", phylum_levels = phylum_levels,
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "union",
  repulsion = 0.85, repulsion2 = 0.7,
  cexNodes = 1.5, cexNodes2 = 0.8,
  cexHubs = 2, cexHubs2 = 1.2,
  cexLabels = 4, cexLabels2 = 5,
  cexHubLabels = 5, cexHubLabels2 = 7,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = sub_dir, prefix = name, showLabel = T
)

############################################# Oral Species and Oral GO
name <- "GO"
sub_dir <- file.path(output_dir, "Dysbiosis/Oral", "Oral_GO", "sub-network")
res_g1 <- readRDS(file.path(output_dir, res_dir, "Oral_GO", "go_g1.rds"))
res_g2 <- readRDS(file.path(output_dir, res_dir, "Oral_GO", "go_g2.rds"))
sig_diff <- read.table(file.path(output_dir, res_dir, "Oral_GO", "go_differential_correlations_sig.txt"),
  sep = "\t", quote = "", check.names = FALSE, comment.char = "", header = T
)
sig_diff <- sig_diff %>%
  group_by(Description) %>%
  mutate(num_description = n()) %>%
  ungroup() %>%
  group_by(Taxa) %>%
  mutate(num_taxa = n())

sig_diff_sub <- sig_diff %>%
  filter(grepl("Haemophilus_parainfluenzae", Taxa)) %>%
  mutate(
    Taxa = paste0("data1 - ", Taxa),
    GO = paste0("data2 - ", GO),
    Description = paste0("data2 - ", Description)
  )
features <- unique(c(sig_diff_sub$Taxa, sig_diff_sub$GO))
g1_mat <- res_g1$corr
g1_mat_cooccur <- res_g1$mat_cooccur
g2_mat <- res_g2$corr
g2_mat_cooccur <- res_g2$mat_cooccur

g1_mat_sub <- g1_mat[rownames(g1_mat) %in% features, colnames(g1_mat) %in% features]
g2_mat_sub <- g2_mat[rownames(g2_mat) %in% features, colnames(g2_mat) %in% features]
g1_mat_cooccur_sub <- g1_mat_cooccur[rownames(g1_mat_cooccur) %in% features, colnames(g1_mat_cooccur) %in% features]
g2_mat_cooccur_sub <- g2_mat_cooccur[rownames(g2_mat_cooccur) %in% features, colnames(g2_mat_cooccur) %in% features]

names_mat <- data.frame(GO = colnames(g1_mat_sub)) %>%
  left_join(sig_diff_sub[, c("GO", "Description")] %>% distinct(GO, .keep_all = T)) %>%
  mutate(Description = ifelse(is.na(Description), GO, Description)) %>%
  arrange(order(match(colnames(g1_mat_sub), GO)))
stopifnot(all.equal(colnames(g1_mat_sub), names_mat$GO))
stopifnot(all.equal(colnames(g2_mat_sub), names_mat$GO))
stopifnot(all.equal(colnames(g1_mat_cooccur_sub), names_mat$GO))
stopifnot(all.equal(colnames(g2_mat_cooccur_sub), names_mat$GO))


colnames(g1_mat_sub) <- names_mat$Description
rownames(g1_mat_sub) <- names_mat$Description
colnames(g2_mat_sub) <- names_mat$Description
rownames(g2_mat_sub) <- names_mat$Description
colnames(g1_mat_cooccur_sub) <- names_mat$Description
rownames(g1_mat_cooccur_sub) <- names_mat$Description
colnames(g2_mat_cooccur_sub) <- names_mat$Description
rownames(g2_mat_cooccur_sub) <- names_mat$Description
sig_diff_sub1 <- data.frame(Taxa = sig_diff_sub$Taxa, name = sig_diff_sub$Description)

create_networks(
  g1_mat = g1_mat_sub, g1_mat_cooccur = g1_mat_cooccur_sub,
  g2_mat = g2_mat_sub, g2_mat_cooccur = g2_mat_cooccur_sub,
  sig_diff = sig_diff_sub1,
  taxonomies = NULL,
  t = NULL, datasets = c("Oral Species", "Oral GO"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "2",
  repulsion = 0.85, repulsion2 = 0.85,
  cexNodes = 1, cexNodes2 = 1,
  cexHubs = 1.5, cexHubs2 = 1.5,
  cexLabels = 1, cexLabels2 = 5,
  cexHubLabels = 2, cexHubLabels2 = 9,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = sub_dir, prefix = "GO", showLabel = T
)

############################################# Oral Species and Oral Metabolites
name <- "Oral_Metabolites"
sub_dir <- file.path(output_dir, "Dysbiosis/Oral", "Oral_MTB_NEG", "sub-network")
res_g1 <- readRDS(file.path(output_dir, res_dir, "Oral_MTB_NEG", "mtb_neg_g1.rds"))
res_g2 <- readRDS(file.path(output_dir, res_dir, "Oral_MTB_NEG", "mtb_neg_g2.rds"))
sig_diff <- read.table(file.path(output_dir, res_dir, "Oral_MTB_NEG", "mtb_neg_differential_correlations_sig.txt"),
  sep = "\t", quote = "", check.names = FALSE, comment.char = "", header = T
)
sig_diff <- sig_diff %>%
  group_by(Metabolites) %>%
  mutate(num_description = n()) %>%
  ungroup() %>%
  group_by(Taxa) %>%
  mutate(num_taxa = n())

#  LysoPC
sig_diff_sub <- sig_diff %>%
  filter(grepl("LysoPC", Metabolites)) %>%
  mutate(
    Taxa = paste0("data1 - ", Taxa),
    metabolite_ID = paste0("data2 - ", metabolite_ID),
    Metabolites = paste0("data2 - ", Metabolites)
  )
features <- unique(c(sig_diff_sub$Taxa, sig_diff_sub$metabolite_ID))
g1_mat <- res_g1$corr
g1_mat_cooccur <- res_g1$mat_cooccur
g2_mat <- res_g2$corr
g2_mat_cooccur <- res_g2$mat_cooccur

g1_mat_sub <- g1_mat[rownames(g1_mat) %in% features, colnames(g1_mat) %in% features]
g2_mat_sub <- g2_mat[rownames(g2_mat) %in% features, colnames(g2_mat) %in% features]
g1_mat_cooccur_sub <- g1_mat_cooccur[rownames(g1_mat_cooccur) %in% features, colnames(g1_mat_cooccur) %in% features]
g2_mat_cooccur_sub <- g2_mat_cooccur[rownames(g2_mat_cooccur) %in% features, colnames(g2_mat_cooccur) %in% features]

names_mat <- data.frame(metabolite_ID = colnames(g1_mat_sub)) %>%
  left_join(sig_diff_sub[, c("metabolite_ID", "Metabolites")] %>% distinct(metabolite_ID, .keep_all = T)) %>%
  mutate(Metabolites = ifelse(is.na(Metabolites), metabolite_ID, Metabolites)) %>%
  arrange(order(match(colnames(g1_mat_sub), metabolite_ID)))
stopifnot(all.equal(colnames(g1_mat_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g1_mat_cooccur_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_cooccur_sub), names_mat$metabolite_ID))


colnames(g1_mat_sub) <- names_mat$Metabolites
rownames(g1_mat_sub) <- names_mat$Metabolites
colnames(g2_mat_sub) <- names_mat$Metabolites
rownames(g2_mat_sub) <- names_mat$Metabolites
colnames(g1_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g1_mat_cooccur_sub) <- names_mat$Metabolites
colnames(g2_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g2_mat_cooccur_sub) <- names_mat$Metabolites

sig_diff_sub1 <- data.frame(Taxa = sig_diff_sub$Taxa, name = sig_diff_sub$Metabolites)

create_networks(
  g1_mat = g1_mat_sub, g1_mat_cooccur = g1_mat_cooccur_sub,
  g2_mat = g2_mat_sub, g2_mat_cooccur = g2_mat_cooccur_sub,
  sig_diff = sig_diff_sub1,
  taxonomies = NULL,
  t = NULL, datasets = c("Oral Species", "Oral MTB-POS"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "2",
  repulsion = 0.85, repulsion2 = 0.7,
  cexNodes = 1, cexNodes2 = 1,
  cexHubs = 1.5, cexHubs2 = 1.5,
  cexLabels = 1, cexLabels2 = 5,
  cexHubLabels = 2, cexHubLabels2 = 3,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = sub_dir, prefix = " LysoPC", showLabel = T
)

############################################# Oral Species and Oral Metabolites
name <- "Oral_Metabolites"
sub_dir <- file.path(output_dir, "Dysbiosis/Oral", "Oral_MTB_POS", "sub-network")
res_g1 <- readRDS(file.path(output_dir, res_dir, "Oral_MTB_POS", "mtb_pos_g1.rds"))
res_g2 <- readRDS(file.path(output_dir, res_dir, "Oral_MTB_POS", "mtb_pos_g2.rds"))
sig_diff <- read.table(file.path(output_dir, res_dir, "Oral_MTB_POS", "mtb_pos_differential_correlations_sig.txt"),
  sep = "\t", quote = "", check.names = FALSE, comment.char = "", header = T
)
sig_diff <- sig_diff %>%
  group_by(Metabolites) %>%
  mutate(num_description = n()) %>%
  ungroup() %>%
  group_by(Taxa) %>%
  mutate(num_taxa = n())

#  Cotinine
sig_diff_sub <- sig_diff %>%
  filter(grepl("Cotinine", Metabolites)) %>%
  mutate(
    Taxa = paste0("data1 - ", Taxa),
    metabolite_ID = paste0("data2 - ", metabolite_ID),
    Metabolites = paste0("data2 - ", Metabolites)
  )
features <- unique(c(sig_diff_sub$Taxa, sig_diff_sub$metabolite_ID))
g1_mat <- res_g1$corr
g1_mat_cooccur <- res_g1$mat_cooccur
g2_mat <- res_g2$corr
g2_mat_cooccur <- res_g2$mat_cooccur

g1_mat_sub <- g1_mat[rownames(g1_mat) %in% features, colnames(g1_mat) %in% features]
g2_mat_sub <- g2_mat[rownames(g2_mat) %in% features, colnames(g2_mat) %in% features]
g1_mat_cooccur_sub <- g1_mat_cooccur[rownames(g1_mat_cooccur) %in% features, colnames(g1_mat_cooccur) %in% features]
g2_mat_cooccur_sub <- g2_mat_cooccur[rownames(g2_mat_cooccur) %in% features, colnames(g2_mat_cooccur) %in% features]

names_mat <- data.frame(metabolite_ID = colnames(g1_mat_sub)) %>%
  left_join(sig_diff_sub[, c("metabolite_ID", "Metabolites")] %>% distinct(metabolite_ID, .keep_all = T)) %>%
  mutate(Metabolites = ifelse(is.na(Metabolites), metabolite_ID, Metabolites)) %>%
  arrange(order(match(colnames(g1_mat_sub), metabolite_ID)))
stopifnot(all.equal(colnames(g1_mat_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g1_mat_cooccur_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_cooccur_sub), names_mat$metabolite_ID))


colnames(g1_mat_sub) <- names_mat$Metabolites
rownames(g1_mat_sub) <- names_mat$Metabolites
colnames(g2_mat_sub) <- names_mat$Metabolites
rownames(g2_mat_sub) <- names_mat$Metabolites
colnames(g1_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g1_mat_cooccur_sub) <- names_mat$Metabolites
colnames(g2_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g2_mat_cooccur_sub) <- names_mat$Metabolites

sig_diff_sub1 <- data.frame(Taxa = sig_diff_sub$Taxa, name = sig_diff_sub$Metabolites)

create_networks(
  g1_mat = g1_mat_sub, g1_mat_cooccur = g1_mat_cooccur_sub,
  g2_mat = g2_mat_sub, g2_mat_cooccur = g2_mat_cooccur_sub,
  sig_diff = sig_diff_sub1,
  taxonomies = NULL,
  t = NULL, datasets = c("Oral Species", "Oral MTB-POS"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "2",
  repulsion = 0.85, repulsion2 = 0.5,
  cexNodes = 1, cexNodes2 = 0.5,
  cexHubs = 1.5, cexHubs2 = 1,
  cexLabels = 1, cexLabels2 = 5,
  cexHubLabels = 2, cexHubLabels2 = 3,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("NC", "SC"), outputdir = sub_dir, prefix = " Cotinine", showLabel = T
)

############################################# Oral Species and Plasma Metabolites
name <- "Plasma_Metabolites"
sub_dir <- file.path(output_dir, "Dysbiosis/Oral", "Plasma_MTB_NEG", "sub-network")
res_g1 <- readRDS(file.path(output_dir, res_dir, "Plasma_MTB_NEG", "plasma_mtb_neg_g1.rds"))
res_g2 <- readRDS(file.path(output_dir, res_dir, "Plasma_MTB_NEG", "plasma_mtb_neg_g2.rds"))

sig_diff <- read.table(file.path(output_dir, res_dir, "Plasma_MTB_NEG", "plasma_mtb_neg_differential_correlations_sig.txt"),
  sep = "\t", quote = "", check.names = FALSE, comment.char = "", header = T
)
sig_diff <- sig_diff %>%
  group_by(Metabolites) %>%
  mutate(num_description = n()) %>%
  ungroup() %>%
  group_by(Taxa) %>%
  mutate(num_taxa = n())

# 4-Ethylphenylsulfate
sig_diff_sub <- sig_diff %>%
  filter(grepl("4-Ethylphenylsulfate", Metabolites)) %>%
  mutate(
    Taxa = paste0("data1 - ", Taxa),
    metabolite_ID = paste0("data2 - ", metabolite_ID),
    Metabolites = paste0("data2 - ", Metabolites)
  )
features <- unique(c(sig_diff_sub$Taxa, sig_diff_sub$metabolite_ID))
g1_mat <- res_g1$corr
g1_mat_cooccur <- res_g1$mat_cooccur
g2_mat <- res_g2$corr
g2_mat_cooccur <- res_g2$mat_cooccur

g1_mat_sub <- g1_mat[rownames(g1_mat) %in% features, colnames(g1_mat) %in% features]
g2_mat_sub <- g2_mat[rownames(g2_mat) %in% features, colnames(g2_mat) %in% features]
g1_mat_cooccur_sub <- g1_mat_cooccur[rownames(g1_mat_cooccur) %in% features, colnames(g1_mat_cooccur) %in% features]
g2_mat_cooccur_sub <- g2_mat_cooccur[rownames(g2_mat_cooccur) %in% features, colnames(g2_mat_cooccur) %in% features]

names_mat <- data.frame(metabolite_ID = colnames(g1_mat_sub)) %>%
  left_join(sig_diff_sub[, c("metabolite_ID", "Metabolites")] %>% distinct(metabolite_ID, .keep_all = T)) %>%
  mutate(Metabolites = ifelse(is.na(Metabolites), metabolite_ID, Metabolites)) %>%
  arrange(order(match(colnames(g1_mat_sub), metabolite_ID)))
stopifnot(all.equal(colnames(g1_mat_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g1_mat_cooccur_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_cooccur_sub), names_mat$metabolite_ID))


colnames(g1_mat_sub) <- names_mat$Metabolites
rownames(g1_mat_sub) <- names_mat$Metabolites
colnames(g2_mat_sub) <- names_mat$Metabolites
rownames(g2_mat_sub) <- names_mat$Metabolites
colnames(g1_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g1_mat_cooccur_sub) <- names_mat$Metabolites
colnames(g2_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g2_mat_cooccur_sub) <- names_mat$Metabolites
sig_diff_sub1 <- data.frame(Taxa = sig_diff_sub$Taxa, name = sig_diff_sub$Metabolites)

create_networks(
  g1_mat = g1_mat_sub, g1_mat_cooccur = g1_mat_cooccur_sub,
  g2_mat = g2_mat_sub, g2_mat_cooccur = g2_mat_cooccur_sub,
  sig_diff = sig_diff_sub1,
  taxonomies = NULL,
  t = NULL, datasets = c("Oral Species", "Plasma MTB-NEG"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "2",
  repulsion = 0.85, repulsion2 = 0.6,
  cexNodes = 1, cexNodes2 = 1,
  cexHubs = 1.5, cexHubs2 = 1.5,
  cexLabels = 1, cexLabels2 = 4,
  cexHubLabels = 2, cexHubLabels2 = 3,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = sub_dir, prefix = "4-Ethylphenylsulfate", showLabel = T
)
#########################
# Positive
name <- "Plasma_Metabolites"
sub_dir <- file.path(output_dir, "Dysbiosis/Oral", "Plasma_MTB_POS", "sub-network")
res_g1 <- readRDS(file.path(output_dir, res_dir, "Plasma_MTB_POS", "plasma_mtb_pos_g1.rds"))
res_g2 <- readRDS(file.path(output_dir, res_dir, "Plasma_MTB_POS", "plasma_mtb_pos_g2.rds"))
sig_diff <- read.table(file.path(output_dir, res_dir, "Plasma_MTB_POS", "plasma_mtb_pos_differential_correlations_sig.txt"),
  sep = "\t", quote = "", check.names = FALSE, comment.char = "", header = T
)
sig_diff <- sig_diff %>%
  group_by(Metabolites) %>%
  mutate(num_description = n()) %>%
  ungroup() %>%
  group_by(Taxa) %>%
  mutate(num_taxa = n())

# Cotinine
sig_diff_sub <- sig_diff %>%
  filter(grepl("Cotinine", Metabolites)) %>%
  mutate(
    Taxa = paste0("data1 - ", Taxa),
    metabolite_ID = paste0("data2 - ", metabolite_ID),
    Metabolites = paste0("data2 - ", Metabolites)
  )
features <- unique(c(sig_diff_sub$Taxa, sig_diff_sub$metabolite_ID))
g1_mat <- res_g1$corr
g1_mat_cooccur <- res_g1$mat_cooccur
g2_mat <- res_g2$corr
g2_mat_cooccur <- res_g2$mat_cooccur

g1_mat_sub <- g1_mat[rownames(g1_mat) %in% features, colnames(g1_mat) %in% features]
g2_mat_sub <- g2_mat[rownames(g2_mat) %in% features, colnames(g2_mat) %in% features]
g1_mat_cooccur_sub <- g1_mat_cooccur[rownames(g1_mat_cooccur) %in% features, colnames(g1_mat_cooccur) %in% features]
g2_mat_cooccur_sub <- g2_mat_cooccur[rownames(g2_mat_cooccur) %in% features, colnames(g2_mat_cooccur) %in% features]

names_mat <- data.frame(metabolite_ID = colnames(g1_mat_sub)) %>%
  left_join(sig_diff_sub[, c("metabolite_ID", "Metabolites")] %>% distinct(metabolite_ID, .keep_all = T)) %>%
  mutate(Metabolites = ifelse(is.na(Metabolites), metabolite_ID, Metabolites)) %>%
  arrange(order(match(colnames(g1_mat_sub), metabolite_ID)))
stopifnot(all.equal(colnames(g1_mat_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g1_mat_cooccur_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_cooccur_sub), names_mat$metabolite_ID))


colnames(g1_mat_sub) <- names_mat$Metabolites
rownames(g1_mat_sub) <- names_mat$Metabolites
colnames(g2_mat_sub) <- names_mat$Metabolites
rownames(g2_mat_sub) <- names_mat$Metabolites
colnames(g1_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g1_mat_cooccur_sub) <- names_mat$Metabolites
colnames(g2_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g2_mat_cooccur_sub) <- names_mat$Metabolites
sig_diff_sub1 <- data.frame(Taxa = sig_diff_sub$Taxa, name = sig_diff_sub$Metabolites)

create_networks(
  g1_mat = g1_mat_sub, g1_mat_cooccur = g1_mat_cooccur_sub,
  g2_mat = g2_mat_sub, g2_mat_cooccur = g2_mat_cooccur_sub,
  sig_diff = sig_diff_sub1,
  taxonomies = NULL,
  t = NULL, datasets = c("Oral Species", "Plasma MTB-NEG"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "2",
  repulsion = 0.85, repulsion2 = 0.6,
  cexNodes = 1, cexNodes2 = 1,
  cexHubs = 1.5, cexHubs2 = 1.5,
  cexLabels = 1, cexLabels2 = 4,
  cexHubLabels = 2, cexHubLabels2 = 2,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = sub_dir, prefix = "Cotinine", showLabel = T
)
# Indole
sig_diff_sub <- sig_diff %>%
  filter(grepl("Indole", Metabolites)) %>%
  mutate(
    Taxa = paste0("data1 - ", Taxa),
    metabolite_ID = paste0("data2 - ", metabolite_ID),
    Metabolites = paste0("data2 - ", Metabolites)
  )
features <- unique(c(sig_diff_sub$Taxa, sig_diff_sub$metabolite_ID))
g1_mat <- res_g1$corr
g1_mat_cooccur <- res_g1$mat_cooccur
g2_mat <- res_g2$corr
g2_mat_cooccur <- res_g2$mat_cooccur

g1_mat_sub <- g1_mat[rownames(g1_mat) %in% features, colnames(g1_mat) %in% features]
g2_mat_sub <- g2_mat[rownames(g2_mat) %in% features, colnames(g2_mat) %in% features]
g1_mat_cooccur_sub <- g1_mat_cooccur[rownames(g1_mat_cooccur) %in% features, colnames(g1_mat_cooccur) %in% features]
g2_mat_cooccur_sub <- g2_mat_cooccur[rownames(g2_mat_cooccur) %in% features, colnames(g2_mat_cooccur) %in% features]

names_mat <- data.frame(metabolite_ID = colnames(g1_mat_sub)) %>%
  left_join(sig_diff_sub[, c("metabolite_ID", "Metabolites")] %>% distinct(metabolite_ID, .keep_all = T)) %>%
  mutate(Metabolites = ifelse(is.na(Metabolites), metabolite_ID, Metabolites)) %>%
  arrange(order(match(colnames(g1_mat_sub), metabolite_ID)))
stopifnot(all.equal(colnames(g1_mat_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g1_mat_cooccur_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_cooccur_sub), names_mat$metabolite_ID))


colnames(g1_mat_sub) <- names_mat$Metabolites
rownames(g1_mat_sub) <- names_mat$Metabolites
colnames(g2_mat_sub) <- names_mat$Metabolites
rownames(g2_mat_sub) <- names_mat$Metabolites
colnames(g1_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g1_mat_cooccur_sub) <- names_mat$Metabolites
colnames(g2_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g2_mat_cooccur_sub) <- names_mat$Metabolites
sig_diff_sub1 <- data.frame(Taxa = sig_diff_sub$Taxa, name = sig_diff_sub$Metabolites)

create_networks(
  g1_mat = g1_mat_sub, g1_mat_cooccur = g1_mat_cooccur_sub,
  g2_mat = g2_mat_sub, g2_mat_cooccur = g2_mat_cooccur_sub,
  sig_diff = sig_diff_sub1,
  taxonomies = NULL,
  t = NULL, datasets = c("Oral Species", "Plasma MTB-NEG"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "2",
  repulsion = 0.85, repulsion2 = 0.6,
  cexNodes = 1, cexNodes2 = 1,
  cexHubs = 1.5, cexHubs2 = 1.5,
  cexLabels = 1, cexLabels2 = 4,
  cexHubLabels = 2, cexHubLabels2 = 4,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = sub_dir, prefix = "Indole", showLabel = T
)


# carnitine
sig_diff_sub <- sig_diff %>%
  filter(grepl("carnitine", Metabolites)) %>%
  mutate(
    Taxa = paste0("data1 - ", Taxa),
    metabolite_ID = paste0("data2 - ", metabolite_ID),
    Metabolites = paste0("data2 - ", Metabolites)
  )
features <- unique(c(sig_diff_sub$Taxa, sig_diff_sub$metabolite_ID))
g1_mat <- res_g1$corr
g1_mat_cooccur <- res_g1$mat_cooccur
g2_mat <- res_g2$corr
g2_mat_cooccur <- res_g2$mat_cooccur

g1_mat_sub <- g1_mat[rownames(g1_mat) %in% features, colnames(g1_mat) %in% features]
g2_mat_sub <- g2_mat[rownames(g2_mat) %in% features, colnames(g2_mat) %in% features]
g1_mat_cooccur_sub <- g1_mat_cooccur[rownames(g1_mat_cooccur) %in% features, colnames(g1_mat_cooccur) %in% features]
g2_mat_cooccur_sub <- g2_mat_cooccur[rownames(g2_mat_cooccur) %in% features, colnames(g2_mat_cooccur) %in% features]

names_mat <- data.frame(metabolite_ID = colnames(g1_mat_sub)) %>%
  left_join(sig_diff_sub[, c("metabolite_ID", "Metabolites")] %>% distinct(metabolite_ID, .keep_all = T)) %>%
  mutate(Metabolites = ifelse(is.na(Metabolites), metabolite_ID, Metabolites)) %>%
  arrange(order(match(colnames(g1_mat_sub), metabolite_ID)))
stopifnot(all.equal(colnames(g1_mat_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g1_mat_cooccur_sub), names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_cooccur_sub), names_mat$metabolite_ID))


colnames(g1_mat_sub) <- names_mat$Metabolites
rownames(g1_mat_sub) <- names_mat$Metabolites
colnames(g2_mat_sub) <- names_mat$Metabolites
rownames(g2_mat_sub) <- names_mat$Metabolites
colnames(g1_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g1_mat_cooccur_sub) <- names_mat$Metabolites
colnames(g2_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g2_mat_cooccur_sub) <- names_mat$Metabolites
sig_diff_sub1 <- data.frame(Taxa = sig_diff_sub$Taxa, name = sig_diff_sub$Metabolites)

create_networks(
  g1_mat = g1_mat_sub, g1_mat_cooccur = g1_mat_cooccur_sub,
  g2_mat = g2_mat_sub, g2_mat_cooccur = g2_mat_cooccur_sub,
  sig_diff = sig_diff_sub1,
  taxonomies = NULL,
  t = NULL, datasets = c("Oral Species", "Plasma MTB-NEG"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "2",
  repulsion = 0.85, repulsion2 = 0.6,
  cexNodes = 1, cexNodes2 = 1,
  cexHubs = 1.5, cexHubs2 = 1.5,
  cexLabels = 1, cexLabels2 = 3,
  cexHubLabels = 2, cexHubLabels2 = 4,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("NC", "SC"), outputdir = sub_dir, prefix = "carnitine", showLabel = T
)

############################ Oral and Gut Species interactions
name <- "Gut_Species"
sub_dir <- file.path(output_dir, "Dysbiosis/Oral", "Gut_Species", "sub-network")
res_g1 <- readRDS(file.path(output_dir, res_dir, "Gut_Species", "species2_g1.rds"))
res_g2 <- readRDS(file.path(output_dir, res_dir, "Gut_Species", "species2_g2.rds"))
sig_diff <- read.table(file.path(output_dir, res_dir, "Gut_Species", "species2_differential_correlations_sig.txt"),
  sep = "\t", quote = "", check.names = FALSE, comment.char = "", header = T
)
sig_diff <- sig_diff %>%
  group_by(name) %>%
  mutate(num_description = n()) %>%
  ungroup() %>%
  group_by(Taxa) %>%
  mutate(num_taxa = n())

# Taxa
sig_diff_sub <- sig_diff %>%
  filter(grepl("Ruminococcus|Bifidobacterium|Roseburia", name)) %>%
  mutate(
    Taxa = paste0("data1 - ", Taxa),
    name = paste0("data2 - ", name)
  )
features <- unique(c(sig_diff_sub$Taxa, sig_diff_sub$name))
g1_mat <- res_g1$corr
g1_mat_cooccur <- res_g1$mat_cooccur
g2_mat <- res_g2$corr
g2_mat_cooccur <- res_g2$mat_cooccur

g1_mat_sub <- g1_mat[rownames(g1_mat) %in% features, colnames(g1_mat) %in% features]
g2_mat_sub <- g2_mat[rownames(g2_mat) %in% features, colnames(g2_mat) %in% features]
g1_mat_cooccur_sub <- g1_mat_cooccur[rownames(g1_mat_cooccur) %in% features, colnames(g1_mat_cooccur) %in% features]
g2_mat_cooccur_sub <- g2_mat_cooccur[rownames(g2_mat_cooccur) %in% features, colnames(g2_mat_cooccur) %in% features]

create_networks(
  g1_mat = g1_mat_sub, g1_mat_cooccur = g1_mat_cooccur_sub,
  g2_mat = g2_mat_sub, g2_mat_cooccur = g2_mat_cooccur_sub,
  sig_diff = sig_diff_sub,
  taxonomies = NULL,
  t = NULL, datasets = c("Oral Species", "Gut Species"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "2",
  repulsion = 0.85, repulsion2 = 0.9,
  cexNodes = 1, cexNodes2 = 0.8,
  cexHubs = 1.5, cexHubs2 = 1.3,
  cexLabels = 1, cexLabels2 = 5,
  cexHubLabels = 2, cexHubLabels2 = 7,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = sub_dir, prefix = "Gut_Speccies", showLabel = T
)
