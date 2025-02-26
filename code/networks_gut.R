# Networks
rm(list = ls())
source("code/load.R")
sample_type <- "Gut"
output_dir <- file.path(output_dir, "Dysbiosis", sample_type)
threshold <- 0.3
pval_t <- 0.01
source("code/network_function.R")
mycols <- c("#D95F02", "#7570B3", "#E7298A", "#1B9E77", "#00FF00", "#E31A1C", "#E6AB02", "#00FFFF", "#FFFF99", "#6A3D9A", "#666666", "#A6761D", "blue", "black", "brown")
cor_names <- c("cor_g1", "cor_g2")
cor_names1 <- c(paste0("cor_", group_names[1]), paste0("cor_", group_names[2]))
cooccurrence_names <- c("cooccurrence_g1", "cooccurrence_g2")
cooccurrence_names1 <- c(paste0("cooccurrence_", group_names[1]), paste0("cooccurrence_", group_names[2]))

# Annotation files for GO terms and metabolites
phy_obj_go <- readRDS("data/processed_data/gut/GO_phy.rds")
ann_go <- microbiome::tax_tibble(phy_obj_go) %>%
  as.data.frame() %>%
  dplyr::rename(GO = FeatureID)
mtb_ann <- read.table(file.path("data/metabolomics_annotation/metabolites_annotation.txt"), sep = "\t", comment.char = "", quote = "", check.names = F, header = T)
mtb_ann <- mtb_ann %>%
  dplyr::rename(metabolite_ID = accession, Metabolites = Compound_name) %>%
  dplyr::select(-c("name"))

######################### Gut species
sub_dir <- file.path(output_dir, "Gut_Species")
name <- "species"
res_nc_gut_species <- readRDS(file.path(sub_dir, paste0(name, "_g1.rds")))
res_sc_gut_species <- readRDS(file.path(sub_dir, paste0(name, "_g2.rds")))
phylum_levels <- c(
  "p__Firmicutes", "p__Bacteroidota", "p__Actinobacteria", "p__Proteobacteria",
  "p__Verrucomicrobia", "p__Spirochaetes", "p__Euryarchaeota", "p__Uroviricota", "p__Apicomplexa", "p__Unclassified", "p__Missing"
)
names(phylum_levels) <- mycols[1:length(phylum_levels)]
taxonomies <- read_rds("data/processed_data/gut/taxonomies.rds")
res_gut_species <- read.csv(file.path(sub_dir, paste0(name, "_differential_correlations.csv")), check.names = F)
sig_diff_gut_species <- res_gut_species[res_gut_species$diff_cor > threshold & res_gut_species$p < pval_t, ]
sig_diff_save <- sig_diff_gut_species %>%
  mutate(across(grep("cor", colnames(.)), ~ round(.x, 3))) %>%
  mutate(sign = ifelse(grepl("g1", sign), gsub("g1", group_names[1], sign), gsub("g2", group_names[2], sign))) %>%
  dplyr::rename_with(~cor_names1, all_of(cor_names)) %>%
  dplyr::rename_with(~cooccurrence_names1, all_of(cooccurrence_names)) %>%
  arrange(p)
write.csv(sig_diff_save, file.path(sub_dir, paste0(name, "_differential_correlations_sig.csv")), row.names = F, quote = F)
create_networks(
  g1_mat = res_nc_gut_species$corr, g1_mat_cooccur = res_nc_gut_species$mat_cooccur,
  g2_mat = res_sc_gut_species$corr, g2_mat_cooccur = res_sc_gut_species$mat_cooccur,
  sig_diff = sig_diff_gut_species,
  # taxonomies = taxonomies,
  datasets = c("Gut Species", "Gut Species"),
  t = "Species", phylum_levels = phylum_levels,
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "union",
  repulsion = 0.85, repulsion2 = 0.7,
  cexNodes = 1.5, cexNodes2 = 1,
  cexHubs = 2, cexHubs2 = 1.5,
  cexLabels = 4, cexLabels2 = 3,
  cexHubLabels = 5, cexHubLabels2 = 5,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = file.path(sub_dir, "network"), prefix = name
)


######################### Gut species and Gut GO
sub_dir <- file.path(output_dir, "Gut_GO")
name <- "go"
res_nc_gut_go <- readRDS(file.path(sub_dir, paste0(name, "_g1.rds")))
res_sc_gut_go <- readRDS(file.path(sub_dir, paste0(name, "_g2.rds")))
res_gut_go <- read.csv(file.path(sub_dir, paste0(name, "_differential_correlations.csv")), check.names = F)
sig_diff_gut_go <- res_gut_go[res_gut_go$diff_cor > threshold & res_gut_go$p < pval_t, ]
sig_diff_save <- sig_diff_gut_go %>%
  mutate(
    Taxa = gsub("data1 - ", "", Taxa),
    name = gsub("data2 - ", "", name)
  ) %>%
  dplyr::rename(GO = name) %>%
  left_join(ann_go) %>%
  relocate(Description, .after = GO) %>%
  mutate(across(grep("cor", colnames(.)), ~ round(.x, 3))) %>%
  mutate(sign = ifelse(grepl("g1", sign), gsub("g1", group_names[1], sign), gsub("g2", group_names[2], sign))) %>%
  dplyr::rename_with(~cor_names1, all_of(cor_names)) %>%
  dplyr::rename_with(~cooccurrence_names1, all_of(cooccurrence_names)) %>%
  arrange(p)
write.csv(sig_diff_save, file.path(sub_dir, paste0(name, "_differential_correlations_sig.csv")), row.names = F, quote = F)
write.table(sig_diff_save, file.path(sub_dir, paste0(name, "_differential_correlations_sig.txt")), sep = "\t", row.names = F, quote = F)

create_networks(
  g1_mat = res_nc_gut_go$corr, g1_mat_cooccur = res_nc_gut_go$mat_cooccur,
  g2_mat = res_sc_gut_go$corr, g2_mat_cooccur = res_sc_gut_go$mat_cooccur,
  sig_diff = sig_diff_gut_go,
  taxonomies = NULL,
  t = NULL, datasets = c("Gut Species", "Gut GO"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "union",
  repulsion = 0.85, repulsion2 = 0.7,
  cexNodes = 1, cexNodes2 = 1,
  cexHubs = 1.5, cexHubs2 = 1,
  cexLabels = 1, cexLabels2 = 2,
  cexHubLabels = 2, cexHubLabels2 = 2,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = file.path(sub_dir, "network"), prefix = name
)


######################### Gut species and Gut MTB Pos
sub_dir <- file.path(output_dir, "Gut_MTB_Pos")
name <- "mtb_pos"
res_nc_gut_mtb_pos <- readRDS(file.path(sub_dir, paste0(name, "_g1.rds")))
res_sc_gut_mtb_pos <- readRDS(file.path(sub_dir, paste0(name, "_g2.rds")))
res_gut_mtb_pos <- read.csv(file.path(sub_dir, paste0(name, "_differential_correlations.csv")), check.names = F)
sig_diff_gut_mtb_pos <- res_gut_mtb_pos[res_gut_mtb_pos$diff_cor > threshold & res_gut_mtb_pos$p < pval_t, ]
sig_diff_save <- sig_diff_gut_mtb_pos %>%
  mutate(
    Taxa = gsub("data1 - ", "", Taxa),
    name = gsub("data2 - ", "", name)
  ) %>%
  dplyr::rename(metabolite_ID = name) %>%
  left_join(mtb_ann) %>%
  relocate(Metabolites, .after = metabolite_ID) %>%
  mutate(across(grep("cor", colnames(.)), ~ round(.x, 3))) %>%
  mutate(sign = ifelse(grepl("g1", sign), gsub("g1", group_names[1], sign), gsub("g2", group_names[2], sign))) %>%
  dplyr::rename_with(~cor_names1, all_of(cor_names)) %>%
  dplyr::rename_with(~cooccurrence_names1, all_of(cooccurrence_names)) %>%
  arrange(p)
write.csv(sig_diff_save, file.path(sub_dir, paste0(name, "_differential_correlations_sig.csv")), row.names = F, quote = F)
write.table(sig_diff_save, file.path(sub_dir, paste0(name, "_differential_correlations_sig.txt")), sep = "\t", row.names = F, quote = F)


create_networks(
  g1_mat = res_nc_gut_mtb_pos$corr, g1_mat_cooccur = res_nc_gut_mtb_pos$mat_cooccur,
  g2_mat = res_sc_gut_mtb_pos$corr, g2_mat_cooccur = res_sc_gut_mtb_pos$mat_cooccur,
  sig_diff = sig_diff_gut_mtb_pos,
  taxonomies = NULL,
  t = NULL, datasets = c("Gut Species", "Gut MTB-POS"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "union",
  repulsion = 0.85, repulsion2 = 0.7,
  cexNodes = 1, cexNodes2 = 1,
  cexHubs = 1.5, cexHubs2 = 2,
  cexLabels = 1, cexLabels2 = 1,
  cexHubLabels = 2, cexHubLabels2 = 2,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = file.path(sub_dir, "network"), prefix = name
)

######################### Gut species and Gut MTB Neg
sub_dir <- file.path(output_dir, "Gut_MTB_Neg")
name <- "mtb_neg"
res_nc_gut_mtb_neg <- readRDS(file.path(sub_dir, paste0(name, "_g1.rds")))
res_sc_gut_mtb_neg <- readRDS(file.path(sub_dir, paste0(name, "_g2.rds")))
res_gut_mtb_neg <- read.csv(file.path(sub_dir, paste0(name, "_differential_correlations.csv")), check.names = F)
sig_diff_gut_mtb_neg <- res_gut_mtb_neg[res_gut_mtb_neg$diff_cor > threshold & res_gut_mtb_neg$p < pval_t, ]
sig_diff_save <- sig_diff_gut_mtb_neg %>%
  mutate(
    Taxa = gsub("data1 - ", "", Taxa),
    name = gsub("data2 - ", "", name)
  ) %>%
  dplyr::rename(metabolite_ID = name) %>%
  left_join(mtb_ann) %>%
  relocate(Metabolites, .after = metabolite_ID) %>%
  mutate(across(grep("cor", colnames(.)), ~ round(.x, 3))) %>%
  mutate(sign = ifelse(grepl("g1", sign), gsub("g1", group_names[1], sign), gsub("g2", group_names[2], sign))) %>%
  dplyr::rename_with(~cor_names1, all_of(cor_names)) %>%
  dplyr::rename_with(~cooccurrence_names1, all_of(cooccurrence_names)) %>%
  arrange(p)
write.csv(sig_diff_save, file.path(sub_dir, paste0(name, "_differential_correlations_sig.csv")), row.names = F, quote = F)
write.table(sig_diff_save, file.path(sub_dir, paste0(name, "_differential_correlations_sig.txt")), sep = "\t", row.names = F, quote = F)

create_networks(
  g1_mat = res_nc_gut_mtb_neg$corr, g1_mat_cooccur = res_nc_gut_mtb_neg$mat_cooccur,
  g2_mat = res_sc_gut_mtb_neg$corr, g2_mat_cooccur = res_sc_gut_mtb_neg$mat_cooccur,
  sig_diff = sig_diff_gut_mtb_neg,
  taxonomies = NULL,
  t = NULL, datasets = c("Gut Species", "Gut MTB-Neg"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "union",
  repulsion = 0.85, repulsion2 = 0.7,
  cexNodes = 1, cexNodes2 = 1,
  cexHubs = 1.5, cexHubs2 = 2,
  cexLabels = 1, cexLabels2 = 1,
  cexHubLabels = 2, cexHubLabels2 = 2,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = file.path(sub_dir, "network"), prefix = name
)

######################### Gut species and Plasma MTB Pos
sub_dir <- file.path(output_dir, "Plasma_MTB_Pos")
name <- "plasma_mtb_pos"
res_nc_plasma_mtb_pos <- readRDS(file.path(sub_dir, paste0(name, "_g1.rds")))
res_sc_plasma_mtb_pos <- readRDS(file.path(sub_dir, paste0(name, "_g2.rds")))
res_plasma_mtb_pos <- read.csv(file.path(sub_dir, paste0(name, "_differential_correlations.csv")), check.names = F)
sig_diff_plasma_mtb_pos <- res_plasma_mtb_pos[res_plasma_mtb_pos$diff_cor > threshold & res_plasma_mtb_pos$p < pval_t, ]
sig_diff_save <- sig_diff_plasma_mtb_pos %>%
  mutate(
    Taxa = gsub("data1 - ", "", Taxa),
    name = gsub("data2 - ", "", name)
  ) %>%
  dplyr::rename(metabolite_ID = name) %>%
  left_join(mtb_ann) %>%
  relocate(Metabolites, .after = metabolite_ID) %>%
  mutate(across(grep("cor", colnames(.)), ~ round(.x, 3))) %>%
  mutate(sign = ifelse(grepl("g1", sign), gsub("g1", group_names[1], sign), gsub("g2", group_names[2], sign))) %>%
  dplyr::rename_with(~cor_names1, all_of(cor_names)) %>%
  dplyr::rename_with(~cooccurrence_names1, all_of(cooccurrence_names)) %>%
  arrange(p)
write.csv(sig_diff_save, file.path(sub_dir, paste0(name, "_differential_correlations_sig.csv")), row.names = F, quote = F)
write.table(sig_diff_save, file.path(sub_dir, paste0(name, "_differential_correlations_sig.txt")), sep = "\t", row.names = F, quote = F)


create_networks(
  g1_mat = res_nc_plasma_mtb_pos$corr, g1_mat_cooccur = res_nc_plasma_mtb_pos$mat_cooccur,
  g2_mat = res_sc_plasma_mtb_pos$corr, g2_mat_cooccur = res_sc_plasma_mtb_pos$mat_cooccur,
  sig_diff = sig_diff_plasma_mtb_pos,
  taxonomies = NULL,
  t = NULL, datasets = c("Gut Species", "Plasma MTB-Pos"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "union",
  repulsion = 0.85, repulsion2 = 0.7,
  cexNodes = 1, cexNodes2 = 1,
  cexHubs = 1.5, cexHubs2 = 2,
  cexLabels = 1, cexLabels2 = 1,
  cexHubLabels = 2, cexHubLabels2 = 2,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = file.path(sub_dir, "network"), prefix = name
)

######################### Gut species and Plasma MTB Neg
sub_dir <- file.path(output_dir, "Plasma_MTB_Neg")
name <- "plasma_mtb_neg"
res_nc_plasma_mtb_neg <- readRDS(file.path(sub_dir, paste0(name, "_g1.rds")))
res_sc_plasma_mtb_neg <- readRDS(file.path(sub_dir, paste0(name, "_g2.rds")))
res_plasma_mtb_neg <- read.csv(file.path(sub_dir, paste0(name, "_differential_correlations.csv")), check.names = F)
sig_diff_plasma_mtb_neg <- res_plasma_mtb_neg[res_plasma_mtb_neg$diff_cor > threshold & res_plasma_mtb_neg$p < pval_t, ]
sig_diff_save <- sig_diff_plasma_mtb_neg %>%
  mutate(
    Taxa = gsub("data1 - ", "", Taxa),
    name = gsub("data2 - ", "", name)
  ) %>%
  dplyr::rename(metabolite_ID = name) %>%
  left_join(mtb_ann) %>%
  relocate(Metabolites, .after = metabolite_ID) %>%
  mutate(across(grep("cor", colnames(.)), ~ round(.x, 3))) %>%
  mutate(sign = ifelse(grepl("g1", sign), gsub("g1", group_names[1], sign), gsub("g2", group_names[2], sign))) %>%
  dplyr::rename_with(~cor_names1, all_of(cor_names)) %>%
  dplyr::rename_with(~cooccurrence_names1, all_of(cooccurrence_names)) %>%
  arrange(p)
write.csv(sig_diff_save, file.path(sub_dir, paste0(name, "_differential_correlations_sig.csv")), row.names = F, quote = F)
write.table(sig_diff_save, file.path(sub_dir, paste0(name, "_differential_correlations_sig.txt")), sep = "\t", row.names = F, quote = F)

create_networks(
  g1_mat = res_nc_plasma_mtb_neg$corr, g1_mat_cooccur = res_nc_plasma_mtb_neg$mat_cooccur,
  g2_mat = res_sc_plasma_mtb_neg$corr, g2_mat_cooccur = res_sc_plasma_mtb_neg$mat_cooccur,
  sig_diff = sig_diff_plasma_mtb_neg,
  taxonomies = NULL,
  t = NULL, datasets = c("Gut Species", "Plasma MTB-Neg"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "union",
  repulsion = 0.85, repulsion2 = 0.7,
  cexNodes = 1, cexNodes2 = 1, # 1,
  cexHubs = 1.5, cexHubs2 = 2, # 2,
  cexLabels = 1, cexLabels2 = 1,
  cexHubLabels = 2, cexHubLabels2 = 2,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = file.path(sub_dir, "network"), prefix = name
)
