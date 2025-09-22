# Validation analysis - External dataset
rm(list = ls())
source("code/load.R")
source("code/bias_correction.R")
source("code/differential_correlations_functions.R")
source("code/create_heatmap.R")

output_dir <- file.path("Outputs", "External_Dataset", "Armstrong_2018")
dir.create(output_dir, recursive = T, showWarnings = F)
method <- "spearman"
cov_cat <- NULL
cov_num <- c("host_age")
assay_name <- "counts"
tax_level <- NULL
prv_cut <- 0.3
threshold <- 0.3
filt_th <- 0.001
groups <- c("MSM_HIV+ART-", "MSM_HIV+ART")

# Load data
meta <- read.table("data/external_dataset/Armstrong_2018/11914_20180904-130826.txt", sep = "\t", header = T, check.names = F, quote = "", comment.char = "")
table(meta$gender, meta$msm_status)
table(meta$msm_status, meta$hiv)
meta_sub <- meta %>% filter(msm_status %in% c("MSM", "MSW"))
table(meta_sub$high_risk, meta_sub$msm_status)
table(meta_sub$high_risk, meta_sub$hiv)

meta_sub <- meta_sub %>%
  mutate(groups = ifelse(msm_status == "MSW" & high_risk == "No" & hiv == "Negative", "Control",
    ifelse(msm_status == "MSM" & high_risk == "Yes" & hiv == "Negative", "MSM_HIV-",
      ifelse(msm_status == "MSM" & high_risk == "HIV_pos" & treatment == "treated", "MSM_HIV+ART",
        ifelse(msm_status == "MSM" & high_risk == "HIV_pos" & treatment == "untreated", "MSM_HIV+ART-", "Others")
      )
    )
  )) %>%
  filter(groups != "Others")


table(meta_sub$groups)


ft_df <- qiime2R::read_qza("data/external_dataset/Armstrong_2018/155590_feature-table.qza")$data

# Create a fasta file seqs
seqs <- rownames(ft_df)
names(seqs) <- paste0("ASV", seq_along(seqs))
rownames(ft_df) <- names(seqs)
fasta_lines <- unlist(lapply(names(seqs), function(nm) {
  c(paste0(">", nm), seqs[nm])
}))
# Write to file (Sequences were taxonomically classified using the SILVA database )
writeLines(fasta_lines, file.path("data/external_dataset/Armstrong_2018", "seqs.fasta"))

tax_df <- qiime2R::read_qza("data/external_dataset/Armstrong_2018/taxonomy.qza")$data
# Clean up white space
tax_df$Taxon <- gsub("\\s+", "", tax_df$Taxon)
df_split <- tax_df %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right")
df_split <- df_split %>%
  rowwise() %>%
  mutate(LKT = {
    last_non_na <- head(na.omit(c(Species, Genus, Family, Order, Class, Phylum, Domain)), 1)
    ifelse(is.null(last_non_na), NA, last_non_na)
  }) %>%
  ungroup() %>%
  column_to_rownames("Feature.ID")
all.equal(rownames(ft_df), rownames(df_split))


sum(colnames(ft_df) %in% meta_sub$sample_name)
asv <- t(ft_df[, meta_sub$sample_name])
rownames(meta_sub) <- meta_sub$sample_name
all.equal(meta_sub$sample_name, rownames(asv))

# Create phyloseq data
OTU <- otu_table(t(asv), taxa_are_rows = TRUE)
META <- sample_data(meta_sub)
TAX <- tax_table(as.matrix(df_split))
phy_taxa <- phyloseq(OTU, TAX, META)


for (group in groups) {
  subdir <- file.path(output_dir, group)
  dir.create(subdir, showWarnings = F, recursive = T)
  name <- paste0("Armstrong_", group)
  # Subset data
  phy_sub <- phy_taxa %>% subset_samples(groups %in% c("MSM_HIV-", group))
  phy_sub <- prune_taxa(taxa_sums(phy_sub) > 0, phy_sub)
  phy_g1 <- phy_sub %>% subset_samples(groups == "MSM_HIV-")
  phy_g2 <- phy_sub %>% subset_samples(groups == group)
  meta1 <- microbiome::meta(phy_g1)
  meta2 <- microbiome::meta(phy_g2)


  g1_bias_corrected <- bias_correction(list(phy_g1),
    tax_level = "LKT", pseudo = 0, prv_cut = prv_cut,
    lib_cut = 0, cov_cat = cov_cat, cov_num = cov_num
  )$residuals

  g2_bias_corrected <- bias_correction(list(phy_g2),
    tax_level = "LKT", pseudo = 0, prv_cut = prv_cut,
    lib_cut = 0, cov_cat = cov_cat, cov_num = cov_num
  )$residuals

  # Gut species
  df.g1 <- t(g1_bias_corrected)
  df.g2 <- t(g2_bias_corrected)
  common_features <- intersect(colnames(df.g1), colnames(df.g2))
  df.g1 <- df.g1[, common_features]
  df.g2 <- df.g2[, common_features]

  species_res <- disco_score(
    data1 = list(df.g1), data2 = list(df.g2),
    output_dir = subdir, prefix = paste0(name, "_gut_species")
  )
}

# Heatmap
no_ART <- read.csv(file.path("Outputs/External_Dataset/Armstrong_2018/MSM_HIV+ART-", paste0("Armstrong_MSM_HIV+ART-", "_gut_species_disco.csv")), check.names = F)
with_ART <- read.csv(file.path("Outputs/External_Dataset/Armstrong_2018/MSM_HIV+ART", paste0("Armstrong_MSM_HIV+ART", "_gut_species_disco.csv")), check.names = F)

res_list <- list(no_ART, with_ART)
names(res_list) <- c("ART-naïve", "ART-treated")

create_heatmap(res_list,
  p_alpha = 0.01, bh_alpha = 0.1,
  col_names = names(res_list), top_n = 10,
  sample_type = "Gut", output_dir = "Outputs/External_Dataset/Armstrong_2018",
  feature = "Taxa"
)
