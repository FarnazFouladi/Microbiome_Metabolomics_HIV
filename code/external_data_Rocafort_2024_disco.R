# Validation analysis - External dataset
rm(list = ls())
source("code/load.R")
source("code/bias_correction.R")
source("code/differential_correlations_functions.R")
source("code/create_heatmap.R")

output_dir <- file.path("Outputs", "External_Dataset", "Rocafort_2024")
dir.create(output_dir, recursive = T, showWarnings = F)
method <- "spearman"
assay_name <- "counts"
tax_level <- "Species"
filt_th <- 0.001
prv_cut <- 0.3
threshold <- 0.3
pval_t <- 0.01
cohorts <- c("uganda_2", "botswana", "boston")
groups <- c("HIV_negative", "HIV_positive")
sexual_orientations <- c("heterosexual", "MSM")
cov_cat <- c("Ethnicity")
cov_num <- c("Age")


# Load data
meta <- read.table(file.path("data/external_dataset", "Rocafort_2024", "41467_2023_44566_MOESM5_ESM.csv"), sep = ";", header = T, check.names = F, quote = "", comment.char = "")
asv <- read.table(file.path("data/external_dataset", "Rocafort_2024", "41467_2023_44566_MOESM4_ESM.csv"), sep = ";", header = T, check.names = F, quote = "", comment.char = "", row.names = 1)
all.equal(meta$SubjectID, rownames(asv))
rownames(meta) <- meta$SubjectID
meta$Sexual_Orientation[is.na(meta$Sexual_Orientation)] <- "heterosexual"

# Create taxonomic table
ann_df <- data.frame(asv = paste0("ASV", 1:ncol(asv)), taxa = colnames(asv), row.names = paste0("ASV", 1:ncol(asv)))
ann_df <- ann_df %>%
  mutate(Genus = sapply(stringi::stri_split(taxa, fixed = "_"), function(x) {
    return(x[6])
  })) %>%
  mutate(Species = sapply(stringi::stri_split(taxa, fixed = "_"), function(x) {
    return(x[7])
  }))
colnames(asv) <- ann_df$asv


ann_df$Species <- sapply(1:nrow(ann_df), function(i) {
  if (!grepl("\\[", ann_df$Species[i])) {
    return(paste(ann_df$Genus[i], ann_df$Species[i]))
  } else {
    return(ann_df$Species[i])
  }
})


meta$hiv_status <- ifelse(meta$HIV_Phenotype == "1_hiv_negative", "HIV_negative", "HIV_positive")

# Create phyloseq data
OTU <- otu_table(t(asv), taxa_are_rows = TRUE)
META <- sample_data(meta)
TAX <- tax_table(as.matrix(ann_df))
phy_taxa <- phyloseq(OTU, TAX, META)


for (cohort in cohorts) {
  for (so in sexual_orientations) {
    if (so == "MSM" & cohort != "boston") {
      next
    }

    if (so == "MSM") {
      cov_cat <- c("Ethnicity")
    } else {
      cov_cat <- c("Ethnicity", "Gender")
    }


    subdir <- file.path(output_dir, cohort, so)
    dir.create(subdir, recursive = T, showWarnings = F)
    name <- paste0("Rocafort_", cohort, "_", so)

    # Subset data
    phy_sub <- subset_samples(phy_taxa, Cohort == cohort)
    phy_sub <- phy_sub %>% subset_samples(Sexual_Orientation == so)

    phy_sub <- prune_taxa(taxa_sums(phy_sub) > 0, phy_sub)

    # Subset to control (g1) and HIV positive
    phy_g1 <- phy_sub %>% subset_samples(hiv_status == groups[1])
    phy_g2 <- phy_sub %>% subset_samples(hiv_status != groups[1])
    meta1 <- microbiome::meta(phy_g1)
    meta2 <- microbiome::meta(phy_g2)


    g1_bias_corrected <- bias_correction(list(phy_g1),
      tax_level = "Species", pseudo = 0, prv_cut = prv_cut,
      lib_cut = 0, cov_cat = cov_cat, cov_num = cov_num
    )$residuals

    g2_bias_corrected <- bias_correction(list(phy_g2),
      tax_level = "Species", pseudo = 0, prv_cut = prv_cut,
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
      output_dir = subdir,prefix = paste0(name,"_gut_species")
    )

  }
}

# Heatmap
boston_msm <- read.csv(file.path("Outputs/External_Dataset/Rocafort_2024/boston/MSM", paste0("Rocafort_boston_MSM", "_gut_species_disco.csv")),check.names = F )
boston_hetr <- read.csv(file.path("Outputs/External_Dataset/Rocafort_2024/boston/heterosexual", paste0("Rocafort_boston_heterosexual", "_gut_species_disco.csv")),check.names = F )
bostwana_hetr <- read.csv(file.path("Outputs/External_Dataset/Rocafort_2024/botswana/heterosexual", paste0("Rocafort_botswana_heterosexual", "_gut_species_disco.csv")),check.names = F )
uganda_hetr <- read.csv(file.path("Outputs/External_Dataset/Rocafort_2024/uganda_2/heterosexual", paste0("Rocafort_uganda_2_heterosexual", "_gut_species_disco.csv")),check.names = F )

res_list <- list(boston_msm,boston_hetr,bostwana_hetr,uganda_hetr)
names(res_list)<-c("Boston MSM","Boston Heterosexual","Botswana","Uganda")

create_heatmap(res_list,
               p_alpha = 0.01,bh_alpha = 0.1,
               col_names = names(res_list),top_n = 10,
               sample_type = "Gut", output_dir = "Outputs/External_Dataset/Rocafort_2024",
               feature = "Taxa"
)
