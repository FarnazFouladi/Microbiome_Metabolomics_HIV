rm(list = ls())
source("code/load.R")
dir.create("data/processed_data/gut", showWarnings = F, recursive = T)
# Load data
expvec <- read_rds("data/metagenomics/gut/SP_HIV_expvec_c.rds")
# Taxonomies
gut_taxonomies <- expvec$LKT@elementMetadata
gut_taxonomies$Species <- gsub("Unclassified","sp.",gut_taxonomies$Species)
gut_taxonomies$LKT <- gsub("Unclassified","sp.",gut_taxonomies$LKT)
saveRDS(gut_taxonomies, file = "data/processed_data/gut/taxonomies.rds")
# Total base sequenced
totalbases <- expvec$LKT@metadata$TotalBasesSequenced
# Get genome completeness
genome_com <- assays(expvec[["LKT"]])$GenomeCompleteness
genome_com <- genome_com[, names(totalbases)]

# Extract feature tables
features <- names(expvec)
features_tables <- lapply(features, function(x) {
  df_tmp <- as.data.frame(assays(expvec[[x]])$BaseCounts)
  df_tmp <- df_tmp[, names(totalbases)]
})
names(features_tables) <- features

# Convert to relative abundance
features_tables_ppm <- lapply(features, function(x) {
  df_tmp <- as.data.frame(assays(expvec[[x]])$BaseCounts)
  df_tmp <- df_tmp[, names(totalbases)]
  stopifnot(all.equal(colnames(df_tmp), names(totalbases)))
  df_relab <- sweep(df_tmp, 2, as.numeric(totalbases), "/") * 1E6
})
names(features_tables_ppm) <- features

# Filter based on relative abundances, prevalence, and genome completeness
features_filtered_ppm <- lapply(features, function(x) {
  df_tmp <- features_tables_ppm[[x]]
  if (x == "LKT") {
    # Step 1: Keep taxa with abundance of 50 ppm or more in more than 10% of samples.
    features_to_keep <- rowMeans(df_tmp >= 50) > 0.1
    df_tmp_filtered <- df_tmp[features_to_keep, ]
    df_tmp_filtered <- df_tmp_filtered[, colSums(df_tmp_filtered) > 0]

    # Step 2 filter based on a minimum genome completeness in at least 10% of samples that the feature is present
    genome_com_filtered <- genome_com[rownames(df_tmp_filtered), ]
    nonzero_samples <- apply(df_tmp_filtered, 1, function(x) which(x > 0))
    features_to_keep2 <- lapply(rownames(genome_com_filtered), function(x) {
      mean(genome_com[x, nonzero_samples[[x]]] >= 0.1) > 0.1
    })
    features_to_keep2 <- unlist(features_to_keep2)
    df_tmp_filtered <- df_tmp_filtered[features_to_keep2, ]
  } else {
    # For functional features, perform filtering only based on 10% prevalence
    features_to_keep <- rowMeans(df_tmp > 0) > 0.1
    df_tmp_filtered <- df_tmp[features_to_keep, ]
    df_tmp_filtered <- df_tmp_filtered[, colSums(df_tmp_filtered) > 0]
  }
  return(df_tmp_filtered)
})

names(features_filtered_ppm) <- features

# Get filtered features for raw counts
features_filtered_tables <- lapply(features, function(x) {
  df_tmp <- features_tables[[x]]
  df_tmp_filtered <- features_filtered_ppm[[x]]
  df_tmp <- df_tmp[rownames(df_tmp_filtered), colnames(df_tmp_filtered)]
})
names(features_filtered_tables) <- features

# Get annotations
annt_list <- lapply(features, function(x) {
  ann_df <- as.data.frame(expvec[[x]]@elementMetadata)
  if (x == "LKT") {
    ann_df$Gram <- NULL
  }
  if ("Accession" %in% colnames(ann_df)) {
    colnames(ann_df)[colnames(ann_df) == "Accession"] <- "Feature"
  }
  if (x == "LKT") {
    colnames(ann_df)[colnames(ann_df) == "LKT"] <- "Feature"
  }

  return(ann_df)
})
names(annt_list) <- features


# Metadata
meta <- read.csv("data/metadata/metadata_c.csv", check.names = F)
meta$Sample <- meta$nci_metagenomic_stool_id
meta <- meta %>% filter(Sample %in% names(totalbases))
meta <- meta[match(names(totalbases), meta$Sample), ]

stopifnot(all.equal(meta$Sample, colnames(features_filtered_tables[[1]])))

# Create new columns for cd4, cd8, and their ratio
meta <- meta %>%
  mutate(cd4 = leu3p, cd8 = leu2p, ratio = log(cd4 / cd8, base = 2)) 

# Create phyloseq objects
for (f in features) {
  print(f)
  count_tmp <- features_filtered_tables[[f]]
  meta_tmp <- meta
  stopifnot(all.equal(meta_tmp$Sample, colnames(count_tmp)))
  ann_tmp <- annt_list[[f]] %>% column_to_rownames("Feature")
  ann_tmp <- ann_tmp[rownames(count_tmp), , drop = FALSE]
  if(f == "LKT"){
    rownames(count_tmp) <- gsub("Unclassified","sp.",rownames(count_tmp))
    rownames(ann_tmp) <- gsub("Unclassified","sp.",rownames(ann_tmp))
    ann_tmp$Species <- gsub("Unclassified","sp.",ann_tmp$Species)
  }
  ann_tmp <- as.matrix(ann_tmp)
  stopifnot(all.equal(rownames(ann_tmp), rownames(count_tmp)))
  OTU <- otu_table(count_tmp, taxa_are_rows = TRUE)
  META <- sample_data(meta_tmp)
  sample_names(META) <- meta_tmp$Sample
  TAX <- tax_table(ann_tmp)
  phy_taxa <- phyloseq(OTU, TAX, META)
  saveRDS(phy_taxa, file.path("data/processed_data/gut", paste0(f, "_phy.rds")))
}
