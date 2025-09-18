rm(list = ls())
source("code/load.R")
output_dir <- "Outputs/Differential_Abundance_Comparisons"
dir.create(output_dir)

# Create a heatmap to show LFC from differential abundance analysis for 
# Pre-HIV versus Non-HIV and trend analysis for sexual activity groups
compare_results <- function(df_sc, df_sa, data_module,sample_type) {
  if (any(grepl("passed_ss", colnames(df_sc)))) {
    # Subset to significant features
    df_sc_sig <- df_sc %>%
      filter(p_statussc < 0.01, passed_ss_statussc) %>%
      filter(!grepl("s__sp", Feature))
    df_sa_sig <- df_sa %>%
      filter(p_val < 0.01, passed_ss) %>%
      filter(!grepl("s__sp", Feature))
  } else {
    df_sc_sig <- df_sc %>%
      filter(p_statussc < 0.01) %>%
      filter(!grepl("s__sp", Feature))
    df_sa_sig <- df_sa %>%
      filter(p_val < 0.01) %>%
      filter(!grepl("s__sp", Feature))
  }
  sig_features <- unique(c(df_sc_sig$Feature, df_sa_sig$Feature))

  if (data_module == "Metabolomics") {
    join <- c("Feature", "mode")
  } else {
    join <- "Feature"
  }

  # Merge data
  df_all <- df_sc %>%
    full_join(df_sa, by = join) %>%
    filter(Feature %in% sig_features)
  if (data_module == "GO") {
    df_all <- df_all %>% mutate(Feature = Description.x)
  }

  # Heatmaps will show the union features with p < 0.01 across the two analyses
  # and asterisks will show the q-values < 0.1. For microbiome data, if sensitivity test 
  # is not passed, we assign 1 to q-value to prevent drawing asterisks if q < 0.1.
  if (data_module != "Metabolomics") {
    df_all <- df_all %>%
      mutate(q_statussc = ifelse(passed_ss_statussc, q_statussc, 1)) %>%
      mutate(q_val = ifelse(passed_ss, q_val, 1))
  }

  df_all <- df_all %>%
    dplyr::select(c(
      "Feature", "lfc_statussc", "lfc_group1g4",
      "q_statussc", "q_val"
    )) %>%
    distinct(Feature, .keep_all = T)

  if (data_module == "Species"){
    df_all <- df_all %>% mutate(Feature = gsub("s__","",Feature)) %>%
      mutate(Feature = gsub("_"," ",Feature))
  }
  
  df_lfc <- df_all %>%
    column_to_rownames("Feature") %>%
    dplyr::select(grep("lfc", colnames(.)))
  df_q <- df_all %>%
    column_to_rownames("Feature") %>%
    dplyr::select(grep("q_", colnames(.)))

  # Save source data
  write.table(df_all %>%
                 column_to_rownames("Feature") %>%
                 dplyr::select(grep("lfc", colnames(.)),grep("q_", colnames(.))) %>%
              dplyr::rename( lfc_preHIVNonHIV = lfc_statussc,
                             q_preHIVNonHIV = q_statussc), 
              file.path(output_dir, paste0(tolower(sample_type),"_",tolower(data_module), "_source_data.txt")), sep = "\t", row.names = T, quote = F)
  
  
  
  # Heatmap
  col_fun <- colorRamp2(c(min(df_lfc), 0, max(df_lfc)), c("darkgreen", "white", "red"))
  hm1 <- Heatmap(as.matrix(df_lfc),
    row_names_gp = gpar(fontsize = 30, fontface = "italic"),
    column_names_gp = gpar(fontsize = 30),
    row_names_side = "left",
    column_title_gp = gpar(fontsize = 50),
    heatmap_height = unit(ifelse(data_module == "GO", 250, 75), "cm"),
    heatmap_width = unit(ifelse(data_module == "GO", 75, 50), "cm"),
    column_labels = c("Pre-HIV vs Non-HIV", "G4 versus G1"),
    border = TRUE,
    row_names_max_width = max_text_width(row.names(df_lfc), gp = gpar(fontsize = 30)),
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (df_q[i, j] < 0.001) {
        grid.text("***", x, y, gp = gpar(fontsize = 35))
      } else if (df_q[i, j] < 0.01) {
        grid.text("**", x, y, gp = gpar(fontsize = 35))
      } else if (df_q[i, j] < 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 35))
      } else if (df_q[i, j] < 0.1) {
        grid.text("+", x, y, gp = gpar(fontsize = 35))
      }
    },
    cluster_columns = FALSE, cluster_rows = TRUE, col = col_fun,
    heatmap_legend_param = list(
      legend_height = unit(5, "cm"),
      legend_width = unit(5, "cm"),
      title_gp = gpar(fontsize = 30, fontface = "bold"),
      labels_gp = gpar(fontsize = 30),
      # direction = "horizontal",
      title = "LFC"
    )
  )

  return(hm1)
}

# Gut Species
df_sc <- read.table(file.path("Outputs/Differential_Abundance_Microbiome/Gut/Species", "Species_ANCOMBC_results.txt"),
  sep = "\t", header = TRUE, check.names = F, quote = "", comment.char = ""
)
df_sa <- read.table(file.path("Outputs/Differential_Abundance_Microbiome_SA/Gut/Species", "Species_ANCOMBC_trend_results.txt"),
  sep = "\t", header = TRUE, check.names = F, quote = "", comment.char = ""
)

hm <- compare_results(df_sc, df_sa, data_module = "Species",sample_type = "Gut")
pdf(file.path(output_dir, paste0("gut_species", "_heatmap.pdf")), 40, 35)
draw(hm)
dev.off()

# Gut GO
df_sc <- read.table(file.path("Outputs/Differential_Abundance_Microbiome/Gut/GO", "GO_ANCOMBC_results.txt"),
  sep = "\t", header = TRUE, check.names = F, quote = "", comment.char = ""
)
df_sa <- read.table(file.path("Outputs/Differential_Abundance_Microbiome_SA/Gut/GO", "GO_ANCOMBC_trend_results.txt"),
  sep = "\t", header = TRUE, check.names = F, quote = "", comment.char = ""
)

hm <- compare_results(df_sc, df_sa, data_module = "GO", sample_type = "Gut")
pdf(file.path(output_dir, paste0("gut_go", "_heatmap.pdf")), 100, 100)
draw(hm)
dev.off()



# Oral Species
df_sc <- read.table(file.path("Outputs/Differential_Abundance_Microbiome/Oral/Species", "Species_ANCOMBC_results.txt"),
  sep = "\t", header = TRUE, check.names = F, quote = "", comment.char = ""
)
df_sa <- read.table(file.path("Outputs/Differential_Abundance_Microbiome_SA/Oral/Species", "Species_ANCOMBC_trend_results.txt"),
  sep = "\t", header = TRUE, check.names = F, quote = "", comment.char = ""
)

hm <- compare_results(df_sc, df_sa, data_module = "Species", sample_type = "Oral")
pdf(file.path(output_dir, paste0("oral_species", "_heatmap.pdf")), 40, 35)
draw(hm)
dev.off()

# Oral GO
df_sc <- read.table(file.path("Outputs/Differential_Abundance_Microbiome/Oral/GO", "GO_ANCOMBC_results.txt"),
  sep = "\t", header = TRUE, check.names = F, quote = "", comment.char = ""
)
df_sa <- read.table(file.path("Outputs/Differential_Abundance_Microbiome_SA/Oral/GO", "GO_ANCOMBC_trend_results.txt"),
  sep = "\t", header = TRUE, check.names = F, quote = "", comment.char = ""
)

hm <- compare_results(df_sc, df_sa, data_module = "GO", sample_type = "Oral")
pdf(file.path(output_dir, paste0("oral_go", "_heatmap.pdf")), 70, 110)
draw(hm)
dev.off()



# Gut Metabolites
df_sc <- read.table(file.path("Outputs/Differential_Abundance_Metabolomics/Gut", "metabolites_ANCOMBC_results.txt"),
  sep = "\t", header = TRUE, check.names = F, quote = "", comment.char = ""
)
df_sa <- read.table(file.path("Outputs/Differential_Abundance_Metabolomics_SA/Gut", "metabolites_ANCOMBC_results.txt"),
  sep = "\t", header = TRUE, check.names = F, quote = "", comment.char = ""
)


hm <- compare_results(df_sc, df_sa, data_module = "Metabolomics", sample_type = "Gut")
pdf(file.path(output_dir, paste0("gut_metabolomics", "_heatmap.pdf")), 40, 40)
print(draw(hm))
dev.off()


# Plasma Metabolites
df_sc <- read.table(file.path("Outputs/Differential_Abundance_Metabolomics/Plasma", "metabolites_ANCOMBC_results.txt"),
  sep = "\t", header = TRUE, check.names = F, quote = "", comment.char = ""
)
df_sa <- read.table(file.path("Outputs/Differential_Abundance_Metabolomics_SA/Plasma", "metabolites_ANCOMBC_results.txt"),
  sep = "\t", header = TRUE, check.names = F, quote = "", comment.char = ""
)

hm <- compare_results(df_sc, df_sa, data_module = "Metabolomics", sample_type = "Plasma")
pdf(file.path(output_dir, paste0("plasma_metabolomics", "_heatmap.pdf")), 40, 40)
draw(hm)
dev.off()


# Oral Metabolites
df_sc <- read.table(file.path("Outputs/Differential_Abundance_Metabolomics/Oral", "metabolites_ANCOMBC_results.txt"),
  sep = "\t", header = TRUE, check.names = F, quote = "", comment.char = ""
)
df_sa <- read.table(file.path("Outputs/Differential_Abundance_Metabolomics_SA/Oral", "metabolites_ANCOMBC_results.txt"),
  sep = "\t", header = TRUE, check.names = F, quote = "", comment.char = ""
)

hm <- compare_results(df_sc, df_sa, data_module = "Metabolomics", sample_type = "Oral")
pdf(file.path(output_dir, paste0("oral_metabolomics", "_heatmap.pdf")), 40, 40)
draw(hm)
dev.off()

