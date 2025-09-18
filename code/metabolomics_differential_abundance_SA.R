# Differential abundance analyses
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
sample_type <- args[1]
dt = "metabolites"
source("code/load.R")
output_dir <- file.path(output_dir, "Differential_Abundance_Metabolomics_SA", sample_type)
sample_id <- paste0(tolower(sample_type), "_sample_id")

# Load metadata
meta <- read_csv("data/metadata/metadata_c.csv")
meta$metab_sample_id <- as.character(meta[[sample_id]])
# Load annotation
mtb_ann <- read.table(file.path("data/metabolomics_annotation/metabolites_annotation.txt"), sep = "\t", comment.char = "", quote = "", check.names = F, header = T)

##################################
# Differential Abundance Testing
##################################

for (mode in c("positive", "negative")) {
  # Load metabolomics data
  df_met <- read.table(file.path("data/metabolomics",tolower(sample_type),paste0(sample_type,"_",mode,"_metabolites.txt")),
                       sep = "\t",header = T,check.names = F,quote = "",comment.char = "")
   
  # Filter metabolites with low confidence detection
  df_met <- df_met[df_met$`Delta(ppm)` < 5, ]
  df_sub <- df_met[, c(which(colnames(df_met) == "Compound_ID"), which(colnames(df_met) == "1"):ncol(df_met))]
  df_sub <- df_sub %>% remove_rownames %>% column_to_rownames("Compound_ID")

  # Remove features with low variation
  sds <- apply(df_sub, 1, sd)
  means <- apply(df_sub, 1, mean)
  rsd <- abs(sds / means)
  rk <- rank(-rsd, ties.method = "random")
  remain.num <- round(nrow(df_sub) * 0.9)
  remain <- rk <= remain.num
  df_filtered <- df_sub[remain, ]
  df_sub <- df_filtered

  meta_filt <- meta %>% filter(metab_sample_id %in% colnames(df_sub))
  df_sub <- df_sub[, meta_filt$metab_sample_id]
  stopifnot(all.equal(meta_filt$metab_sample_id, colnames(df_sub)))

  meta_filt$loc <- as.factor(meta_filt$loc)

  # Create a metabolite dataframe
  met <- data.frame(accession = rownames(df_sub), row.names = rownames(df_sub)) %>%
    left_join(mtb_ann) %>%
    column_to_rownames("accession") %>%
    as.matrix()

  # Create phyloseq object
  OTU <- otu_table(df_sub, taxa_are_rows = TRUE)
  META <- sample_data(meta_filt)
  sample_names(META) <- meta_filt$metab_sample_id
  TAX <- tax_table(met)
  phy_obj <- phyloseq(OTU, TAX, META)
  if (!dir.exists(file.path("data/processed_data", sample_type))) {
    dir.create(file.path("data/processed_data", sample_type))
  }
  saveRDS(phy_obj, file.path("data/processed_data", sample_type, paste0(dt, "_", mode, "_phy.rds")))

  dir.create(file.path(output_dir, mode, dt), recursive = T)

  # Differential abundance testing
  diff_out <- ANCOMBC::ancombc2(
    data = phy_obj, assay_name = "counts",
    fix_formula = "group1 + age + loc + abx_use", rand_formula = NULL,
    p_adj_method = "BH", pseudo_sens = FALSE,
    prv_cut = 0.2, lib_cut = 0,
    group = "group1", struc_zero = FALSE,
    alpha = 0.05, n_cl = 5, verbose = TRUE,trend = T, trend_control = list(       
      contrast = list(matrix(c( 1, 0, 0, 
                                -1, 1, 0,
                                0, -1, 1), # Increasing trend 
                             nrow = 3,
                             byrow = TRUE
      ),
      matrix(c(-1, 0,0, 
               1, -1,0,
               0, 1, -1), # decreasing trend
             nrow = 3,
             byrow = TRUE
      )),
      node = list(3, 3),
      solver = "ECOS",
      B = 1000  
    )
  )

  # Bias corrected abundances
  bias_corrected_df <- diff_out$bias_correct_log_table
  # Extract results
  res <- diff_out$res %>%
    dplyr::rename(accession = taxon) %>%
    left_join(met %>% as.data.frame() %>% rownames_to_column("accession")) %>%
    relocate(Compound_name, .after = accession) %>%
    dplyr::rename(Feature = Compound_name)
  
  res_dunnett <- diff_out$res_dunn %>%
    dplyr::rename(accession = taxon) %>%
    left_join(met %>% as.data.frame() %>% rownames_to_column("accession")) %>%
    relocate(Compound_name, .after = accession) %>%
    dplyr::rename(Feature = Compound_name)
  
  res_trend <- diff_out$res_trend %>% 
    dplyr::rename(accession = taxon) %>%
    left_join(met %>% as.data.frame() %>% rownames_to_column("accession")) %>%
    relocate(Compound_name, .after = accession) %>%
    dplyr::rename(Feature = Compound_name)

  # Save results
  write.table(res, file.path(output_dir, mode, dt, paste0(mode, "_", dt, "_ANCOMBC_results.txt")), row.names = F, quote = F, sep = "\t")
  write.table(res_dunnett, file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_dunnett_results.txt")), row.names = F, quote = F, sep = "\t")
  write.table(res_trend, file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_trend_results.txt")), row.names = F, quote = F, sep = "\t")
  write.table(bias_corrected_df, file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_bias_corrected_data.txt")), row.names = T, quote = F, sep = "\t")
  
}

##################################
# Plotting
##################################

output_dir_sub <- file.path(output_dir)

# Upload the results
ancombc_negative <- read.table(file.path(output_dir_sub, "negative/metabolites", paste0("negative", "_", dt, "_ANCOMBC_trend_results.txt")),
  sep = "\t", check.names = F, quote = "", comment.char = "", header = T
) %>%
  mutate(mode = "negative")
df_negative <- read.table(file.path(output_dir_sub, "negative/metabolites", paste0("negative", "_", dt, "_ANCOMBC_bias_corrected_data.txt")), sep = "\t", check.names = F, quote = "", comment.char = "", header = T)

ancombc_positive <- read.table(file.path(output_dir_sub, "positive/metabolites", paste0("positive", "_", dt, "_ANCOMBC_trend_results.txt")),
  sep = "\t", check.names = F, quote = "", comment.char = "", header = T
) %>%
  mutate(mode = "positive")
df_positive <- read.table(file.path(output_dir_sub, "positive/metabolites", paste0("positive", "_", dt, "_ANCOMBC_bias_corrected_data.txt")), sep = "\t", check.names = F, quote = "", comment.char = "", header = T)

ancombc <- rbind(ancombc_negative, ancombc_positive)
write.table(ancombc, file.path(output_dir_sub, paste0("metabolites", "_ANCOMBC_results.txt")), row.names = F, quote = F, sep = "\t")



# Waterfall plot for significant results
ancombc_sig <- ancombc %>%
  filter(p_val < 0.01) %>%
  arrange(lfc_group1g4) %>%
  mutate(
    Feature = factor(Feature, levels = unique(Feature)),
    change = ifelse(lfc_group1g4 < 0, paste("Lower in", "G4"), paste("Higher in", "G4")),
    change = factor(change, levels = c(paste("Lower in", "G4"), paste("Higher in", "G4"))),
    p_sig = ifelse(q_val < 0.001, "***", ifelse(q_val < 0.01, "**", ifelse(q_val < 0.05, "*" ,ifelse(q_val < 0.1, "+", "")))),
    q = round(q_val, 3),
    y_p = ifelse(lfc_group1g4 > 0, lfc_group1g4 + se_group1g4 + 0.2, lfc_group1g4 - se_group1g4 - 0.2)
  ) %>%
  distinct(Feature,.keep_all = T)

if (nrow(ancombc_sig) > 0) {
  # Save the significant results
  ancombc_save <- ancombc_sig %>%
    dplyr::select(c(
      "Feature","lfc_group1g2","lfc_group1g3" ,"lfc_group1g4", "p_val", "q_val", "smiles", "kingdom", "super_class",
      "class", "sub_class", "direct_parent", "mode"
    )) %>%
    mutate(across(grep("lfc", colnames(.)), ~ round(.x, 3))) %>%
    arrange(q_val) %>%
    dplyr::rename(p = p_val, q = q_val)

  write.table(ancombc_save, file.path(output_dir_sub, paste0(dt, "_ANCOMBC_results_sig.txt")), row.names = F, quote = F, sep = "\t")

  # if (nrow(ancombc_sig) > 25) {
  #   ancombc_sig <- ancombc_sig %>%
  #     arrange(q_val) %>%
  #     slice_head(n = 25)
  # }

  if (length(unique(ancombc_sig$change)) == 1 && unique(ancombc_sig$change) == paste("Higher in", "G4")) {
    cols_sub <-"#bc3c29ff" 
  } else if (length(unique(ancombc_sig$change)) == 1 && unique(ancombc_sig$change) == paste("Lower in", "G4")) {
    cols_sub <- "#0072b5ff"
  } else {
    cols_sub <- c("#0072b5ff","#bc3c29ff")
  }

  # Waterfall plot
  p1 <- ggplot(ancombc_sig, aes(y = .data[["Feature"]], x = .data[["lfc_group1g4"]], fill = .data[["change"]])) +
    geom_bar(stat = "identity") +
    geom_errorbar(
      aes(
        xmin = lfc_group1g4 - se_group1g4,
        xmax = lfc_group1g4 + se_group1g4
      ),
      width = 0.2,
      position = position_dodge(0.05), color = "black"
    ) +
    theme_bw(base_size = 16) +
    scale_fill_manual(values = as.character(cols_sub)) +
    geom_text(aes(x = y_p, label = p_sig), size = 5) +
    labs(y = "", title = paste(sample_type, "Metabolites"), x = "Log Fold Change", fill = "")

  pdf(file.path(output_dir_sub, paste0(dt, "_wf_fig.pdf")), 10, ifelse(sample_type=="Oral",5,10))
  print(p1)
  dev.off()
  
  
  # Heatmap
  trend_sig <- ancombc %>%
    filter(p_val < 0.01 ) %>%
    mutate(trend = ifelse(lfc_group1g4 > 0, "increasing","decreasing")) %>%
    distinct(Feature,.keep_all = T) %>%
    pivot_longer(
      cols = lfc_group1g2:se_group1g4,
      cols_vary = "slowest",
      names_to = c(".value", "study_condition"),
      names_pattern = "(.*)_(.*)"
    ) %>%
    mutate_if(is.numeric, function(x) round(x, 2)) %>%
    mutate(study_condition = factor(gsub("group1", "", study_condition), levels = c("g2","g3","g4"))) %>%
    mutate(study_condition = toupper(study_condition) ) %>%
    mutate(p_sig = ifelse(q_val < 0.001, "***", ifelse(q_val < 0.01, "**", ifelse(q_val < 0.05,"*",ifelse(q_val < 0.1, "+", ""))))) %>%
    mutate(lb = paste(lfc,p_sig)) %>%
    arrange(trend,lfc,Feature) %>%
    mutate(Feature = factor(Feature,levels = unique(Feature)))
  
  # Save source data
  write.table(trend_sig[,c("Feature","study_condition","lfc","p_val","q_val","p_sig")] %>%
                arrange(Feature,study_condition),
              file.path(output_dir_sub, paste0(dt, "_source_data.txt")), sep = "\t", row.names = F, quote = F)
  
  trend_fig <- ggplot(trend_sig, aes(y = Feature, x = study_condition, fill = lfc)) +
    geom_tile(color = "black") +
    geom_text(aes(label = lb), color = "black", size = 4) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 16) +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white",
      na.value = "white", midpoint = 0, limit = c(min(trend_sig$lfc), 
                                                  max(trend_sig$lfc)),
      name = "LFC"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) 


  pdf(file.path(output_dir_sub, paste0(dt, "_trend_fig.pdf")), 10, ifelse(sample_type=="Oral",5,15))
  print(trend_fig)
  dev.off()
  
  
  
}

