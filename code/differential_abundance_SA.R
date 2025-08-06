# Differential abundance testing
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
sample_type <- args[1]
source("code/load.R")
dir_diff <- file.path(output_dir, "Differential_Abundance_Microbiome_SA", sample_type)
dir.create(dir_diff, showWarnings = F, recursive = T)
features_all <- c("Species", "GO")


##################################
# Differential Abundance Testing
##################################
for (feature in features_all) {
  phy_obj <- readRDS(file.path(
    "data/processed_data", tolower(sample_type),
    paste0(ifelse(feature == "Species", "LKT", feature), "_phy.rds")
  ))

  dir_diff_sub <- file.path(dir_diff, feature)
  dir.create(dir_diff_sub, recursive = T, showWarnings = F)

  sample_data(phy_obj)$loc <- as.factor(sample_data(phy_obj)$loc)

  # Get annotation
  ann_tmp <- tax_table(phy_obj) %>%
    as.data.frame() %>%
    rownames_to_column("Feature")


  if (feature == "Species") {
    tax <- feature
  } else {
    tax <- NULL
  }

  # Perform ANCOMBC2
  diff_out <- ANCOMBC::ancombc2(
    data = phy_obj, assay_name = "counts", tax_level = tax,
    fix_formula = "group1 + age + loc + abx_use", rand_formula = NULL,
    p_adj_method = "BH", pseudo_sens = TRUE,
    prv_cut = 0.2, s0_perc = 0.05,
    group = "group1", struc_zero = TRUE, neg_lb = TRUE, dunnet = TRUE,
    alpha = 0.05, n_cl = 20, verbose = TRUE, trend = T, trend_control = list(       
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
  if (feature == "Species") {
    res <- diff_out$res %>% dplyr::rename(Feature = taxon)
    res_dunnett <- diff_out$res_dunn %>% dplyr::rename(Feature = taxon)
    res_trend <- diff_out$res_trend %>% dplyr::rename(Feature = taxon)
  } else {
    res <- diff_out$res %>%
      dplyr::rename(Feature = taxon) %>%
      left_join(ann_tmp, by = c("Feature"))
    
    res_dunnett <- diff_out$res_dunn %>%
      dplyr::rename(Feature = taxon) %>%
      left_join(ann_tmp, by = c("Feature"))
    
    res_trend <- diff_out$res_trend %>%
      dplyr::rename(Feature = taxon) %>%
      left_join(ann_tmp, by = c("Feature"))
  }

  # Save results
  write.table(res, file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_results.txt")), row.names = F, quote = F, sep = "\t")
  write.table(res_dunnett, file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_dunnett_results.txt")), row.names = F, quote = F, sep = "\t")
  write.table(res_trend, file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_trend_results.txt")), row.names = F, quote = F, sep = "\t")
  write.table(bias_corrected_df, file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_bias_corrected_data.txt")), row.names = T, quote = F, sep = "\t")
}


##################################
# Plotting
##################################

for (feature in features_all) {
 dir_diff_sub <- file.path(dir_diff, feature)
  
  # Upload the results
  ancombc <- read.csv(file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_trend_results.txt")), sep = "\t", check.names = F, quote = "", comment.char = "")
  df <- read.csv(file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_bias_corrected_data.txt")), sep = "\t", check.names = F, quote = "", comment.char = "")
  phy_obj <- readRDS(file.path(
    "data/processed_data", tolower(sample_type),
    paste0(ifelse(feature == "Species", "LKT", feature), "_phy.rds")
  ))
  meta <- microbiome::meta(phy_obj)
  
  # Remove any missing or unclassified species
  ancombc <- ancombc %>% filter(!grepl("__sp|__Missing", Feature))
  df <- df[!grepl("__sp|__Missing", rownames(df)), ]
  
  if(feature == "GO")
    ancombc$Description[ancombc$Description=="oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of two atoms of oxygen"] <-
    "oxidoreductase activity"
  
  
  # Waterfall plot for significant results
  ancombc_sig <- ancombc %>%
    filter(p_val < 0.01, passed_ss) %>%
    arrange(lfc_group1g4) %>%
    mutate(
      Feature = gsub("^[a-z]__", "", Feature),
      Feature = factor(Feature, levels = unique(Feature)),
      change = ifelse(lfc_group1g4 < 0, paste("Lower in", "G4"), paste("Higher in", "G4")),
      change = factor(change, levels = c(paste("Lower in", "G4"), paste("Higher in", "G4"))),
      p_sig = ifelse(q_val < 0.001, "***", ifelse(q_val < 0.01, "**", ifelse(q_val < 0.1, "*", ""))),
      q = round(q_val, 3),
      y_p = ifelse(lfc_group1g4 > 0, lfc_group1g4 + se_group1g4 + 0.2, lfc_group1g4 - se_group1g4 - 0.2)
    )
  
  if (feature == "Species") {
    ancombc_sig <- ancombc_sig %>%
      mutate(
        Feature = gsub("_", " ", Feature),
        Feature = factor(Feature, levels = unique(Feature))
      )
  }
  
  
  if (nrow(ancombc_sig) > 0) {
    # Save the significant results
    if (feature == "Species") {
      ancombc_save <- ancombc_sig %>%
        dplyr::select(c("Feature","lfc_group1g2","lfc_group1g3", "lfc_group1g4", "p_val", "q_val")) %>%
        dplyr::rename(p = p_val, q = q_val) %>%
        mutate(across(grep("lfc", colnames(.)), ~ round(.x, 3))) %>%
        arrange(q)
      variable <- "Feature"
    } else {
      ancombc_save <- ancombc_sig %>%
        dplyr::select(c("Feature", "Description","lfc_group1g2","lfc_group1g3", "lfc_group1g4", "p_val", "q_val", "Ontology", "Definition")) %>%
        dplyr::rename( p = p_val, q = q_val) %>%
        mutate(across(grep("lfc", colnames(.)), ~ round(.x, 3))) %>%
        arrange(q)
      variable <- "Description"
    }
    
    write.table(ancombc_save, file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_results_sig.txt")), sep = "\t", row.names = F, quote = F)
    
    # If there are more than 25 features, only show the first 25 most significant ones
    if (nrow(ancombc_sig) > 30 & feature == "GO") {
      ancombc_sig <- ancombc_sig %>%
        # slice_head(n = 25) %>%
        arrange(lfc_group1g4) %>%
        dplyr::slice(1:15, (n() - 14):n()) %>%
        mutate(
          Description = factor(Description, levels = unique(Description))
        )
    }
    
    # Set the colors
    if (length(unique(ancombc_sig$change)) == 1 && unique(ancombc_sig$change) == paste("Higher in","G4")) {
      cols_sub <- "#bc3c29ff"
    } else if (length(unique(ancombc_sig$change)) == 1 && unique(ancombc_sig$change) == paste("Lower in","G4")) {
      cols_sub <- "#0072b5ff"
    } else {
      cols_sub <- c("#0072b5ff","#bc3c29ff" )
    }
    
     # Waterfall plot
    pdf(file.path(dir_diff_sub, paste0(feature, "_wf_fig.pdf")), 12, ifelse(nrow(ancombc_sig) > 4, 8, 4))
    
    p1 <- ggplot(ancombc_sig, aes(y = .data[[variable]], x = .data[["lfc_group1g4"]], fill = .data[["change"]])) +
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
      {
        if (feature == "Species") theme(axis.text.y = element_text(face = "italic"))
      } +
      geom_text(aes(x = y_p, label = p_sig), size = 5) +
      labs(y = "", title = paste(sample_type, feature), x = "Log Fold Change", fill = "")
    
    print(p1)
    dev.off()
  
  
  # Heatmap
  trend_sig <- ancombc %>%
    filter(p_val < 0.01, passed_ss)
  
  # For GO terms, select top 30
  if (nrow(trend_sig) > 30 & feature == "GO") {
    trend_sig <- trend_sig %>%
      # slice_head(n = 25) %>%
      arrange(lfc_group1g4) %>%
      dplyr::slice(1:15, (n() - 14):n()) %>%
      mutate(
        Feature = factor(Description, levels = unique(Description))
      )
  }  
  
  # trim the species name
  if (feature == "Species") {
    trend_sig <- trend_sig %>%
      mutate(
        Feature = gsub("s__", "", Feature),
        Feature = gsub("_", " ", Feature),
        Feature = factor(Feature, levels = unique(Feature))
      )
  }
  
  # Convert to long format
  trend_sig <- trend_sig  %>%
    mutate(trend = ifelse(lfc_group1g4 > 0, "increasing","decreasing")) %>%
    pivot_longer(
      cols = lfc_group1g2:se_group1g4,
      cols_vary = "slowest",
      names_to = c(".value", "study_condition"),
      names_pattern = "(.*)_(.*)"
    ) %>%
    mutate_if(is.numeric, function(x) round(x, 2)) %>%
    mutate(study_condition = factor(gsub("group1", "", study_condition), levels = c("g2","g3","g4"))) %>%
    mutate(study_condition = toupper(study_condition) ) %>%
    mutate(p_sig = ifelse(q_val < 0.001, "***", ifelse(q_val < 0.01, "**", ifelse(q_val < 0.1, "*", "")))) %>%
    mutate(lb = paste(lfc,p_sig)) %>%
    arrange(trend,lfc,Feature) %>%
    mutate(Feature = factor(Feature,levels = unique(Feature)))
  
  
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
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    {
      if (feature == "Species") theme(axis.text.y = element_text(face = "italic"))
    }
  
  
  pdf(file.path(dir_diff_sub, paste0(feature, "_trend_fig.pdf")), 10, ifelse(sample_type == "Gut",10,8))
  print(trend_fig)
  dev.off()
  }
}
