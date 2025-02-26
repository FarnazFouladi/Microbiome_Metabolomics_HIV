# Differential abundance testing
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
sample_type <- args[1]
source("code/load.R")
dir_diff <- file.path(output_dir, "Differential_Abundance_Microbiome", sample_type)
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
    fix_formula = "status + age + loc + abx_use", rand_formula = NULL,
    p_adj_method = "BH", pseudo_sens = TRUE,
    prv_cut = 0.2, s0_perc = 0.05,
    group = "status", struc_zero = TRUE, neg_lb = TRUE,
    alpha = 0.05, n_cl = 5, verbose = TRUE
  )
  # Bias corrected abundances
  bias_corrected_df <- diff_out$bias_correct_log_table
  # Extract results
  if (feature == "Species") {
    res <- diff_out$res %>% dplyr::rename(Feature = taxon)
  } else {
    res <- diff_out$res %>%
      dplyr::rename(Feature = taxon) %>%
      left_join(ann_tmp, by = c("Feature"))
  }

  # Save results
  write.table(res, file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_results.txt")), row.names = F, quote = F, sep = "\t")
  write.table(bias_corrected_df, file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_bias_corrected_data.txt")), row.names = T, quote = F, sep = "\t")
}

##################################
# Plotting
##################################

for (feature in features_all) {
  dir_diff_sub <- file.path(dir_diff, feature)

  # Upload the results
  ancombc <- read.csv(file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_results.txt")), sep = "\t", check.names = F, quote = "", comment.char = "")
  df <- read.csv(file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_bias_corrected_data.txt")), sep = "\t", check.names = F, quote = "", comment.char = "")
  phy_obj <- readRDS(file.path(
    "data/processed_data", tolower(sample_type),
    paste0(ifelse(feature == "Species", "LKT", feature), "_phy.rds")
  ))
  meta <- microbiome::meta(phy_obj)

  # Remove any missing or unclassified species
  ancombc <- ancombc %>% filter(!grepl("__sp|__Missing", Feature))
  df <- df[!grepl("__sp|__Missing", rownames(df)), ]

  # Waterfall plot for significant results
  ancombc_sig <- ancombc %>%
    filter(p_statussc < 0.01, passed_ss_statussc) %>%
    arrange(lfc_statussc) %>%
    mutate(
      Feature = gsub("^[a-z]__", "", Feature),
      Feature = factor(Feature, levels = unique(Feature)),
      change = ifelse(lfc_statussc < 0, paste("Lower in", group_names[2]), paste("Higher in", group_names[2])),
      change = factor(change, levels = c(paste("Lower in", group_names[2]), paste("Higher in", group_names[2]))),
      p_sig = ifelse(q_statussc < 0.001, "***", ifelse(q_statussc < 0.01, "**", ifelse(q_statussc < 0.1, "*", ""))),
      q = round(q_statussc, 3),
      y_p = ifelse(lfc_statussc > 0, lfc_statussc + se_statussc + 0.2, lfc_statussc - se_statussc - 0.2)
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
        dplyr::select(c("Feature", "lfc_statussc", "p_statussc", "q_statussc")) %>%
        dplyr::rename(lfc = lfc_statussc, p = p_statussc, q = q_statussc) %>%
        mutate(across(grep("lfc", colnames(.)), ~ round(.x, 3))) %>%
        arrange(q)
      variable <- "Feature"
    } else {
      ancombc_save <- ancombc_sig %>%
        dplyr::select(c("Feature", "Description", "lfc_statussc", "p_statussc", "q_statussc", "Ontology", "Definition")) %>%
        dplyr::rename(lfc = lfc_statussc, p = p_statussc, q = q_statussc) %>%
        mutate(across(grep("lfc", colnames(.)), ~ round(.x, 3))) %>%
        arrange(q)
      variable <- "Description"
    }

    write.table(ancombc_save, file.path(dir_diff_sub, paste0(feature, "_ANCOMBC_results_sig.txt")), sep = "\t", row.names = F, quote = F)

    # If there are more than 25 features, only show the first 25 most significant ones
    if (nrow(ancombc_sig) > 30) {
      ancombc_sig <- ancombc_sig %>%
        # slice_head(n = 25) %>%
        arrange(lfc_statussc) %>%
        dplyr::slice(1:15, (n() - 14):n()) %>%
        mutate(
          Description = factor(Description, levels = unique(Description))
        )
    }

    # Set the colors
    if (length(unique(ancombc_sig$change)) == 1 && unique(ancombc_sig$change) == paste("Higher in", group_names[2])) {
      cols_sub <- status_cols["sc"]
    } else if (length(unique(ancombc_sig$change)) == 1 && unique(ancombc_sig$change) == paste("Lower in", group_names[2])) {
      cols_sub <- status_cols["nc"]
    } else {
      cols_sub <- status_cols
    }

    # Waterfall plot
    pdf(file.path(dir_diff_sub, paste0(feature, "_wf_fig.pdf")), 12, ifelse(nrow(ancombc_sig) > 4, 8, 4))

    p1 <- ggplot(ancombc_sig, aes(y = .data[[variable]], x = .data[["lfc_statussc"]], fill = .data[["change"]])) +
      geom_bar(stat = "identity") +
      geom_errorbar(
        aes(
          xmin = lfc_statussc - se_statussc,
          xmax = lfc_statussc + se_statussc
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
  }

  # Boxplots for significant species
  # Trim species names
  rownames(df) <- gsub("^[a-z]__", "", rownames(df))
  if (feature == "Species") {
    rownames(df) <- gsub("_", " ", rownames(df))
  }
  df_merged <- t(df) %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    left_join(meta)

  if (nrow(ancombc_sig) > 0) {
    boxplots <- create_boxplots(
      data = df_merged, ancomcs_results = ancombc_sig,
      group = "status", cols = status_cols, subt = "Visit 1",
      titley = "Bias corrected abudance", pval = "p_statussc", group_names = group_names,
      qval = "q_statussc", legend = F
    )

    pdf(file.path(dir_diff_sub, paste0(feature, "_boxplots.pdf")), 10, 10)
    print(ggarrange(plotlist = boxplots, ncol = 2, nrow = 2))
    dev.off()
  }
}
