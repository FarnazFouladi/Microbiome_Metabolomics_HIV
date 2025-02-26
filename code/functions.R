# Linear model accounting for age, location, antibiotic use
lm_test <- function(df = scf_tmp, x = "status", y = "value", plot_title = "Plasma Acetic Acid", y_label = "Log transformed", isPaired = FALSE) {
  # Remove NAs
  df <- df %>% filter(!is.na(!!sym(x)))

  fit <- summary(lm(as.formula(paste0(y, "~", x, "+ age + loc + abx_use")), data = df))

  fig_annotate <- df %>%
    rstatix::wilcox_test(as.formula(paste0(y, "~", x))) %>%
    add_xy_position() %>%
    mutate(p = fit$coefficients["statussc", "Pr(>|t|)"]) %>%
    mutate(p.format = paste0("p = ", round(p, 3)))

  plot <- ggplot(data = df, aes(x = .data[[x]], .data[[y]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha = 0.5, size = 2, position = position_jitter(width = 0.1), aes(color = .data[[x]])) +
    theme_bw(base_size = 16) +
    labs(title = plot_title, y = y_label, x = "") +
    lims(y = c(min(df[, y], na.rm = T), fig_annotate$y.position[1] + 0.1)) +
    theme(plot.title = element_text(face = "bold", size = 16), legend.position = "none") +
    scale_color_manual(values = c("#0072b5ff", "#bc3c29ff")) +
    stat_pvalue_manual(fig_annotate, label = "p.format")

  return(list(fit, plot))
}


# Functions
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  diag(cormat) <- NA
  return(cormat)
}

# Spearman correlation test
my_cor_test <- function(features = all_scfa_names, taxa = rownames(bias_corrected), df = scfa_taxa_log, fdr_method = "BH") {
  mylist <- list()
  index <- 1
  for (s in features) {
    for (t in taxa) {
      df_tmp <- df[!is.na(df[, s]) & !is.na(df[, t]), ]
      if (nrow(df_tmp) > 2) {
        fit <- cor.test(df_tmp[, s], df_tmp[, t], method = "spearman")
        mylist[[index]] <- data.frame(var1 = s, var2 = t, cor = fit$estimate, p = fit$p.value, num_obs = nrow(df_tmp))
        index <- index + 1
      } else {
        mylist[[index]] <- data.frame(var1 = s, var2 = t, cor = 0, p = 1, num_obs = nrow(df_tmp))
      }
    }
  }
  
  if (fdr_method == "BH") {
    cor_df <- mylist %>%
      bind_rows() %>%
      remove_rownames() %>%
      mutate(adjp = p.adjust(p, method = "BH"))
  } else if (fdr_method == "lfdr") {
    cor_df <- mylist %>%
      bind_rows() %>%
      remove_rownames() %>%
      mutate(adjp = fdrtool(p, statistic = "pvalue")$lfdr)
  } else {
    stop("FDR method should be either BH or lfdr")
  }
  
  return(cor_df)
}


waterfall_plot <- function(data = ancombc_sig, group,t ,subt, cols) {
  wf_fig <- ggplot(data, aes(x = .data[["Feature"]], y = .data[[group]], fill = .data[["change"]])) +
    geom_bar(stat = "identity") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10)) +
    scale_fill_manual(values = as.character(cols)) +
    geom_text(aes(y = y_p, label = p_sig)) +
    # geom_text(aes(y = y_q, label = q),size = 2) +
    labs(x = "", title = t, subtitle = subt, y = "Log Fold Change", fill = "")
}

create_boxplots <- function(data = df_merged, ancomcs_results, group = "status", 
                            cols = status_cols, subt = "v1", titley ="Normalized peak intesity",
                            group_names,
                            pval = "p_statussc", 
                            qval = "q_statussc", legend = TRUE) {
  boxplots_list <- lapply(as.character(ancomcs_results$Feature), function(a) {
    print(a)
    ancomcs_results_sub <- ancomcs_results %>% filter(Feature == a)
    stat_summary <- paste(
       p_format(ancomcs_results_sub[, pval], 3,space = TRUE,add.p = T,type = "p"), "\n",
       p_format(ancomcs_results_sub[, qval], 3, space = TRUE, add.p = T, type = "p.adj")
    )
    stat_summary<- gsub("\n p","\n q",stat_summary)

    df_tmp <- data.frame(feature =  data[,as.character(a)],group = data[,group])
    plot <- df_tmp %>% filter(!is.na(feature)) %>% ggplot(aes(y = feature, x = group )) +
      geom_boxplot(outlier.colour = NA) +
      geom_jitter(position = position_jitter(width = 0.1, height = 0), size = 2, alpha = 0.5, aes(color = group)) +
      scale_color_manual(values = cols) +
      theme_bw(base_size = 16) +
      scale_x_discrete(labels=group_names)+
      labs(title = a, caption = stat_summary, y = titley, subtitle = subt, x = "")+{if(!legend) theme(legend.position = "none")}
  })

  return(boxplots_list)
}


create_scatterplots <- function(data = df, results_sig, group, cols, path, file_name) {
  scatterplots <- lapply(1:nrow(results_sig), function(p) {
    var_p <- as.vector(c(results_sig$var1[p], results_sig$var2[p]))
    results_sig_sub <- results_sig[p, ]

    plot <- ggplot(data, aes(x = .data[[var_p[1]]], y = .data[[var_p[2]]])) +
      geom_point(alpha = 0.7, aes(color = .data[[group]])) +
      theme_classic() +
      # theme(legend.position = "none")+
      scale_color_manual(values = cols) +
      labs(title = paste("R2 =", round(results_sig_sub$cor, 3), ";", "p =", round(results_sig_sub$p, digits = 4)))
    # "; adj.p =",round(results_sig_sub$adjp,digits = 4)) )
  })

  pdf(file.path(path, paste0(file_name, "_scatterplots.pdf")), 6, 5)
  print(ggarrange(plotlist = scatterplots, ncol = 1, nrow = 1))
  dev.off()

  # return(scatterplots)
}

create_heatmap <- function(cor_results = cor_nc, row_annotation = NULL, row_ann_name = NULL, row_cluster = FALSE, row.cex = 18, col_annotation = NULL, col_ann_name = NULL, col_cluster = FALSE, row_ann_cols = NULL,
                           col_ann_cols = NULL, show_col_legend = TRUE, col.cex = 18, show_row_legend = TRUE, heatmap_width, heatmap_height, creat_annottaion = TRUE, path, file_name) {
  # Correlation matrix
  cor_mat <- cor_results

  if (creat_annottaion) {
    # order based on annotation
    col_annotation_sub <- col_annotation %>%
      filter(Feature %in% colnames(cor_mat)) %>%
      arrange(!!sym(col_ann_name))
    cor_mat <- cor_mat[, col_annotation_sub$Feature]
    row_annotation_sub <- row_annotation %>%
      filter(Feature %in% rownames(cor_mat)) %>%
      arrange(!!sym(col_ann_name))
    cor_mat <- cor_mat[row_annotation_sub$Feature, ]

    # Determine the annotation colors
    row_ann_cols <- row_ann_cols[1:length(unique(row_annotation_sub[, row_ann_name]))]
    names(row_ann_cols) <- unique(row_annotation_sub[, row_ann_name])
    col_ann_cols <- col_ann_cols[1:length(unique(col_annotation_sub[, col_ann_name]))]
    names(col_ann_cols) <- unique(col_annotation_sub[, col_ann_name])

    # Create annotation
    annotation_row <- rowAnnotation(
      row_Feature = row_annotation_sub[, row_ann_name],
      col = list(
        row_Feature = row_ann_cols
      ),
      annotation_legend_param = list(
        title_gp = gpar(
          fontsize = 20,
          fontface = "bold"
        ),
        title = row_ann_name,
        labels_gp = gpar(fontsize = 20)
      ),
      show_legend = show_col_legend
    )
    annotation_col <- HeatmapAnnotation(
      col_Feature = col_annotation_sub[, col_ann_name],
      col = list(
        col_Feature = col_ann_cols
      ),
      annotation_legend_param = list(
        title_gp = gpar(
          fontsize = 20,
          fontface = "bold"
        ),
        title = col_ann_name,
        labels_gp = gpar(fontsize = 20)
      ),
      show_legend = show_row_legend
    )

    if (col_cluster) {
      g_col <- col_annotation_sub[, col_ann_name]
      dend_col <- ComplexHeatmap::cluster_within_group(cor_mat, factor = g_col)
    } else {
      dend_col <- FALSE
    }
    if (row_cluster) {
      g_row <- row_annotation_sub[, row_ann_name]
      dend_row <- ComplexHeatmap::cluster_within_group(t(cor_mat), factor = g_row)
    } else {
      dend_row <- FALSE
    }


    cor_mat <- as.matrix(cor_mat)
    col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    hm1 <- Heatmap(cor_mat,
      row_names_gp = gpar(fontsize = row.cex),
      column_names_gp = gpar(fontsize = col.cex),
      right_annotation = annotation_row, bottom_annotation = annotation_col,
      heatmap_width = unit(heatmap_width, "cm"),
      heatmap_height = unit(heatmap_height, "cm"),
      row_names_max_width = max_text_width(row.names(cor_mat), gp = gpar(fontsize = 16)),
      cluster_columns = dend_col, cluster_rows = dend_row, col = col_fun,
      show_column_names = T, name = "Spearman Coefficients", na_col = "black",
      heatmap_legend_param = list(legend_height = unit(4, "cm"), legend_width = unit(4, "cm"), title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20))
    )
  } else {
    cor_mat <- as.matrix(cor_mat)
    col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    hm1 <- Heatmap(cor_mat,
      row_names_gp = gpar(fontsize = row.cex),
      column_names_gp = gpar(fontsize = col.cex),
      heatmap_width = unit(heatmap_width, "cm"),
      heatmap_height = unit(heatmap_height, "cm"),
      row_names_max_width = max_text_width(row.names(cor_mat), gp = gpar(fontsize = 16)),
      cluster_columns = col_cluster, cluster_rows = row_cluster, col = col_fun,
      show_column_names = T, name = "Spearman Coefficients", na_col = "black",
      heatmap_legend_param = list(legend_height = unit(4, "cm"), legend_width = unit(4, "cm"), title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20))
    )
  }


  pdf(file.path(path, paste0(file_name, "_heatmap.pdf")),
    width = convertX(ComplexHeatmap:::width(draw(hm1)), "inch", valueOnly = TRUE),
    height = 1.2 * convertY(ComplexHeatmap:::height(draw(hm1)), "inch", valueOnly = TRUE)
  )
  if (creat_annottaion) {
    draw(hm1, heatmap_legend_side = "left", annotation_legend_side = "left")
  } else {
    draw(hm1)
  }

  dev.off()

  return(draw(hm1))
}


create_differential_scatterplot <- function(sig_diff, df_1, df_2, name_1, name_2, cols, output_dir, prefix, ycol = NULL, xcol = NULL) {
  plots <- lapply(1:nrow(sig_diff), function(x) {
    print(x)
    t1 <- sig_diff$Taxa[x]
    t2 <- sig_diff$name[x]


    taxa_merged_1 <- data.frame(tax1 = as.numeric(df_1[t1, ]), tax2 = as.numeric(df_1[t2, ]))
    taxa_merged_2 <- data.frame(tax1 = as.numeric(df_2[t1, ]), tax2 = as.numeric(df_2[t2, ]))

    p1 <- ggplot(taxa_merged_1, aes(x = tax2, y = tax1)) +
      geom_point(color = cols[1]) +
      geom_smooth(method = "lm", color = "darkgray", alpha = 0.1) +
      theme_classic(base_size = 14) +
      labs(
        caption = paste0("r = ", signif(sig_diff$cor_g1[x], 3)), x = gsub("s__|g__", "", t2), y = gsub("s__|g__", "", t1),
        title = name_1
      )

    p2 <- ggplot(taxa_merged_2, aes(x = tax2, y = tax1)) +
      geom_point(color = cols[2]) +
      geom_smooth(method = "lm", color = "darkgray", alpha = 0.1) +
      theme_classic(base_size = 14) +
      labs(
        caption = paste0("r = ", signif(sig_diff$cor_g2[x], 3)), x = gsub("s__|g__", "", t2), y = gsub("s__|g__", "", t1),
        title = name_2
      )


    if (!is.null(ycol)) {
      ylab <- sig_diff[, ycol][x]
    } else {
      ylab <- t1
    }


    if (!is.null(xcol)) {
      xlab <- sig_diff[, xcol][x]
    } else {
      xlab <- t2
    }

    p1 <- p1 + labs(x = xlab, y = ylab)
    p2 <- p2 + labs(x = xlab, y = ylab)


    p <- ggarrange(p1, p2, nrow = 1)
    p <- annotate_figure(p, top = text_grob(paste0("Fisher's z-Test\np = ", signif(sig_diff$p[x], 3)),
      color = "darkgreen", face = "bold", size = 14
    ))
  })

  # Save the plots
  pdf(file.path(output_dir, paste0(prefix, "_differential_scatterplots.pdf")), 11, 5)
  print(ggarrange(plotlist = plots, nrow = 1, ncol = 1))
  dev.off()
}




# Bias correction normalization (Extracted from ANCOMBC)
bias_correction <- function(data, assay.type = assay_name, assay_name = "counts",
                            rank = tax_level, tax_level ="Species", pseudo = 0, prv_cut = 0.5,
                            lib_cut = 1000){

  source("code/ANCOMBC_mod_covariates/R/ancombc_prep.R")
  source("code/ANCOMBC_mod_covariates/R/secom_prep.R")
  source("code/ANCOMBC_mod_covariates/R/utils.R")

  if (!is.null(assay.type)) {
    assay_name = assay.type
  }
  if (!is.null(rank)) {
    tax_level = rank
  }
  if (length(data) == 1) {
    tse_obj = .tse_construct(data = data[[1]], assay_name = assay_name[1],
                             tax_level = tax_level[1], phyloseq = NULL)
    abn_list = .abn_est(tse = tse_obj$tse, tax_level = tse_obj$tax_level,
                        assay_name = tse_obj$assay_name, pseudo = pseudo,
                        prv_cut = prv_cut, lib_cut = lib_cut)
    s_diff_hat = abn_list$s_diff_hat
    y_hat = abn_list$y_hat
  }
  else {
    if (is.null(names(data)))
      names(data) = paste0("data", seq_along(data))
    samp_names = lapply(data, function(x) colnames(x))
    samp_common = Reduce(intersect, samp_names)
    samp_txt = sprintf(paste0("Number of common samples ",
                              "across datasets: ", length(samp_common)))
    message(samp_txt)
    if (length(samp_common) < 10) {
      stop_txt = paste0("Insufficient common samples: ",
                        "Multi-dataset computation not recommended")
      stop(stop_txt)
    }
    tse_list = lapply(seq_along(data), function(i) {
      tse_obj = .tse_construct(data = data[[i]], assay_name = assay_name[i],
                               tax_level = tax_level[i], phyloseq = NULL)
      return(tse_obj)
    })
    for (i in seq_along(tse_list)) {
      rownames(SingleCellExperiment::altExp(tse_list[[i]]$tse,
                                            tse_list[[i]]$tax_level)) = paste(names(data)[[i]],
                                                                              rownames(SingleCellExperiment::altExp(tse_list[[i]]$tse,
                                                                                                                    tse_list[[i]]$tax_level)), sep = " - ")
    }
    abn_list = lapply(seq_along(tse_list), function(i) {
      .abn_est(tse = tse_list[[i]]$tse, tax_level = tse_list[[i]]$tax_level,
               assay_name = assay_name[i], pseudo = pseudo,
               prv_cut = prv_cut, lib_cut = lib_cut)
    })
    s_diff_hat = lapply(abn_list, function(x) x$s_diff_hat)
    y_hat = do.call(gtools::smartbind, lapply(abn_list, function(x) as.data.frame(x$y_hat)))
    y_hat_rownames = do.call(c, lapply(abn_list, function(x) rownames(x$y_hat)))
    y_hat = as.matrix(y_hat)
    rownames(y_hat) = y_hat_rownames
  }

  res = c(list(s_diff_hat = s_diff_hat, y_hat = y_hat))
  return(res)
}

