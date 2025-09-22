create_heatmap <- function(list_res, p_alpha = 0.01, bh_alpha = 0.1, top_n, min_corr = 0.2, col_names, sample_type, output_dir, feature = "Species") {
  # Extract p-values
  list_res_p <- lapply(1:length(list_res), function(x) {
    list_tmp_p <- list_res[[x]] %>%
      dplyr::select(c("species", , "p"))
    colnames(list_tmp_p)[2] <- names(list_res)[x]

    return(list_tmp_p)
  })

  # Extract bh
  list_res_bh <- lapply(1:length(list_res), function(x) {
    list_tmp_p <- list_res[[x]] %>%
      dplyr::select(c("species", , "BH_p"))
    colnames(list_tmp_p)[2] <- names(list_res)[x]

    return(list_tmp_p)
  })

  # Extract number of correlation
  list_res_cor <- lapply(1:length(list_res), function(x) {
    list_tmp_p <- list_res[[x]] %>%
      dplyr::select(c("species", "num_corr"))
    colnames(list_tmp_p)[2] <- names(list_res)[x]

    return(list_tmp_p)
  })


  df_p <- Reduce(function(x, y) full_join(x, y, by = c("species")), list_res_p) %>% column_to_rownames("species")
  df_bh <- Reduce(function(x, y) full_join(x, y, by = c("species")), list_res_bh) %>% column_to_rownames("species")
  df_cor <- Reduce(function(x, y) full_join(x, y, by = c("species")), list_res_cor) %>% column_to_rownames("species")
  df_p[is.na(df_p)] <- 1
  df_bh[is.na(df_bh)] <- 1
  df_cor[is.na(df_cor)] <- 0

  # Filter out those correlations less than colMaxs(as.matrix(df_cor)) * min_corr
  df_cor_keep <- sweep(df_cor, 2, FUN = ">", colMaxs(as.matrix(df_cor)) * min_corr)
  # Filter out those correlations with p >= p_alpha
  df_p_keep <- df_p < p_alpha
  # Filter out those correlations with bh_p >= bh_alpha
  df_bh_keep <- df_bh < bh_alpha
  taxa_keep <- df_cor_keep & df_p_keep & df_bh_keep
  df_p_sig <- df_p[rowSums(taxa_keep) > 0, ]
  df_bh_sig <- df_bh[rownames(df_p_sig), ]
  df_cor_keep <- df_cor_keep[rownames(df_p_sig), ]

  # For each data modality select top [top_n]
  top_taxa <- lapply(1:ncol(df_p_sig), function(x) {
    df_p_sig_sub <- df_p_sig[df_p_sig[, x] < p_alpha & df_bh_sig[, x] < bh_alpha & df_cor_keep[, x], ]
    na.omit(rownames(df_p_sig_sub)[order(df_p_sig_sub[, x])][1:top_n])
  })

  df_p_sig <- df_p_sig[Reduce(union, top_taxa), ]
  df_bh_sig <- df_bh_sig[Reduce(union, top_taxa), ]
  df_cor_keep <- df_cor_keep[Reduce(union, top_taxa), ]


  df_p_sig_log <- -log10(df_p_sig)
  rownames(df_p_sig_log) <- gsub("s__|g__|f__|c__", "", rownames(df_p_sig_log))
  rownames(df_p_sig_log) <- gsub("_", " ", rownames(df_p_sig_log))


  library(circlize)
  library(ComplexHeatmap)
  col_fun <- colorRamp2(c(min(df_p_sig_log), max(df_p_sig_log)), c("white", "red"))
  # Save source data
  write.table(df_p_sig_log, 
              file.path(output_dir, paste0(sample_type, "_source_data.txt")), sep = "\t", row.names = T, quote = F)
  
  hm <- Heatmap(as.matrix(df_p_sig_log),
    row_names_gp = gpar(fontsize = 40, fontface = "italic"),
    column_names_gp = gpar(fontsize = 40),
    row_names_side = "left",
    row_dend_side = "right",
    column_title = paste(sample_type, feature),
    column_title_gp = gpar(fontsize = 90),
    heatmap_height = unit(100, "cm"),
    heatmap_width = unit(45, "cm"),
    column_labels = col_names,
    border = TRUE,
    row_names_max_width = max_text_width(row.names(df_p_sig_log), gp = gpar(fontsize = 20)),
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (df_p_sig[i, j] < p_alpha & df_cor_keep[i, j] & df_bh_sig[i, j] < bh_alpha / 100) {
        grid.text("***", x, y, gp = gpar(fontsize = 35))
      } else if (df_p_sig[i, j] < p_alpha & df_cor_keep[i, j] & df_bh_sig[i, j] < bh_alpha / 10) {
        grid.text("**", x, y, gp = gpar(fontsize = 35))
      } else if (df_p_sig[i, j] < p_alpha & df_cor_keep[i, j] & df_bh_sig[i, j] < bh_alpha/2) {
        grid.text("*", x, y, gp = gpar(fontsize = 35))
      } else if (df_p_sig[i, j] < p_alpha & df_cor_keep[i, j] & df_bh_sig[i, j] < bh_alpha) {
        grid.text("+", x, y, gp = gpar(fontsize = 35))
      }
    },
    cluster_columns = FALSE, cluster_rows = TRUE, col = col_fun,
    heatmap_legend_param = list(
      legend_height = unit(5, "cm"),
      legend_width = unit(10, "cm"),
      title_gp = gpar(fontsize = 40, fontface = "bold"),
      labels_gp = gpar(fontsize = 40),
      direction = "vertical",
      title = "-log10 p-value"
    )
  )


 
  pdf(file.path(output_dir, paste0(sample_type, "_disco_heatmap_top.pdf")), 33, 50)
  draw(hm, heatmap_legend_side = "left")
  dev.off()
}
