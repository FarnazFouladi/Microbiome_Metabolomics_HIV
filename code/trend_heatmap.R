trend_heatmap <- function(mat, annot, output_dir, prefix) {
  
  col_fun <- colorRamp2(c(min(mat, na.rm = T), 0, max(mat, na.rm = T)), c("darkorange", "white", "purple"))
  
  # Adjusted p-values
  ha <- rowAnnotation(
    p_values = anno_text(formatC(annot$adj.p.trend, format = "f", digits = 2),
      just = "right",
      location = 0.9,
      gp = gpar(fontsize = 8)
    )
  )

  hm1 <- Heatmap(as.matrix(mat),
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 15),
    # bottom_annotation = annotation_col,
    left_annotation = ha,
    heatmap_width = unit(ifelse(prefix == "oral_go", 30, 20), "cm"),
    heatmap_height = unit(ifelse(prefix == "oral_go", 35, 20), "cm"),
    column_labels = c("G1","G2","G3","G4"),
    border = TRUE,
    row_names_max_width = max_text_width(row.names(mat), gp = gpar(fontsize = 16)),
    cluster_columns = FALSE, cluster_rows = FALSE, col = col_fun,
    heatmap_legend_param = list(
      legend_height = unit(4, "cm"),
      legend_width = unit(4, "cm"),
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 10),
      title = "Fisher Z transformed \ncorrelation coefficient"
    )
  )

  pdf(file.path(output_dir, paste0(prefix, "_heatmap.pdf")), 14, ifelse(prefix == "oral_go", 20, 10))
  draw(hm1, heatmap_legend_side = "left", annotation_legend_side = "left")
  dev.off()
}
