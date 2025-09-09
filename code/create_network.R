create_network <- function(mat,pval_mat,
                           taxonomies = NULL, t = "Species", datasets = NULL, phylum_levels = NULL,
                           repulsion = 0.85,
                           cexNodes = 1.5,
                           cexHubs = 2,
                           cexLabels = 4,
                           cexHubLabels = 4,
                           cexTitle = 3.8,
                           legend_cex = 3,
                           association_cex = 3,
                           outputdir = output_dir,
                           prefix = "status",
                           showLabel = TRUE) {
  diag(mat) <- 1
  net <- netConstruct(
    data = as.matrix(mat),
    dataType = "correlation",
    dissFunc = "unsigned",
    sparsMethod = "none"
  )
  # Analyze the network
  props <- netAnalyze(net, gcmHeat = FALSE, hubPar = c("degree"))


  # Set the color of nodes depending on datasets and whether taxonomies are provided
  if (!is.null(taxonomies)) {
    # Extract phyla for colors in the plot
    taxonomies_tmp <- taxonomies %>%
      as.data.frame() %>%
      distinct(!!sym(t), .keep_all = T) %>%
      filter(!!sym(t) %in% rownames(mat))
    phylum_levels_sub <- phylum_levels[phylum_levels %in% unique(taxonomies_tmp$Phylum)]
    cols <- names(phylum_levels_sub)
    nodeCol <- factor(taxonomies_tmp$Phylum, levels = unique(phylum_levels_sub))
    names(nodeCol) <- taxonomies_tmp[, t]
    p <- "taxa "
  } else if (is.null(taxonomies) & length(unique(datasets)) > 1) {
    nodeCol <- gsub(" .+", "", rownames(mat))
    names(nodeCol) <- rownames(mat)
    nodeCol[grepl("taxa ", rownames(mat))] <- datasets[1]
    nodeCol[!grepl("taxa ", rownames(mat))] <- datasets[2]
    nodeCol <- factor(nodeCol, levels = c(datasets[1], datasets[2]))
    cols <- c("#D95F02", "#7570B3")
    p <- "taxa |s__ |gut s "
  } else {
    cols <- "#D95F02"
    nodeCol <- rep(datasets[1], nrow(mat))
    names(nodeCol) <- rownames(mat)
    p <- "taxa "
  }

  # Edges (if adjusted p is between 0.1 and 0.05, replace the solid line with dashed line)
  edge_lty <- ifelse(pval_mat < 0.05, 1, 2)

  pdf(file.path(outputdir, paste0(prefix, "_network_diff.pdf")), 11, 10)
  plot(props,
    repulsion = repulsion,
    shortenLabels = "simple",
    rmSingles = "all",
    labelLength = 70,
    charToRm = p,
    nodeSizeSpread = 1,
    nodeColor = "feature",
    featVecCol = nodeCol,
    colorVec = cols,
    nodeTransp = 30,
    lty = edge_lty,
    sameClustCol = TRUE,
    labelScale = TRUE,
    nodeSize = "degree",
    cexNodes = cexNodes,
    cexHubs = cexHubs,
    cexLabels = cexLabels,
    cexHubLabels = cexHubLabels,
    negDiffCol = TRUE,
    labels = showLabel,
    # highlightHubs = FALSE,
    hubBorderCol = "gray40",
    esize = 1,
    asize = 1,
    edgeWidth = 1.5
  )

  dev.off()
}
