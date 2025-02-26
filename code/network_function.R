get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  diag(cormat) <- NA
  return(cormat)
}


create_networks <- function(g1_mat, g1_mat_cooccur, g2_mat, g2_mat_cooccur, sig_diff,
                            taxonomies = NULL, t = "Species", datasets = NULL, phylum_levels = NULL,
                            threshold = 0.3, pval_th = 0.05,
                            layoutGroup = "union", layoutGroup2 = "union",
                            repulsion = 0.85, repulsion2 = 0.85,
                            cexNodes = 1.5, cexNodes2 = 1.5,
                            cexHubs = 2, cexHubs2 = 2,
                            cexLabels = 4, cexLabels2 = 4,
                            cexHubLabels = 4, cexHubLabels2 = 4,
                            cexTitle = 3.8, cexTitle2 = 3.8,
                            legend_cex = 3, legend_cex2 = 3,
                            association_cex = 3, association_cex2 = 3,
                            group_names = c("NC", "SC"), outputdir = output_dir, prefix = "status",
                            showLabel = FALSE) {
  # For network only show common taxa
  comm <- intersect(colnames(g2_mat), colnames(g1_mat))
  # Remove missing or unclassified
  comm <- comm[!grepl("__Missing|__Unclassified", comm)]
  g2_mat <- g2_mat[comm, comm]
  g1_mat <- g1_mat[comm, comm]

  stopifnot(all.equal(colnames(g2_mat), colnames(g1_mat)))

  # Replace correlations within dataset1 and dataset2  to zero 
  # We want to visualize correlations between dataset1 and dataset2 only.
  g1_mat[grepl("data2 ", rownames(g1_mat)), grepl("data2 ", colnames(g1_mat))] <- 0
  g1_mat[grepl("data1 ", rownames(g1_mat)), grepl("data1 ", colnames(g1_mat))] <- 0
  g2_mat[grepl("data2 ", rownames(g2_mat)), grepl("data2 ", colnames(g2_mat))] <- 0
  g2_mat[grepl("data1 ", rownames(g2_mat)), grepl("data1 ", colnames(g2_mat))] <- 0

  # Keep only pairs that at least in one group the coefficient > threshold.
  g1_mat_org <- g1_mat
  g2_mat_org <- g2_mat
  g1_mat[abs(g1_mat_org) < threshold & abs(g2_mat_org) < threshold] <- 0
  g2_mat[abs(g1_mat_org) < threshold & abs(g2_mat_org) < threshold] <- 0

  g1_mat_cooccur <- g1_mat_cooccur[colnames(g1_mat), colnames(g1_mat)]
  g2_mat_cooccur <- g2_mat_cooccur[colnames(g2_mat), colnames(g2_mat)]


  # Set the color of nodes depending on datasets and whether taxonomies are provided
  if (!is.null(taxonomies)) {
    # Extract phyla for colors in the plot
    taxonomies_tmp <- taxonomies %>%
      as.data.frame() %>%
      distinct(!!sym(t), .keep_all = T) %>%
      filter(!!sym(t) %in% rownames(g1_mat))
    phylum_levels_sub <- phylum_levels[phylum_levels %in% unique(taxonomies_tmp$Phylum)]
    cols <- names(phylum_levels_sub)
    nodeCol <- factor(taxonomies_tmp$Phylum, levels = unique(phylum_levels_sub))
    names(nodeCol) <- taxonomies_tmp[, t]
    p <- "s__"
  } else if (is.null(taxonomies) & length(unique(datasets))>1) {
    nodeCol <- gsub(" .+", "", rownames(g2_mat))
    names(nodeCol) <- rownames(g2_mat)
    nodeCol[grepl("data1",rownames(g2_mat))] <- datasets[1]
    nodeCol[grepl("data2",rownames(g2_mat))] <- datasets[2]
    nodeCol <- factor(nodeCol, levels = c(datasets[1], datasets[2]))
    cols <- c("#D95F02", "#7570B3")
    p <-  "data[1-2] - [a-z]__|data[1-2] - "
  } else{
    cols <- "#D95F02"
    nodeCol <- rep(datasets[1],nrow(g2_mat))
    names(nodeCol) <- rownames(g2_mat)
    
    p <- "s__"
    
  }

    # Create output directory
    if (!dir.exists(outputdir)) {
      dir.create(outputdir, showWarnings = F, recursive = T)
    }

  ########################
  # Networks of correlations
  ########################

    # Construct a network
    # net <- netConstruct(
    #   data = g1_mat,
    #   data2 = g2_mat,
    #   dataType = "correlation",
    #   dissFunc = "unsigned", sparsMethod = "none"
    # )
    # # Analyze the network
    # props <- netAnalyze(net,
    #   centrLCC = TRUE,
    #   normDeg = TRUE,
    #   clustMethod = "cluster_fast_greedy",
    #   hubPar = c("degree"), gcmHeat = FALSE
    # )

    # pdf(file.path(outputdir, paste0(prefix, "_network.pdf")), 60, 40)
    # 
    # plot(props,
    #   sameLayout = TRUE,
    #   repulsion = repulsion,
    #   shortenLabels = "simple",
    #   labelLength = 45,
    #   charToRm = "[a-z]__",
    #   layoutGroup = "union",
    #   rmSingles = "inboth",
    #   nodeSizeSpread = 1,
    #   nodeColor = "feature",
    #   featVecCol = nodeCol,
    #   colorVec = cols,
    #   nodeTransp = 30,
    #   sameClustCol = TRUE,
    #   labelScale = TRUE,
    #   # labelFont=0.8,
    #   nodeSize = "degree",
    #   cexNodes = cexNodes,
    #   cexHubs = cexHubs,
    #   cexLabels = cexLabels,
    #   cexHubLabels = cexHubLabels,
    #   cexTitle = cexTitle,
    #   edgeFilter = "threshold",
    #   edgeFilterPar = threshold,
    #   negDiffCol = TRUE,
    #   edgeTranspLow = 50,
    #   edgeTranspHigh = 50,
    #   groupNames = group_names,
    #   hubBorderCol = "gray40"
    # )
    # 
    # # Colors used in the legend should be equally transparent as in the plot
    # phylcol_transp <- colToTransp(cols, 30)
    # 
    # legend(-1.2, 1.2,
    #   cex = legend_cex, pt.cex = 3, title = "Phylum:",
    #   legend = unique(nodeCol), col = phylcol_transp, bty = "n", pch = 16
    # )
    # 
    # legend("bottom",
    #   title = "estimated association:", legend = c("+", "-"),
    #   col = c("#009900", "red"), inset = 0.02, cex = association_cex, lty = 1, lwd = 4,
    #   bty = "n", horiz = TRUE
    # )
    # 
    # 
    # dev.off()
    # 

  #######################################
  # Networks of differential correlations
  #######################################
    
    if (nrow(sig_diff) > 0) {
      g1_mat_n <- matrix(0, nrow(g1_mat), ncol(g1_mat), dimnames = list(rownames(g1_mat), colnames(g1_mat)))
      g2_mat_n <- matrix(0, nrow(g2_mat), ncol(g2_mat), dimnames = list(rownames(g2_mat), colnames(g2_mat)))

      # Only keep pairs with significant differential correlations
      for (i in 1:nrow(sig_diff)) {
        p1 <- sig_diff[i, ]$Taxa
        p2 <- sig_diff[i, ]$name

        g1_mat_n[p1, p2] <- g1_mat[p1, p2]
        g2_mat_n[p1, p2] <- g2_mat[p1, p2]

        g1_mat_n[p2, p1] <- g1_mat[p2, p1]
        g2_mat_n[p2, p1] <- g2_mat[p2, p1]
      }

      # Construct a network
      net <- netConstruct(
        data = g1_mat_n,
        data2 = g2_mat_n,
        dataType = "correlation",
        dissFunc = "unsigned", sparsMethod = "none"
      )
      # Analyze the network
      props <- netAnalyze(net,
        centrLCC = FALSE,
        clustMethod = "cluster_fast_greedy",
        normDeg = TRUE,
        hubPar = c("degree"), gcmHeat = FALSE
      )

      pdf(file.path(outputdir, paste0(prefix, "_network_diff.pdf")), 70, 50)

      plot(props,
        sameLayout = TRUE,
        labels = showLabel,
        repulsion = repulsion2,
        shortenLabels = "simple",
        labelLength = 40,
        charToRm = p,
        layoutGroup = layoutGroup2,
        rmSingles = "inboth",
        nodeSizeSpread = 1,
        nodeColor = "feature",
        featVecCol = nodeCol,
        colorVec = cols,
        labelScale = TRUE,
        nodeTransp = 30,
        # labelFont=0.8,
        # nodeSize = "degree",
        cexNodes = cexNodes2,
        cexHubs = cexHubs2,
        cexLabels = cexLabels2,
        cexHubLabels = cexHubLabels2,
        cexTitle = cexTitle2,
        negDiffCol = TRUE,
        edgeTranspLow = 50,
        edgeTranspHigh = 50,
        groupNames = group_names,
        hubBorderCol = "gray40"
      )

      # # Colors used in the legend should be equally transparent as in the plot
      # phylcol_transp <- colToTransp(cols, 30)
      # 
      # legend(-1.2, 1.2,
      #   cex = legend_cex2, pt.cex = 4, title = "Phylum:",
      #   legend = unique(nodeCol), col = phylcol_transp, bty = "n", pch = 16
      # )
      # 
      # legend("bottom",
      #   title = "estimated association:", legend = c("+", "-"),
      #   col = c("#009900", "red"), inset = 0.02, cex = association_cex2, lty = 1, lwd = 4,
      #   bty = "n", horiz = TRUE
      # )


      dev.off()
    }
}
