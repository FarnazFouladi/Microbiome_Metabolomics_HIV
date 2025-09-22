# Fisher's Z test for one sample
fisher_z_p <- function(R, N) {
  Z <- .5 * log((1 + R) / (1 - R))
  # Correct for the bias related to sample size
  Z_bc <- Z - R / (2 * (N - 1))
  se <- 1 / sqrt(pmax(N - 3, 1))
  2 * pnorm(-abs(Z_bc / se))
}

# Fisher's Z test for two samples
fisher_z_p_2sample <- function(R1, N1, R2, N2) {
  Z1 <- 0.5 * log((1 + R1) / (1 - R1))
  Z2 <- 0.5 * log((1 + R2) / (1 - R2))

  # Correct for different sample sizes between the two groups
  Z1_bc <- Z1 - R1 / (2 * (N1 - 1))
  Z2_bc <- Z2 - R2 / (2 * (N2 - 1))
  se <- sqrt(1 / pmax(N1 - 3, 1) + 1 / pmax(N2 - 3, 1))
  Z <- (Z1_bc - Z2_bc) / se
  pvals <- 2 * pnorm(-abs(Z))
  return(list(pvals, Z2_bc-Z1_bc))
  
}

get_corr <- function(D, C) Hmisc::rcorr(as.matrix(D), type = "spearman")$r
get_npairs <- function(D) {
  outer(1:ncol(D), 1:ncol(D), Vectorize(function(i, j) {
    sum(complete.cases(D[, c(i, j)]))
  }))
}


disco_score <- function(data1, # [List] Data for group 1 
                        data2, # [List] Data for group 2
                        cooccur = 10, 
                        filt_threshold = 0.001, 
                        diff_threshold = 0.3, 
                        p_alpha = 0.01, 
                        bh_alpha = 0.1, 
                        output_dir, 
                        prefix) {
  
  taxa <- colnames(data1[[1]])

  if (length(data1) == 2 & length(data2) == 2) {
    two_modal <- T
    p <- ncol(data1[[1]])
    q <- ncol(data1[[2]])
 
    # Concatenating data modality 1 and data modality 2 for group 1 
    d_merged1 <- merge(data1[[1]], data1[[2]], by = "row.names", all = T) %>%
      column_to_rownames("Row.names")
    # Concatenating data modality 1 and data modality 2 for group 2
    d_merged2 <- merge(data2[[1]], data2[[2]], by = "row.names", all = T) %>%
      column_to_rownames("Row.names")
    
  } else if (length(data1) == 1 & length(data2) == 1) {
    two_modal <- F
    p <- ncol(data1[[1]])
    q <- 0
    d_merged1 <- data1[[1]]
    d_merged2 <- data2[[1]]
  } else {
    stop("Number of data modalities should be one or two for each group")
  }

  # Estimate correlations and cooccurrences
  C1 <- get_npairs(d_merged1)
  C2 <- get_npairs(d_merged2)
  R1 <- get_corr(d_merged1)
  R2 <- get_corr(d_merged2)
  R1[C1 < cooccur] <- 0
  R2[C2 < cooccur] <- 0
  
  # Fisher's Z test
  P1 <- fisher_z_p(R1, C1)
  P2 <- fisher_z_p(R2, C2)

  # Filter correlations with p >= threshold in both groups
  filt <- P1 >= filt_threshold & P2 >= filt_threshold
  diag(filt) <- TRUE

  # Fisher Z test - two sample test
  fisher_res <- fisher_z_p_2sample(R1, C1, R2, C2)
  pvals <- fisher_res[[1]]
  delta_R <- fisher_res[[2]]
  delta_R[is.na(delta_R)] <-0
  
  # Replace p-values with NAs for correlations that are filtered out
  pvals_tmp <- pvals
  pvals_tmp[filt] <- NA

  # Create a long dataframe
  if (two_modal) {
    inds <- as.matrix(expand.grid(row = 1:p, col = (p + 1):(p + q)))
  } else {
    inds <- which(upper.tri(R1), arr.ind = TRUE)
  }

  t1 <- rownames(R1)[inds[, 1]]
  t2 <- rownames(R1)[inds[, 2]]
  r1 <- R1[inds]
  r2 <- R2[inds]
  n1 <- C1[inds]
  n2 <- C2[inds]
  pvals <- pvals_tmp[inds]
  adjp <- p.adjust(pvals,method = "BH")
  
  cat("Number of all pairwise correlations:",length(pvals),"\n")
  cat("Number of retained pairwise correlations:",sum(!is.na(pvals)),"\n")
  
  # Create a matrix for adjusted p-values
  adj_pvals <- matrix(NA, nrow = nrow(pvals_tmp), ncol = ncol(pvals_tmp))
  adj_pvals[inds] <- adjp

  df_tmp <- data.frame(t1, t2, r1, r2, n1, n2, p = pvals, BH_p = adjp) %>%
    filter(!is.na(pvals))

  # For each feature in t1 and t2 columns, find the number of significant correlations
  df_tmp_sig <- df_tmp %>%
    mutate(abs_diff_cor = abs(r1 - r2)) %>%
    filter(p < p_alpha & abs_diff_cor > diff_threshold & BH_p < bh_alpha) %>%
    group_by(t1) %>%
    mutate(num_sig_corr_t1 = n()) %>%
    ungroup() %>%
    group_by(t2) %>%
    mutate(num_sig_corr_t2 = n())

  cat("Number of significant pairwise correlations:",dim(df_tmp_sig),"\n")
  
  # Calculating disco score
  if (two_modal) {
    disco <- bind_rows(sapply(1:p, function(i) {
      p_vect <- as.vector(na.omit(pvals_tmp[(p + 1):(p + q), i]))
      t <- 1 / length(p_vect) * sum(tan((0.5 - p_vect) * pi))
      pval <- pcauchy(t, location = 0, scale = 1, lower.tail = FALSE)
      return(data.frame(species = colnames(pvals_tmp)[i], t, p = pval, num_corr = length(p_vect)))
    }, simplify = F))
  } else {
    disco <- bind_rows(sapply(1:p, function(i) {
      p_vect <- as.vector(na.omit(pvals_tmp[-i, i]))
      t <- 1 / length(p_vect) * sum(tan((0.5 - p_vect) * pi))
      pval <- pcauchy(t, location = 0, scale = 1, lower.tail = FALSE)
      return(data.frame(species = colnames(pvals_tmp)[i], t, p = pval, num_corr = length(p_vect)))
    }, simplify = F))
  }

  # Adjustment of p-values
  disco <- disco %>%
    mutate(BH_p = p.adjust(p,method = "BH")) %>%
    arrange(BH_p)
  
  # Save the results
  write.csv(df_tmp, file.path(output_dir, paste0(prefix, "_cor_diff.csv")), quote = T, row.names = F)
  write.csv(df_tmp_sig, file.path(output_dir, paste0(prefix, "_cor_diff_sig.csv")), quote = T, row.names = F)
  write.csv(disco, file.path(output_dir, paste0(prefix, "_disco.csv")), quote = T, row.names = F)
  write.csv(delta_R, file.path(output_dir, paste0(prefix, "_delta_cor.csv")), quote = T, row.names = T)
  write.csv(pvals_tmp, file.path(output_dir, paste0(prefix, "_p.csv")), quote = T, row.names = T)
  write.csv(adj_pvals, file.path(output_dir, paste0(prefix, "_bh.csv")), quote = T, row.names = T)
  
  return(list(df_tmp,df_tmp_sig, disco, delta_R, pvals_tmp, adj_pvals))
}
