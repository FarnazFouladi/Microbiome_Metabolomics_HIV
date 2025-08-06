library(isotone)
pava_estimator <- function(x, v) { # x is the vector of correlation estimates and v is their variances 
  E1 <- gpava(1:length(x), x, weights = 1 / v)$x # Increasing trend
  E2 <- gpava(1:length(rev(x)), rev(x), weights = rev(1 / v))$x # Decreasing trend
  Z1 <- (E1[length(x)] - E1[1]) / sqrt(v[length(x)] + v[1])
  Z2 <- (E2[length(x)] - E2[1]) / sqrt(v[length(x)] + v[1])
  Z <- max(Z1, Z2)
  return(list(Z, max = ifelse(Z == Z1, TRUE, FALSE)))
}

# Null distribution
make_null_distribution <- function(sample_num,seed){
  set.seed(seed)
  null <- sapply(1:length(sample_num), function(d) {
    rnorm(1000000, mean = 0, sd = 1)
  })
  
  null_s <- sweep(null,MARGIN = 2,FUN = "*",c(1/sqrt(sample_num[1]-3), 1/sqrt(sample_num[2]-3), 1/sqrt(sample_num[3]-3), 1/sqrt(sample_num[4]-3)))
  
  pava_null <- sapply(1:nrow(null_s), function(i) {
    pava_estimator(null_s[i, ], c(1/(sample_num[1]-3), 1/(sample_num[2]-3), 1/(sample_num[3]-3), 1/(sample_num[4]-3)))[[1]]
  })
  return(pava_null)
}



# Fisher's Z test for one sample
fisher_z_p <- function(R, N) {
  Z <- .5 * log((1 + R) / (1 - R))
  Z_bc <- Z - R / (2 * (N - 1))
  var <- 1 / pmax(N - 3, 1)
  colnames(var) <- colnames(R)
  rownames(var) <- rownames(R)
  se <- 1 / sqrt(pmax(N - 3, 1))
  p <- 2 * pnorm(-abs(Z_bc / se))
  return(list(Z_bc, var, p))
}


get_corr <- function(D) Hmisc::rcorr(as.matrix(D), type = "spearman")$r
get_npairs <- function(D) {
  outer(1:ncol(D), 1:ncol(D), Vectorize(function(i, j) {
    sum(complete.cases(D[, c(i, j)]))
  }))
}


trend_test <- function(data_modality1 = df.list, # List of data modality 1 for g1, g2, g3, g4
                       data_modality2 = df.go.list, # List of data modality 2 for g1, g2, g3, g4
                       cooccur = 10, 
                       output_dir, 
                       prefix) {
  # Reduce to common features across groups
  taxa <- Reduce(intersect, list(rownames(data_modality1[[1]]), rownames(data_modality1[[2]]), rownames(data_modality1[[3]]), rownames(data_modality1[[4]])))
  taxa <- taxa[!grepl("s__sp.", taxa)]
  data_modality1 <- lapply(data_modality1, function(x) t(x[taxa, ]))
  # Number of groups
  g <- length(data_modality1)
  # Number of taxa
  p <- ncol(data_modality1[[1]])

  if (!is.null(data_modality2)) {
    # Reduce to common features across groups
    features <- Reduce(intersect, list(rownames(data_modality2[[1]]), rownames(data_modality2[[2]]), rownames(data_modality2[[3]]), rownames(data_modality2[[4]])))
    data_modality2 <- lapply(data_modality2, function(x) t(x[features, ]))
    if (length(data_modality2) != g) {
      stop("The number of groups in twp data modalities should be the same!")
    }
    # Number of features
    q <- ncol(data_modality2[[1]])

    # Merge data modality 2 with taxa
    df.all <- lapply(1:g, function(x) {
      df_tmp <- merge(data_modality1[[x]], data_modality2[[x]], by = "row.names", all = TRUE) %>% column_to_rownames("Row.names")
    })
  } else {
    df.all <- data_modality1
    q <- 0
  }

  # Correlation analysis
  C <- lapply(df.all, function(d) get_npairs(d))
  R <- lapply(df.all, function(d) get_corr(d))
  fisher_test <- lapply(1:g, function(i) {
    test <- fisher_z_p(R[[i]], C[[i]])
    bias_corrected_r <- test[[1]]
    var <- test[[2]]
    bias_corrected_r[C[[i]] < cooccur] <- 0
    return(list(bias_corrected_r, var))
  })
  
  
  # Create Null distribution
  pava_null <- make_null_distribution(sample_num = sapply(df.all,nrow),seed = 123)

  # Stack into arrays: [group, row, column]
  cors <- array(NA, dim = c(g, p + q, p + q))
  vars <- array(NA, dim = c(g, p + q, p + q))

  for (i in 1:g) {
    cors[i, , ] <- as.matrix(fisher_test[[i]][[1]]) # Bias corrected coefficients
    vars[i, , ] <- as.matrix(fisher_test[[i]][[2]]) # variances
  }

  # Results will be stored in matrices
  mat_names <- list(colnames(R[[1]]), colnames(R[[1]]))
  
  # Pava Estimate
  pava_e <- matrix(NA, ncol = p + q, nrow = p + q, dimnames = mat_names)
  # Pava p-values
  pava_p <- matrix(NA, ncol = p + q, nrow = p + q, dimnames = mat_names)
  # Pava trend (increasing or decreasing)
  pava_trend <- matrix(NA, ncol = p + q, nrow = p + q, dimnames = mat_names)

  if (is.null(data_modality2)) {
    inds <- which(upper.tri(pava_p), arr.ind = TRUE)
  } else {
    inds <- as.matrix(expand.grid(row = 1:p, col = (p + 1):(p + q)))
  }


  for (i in 1:nrow(inds)) {
    D <- cors[, inds[i, 1], inds[i, 2]]
    V <- vars[, inds[i, 1], inds[i, 2]]
    E <- pava_estimator(D, V)
    pava_e[inds[i, 1], inds[i, 2]] <- E[[1]]
    pava_p[inds[i, 1], inds[i, 2]] <- sum(pava_null > E[[1]]) / length(pava_null)
    pava_trend[inds[i, 1], inds[i, 2]] <- E[[2]]
  }

  # Create a dataframe
  df_tmp <- data.frame(
    t1 = rownames(pava_p)[inds[, 1]],
    t2 = rownames(pava_p)[inds[, 2]],
    r1 = R[[1]][inds],
    r2 = R[[2]][inds],
    r3 = R[[3]][inds],
    r4 = R[[4]][inds],
    E = pava_e[inds],
    pvals = pava_p[inds],
    trend = ifelse(pava_trend[inds],"increase","decrease")
  ) %>% arrange(pvals)


  write.csv(df_tmp, file.path(output_dir, paste0(prefix, "_pava_results.csv")), quote = F, row.names = F)
  return(df_tmp)
}
