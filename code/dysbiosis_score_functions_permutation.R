get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  diag(cormat) <- NA
  return(cormat)
}


dysbiosis_score <- function(data_g1 = data_list_nc_tse[c(1, 3)],
                            data_g2 = data_list_sc_tse[c(1, 3)],
                            assay_name, tax_level = c(NULL, NULL), prv_cut, method, cov_cat, cov_num, threshold = 0.3, pval_t, name = "gut_mtb_pos", two_datasets = TRUE, output_dir, perm) {
  # Create directory
  dir.create(output_dir, showWarnings = F, recursive = T)

  # Correlation
  res_secom_g1 <- secom_linear(
    data = data_g1,
    assay_name = assay_name,
    tax_level = tax_level,
    prv_cut = prv_cut,
    method = method, cov_cat = cov_cat, cov_num = cov_num
  )


  res_secom_g2 <- secom_linear(
    data = data_g2,
    assay_name = assay_name,
    tax_level = tax_level,
    prv_cut = prv_cut,
    method = method, cov_cat = cov_cat, cov_num = cov_num
  )

  # Differential correlations
  g1_mat <- res_secom_g1$corr
  g1_cooccur <- res_secom_g1$mat_cooccur

  g2_mat <- res_secom_g2$corr
  g2_cooccur <- res_secom_g2$mat_cooccur


  # Subset to common taxa
  comm <- intersect(colnames(g2_mat), colnames(g1_mat))
  # Remove missing or unclassified
  comm <- comm[!grepl("__Missing|__sp.", comm)]
  g2_mat <- g2_mat[comm, comm]
  g1_mat <- g1_mat[comm, comm]

  stopifnot(all.equal(colnames(g2_mat), colnames(g1_mat)))


  g1_cooccur <- g1_cooccur[rownames(g1_mat), colnames(g1_mat)]
  g2_cooccur <- g2_cooccur[rownames(g2_mat), colnames(g2_mat)]

  # Replace correlations within data1 features and  within data2 features to zero
  g1_mat[grepl("data2 ", rownames(g1_mat)), grepl("data2 ", colnames(g1_mat))] <- 0
  g1_mat[grepl("data1 ", rownames(g1_mat)), grepl("data1 ", colnames(g1_mat))] <- 0
  g2_mat[grepl("data2 ", rownames(g2_mat)), grepl("data2 ", colnames(g2_mat))] <- 0
  g2_mat[grepl("data1 ", rownames(g2_mat)), grepl("data1 ", colnames(g2_mat))] <- 0

  # Keep only pairs that at least in one group the coefficient > threshold.
  g1_mat_org <- g1_mat
  g2_mat_org <- g2_mat
  g1_mat[abs(g1_mat_org) < threshold & abs(g2_mat_org) < threshold] <- 0
  g2_mat[abs(g1_mat_org) < threshold & abs(g2_mat_org) < threshold] <- 0

  # Get all possible pairwise correlations
  df_num <- get_upper_tri(g1_mat) %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    pivot_longer(cols = -Taxa, values_to = "cor") %>%
    filter(!is.na(cor))

  # Differential association
  res <- fisher_test(mat1 = g1_mat, mat2 = g2_mat, mat1_cooccur = g1_cooccur, mat2_cooccur = g2_cooccur)
  res_sig <- res %>% filter(p < pval_t & diff_cor > threshold, !is.na(p))
  write.csv(res, file.path(output_dir, paste0(name, "_differential_correlations_", perm, ".csv")), row.names = F, quote = F)

  n <- nrow(res_sig)

  if (two_datasets) {
    # Total number of features (fm)
    res_fm <- df_num %>%
      filter(grepl("data1", Taxa) & grepl("data2", name)) %>%
      group_by(Taxa) %>%
      summarise(num = n()) %>%
      mutate(Taxa = gsub("data1 - |data2 - ", "", Taxa))
    colnames(res_fm)[2] <- name


    # Number of correlations that are higher in SC
    res_sc_p <- res %>%
      filter(p < pval_t & diff_cor > threshold, !is.na(p)) %>%
      filter(sign == paste0("SC", "+") | sign == paste0("SC", "-")) %>%
      group_by(Taxa) %>%
      summarise(cor2 = n()) %>%
      mutate(Taxa = gsub("data1 - |data2 - ", "", Taxa))

    # Number of correlations that are higher in NC
    res_nc_p <- res %>%
      filter(p < pval_t & diff_cor > threshold, !is.na(p)) %>%
      filter(sign == paste0("NC", "+") | sign == paste0("NC", "-")) %>%
      group_by(Taxa) %>%
      summarise(cor1 = n()) %>%
      mutate(Taxa = gsub("data1 - |data2 - ", "", Taxa))

    res_all_p <- res_sc_p %>%
      full_join(res_nc_p) %>%
      pivot_longer(cols = -Taxa, names_to = "group", values_to = name)
  } else {
    taxa <- unique(c(res$Taxa, res$name))

    res_fm <- lapply(taxa, function(t) {
      # Total number of features (fm)
      fm <- df_num %>%
        filter(grepl(t, Taxa) | grepl(t, name)) %>%
        summarise(num = n())

      return(data.frame(Taxa = t, fm = fm$num))
    })
    res_fm <- res_fm %>% bind_rows()
    colnames(res_fm)[2] <- name

    res_all_p <- lapply(taxa, function(t) {
      # Number of correlations that are higher in SC
      res_sc <- res %>%
        filter(p < pval_t & diff_cor > threshold, !is.na(p)) %>%
        filter(sign == paste0("SC", "+") | sign == paste0("SC", "-")) %>%
        filter(grepl(t, Taxa) | grepl(t, name)) %>%
        summarise(cor2 = n())
      # Number of correlations that are higher in NC
      res_nc <- res %>%
        filter(p < pval_t & diff_cor > threshold, !is.na(p)) %>%
        filter(sign == paste0("NC", "+") | sign == paste0("NC", "-")) %>%
        filter(grepl(t, Taxa) | grepl(t, name)) %>%
        summarise(cor1 = n())

      return(data.frame(Taxa = t, group = c("cor2", "cor1"), value = c(as.numeric(res_sc), as.numeric(res_nc))))
    })

    res_all_p <- res_all_p %>% bind_rows()
    colnames(res_all_p)[3] <- name
  }

  res_all_p[is.na(res_all_p)] <- 0

  return(list(res_all_p, res_fm))
}




fisher_test <- function(mat1 = g1_mat, mat2 = g2_mat, mat1_cooccur = g1_mat_cooccur, mat2_cooccur = g2_mat_cooccur) {
  stopifnot(all.equal(colnames(mat1), colnames(mat2_cooccur)))
  stopifnot(all.equal(colnames(mat1), colnames(mat1_cooccur)))


  df1 <- get_upper_tri(mat1) %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    pivot_longer(cols = -Taxa, values_to = "cor_1") %>%
    filter(!is.na(cor_1))
  df2 <- get_upper_tri(mat2) %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    pivot_longer(cols = -Taxa, values_to = "cor_2") %>%
    filter(!is.na(cor_2))


  co1 <- get_upper_tri(mat1_cooccur) %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    pivot_longer(cols = -Taxa, values_to = "cooccur_1") %>%
    filter(!is.na(cooccur_1))
  co2 <- get_upper_tri(mat2_cooccur) %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    pivot_longer(cols = -Taxa, values_to = "cooccur_2") %>%
    filter(!is.na(cooccur_2))

  df <- df1 %>%
    full_join(df2) %>%
    full_join(co1) %>%
    full_join(co2)

  df <- df %>% filter(!(cor_1 == 0 & cor_2 == 0))


  res <- list()

  res <- lapply(1:nrow(df), function(i) {
    z1 <- atanh(df[i, ]$cor_1)
    z2 <- atanh(df[i, ]$cor_2)
    n1 <- df[i, ]$cooccur_1
    n2 <- df[i, ]$cooccur_2

    diff_z <- (z1 - z2) / sqrt(1 / (n1 - 3) + (1 / (n2 - 3)))
    if (is.na(diff_z)) {
      p_val <- 1
    } else {
      p_val <- 2 * (1 - pnorm(abs(diff_z)))
    }

    df_tmp <- data.frame(
      Taxa = df[i, ]$Taxa, name = df[i, ]$name,
      cooccurrence_g1 = df$cooccur_1[i],
      cooccurrence_g2 = df$cooccur_2[i],
      cor_g1 = df$cor_1[i],
      cor_g2 = df$cor_2[i],
      diff_cor = abs(df[i, ]$cor_1 - df[i, ]$cor_2), p = p_val,
      sign = ifelse(abs(df$cor_1[i]) > abs(df$cor_2[i]) & df$cor_1[i] > 0, "NC+", ifelse(
        abs(df$cor_1[i]) > abs(df$cor_2[i]) & df$cor_1[i] < 0, "NC-", ifelse(
          abs(df$cor_1[i]) < abs(df$cor_2[i]) & df$cor_2[i] > 0, "SC+", "SC-"
        )
      ))
    )
  })

  res_all <- bind_rows(res)
  return(res_all)
}
