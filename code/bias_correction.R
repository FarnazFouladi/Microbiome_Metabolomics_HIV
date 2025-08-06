# Bias correction normalization
# From SECOM:
# Lin H, Eggesbø M, Peddada SD. Linear and nonlinear correlation estimators unveil undescribed taxa interactions in microbiome data. Nat Commun. 2022 Aug 23;13(1):4946. doi: 10.1038/s41467-022-32243-x. PMID: 35999204; PMCID: PMC9399263.

bias_correction <- function(data, assay.type = assay_name, assay_name = "counts",
                            rank = tax_level, tax_level = "Species", pseudo = 0, prv_cut = 0.5,
                            lib_cut = 0,cov_cat,cov_num) {
  source(file.path("code", "ancombc_prep.R"))
  source(file.path("code", "secom_prep.R"))
  source(file.path("code", "utils.R"))
  
  if (!is.null(assay.type)) {
    assay_name <- assay.type
  }
  if (!is.null(rank)) {
    tax_level <- rank
  }
  if (length(data) == 1) {
    tse_obj <- .tse_construct(
      data = data[[1]], assay_name = assay_name[1],
      tax_level = tax_level[1], phyloseq = NULL
    )
    abn_list <- .abn_est(
      tse = tse_obj$tse, tax_level = tse_obj$tax_level,
      assay_name = tse_obj$assay_name, pseudo = pseudo,
      prv_cut = prv_cut, lib_cut = lib_cut
    )
    s_diff_hat <- abn_list$s_diff_hat
    y_hat <- abn_list$y_hat
  } else {
    if (is.null(names(data))) {
      names(data) <- paste0("data", seq_along(data))
    }
    samp_names <- lapply(data, function(x) colnames(x))
    samp_common <- Reduce(intersect, samp_names)
    samp_txt <- sprintf(paste0(
      "Number of common samples ",
      "across datasets: ", length(samp_common)
    ))
    message(samp_txt)
    if (length(samp_common) < 10) {
      stop_txt <- paste0(
        "Insufficient common samples: ",
        "Multi-dataset computation not recommended"
      )
      stop(stop_txt)
    }
    tse_list <- lapply(seq_along(data), function(i) {
      tse_obj <- .tse_construct(
        data = data[[i]], assay_name = assay_name[i],
        tax_level = tax_level[i], phyloseq = NULL
      )
      return(tse_obj)
    })
    for (i in seq_along(tse_list)) {
      rownames(SingleCellExperiment::altExp(
        tse_list[[i]]$tse,
        tse_list[[i]]$tax_level
      )) <- paste(names(data)[[i]],
                  rownames(SingleCellExperiment::altExp(
                    tse_list[[i]]$tse,
                    tse_list[[i]]$tax_level
                  )),
                  sep = " - "
      )
    }
    abn_list <- lapply(seq_along(tse_list), function(i) {
      .abn_est(
        tse = tse_list[[i]]$tse, tax_level = tse_list[[i]]$tax_level,
        assay_name = assay_name[i], pseudo = pseudo,
        prv_cut = prv_cut, lib_cut = lib_cut
      )
    })
    s_diff_hat <- lapply(abn_list, function(x) x$s_diff_hat)
    y_hat <- do.call(gtools::smartbind, lapply(abn_list, function(x) as.data.frame(x$y_hat)))
    y_hat_rownames <- do.call(c, lapply(abn_list, function(x) rownames(x$y_hat)))
    y_hat <- as.matrix(y_hat)
    rownames(y_hat) <- y_hat_rownames
  }
  
  ##### Adjust for covariates
  if(!is.null(cov_cat) & !is.null(cov_num)){
    
    # Extract the meta data including the covariates
    if(class(data[[1]]) == "TreeSummarizedExperiment"){
      
      meta_list <- lapply(seq_along(data), function(i) colData(data[[i]]) %>% 
                            as.data.frame() %>% 
                            rownames_to_column("sample_names") %>%
                            dplyr::select(c("sample_names",cov_cat,cov_num)))
    }else{
      meta_list <- lapply(seq_along(data), function(i) microbiome::meta(data[[i]]) %>% 
                            rownames_to_column("sample_names") %>%
                            dplyr::select(c("sample_names",cov_cat,cov_num)))
    }
    
    # Convert covariates to character and numeric 
    meta_list <- lapply(seq_along(meta_list), function(i) {
      meta_tmp <- meta_list[[i]]
      meta_tmp[,cov_cat] <- as.data.frame(apply(meta_tmp[,cov_cat, drop = F],2,as.character))
      meta_tmp[,cov_num] <- as.data.frame(apply(meta_tmp[,cov_num,drop = F],2,as.numeric))
      stopifnot(!any(sapply(meta_tmp[,cov_cat,drop = FALSE],class) == "numeric"))
      stopifnot(all(sapply(meta_tmp[,cov_num,drop = FALSE],class) == "numeric"))
      return(meta_tmp)
    })
    
    # If there are more than one metadata, join them
    if(length(meta_list)>1){
      meta = Reduce(full_join, meta_list)
    }else{
      meta <- meta_list[[1]]
    }
    
    meta <- meta %>% dplyr::filter(sample_names %in% colnames(y_hat)) %>% arrange(order(match(colnames(y_hat), sample_names)))
    stopifnot(all.equal(meta$sample_names,colnames(y_hat)))
    
    cov_all <- c(cov_cat,cov_num)
    
    # Get the residuals
    res_list = lapply(seq_len(nrow(y_hat)), function(i) {
      df = data.frame(y = y_hat[i, ], meta)
      cov_all_pass <- cov_all[apply(na.omit(df)[,cov_all],2,function(x) length(unique(x)))!=1]
      fit <- stats::lm(as.formula(paste0("y ~",paste(cov_all_pass,collapse = "+"))), data = df)
      er <- resid(fit)
      return(er)
    })
    
    res_df <- bind_rows(res_list) %>% as.data.frame()
    rownames(res_df) <- rownames(y_hat)
    
  } else{
    res_df <- y_hat
  }
  
  
  res <- c(list(s_diff_hat = s_diff_hat, y_hat = y_hat, residuals = res_df))
  return(res)
}
