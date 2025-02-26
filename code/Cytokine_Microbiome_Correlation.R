rm(list = ls())
source("code/load.R")
output_dir <- file.path(output_dir, "Cytokine/Correlation_Analysis")
source("code/dysbiosis_score_functions.R")
dir.create(output_dir, recursive = T, showWarnings = F)
phy_gut <- readRDS(paste0("data/processed_data/", "gut", "/LKT_phy.rds"))
phy_oral <- readRDS(paste0("data/processed_data/", "oral", "/LKT_phy.rds"))
cytokines <- c("cd163", "il6", "ip10", "crp", "lbp", "cd14", "ratio")

# Set sample names to subjid
sample_names(phy_gut) <- sample_data(phy_gut)$subjid
sample_names(phy_oral) <- sample_data(phy_oral)$subjid

meta <- read_csv("data/metadata/metadata_c.csv")
meta$loc <- factor(meta$loc)
meta <- meta %>% mutate(cd4 = leu3p, cd8 = leu2p, ratio = log(cd4 / cd8, base = 2)) 

#####################################
# Correlation analysis
#####################################

for (s in status) {
  # Output directory
  output_dir_sub <- file.path(output_dir, s)
  dir.create(output_dir_sub, showWarnings = F, recursive = T)
  prefix <- s
  
  # Subset to s
  phy_gut_sub <- subset_samples(phy_gut, status == s)
  phy_oral_sub <- subset_samples(phy_oral, status == s)
  meta_sub <- meta %>% filter( status == s)

  # Bias correction
  gut_bias_corrected <- bias_correction(list(phy_gut_sub), prv_cut = 0.5, tax_level = "Species")$y_hat
  oral_bias_corrected <- bias_correction(list(phy_oral_sub), prv_cut = 0.5, tax_level = "Species")$y_hat

  # Merge oral and gut
  rownames(gut_bias_corrected) <- paste("gut", rownames(gut_bias_corrected))
  rownames(oral_bias_corrected) <- paste("oral", rownames(oral_bias_corrected))
  gut_bias_corrected_t <- gut_bias_corrected %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("subjid")
  oral_bias_corrected_t <- oral_bias_corrected %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("subjid")
  bias_corrected_merged <- gut_bias_corrected_t %>% full_join(oral_bias_corrected_t)

  # Merge with metadata
  bias_corrected_meta <- bias_corrected_merged %>% left_join(meta_sub)
  bias_corrected_merged <- bias_corrected_merged %>% column_to_rownames("subjid")

  # Get residuals for taxa
  res_list <- lapply(colnames(bias_corrected_merged), function(i) {
    df <- data.frame(y = bias_corrected_meta[, i], bias_corrected_meta[, c("age", "loc", "abx_use")])
    rownames(df) <- bias_corrected_meta$subjid
    fit <- stats::lm(as.formula(paste0("y ~", paste(c("age", "loc", "abx_use"), collapse = "+"))), data = df)
    er <- resid(fit)
    return(er)
  })

  res_df <- bind_rows(res_list) %>%
    t() %>%
    as.data.frame()
  colnames(res_df) <- colnames(bias_corrected_merged)

  # Get residuals for cytokines
  res_list_cyt <- lapply(cytokines, function(i) {
    df <- data.frame(y = bias_corrected_meta[, i], bias_corrected_meta[, c("age", "loc", "abx_use")])
    rownames(df) <- bias_corrected_meta$subjid
    fit <- stats::lm(as.formula(paste0("y ~", paste(c("age", "loc", "abx_use"), collapse = "+"))), data = df)
    er <- resid(fit)
    return(er)
  })

  res_df_cyt <- bind_rows(res_list_cyt) %>%
    t() %>%
    as.data.frame()
  colnames(res_df_cyt) <- cytokines

  # Merge taxa and cytokines data
  taxa_merged <- merge(res_df, res_df_cyt, by = "row.names", all = T)
  taxa_merged$status <- s

  # Spearman correlation
  res_cor <- my_cor_test(features = cytokines, taxa = colnames(bias_corrected_merged), df = taxa_merged, fdr_method = "BH")
  # Remove unclassified species
  res_cor <- res_cor %>% filter(!grepl("s__sp.", var2))

  write.table(res_cor, file.path(output_dir_sub, paste0(prefix, "_cytokines_microbiome.txt")), sep = "\t", quote = F, row.names = F)
  write.csv(taxa_merged, file.path(output_dir_sub, paste0(prefix, "_cytokines_microbiome_merged.csv")), row.names = F, quote = F)
}

#####################################
# Networks
#####################################
v <- "v1"
t <- "Species"
output_dir_sub <- file.path(output_dir, "network")
dir.create(output_dir_sub, showWarnings = F)
df_sc <- read.table(file.path(output_dir, "sc", paste0("sc_cytokines_microbiome.txt")), sep = "\t", header = T, check.names = F)
df_nc <- read.table(file.path(output_dir, "nc", paste0("nc_cytokines_microbiome.txt")), header = T, check.names = F, sep = "\t")

# Convert to wide format - SC
mat_sc <- df_sc %>%
  dplyr::select(all_of(c("var1", "var2", "cor"))) %>%
  pivot_wider(names_from = "var2", values_from = "cor") %>%
  column_to_rownames("var1") %>%
  as.matrix()
co_sc <- df_sc %>%
  dplyr::select(all_of(c("var1", "var2", "num_obs"))) %>%
  pivot_wider(names_from = "var2", values_from = "num_obs") %>%
  column_to_rownames("var1") %>%
  as.matrix()

# Create matrix correlations and co-occurrences including cytokines and microbiome - SC
mat_sp <- matrix(0, nrow = ncol(mat_sc), ncol = ncol(mat_sc), dimnames = list(colnames(mat_sc), colnames(mat_sc)))
diag(mat_sp) <- 1
mat_cy <- matrix(0, nrow = nrow(mat_sc), ncol = nrow(mat_sc), dimnames = list(rownames(mat_sc), rownames(mat_sc)))
diag(mat_cy) <- 1

mat_tmp <- rbind(mat_sp, mat_sc)
co_tmp <- rbind(mat_sp, co_sc)

g2_mat <- cbind(mat_tmp, rbind(t(mat_sc), mat_cy))
g2_co <- cbind(co_tmp, rbind(t(co_sc), mat_cy))

# Do the same for NC
mat_nc <- df_nc %>%
  dplyr::select(all_of(c("var1", "var2", "cor"))) %>%
  pivot_wider(names_from = "var2", values_from = "cor") %>%
  column_to_rownames("var1") %>%
  as.matrix()
co_nc <- df_nc %>%
  dplyr::select(all_of(c("var1", "var2", "num_obs"))) %>%
  pivot_wider(names_from = "var2", values_from = "num_obs") %>%
  column_to_rownames("var1") %>%
  as.matrix()


mat_sp <- matrix(0, nrow = ncol(mat_nc), ncol = ncol(mat_nc), dimnames = list(colnames(mat_nc), colnames(mat_nc)))
diag(mat_sp) <- 1
mat_cy <- matrix(0, nrow = nrow(mat_nc), ncol = nrow(mat_nc), dimnames = list(rownames(mat_nc), rownames(mat_nc)))
diag(mat_cy) <- 1

mat_tmp <- rbind(mat_sp, mat_nc)
co_tmp <- rbind(mat_sp, co_nc)

g1_mat <- cbind(mat_tmp, rbind(t(mat_nc), mat_cy))
g1_co <- cbind(co_tmp, rbind(t(co_nc), mat_cy))

comm <- intersect(colnames(g2_mat), colnames(g1_mat))
g1_mat <- g1_mat[comm, comm]
g2_mat <- g2_mat[comm, comm]
g1_co <- g1_co[comm, comm]
g2_co <- g2_co[comm, comm]
stopifnot(all.equal(colnames(g2_mat), colnames(g1_mat)))
stopifnot(all.equal(colnames(g2_co), colnames(g1_co)))
stopifnot(all.equal(colnames(g2_co), colnames(g2_mat)))
stopifnot(all.equal(colnames(g1_co), colnames(g1_mat)))

g1_mat_org <- g1_mat
g2_mat_org <- g2_mat

# Keep only pairs that at least in one group the coefficient > 0.3.
g1_mat[abs(g1_mat_org) < 0.3 & abs(g2_mat_org) < 0.3] <- 0
g2_mat[abs(g1_mat_org) < 0.3 & abs(g2_mat_org) < 0.3] <- 0

g1_co <- g1_co[colnames(g1_mat), colnames(g1_mat)]
g2_co <- g2_co[colnames(g2_mat), colnames(g2_mat)]

# Set the colors
env <- ifelse(grepl("gut s__", rownames(g1_mat)), "Gut Species", ifelse(grepl("oral s__", rownames(g1_mat)), "Oral Species", "Cytokines"))
env <- factor(env, levels = c("Cytokines", "Gut Species", "Oral Species"))
names(env) <- rownames(g1_mat)

unique(diag(g1_mat))

# Construct a network
net <- netConstruct(
  data = g1_mat,
  data2 = g2_mat,
  dataType = "correlation",
  dissFunc = "unsigned", sparsMethod = "none"
)
# Analyze the network
props <- netAnalyze(net,
  centrLCC = TRUE,
  normDeg = TRUE,
  clustMethod = "cluster_fast_greedy", gcmHeat = FALSE,
  hubPar = c("degree")
)

pdf(file.path(output_dir_sub, paste0("network.pdf")), 60, 40)

plot(props,
  sameLayout = TRUE,
  # layout = "spring",
  repulsion = 0.3,
  shortenLabels = "simple",
  labelLength = 35,
  charToRm = "gut [a-z]__|oral [a-z]__",
  layoutGroup = "union",
  rmSingles = "inboth",
  nodeSizeSpread = 1,
  nodeColor = "feature",
  featVecCol = env,
  colorVec = mycols,
  nodeShape = c("triangle", "circle", "square"),
  featVecShape = env,
  nodeTransp = 30,
  sameClustCol = TRUE,
  labelScale = TRUE,
  # labelFont=0.8,
  nodeSize = "degree",
  cexNodes = 4,
  # cexLabels = 5,
  cexHubLabels = 3,
  cexTitle = 3.8,
  edgeFilter = "threshold",
  edgeFilterPar = 0.3,
  negDiffCol = TRUE,
  edgeTranspLow = 50,
  edgeTranspHigh = 50,
  groupNames = group_names,
  hubBorderCol = "gray40"
)

# Colors used in the legend should be equally transparent as in the plot
phylcol_transp <- colToTransp(mycols, 30)

legend(-1.2, 1.2,
  cex = 4, pt.cex = 6, title = "Feature:",
  legend = toupper(levels(env)), col = phylcol_transp, bty = "n", pch = c(17, 16,15)
)

legend("bottom",
  title = "estimated association:", legend = c("+", "-"),
  col = c("#009900", "red"), inset = 0.02, cex = 4, lty = 1, lwd = 4,
  bty = "n", horiz = TRUE
)


dev.off()

# Differential correlations
res <- fisher_test(mat1 = g1_mat, mat2 = g2_mat, mat1_cooccur = g1_co, mat2_cooccur = g2_co)
sig_diff <- res[res$diff_cor > 0.3 & res$p < 0.01, ]

write.csv(res, file.path(output_dir_sub, paste0("differential_correlation.csv")), row.names = F)

res_filtered <- sig_diff %>%
  dplyr::select(c(
    "Taxa", "name", "cooccurrence_g1", "cooccurrence_g2",
    "cor_g1", "cor_g2", "diff_cor", "sign", "p"
  )) %>%
  dplyr::rename(
    Cytokines = name, `cooccurrence_Non-HIV` = cooccurrence_g1, `cooccurrence_Pre-HIV` = cooccurrence_g2,
    `cor_Non-HIV` = cor_g1, `cor_Pre-HIV` = cor_g2
  ) %>%
  mutate(across(grep("cor", colnames(.)), ~ round(.x, 3))) %>%
  mutate(sign = ifelse(grepl("g1",sign),gsub("g1",group_names[1],sign),gsub("g2",group_names[2],sign))) 

write.csv(res_filtered, file.path(output_dir_sub, paste0("Cytokines_differential_correlation_sig.csv")), row.names = F)


if (nrow(sig_diff) > 0) {
  g1_mat_n <- matrix(0, nrow(g1_mat), ncol(g1_mat), dimnames = list(rownames(g1_mat), colnames(g1_mat)))
  g2_mat_n <- matrix(0, nrow(g2_mat), ncol(g2_mat), dimnames = list(rownames(g2_mat), colnames(g2_mat)))

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


  pdf(file.path(output_dir_sub, paste0("network_diff.pdf")), 70, 50)

  plot(props,
    sameLayout = TRUE,
    # repulsion = 0.1,
    shortenLabels = "simple",
    labelLength = 35,
    charToRm = "gut [a-z]__|oral [a-z]__",
    layoutGroup = "2",
    rmSingles = "inboth",
    nodeSizeSpread = 0.5,
    nodeColor = "feature",
    featVecCol = env,
    nodeShape = c("triangle", "circle", "square"),
    featVecShape = env,
    colorVec = mycols,
    labelScale = FALSE,
    nodeTransp = 30,
    # labelFont=0.8,
    nodeSize = "degree",
    # scalelabel = FALSE,
    cexNodes = 3,
    cexLabels = 2,
    cexHubLabels = 4,
    cexTitle = 3.8,
    negDiffCol = TRUE,
    edgeTranspLow = 50,
    edgeTranspHigh = 50,
    groupNames = group_names,
    hubBorderCol = "gray40"
  )

  # Colors used in the legend should be equally transparent as in the plot
  phylcol_transp <- colToTransp(mycols, 30)

  legend(-1.2, 1.2,
    cex = 6, pt.cex = 6, title = "Feature:",
    legend = toupper(levels(env)), col = phylcol_transp, bty = "n", pch = c(17, 16,15)
  )

  legend("bottom",
    title = "estimated association:", legend = c("+", "-"),
    col = c("#009900", "red"), inset = 0.02, cex = 4, lty = 1, lwd = 4,
    bty = "n", horiz = TRUE
  )



  dev.off()
}


##############################
# Scatterplts
##############################

fisher_res <- read.csv(file.path(output_dir, "Network", paste0("differential_correlation.csv")))
taxa_merged_sc <- read.csv(file.path(output_dir, "sc", paste0("sc_Cytokines_microbiome_merged.csv")), check.names = F)
taxa_merged_nc <- read.csv(file.path(output_dir, "nc", paste0("nc_Cytokines_microbiome_merged.csv")), check.names = F)


sig_diff <- fisher_res[fisher_res$diff_cor > 0.3 & fisher_res$p < 0.01, ]


plots <- lapply(1:nrow(sig_diff), function(x) {
  print(x)
  cyt <- sig_diff$name[x]
  taxa <- sig_diff$Taxa[x]

  p1 <- ggplot(taxa_merged_nc, aes(x = .data[[cyt]], y = .data[[taxa]])) +
    geom_point(aes(color = .data[["status"]])) +
    scale_color_manual(values = status_cols[1]) +
    geom_smooth(method = "lm", color = "darkgray", alpha = 0.1) +
    theme_classic() +
    labs(caption = paste0(
      "r = ", signif(sig_diff$cor_g1[x], 3)
    )) +
    theme(legend.position = "none") +
    labs(title = "Non-HIV")


  p2 <- ggplot(taxa_merged_sc, aes(x = .data[[cyt]], y = .data[[taxa]])) +
    geom_point(aes(color = .data[["status"]])) +
    scale_color_manual(values = status_cols[2]) +
    geom_smooth(method = "lm", color = "darkgray", alpha = 0.1) +
    theme_classic() +
    labs(caption = paste0(
      "r = ", signif(sig_diff$cor_g2[x], 3)
    )) +
    theme(legend.position = "none") +
    labs(title = "Pre-HIV")

  p <- ggarrange(p1, p2, nrow = 1)
  p <- annotate_figure(p, top = text_grob(paste0("Fisher'z test\np = ", signif(sig_diff$p[x], 3)),
    color = "darkgreen", face = "bold", size = 14
  ))
})

# Save the plots
pdf(file.path(output_dir, "Network", paste0("Cytokines_Fisher_significant.pdf")), 11, 5)
print(ggarrange(plotlist = plots, nrow = 1, ncol = 1))
dev.off()


######################
# Heatmap
######################
cor_sc <- read.table(file.path(output_dir, "sc", "sc_cytokines_microbiome.txt"), header = T, check.names = F, sep = "\t")
cor_nc <- read.table(file.path(output_dir, "nc", "nc_cytokines_microbiome.txt"), header = T, check.names = F, sep = "\t")
cor_sc$status <- "Pre-HIV"
cor_nc$status <- "Non-HIV"
cor_all <- rbind(cor_sc, cor_nc)
sig_diff$taxa_cyt <- paste0(sig_diff$Taxa, "_", sig_diff$name)

cor_all <- cor_all %>%
  filter(var2 %in% c(sig_diff$Taxa)) %>%
  mutate(taxa_cyt = paste0(var2, "_", var1)) %>%
  mutate(lab = ifelse(taxa_cyt %in% sig_diff$taxa_cyt, "*", ""))%>%
  mutate(var2 = gsub("s__","",var2)) %>%
  mutate(var2 = gsub("_"," ",var2))


p <- ggplot(cor_all, aes(y = var2, x = var1, fill = cor)) +
  geom_tile(color = "black") +
  geom_text(aes(label = lab), color = "black", size = 8) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    na.value = "white", midpoint = 0, limit = c(-max(abs(cor_all$cor)), max(abs(cor_all$cor))),
    name = "Correlation Coefficients"
  ) +
  scale_x_discrete(labels = c("sCD14", "sCD163", "CRP", "IL-6", "IP-10", "LBP", "CD4+/CD8+")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0.9)) +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(strip.text = element_text(size = 25)) +
  facet_wrap(~status, scales = "free")



pdf(file.path(output_dir, "heatmap_microbiome_differential_correlations.pdf"), 20, 5)
print(p)
dev.off()

