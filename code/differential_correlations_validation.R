source("code/load.R")
source("code/differential_correlations_functions.R")
library(doParallel)
library(doRNG)

set.seed(1654)

# Load data
sc_bias_corrected <- read_rds(file.path("Outputs/Processed_data", "processed_pre_hiv.rds"))
nc_bias_corrected <- read_rds(file.path("Outputs/Processed_data", "processed_non_hiv.rds"))

output_dir <- file.path("Outputs", "Differential_Correlations_Validations", "Tables")
dir.create(output_dir, showWarnings = F, recursive = T)

# Gut species
df.nc <- t(nc_bias_corrected[[1]])
df.sc <- t(sc_bias_corrected[[1]])
common_features <- intersect(colnames(df.nc), colnames(df.sc))
common_features <- common_features[!grepl("s__sp.", common_features)]
df.nc <- df.nc[, common_features]
df.sc <- df.sc[, common_features]


####################################################
# Re-sampling for real data
####################################################

B <- 1000
cl <- makeCluster(5)
registerDoParallel(cl = cl)
res_sc <- foreach(i = 1:B, .verbose = T, .packages = c("tidyverse")) %dorng% {
  df.nc1 <- df.nc[sample(1:nrow(df.nc), nrow(df.nc), replace = T), ]
  df.sc2 <- df.sc[sample(1:nrow(df.sc), nrow(df.sc), replace = T), ]

  species_res <- dysco_score(
    data1 = list(df.nc1), data2 = list(df.sc2),
    output_dir = output_dir, prefix = "gut_species"
  )[[3]]
}
stopCluster(cl)


# Median of test statistic
df_t_real <- lapply(1:B, function(i) {
  df_tmp <- res_sc[[i]]
  df_tmp <- df_tmp[, c("species", "t"), drop = F]
  colnames(df_tmp)[2] <- paste0("t_", i)
  return(df_tmp)
})
df_t_merged_real <- Reduce(full_join, df_t_real)
df_t_merged_real <- df_t_merged_real %>% column_to_rownames("species")
mean_t_sc <- data.frame(taxa = rownames(df_t_merged_real), median_t = rowMedians(as.matrix(df_t_merged_real), na.rm = T))
mean_t_sc <- mean_t_sc %>% arrange(desc(median_t))
write.csv(mean_t_sc, file.path("Outputs/Differential_Correlations_Validations", "median_t_real_data_resampling.csv"), quote = F, row.names = F)


####################################################
# Re-sampling for pooled data
####################################################

B <- 1000
cl <- makeCluster(5)
registerDoParallel(cl = cl)

df.pooled <- rbind(df.nc, df.sc)
res_pool <- foreach(i = 1:B, .verbose = T, .packages = c("tidyverse")) %dorng% {
  df.nc1 <- df.pooled[sample(1:nrow(df.pooled), nrow(df.nc), replace = T), ]
  df.sc2 <- df.pooled[sample(1:nrow(df.pooled), nrow(df.sc), replace = T), ]

  species_res <- dysco_score(
    data1 = list(df.nc1), data2 = list(df.sc2),
    output_dir = output_dir, prefix = "gut_species"
  )[[3]]
}
stopCluster(cl)


# Median of test statistic
df_t_pool <- lapply(1:B, function(i) {
  df_tmp <- res_pool[[i]]
  df_tmp <- df_tmp[, c("species", "t"), drop = F]
  colnames(df_tmp)[2] <- paste0("t_", i)
  return(df_tmp)
})
df_t_merged_pool <- Reduce(full_join, df_t_pool)
df_t_merged_pool <- df_t_merged_pool %>% column_to_rownames("species")
mean_t_pool <- data.frame(taxa = rownames(df_t_merged_pool), median_t = rowMedians(as.matrix(df_t_merged_pool), na.rm = T))
mean_t_pool <- mean_t_pool %>% arrange(desc(median_t))
write.csv(mean_t_pool, file.path("Outputs/Differential_Correlations_Validations", "median_t_pooled_data_resampling.csv"), quote = F, row.names = F)

######## Plotting:
df_m <- merge(mean_t_sc, mean_t_pool, by = "taxa", suffixes = c("_real", "_pooled"))
df_long <- df_m %>%
  pivot_longer(cols = -c("taxa")) %>%
  mutate(name = ifelse(name == "median_t_real", "Real data", "Null data"))

p <- ggplot(df_long, aes(x = name, y = value, color = name)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.1), alpha = 0.5) +
  ggrepel::geom_text_repel(aes(label = taxa),
    size = 2,
    box.padding = 0.15,
    point.padding = 0.5,
    segment.color = "grey50", color = "black"
  ) +
  theme_bw(base_size = 12) +
  labs(y = "Median DISCO Score by Species", x = "", title = "Comparison of Median Scores between Null and Real Data", color = "Data") +
  theme(legend.position = "none") +
  theme(title = element_text(face = "bold"))

pdf(file.path("Outputs/Differential_Correlations_Validations", "null_vs_real_boxplot.pdf"), 6.5, 5)
print(p)
dev.off()
