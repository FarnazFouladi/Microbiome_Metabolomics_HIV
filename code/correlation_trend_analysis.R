# Trend analysis
output_dir <- "Outputs/Correlations_trend_analysis"
dir.create(output_dir, showWarnings = F, recursive = T)
source("code/pava.R")

g1_bias_corrected <- read_rds(file.path("Outputs/Processed_data", "processed_g1.rds"))
g2_bias_corrected <- read_rds(file.path("Outputs/Processed_data", "processed_g2.rds"))
g3_bias_corrected <- read_rds(file.path("Outputs/Processed_data", "processed_g3.rds"))
g4_bias_corrected <- read_rds(file.path("Outputs/Processed_data", "processed_g4.rds"))

# Gut species
gut.species.list <- list(g1_bias_corrected[[1]], g2_bias_corrected[[1]], g3_bias_corrected[[1]], g4_bias_corrected[[1]])
gut_species_trend <- trend_test(data_modality1 = gut.species.list, data_modality2 = NULL, cooccur = 10, output_dir = output_dir, prefix = "gut_species")

# Gut GO
gut.go.list <- list(g1_bias_corrected[[2]], g2_bias_corrected[[2]], g3_bias_corrected[[2]], g4_bias_corrected[[2]])
gut_go_trend <- trend_test(data_modality1 = gut.species.list, data_modality2 = gut.go.list, cooccur = 10, output_dir, "gut_go")

# Gut metabolites
gut.meta.p.list <- list(g1_bias_corrected[[3]], g2_bias_corrected[[3]], g3_bias_corrected[[3]], g4_bias_corrected[[3]])
gut_meta.p_trend <- trend_test(data_modality1 = gut.species.list, data_modality2 = gut.meta.p.list, cooccur = 10, output_dir, "gut_meta.p")

gut.meta.n.list <- list(g1_bias_corrected[[4]], g2_bias_corrected[[4]], g3_bias_corrected[[4]], g4_bias_corrected[[4]])
gut_meta.n_trend <- trend_test(data_modality1 = gut.species.list, data_modality2 = gut.meta.n.list, cooccur = 10, output_dir, "gut_meta.n")

# Plasma metabolites
plasma.meta.p.list <- list(g1_bias_corrected[[5]], g2_bias_corrected[[5]], g3_bias_corrected[[5]], g4_bias_corrected[[5]])
plasma_meta.p_trend <- trend_test(data_modality1 = gut.species.list, data_modality2 = plasma.meta.p.list, cooccur = 10, output_dir, "gut_plasma.meta.p")

plasma.meta.n.list <- list(g1_bias_corrected[[6]], g2_bias_corrected[[6]], g3_bias_corrected[[6]], g4_bias_corrected[[6]])
plasma_meta.n_trend <- trend_test(data_modality1 = gut.species.list, data_modality2 = plasma.meta.n.list, cooccur = 10, output_dir, "gut_plasma.meta.n")

# oral species
oral.species.list <- list(g1_bias_corrected[[7]], g2_bias_corrected[[7]], g3_bias_corrected[[7]], g4_bias_corrected[[7]])
oral_species_trend <- trend_test(data_modality1 = oral.species.list, data_modality2 = NULL, cooccur = 10, output_dir, "oral_species")

# oral GO
oral.go.list <- list(g1_bias_corrected[[8]], g2_bias_corrected[[8]], g3_bias_corrected[[8]], g4_bias_corrected[[8]])
oral_go_trend <- trend_test(data_modality1 = oral.species.list, data_modality2 = oral.go.list, cooccur = 10, output_dir, "oral_go")

# oral metabolites
oral.meta.p.list <- list(g1_bias_corrected[[9]], g2_bias_corrected[[9]], g3_bias_corrected[[9]], g4_bias_corrected[[9]])
oral_meta.p_trend <- trend_test(data_modality1 = oral.species.list, data_modality2 = oral.meta.p.list, cooccur = 10, output_dir, "oral_meta.p")

oral.meta.n.list <- list(g1_bias_corrected[[10]], g2_bias_corrected[[10]], g3_bias_corrected[[10]], g4_bias_corrected[[10]])
oral_meta.n_trend <- trend_test(data_modality1 = oral.species.list, data_modality2 = oral.meta.n.list, cooccur = 10, output_dir, "oral_meta.n")

# Plasma metabolites
oral.plasma_meta.p_trend <- trend_test(data_modality1 = oral.species.list, data_modality2 = plasma.meta.p.list, cooccur = 10, output_dir, "oral_plasma.meta.p")
oral.plasma_meta.n_trend <- trend_test(data_modality1 = oral.species.list, data_modality2 = plasma.meta.n.list, cooccur = 10, output_dir, "oral_plasma.meta.n")


# oral gut species
gut.species.list <- lapply(gut.species.list, function(x) {
  rownames(x) <- paste0("gut_", rownames(x))
  return(x)
})
oral_gut_species_trend <- trend_test(data_modality1 = oral.species.list, data_modality2 = gut.species.list, cooccur = 10, output_dir, "oral_gut_species")



########## Add the trend analysis to the significant differential correlation
gut_species_diff <- read.csv(file.path("Outputs/Differential_Correlations/Tables", "gut_species_cor_diff_sig.csv"), check.names = F)
gut_species_diff <- gut_species_diff %>% mutate(
  t2_t1 = paste0(t2, "_", t1),
  t1_t2 = paste0(t1, "_", t2)
)

gut_species_trend <- read.csv(file.path("Outputs/Correlations_trend_analysis", paste0("gut_species", "_pava_results.csv")), check.names = F)

gut_species_trend <- gut_species_trend %>% mutate(
  t2_t1 = paste0(t2, "_", t1),
  t1_t2 = paste0(t1, "_", t2)
)
pair_correlations <- unique(gut_species_trend$t2_t1, gut_species_trend$t1_t2)

gut_species_res <- list()
for (i in 1:nrow(gut_species_diff)) {
  p1 <- gut_species_diff$t1_t2[i]
  p2 <- gut_species_diff$t2_t1[i]

  if (p1 %in% pair_correlations | p2 %in% pair_correlations) {
    df_sub <- gut_species_trend %>% filter(`t1_t2` == p1 | `t2_t1` == p1 | `t1_t2` == p2 | `t2_t1` == p2)
    gut_species_res[[i]] <- df_sub %>% select(colnames(.)[3:9])
  } else {
    gut_species_res[[i]] <- data.frame(r1 = NA, r2 = NA, r3 = NA, r4 = NA, E = NA, pvals = NA, trend = NA)
  }
}

gut_species_res <- bind_rows(gut_species_res) %>%
  dplyr::rename(r1.g1 = r1, r2.g2 = r2, r3.g3 = r3, r4.g4 = r4, p.trend = pvals)
gut_species_diff <- cbind(gut_species_diff, gut_species_res) %>%
  as.data.frame() %>%
  select(-c("t2_t1", "t1_t2")) %>%
  arrange(`p.trend`) %>%
  mutate(`adj.p.trend` = p.adjust(`p.trend`, method = "BH"))

write.csv(gut_species_diff, file.path(output_dir, "gut_species_cor_diff_sig_with_trend_test.csv"), row.names = F, quote = F)



# Other datasets
for (d in c("gut_go", "gut_meta.p", "gut_meta.n", "gut_plasma.meta.p", "gut_plasma.meta.n", "oral_go", "oral_meta.n", "oral_meta.p", "oral_plasma.meta.n", "oral_plasma.meta.p", "oral_gut_species")) {
  diff <- read.csv(file.path("Outputs/Differential_Correlations/Tables", paste0(d, "_cor_diff_sig.csv")), check.names = F)
  trend <- read.csv(file.path("Outputs/Correlations_trend_analysis", paste0(d, "_pava_results.csv")), check.names = F)
  diff_with_trend <- diff %>%
    left_join(
      trend %>%
        dplyr::rename(r1.g1 = r1, r2.g2 = r2, r3.g3 = r3, r4.g4 = r4, p.trend = pvals),
      by = c("t1", "t2")
    ) %>%
    arrange(`p.trend`) %>%
    mutate(`adj.p.trend` = p.adjust(`p.trend`, method = "BH"))

  write.csv(diff_with_trend, file.path("Outputs/Correlations_trend_analysis", paste0(d, "_cor_diff_sig_with_trend_test.csv")), row.names = F, quote = T)
}
