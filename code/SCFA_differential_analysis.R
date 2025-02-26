rm(list = ls())
source("code/load.R")
# Create output directory
outdir <- file.path(output_dir, "Differential_Abundance_SCFA")
dir.create(outdir, recursive = T, showWarnings = T)

scfa <- c(
  "acetic_acid", "propionic_acid", "butyric_acid", "valeric_acid",
  "isobutyric_acid", "isovaleric_acid", "hexanoic_acid"
)
all_scfa_names <- c(paste0("oral_", scfa), paste0("fecal_", scfa), paste0("plasma_", scfa))

# Load SCFA data
df <- read.csv("data/scfa/scfa_combined.csv")

# Get the number of samples
for (i in c("plasma","oral","fecal")){
  print(i)
  df_sub <- df[,c("status",all_scfa_names[grepl(i,all_scfa_names)])]
  df_sub$group <-  ifelse(rowSums(is.na(df_sub[,-1]))==7,"not_measure","measured")
  print(table(df_sub$status,df_sub$group))
}

# data transformation
scfa_log <- df %>%
  column_to_rownames("subjid") %>%
  dplyr::select(all_of(c(all_scfa_names))) %>%
  mutate_all(log)
meta_df <- df %>% dplyr::select(colnames(.)[!colnames(.) %in% colnames(scfa_log)])

# Merge with metadata and convert to long format
scfa_log <- scfa_log %>%
  rownames_to_column("subjid") %>%
  left_join(meta_df)

# Save the results
dir.create("data/processed_data/SCFA", showWarnings = F, recursive = T)
write.csv(scfa_log, file.path("data/processed_data/SCFA", "SCFA_log.csv"), row.names = F, quote = F)

# Convert to long format
scfa_long <- scfa_log %>% pivot_longer(cols = all_of(all_scfa_names), names_to = "scfa")

# Loc should be a factor
scfa_long$loc <- as.factor(scfa_long$loc)

##### NC versus SC at baseline
df_p <- data.frame(SCFA = character(), Estimate = numeric(), p = numeric())
pdf(file.path(outdir, "boxplots_status.pdf"), 5, 5)
for (s in all_scfa_names) {
  scf_sub <- scfa_long %>% filter(scfa == s)
  s1 <- gsub("_", " ", stringr::str_to_title(s))
  fit <- lm_test(scf_sub, x = "status", y = "value", plot_title = s1, y_label = "Log transformed")
  df_tmp <- data.frame(SCFA = s, Estimate = fit[[1]]$coefficients[2, 1], p = fit[[1]]$coefficients[2, "Pr(>|t|)"])
  df_p <- rbind(df_p, df_tmp)
  plot <- fit[[2]] + scale_x_discrete(labels = c("nc" = group_names[1], "sc" = group_names[2]))
  print(plot)
}
dev.off()

df_p_sig <- df_p %>%
  filter(p < 0.05) %>%
  arrange()

write.csv(df_p, file.path(outdir, "SCFA_Diff.csv"), row.names = F, quote = F)
write.csv(df_p_sig, file.path(outdir, "SCFA_Diff_sig.csv"), row.names = F, quote = F)

# Plot only significant SCFAs
scf_sub <- scfa_long %>%
  filter(scfa %in% df_p_sig$SCFA) %>%
  mutate(scfa = factor(scfa, levels = df_p_sig$SCFA))

scf_sub$env <- stringr::str_to_title(str_extract(scf_sub$scfa, "oral|fecal|plasma"))

plots <- lapply(unique(scf_sub$env), function(x) {
  scf_sub1 <- scf_sub %>% filter(env == x)
  df_p_sig1 <- df_p_sig %>% filter(grepl(tolower(x), SCFA))

  # For plotting, use rstatix::wilcox_test to get the table with x and y positions but then
  # replace the p-values from the linear model.
  stat <- scf_sub1 %>%
    group_by(env, scfa) %>%
    rstatix::wilcox_test(value ~ status) %>%
    add_xy_position() %>%
    arrange(order(match(df_p_sig1$SCFA, scfa))) %>%
    mutate(num = 1:n()) %>%
    mutate(xmin = num - 0.1, xmax = num + 0.1) %>%
    ungroup()

  stopifnot(all.equal(as.character(stat$scfa), df_p_sig1$SCFA))
  stat$p <- round(df_p_sig1$p, 3)


  stat$scfa <- gsub("oral_|fecal_|plasma_", "", stat$scfa)
  scf_sub1$scfa <- gsub("oral_|fecal_|plasma_", "", scf_sub1$scfa)
  scf_sub1$scfa <- gsub("_", " ", scf_sub1$scfa)
  scf_sub1$scfa <- stringr::str_to_sentence(scf_sub1$scfa)
  scf_sub1$scfa <- factor(scf_sub1$scfa, levels = stringr::str_to_sentence(scfa))


  p <- ggplot(scf_sub1, aes(x = scfa, y = value, color = status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5, size = 1) +
    scale_color_manual(values = status_cols, labels = group_names) +
    theme_bw(base_size = 16) +
    stat_pvalue_manual(stat, label = "p = {p}") +
    theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
    labs(y = "log transformed", x = "", title = paste0(ifelse(x == "Fecal", "Gut", x), " SCFA"))
})

pdf(file.path(outdir, "boxplots_status_oral.pdf"), 7, 5.5)
print(plots[[1]])
dev.off()
pdf(file.path(outdir, "boxplots_status_gut.pdf"), 7, 5.5)
print(plots[[2]])
dev.off()
pdf(file.path(outdir, "boxplots_status_plasma.pdf"), 7, 5.5)
print(plots[[3]])
dev.off()
