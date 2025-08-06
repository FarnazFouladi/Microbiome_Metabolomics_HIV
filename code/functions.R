# Linear model accounting for age, location, antibiotic use
lm_test <- function(df = scf_tmp, x = "status", y = "value", plot_title = "Plasma Acetic Acid", y_label = "Log transformed", isPaired = FALSE) {
  # Remove NAs
  df <- df %>% filter(!is.na(!!sym(x)))
  
  fit <- summary(lm(as.formula(paste0(y, "~", x, "+ age + loc + abx_use")), data = df))
  
  fig_annotate <- df %>%
    rstatix::wilcox_test(as.formula(paste0(y, "~", x))) %>%
    add_xy_position() %>%
    mutate(p = fit$coefficients["statussc", "Pr(>|t|)"]) %>%
    mutate(p.format = paste0("p = ", round(p, 3)))
  
  plot <- ggplot(data = df, aes(x = .data[[x]], .data[[y]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha = 0.5, size = 2, position = position_jitter(width = 0.1), aes(color = .data[[x]])) +
    theme_bw(base_size = 16) +
    labs(title = plot_title, y = y_label, x = "") +
    lims(y = c(min(df[, y], na.rm = T), fig_annotate$y.position[1] + 0.1)) +
    theme(plot.title = element_text(face = "bold", size = 16), legend.position = "none") +
    scale_color_manual(values = c("#0072b5ff", "#bc3c29ff")) +
    stat_pvalue_manual(fig_annotate, label = "p.format")
  
  return(list(fit, plot))
}


create_boxplots <- function(data = df_merged, ancomcs_results, group = "status", 
                            cols = status_cols, subt = "v1", titley ="Normalized peak intesity",
                            group_names,
                            pval = "p_statussc", 
                            qval = "q_statussc", legend = TRUE) {
  boxplots_list <- lapply(as.character(ancomcs_results$Feature), function(a) {
    print(a)
    ancomcs_results_sub <- ancomcs_results %>% filter(Feature == a)
    stat_summary <- paste(
      p_format(ancomcs_results_sub[, pval], 3,space = TRUE,add.p = T,type = "p"), "\n",
      p_format(ancomcs_results_sub[, qval], 3, space = TRUE, add.p = T, type = "p.adj")
    )
    stat_summary<- gsub("\n p","\n q",stat_summary)
    
    df_tmp <- data.frame(feature =  data[,as.character(a)],group = data[,group])
    plot <- df_tmp %>% filter(!is.na(feature)) %>% ggplot(aes(y = feature, x = group )) +
      geom_boxplot(outlier.colour = NA) +
      geom_jitter(position = position_jitter(width = 0.1, height = 0), size = 2, alpha = 0.5, aes(color = group)) +
      scale_color_manual(values = cols) +
      theme_bw(base_size = 16) +
      scale_x_discrete(labels=group_names)+
      labs(title = a, caption = stat_summary, y = titley, subtitle = subt, x = "")+{if(!legend) theme(legend.position = "none")}
  })
  
  return(boxplots_list)
}
