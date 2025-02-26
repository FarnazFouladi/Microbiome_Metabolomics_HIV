library(tidyverse)
permutation = 1000
args <- commandArgs(TRUE)
# Rscript dysbiosis_score_permutation_pvals.R [Sample_Type] [Directory_of_Permutaion_Results] [Directory_of_Original_Results]

score_dir <- file.path("Outputs",args[2],args[1],"Scores")
score_dir1 <- file.path("Outputs",args[3],args[1],"Scores")
  
# Scores for 1000 permutations
perm_list <- lapply(1:permutation, function(i){
  df_tmp <- read.csv(file.path(score_dir, paste0("scores_wide_diff_p_", i, ".csv")),header = T,check.names = F)
  colnames(df_tmp)[1] <- "Taxa"
  df_tmp <- df_tmp %>% column_to_rownames("Taxa")
})

# Real scores
res <- read.csv(file.path(score_dir1, paste0("scores_wide_diff_p.csv")),header = T,check.names = F)
colnames(res)[1] <- "Taxa"
res <- res  %>% column_to_rownames("Taxa")

# p-values
perm_compare <- lapply(perm_list, function(s) {
  
  stopifnot(all(rownames(res) %in% rownames(s)))
  s1 <-s[rownames(res),colnames(res)]
  return(abs(s1) >= abs(res) )
  
})

pvals <- Reduce(`+`, perm_compare) / permutation

write.csv(pvals, file.path(score_dir1,"pvals_p_diff.csv"),row.names = T,quote = F)
