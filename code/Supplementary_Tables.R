rm(list = ls())
output_dir <- file.path("Supplementary_Tables")
dir.create(output_dir)
library(tidyverse)
library(openxlsx)

make_supplementary_tables <- function(name,sheet_name,data,output_dir){
  wb <- createWorkbook(name)
  
  for (i in 1:length(data)){
    print(i)
    addWorksheet(wb, sheet_name[i])
    writeData(wb,sheet_name[i] ,data[[i]])
  }
  saveWorkbook(wb, file.path(output_dir,paste0(name,".xlsx")), overwrite = TRUE)
  
}

# Supplemental Table 1: Differential abundance - Seroconversion
Gut_Species_Diff <- "Outputs/Differential_Abundance_Microbiome/Gut/Species/Species_ANCOMBC_results_sig.txt"
Gut_GO_Diff <- "Outputs/Differential_Abundance_Microbiome/Gut/GO/GO_ANCOMBC_results_sig.txt"
Gut_MTB_Diff <- "Outputs/Differential_Abundance_Metabolomics/Gut/metabolites_ANCOMBC_results_sig.txt"
Plasma_MTB_Diff <- "Outputs/Differential_Abundance_Metabolomics/Plasma/metabolites_ANCOMBC_results_sig.txt"
SCFA_Diff <- "Outputs/Differential_Abundance_SCFA/SCFA_Diff_sig.csv"


df1 <- read.table(Gut_Species_Diff, sep = "\t",header = T,check.names = F,quote = "",comment.char = "") %>%
  arrange(q) 
df2 <- read.table(Gut_GO_Diff, sep = "\t",header = T,check.names = F,quote = "",comment.char = "") %>%
  arrange(q) 
df3 <- read.table(Gut_MTB_Diff, sep = "\t",header = T,check.names = F,quote = "",comment.char = "") %>%
  arrange(q) 
df4 <- read.table(Plasma_MTB_Diff, sep = "\t",header = T,check.names = F,quote = "",comment.char = "") %>%
  arrange(q) 
df5 <- read.table(SCFA_Diff, sep = ",",header = T,check.names = F,comment.char = "") %>%
  mutate(Estimate = round(Estimate ,digits = 3))%>%
  filter(!grepl("Oral",SCFA))


make_supplementary_tables(name = "Supplementary_Table1",
                          sheet_name = c("A.Gut Species","B.Gut GO","C.Gut Metabolites","D.Plasma Metabolites","E.SCFA"),
                          data = list(df1,df2,df3,df4,df5),output_dir)



# Supplemental Table 2 Differential abundance - Trend analysis for sexual activity
Gut_Species_Diff <- "Outputs/Differential_Abundance_Microbiome_SA/Gut/Species/Species_ANCOMBC_results_sig.txt"
Gut_GO_Diff <- "Outputs/Differential_Abundance_Microbiome_SA/Gut/GO/GO_ANCOMBC_results_sig.txt"
Gut_MTB_Diff <- "Outputs/Differential_Abundance_Metabolomics_SA/Gut/metabolites_ANCOMBC_results_sig.txt"
Plasma_MTB_Diff <- "Outputs/Differential_Abundance_Metabolomics_SA/Plasma/metabolites_ANCOMBC_results_sig.txt"

df1 <- read.table(Gut_Species_Diff, sep = "\t",header = T,check.names = F,quote = "",comment.char = "") %>%
  dplyr::rename(`lfc_G2:G1` = lfc_group1g2, `lfc_G3:G1` = lfc_group1g3, `lfc_G4:G1` = lfc_group1g4 ) %>% arrange(q) 
df2 <- read.table(Gut_GO_Diff, sep = "\t",header = T,check.names = F,quote = "",comment.char = "") %>%
  dplyr::rename(`lfc_G2:G1` = lfc_group1g2, `lfc_G3:G1` = lfc_group1g3, `lfc_G4:G1` = lfc_group1g4 ) %>% arrange(q) 
df3 <- read.table(Gut_MTB_Diff, sep = "\t",header = T,check.names = F,quote = "",comment.char = "") %>%
  dplyr::rename(`lfc_G2:G1` = lfc_group1g2, `lfc_G3:G1` = lfc_group1g3, `lfc_G4:G1` = lfc_group1g4 ) %>% arrange(q) 
df4 <- read.table(Plasma_MTB_Diff, sep = "\t",header = T,check.names = F,quote = "",comment.char = "") %>%
  dplyr::rename(`lfc_G2:G1` = lfc_group1g2, `lfc_G3:G1` = lfc_group1g3, `lfc_G4:G1` = lfc_group1g4 ) %>% arrange(q) 


make_supplementary_tables(name = "Supplementary_Table2",
                          sheet_name = c("A.Gut Species","B.Gut GO","C.Gut Metabolites","D.Plasma Metabolites"),
                          data = list(df1,df2,df3,df4),output_dir)



# Supplemental Table 3: Mixomics

####### Supplemental Table 4: Differential Correlations for gut species
Gut_Species <- "Outputs/Differential_Correlations/Tables/gut_species_cor_diff.csv"
Gut_GO <- "Outputs/Differential_Correlations/Tables/gut_go_cor_diff.csv"
Gut_MTB_NEG <- "Outputs/Differential_Correlations/Tables/gut_meta.n_cor_diff.csv"
Gut_MTB_POS <- "Outputs/Differential_Correlations/Tables/gut_meta.p_cor_diff.csv"
Plasma_MTB_NEG <- "Outputs/Differential_Correlations/Tables/gut_plasma.meta.n_cor_diff.csv"
Plasma_MTB_POS <- "Outputs/Differential_Correlations/Tables/gut_plasma.meta.p_cor_diff.csv"


df1 <- read.csv(Gut_Species, check.names = F) %>% mutate(across(c(r1, r2), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2) %>% arrange(BH_p)
df2 <- read.csv(Gut_GO, check.names = F) %>% mutate(across(c(r1, r2), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2) %>% arrange(BH_p)
df3 <- read.csv(Gut_MTB_NEG,check.names = F) %>% mutate(across(c(r1, r2), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2) %>% arrange(BH_p)
df4 <- read.csv(Gut_MTB_POS, check.names = F) %>% mutate(across(c(r1, r2), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2) %>% arrange(BH_p)
df5 <- read.csv(Plasma_MTB_NEG, check.names = F) %>% mutate(across(c(r1, r2), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2) %>% arrange(BH_p)
df6 <- read.csv(Plasma_MTB_POS, check.names = F) %>% mutate(across(c(r1, r2), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2) %>% arrange(BH_p)


make_supplementary_tables(name = "Supplementary_Table4",
                          sheet_name = c("A.Gut Species","B.Gut GO","C.Gut Metabolome-NEG","D.Gut Metabolome-POS","E.Plasma Metabolome-NEG","F.Plasma Metabolome-POS"),
                          data = list(df1,df2,df3,df4,df5,df6),output_dir)


# Supplemental Table 5: Significant differential Correlations for gut species with trend test
Gut_Species <- "Outputs/Correlations_trend_analysis/gut_species_cor_diff_sig_with_trend_test.csv"
Gut_GO <- "Outputs/Correlations_trend_analysis/gut_go_cor_diff_sig_with_trend_test.csv"
Gut_MTB_NEG <- "Outputs/Correlations_trend_analysis/gut_meta.n_cor_diff_sig_with_trend_test.csv"
Gut_MTB_POS <- "Outputs/Correlations_trend_analysis/gut_meta.p_cor_diff_sig_with_trend_test.csv"
Plasma_MTB_NEG <- "Outputs/Correlations_trend_analysis/gut_plasma.meta.n_cor_diff_sig_with_trend_test.csv"
Plasma_MTB_POS <- "Outputs/Correlations_trend_analysis/gut_plasma.meta.p_cor_diff_sig_with_trend_test.csv"


df1 <- read.csv(Gut_Species, check.names = F) %>% select(-c("num_sig_corr_t1","num_sig_corr_t2")) %>% mutate(across(c(r1, r2,`r1.g1`,`r2.g2`,`r3.g3`,`r4.g4`), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2,`ρG1`=`r1.g1`,`ρG2`=`r2.g2`,`ρG3`=`r3.g3`,`ρG4`=`r4.g4`) %>% arrange(`adj.p.trend`)
df2 <- read.csv(Gut_GO, check.names = F)%>% select(-c("num_sig_corr_t1","num_sig_corr_t2")) %>% mutate(across(c(r1, r2,`r1.g1`,`r2.g2`,`r3.g3`,`r4.g4`), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2,`ρG1`=`r1.g1`,`ρG2`=`r2.g2`,`ρG3`=`r3.g3`,`ρG4`=`r4.g4`) %>% arrange(`adj.p.trend`)
df3 <- read.csv(Gut_MTB_NEG,check.names = F)%>% select(-c("num_sig_corr_t1","num_sig_corr_t2")) %>% mutate(across(c(r1, r2,`r1.g1`,`r2.g2`,`r3.g3`,`r4.g4`), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2,`ρG1`=`r1.g1`,`ρG2`=`r2.g2`,`ρG3`=`r3.g3`,`ρG4`=`r4.g4`) %>% arrange(`adj.p.trend`)
df4 <- read.csv(Gut_MTB_POS, check.names = F)%>% select(-c("num_sig_corr_t1","num_sig_corr_t2")) %>% mutate(across(c(r1, r2,`r1.g1`,`r2.g2`,`r3.g3`,`r4.g4`), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2,`ρG1`=`r1.g1`,`ρG2`=`r2.g2`,`ρG3`=`r3.g3`,`ρG4`=`r4.g4`) %>% arrange(`adj.p.trend`)
df5 <- read.csv(Plasma_MTB_NEG, check.names = F)%>% select(-c("num_sig_corr_t1","num_sig_corr_t2")) %>% mutate(across(c(r1, r2,`r1.g1`,`r2.g2`,`r3.g3`,`r4.g4`), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2,`ρG1`=`r1.g1`,`ρG2`=`r2.g2`,`ρG3`=`r3.g3`,`ρG4`=`r4.g4`) %>% arrange(`adj.p.trend`)
df6 <- read.csv(Plasma_MTB_POS, check.names = F)%>% select(-c("num_sig_corr_t1","num_sig_corr_t2")) %>% mutate(across(c(r1, r2,`r1.g1`,`r2.g2`,`r3.g3`,`r4.g4`), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2,`ρG1`=`r1.g1`,`ρG2`=`r2.g2`,`ρG3`=`r3.g3`,`ρG4`=`r4.g4`) %>% arrange(`adj.p.trend`)

make_supplementary_tables(name = "Supplementary_Table5",
                          sheet_name = c("A.Gut Species","B.Gut GO","C.Gut Metabolome-NEG","D.Gut Metabolome-POS","E.Plasma Metabolome-NEG","F.Plasma Metabolome-POS"),
                          data = list(df1,df2,df3,df4,df5,df6),output_dir)


# Supplemental Table 6: Differential abundance - Oral Species
Oral_Species_Diff <- "Outputs/Differential_Abundance_Microbiome/Oral/Species/Species_ANCOMBC_results_sig.txt"
Oral_GO_Diff <- "Outputs/Differential_Abundance_Microbiome/Oral/GO/GO_ANCOMBC_results_sig.txt"
SCFA_Diff <- "Outputs/Differential_Abundance_SCFA/SCFA_Diff_sig.csv"

df1 <- read.table(Oral_Species_Diff, sep = "\t",header = T,check.names = F,quote = "",comment.char = "") %>%
  arrange(q) 
df2 <- read.table(Oral_GO_Diff, sep = "\t",header = T,check.names = F,quote = "",comment.char = "") %>%
  arrange(q) 
df3 <- read.table(SCFA_Diff, sep = ",",header = T,check.names = F,comment.char = "") %>%
  mutate(Estimate = round(Estimate ,digits = 3)) %>%
  filter(grepl("oral",SCFA))
  
make_supplementary_tables(name = "Supplementary_Table6",
                          sheet_name = c("A.Oral Species","B.Oral GO","C.Oral SCFA"),
                          data = list(df1,df2,df3),output_dir)


# Supplemental Table 7 Differential abundance - Trend analysis for sexual activity - Oral
Oral_Species_Diff <- "Outputs/Differential_Abundance_Microbiome_SA/Oral/Species/Species_ANCOMBC_results_sig.txt"
Oral_GO_Diff <- "Outputs/Differential_Abundance_Microbiome_SA/Oral/GO/GO_ANCOMBC_results_sig.txt"
Oral_MTB_Diff <- "Outputs/Differential_Abundance_Metabolomics_SA/Oral/metabolites_ANCOMBC_results_sig.txt"

df1 <- read.table(Oral_Species_Diff, sep = "\t",header = T,check.names = F,quote = "",comment.char = "") %>%
  dplyr::rename(`lfc_G2:G1` = lfc_group1g2, `lfc_G3:G1` = lfc_group1g3, `lfc_G4:G1` = lfc_group1g4 ) %>% arrange(q) 
df2 <- read.table(Oral_GO_Diff, sep = "\t",header = T,check.names = F,quote = "",comment.char = "") %>%
  dplyr::rename(`lfc_G2:G1` = lfc_group1g2, `lfc_G3:G1` = lfc_group1g3, `lfc_G4:G1` = lfc_group1g4 ) %>% arrange(q) 
df3 <- read.table(Oral_MTB_Diff, sep = "\t",header = T,check.names = F,quote = "",comment.char = "") %>%
  dplyr::rename(`lfc_G2:G1` = lfc_group1g2, `lfc_G3:G1` = lfc_group1g3, `lfc_G4:G1` = lfc_group1g4 ) %>% arrange(q) 

make_supplementary_tables(name = "Supplementary_Table7",
                          sheet_name = c("A.Oral Species","B.Oral GO","C.Oral Metabolites"),
                          data = list(df1,df2,df3),output_dir)


# Supplemental Table 8, differential correlation for oral species
Oral_Species <- "Outputs/Differential_Correlations/Tables/oral_species_cor_diff.csv"
Oral_GO <- "Outputs/Differential_Correlations/Tables/oral_go_cor_diff.csv"
Oral_MTB_NEG <- "Outputs/Differential_Correlations/Tables/oral_meta.n_cor_diff.csv"
Oral_MTB_POS <- "Outputs/Differential_Correlations/Tables/oral_meta.p_cor_diff.csv"
Plasma_MTB_NEG <- "Outputs/Differential_Correlations/Tables/oral_plasma.meta.n_cor_diff.csv"
Plasma_MTB_POS <- "Outputs/Differential_Correlations/Tables/oral_plasma.meta.p_cor_diff.csv"
Gut_Species <- "Outputs/Differential_Correlations/Tables/oral_gut_species_cor_diff.csv"

df1 <- read.csv(Oral_Species, check.names = F) %>% mutate(across(c(r1, r2), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2) %>% arrange(BH_p)
df2 <- read.csv(Oral_GO, check.names = F) %>% mutate(across(c(r1, r2), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2) %>% arrange(BH_p)
df3 <- read.csv(Oral_MTB_NEG, check.names = F) %>% mutate(across(c(r1, r2), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2) %>% arrange(BH_p)
df4 <- read.csv(Oral_MTB_POS, check.names = F) %>% mutate(across(c(r1, r2), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2) %>% arrange(BH_p)
df5 <- read.csv(Plasma_MTB_NEG, check.names = F) %>% mutate(across(c(r1, r2), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2) %>% arrange(BH_p)
df6 <- read.csv(Plasma_MTB_POS, check.names = F) %>% mutate(across(c(r1, r2), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2) %>% arrange(BH_p)
df7 <- read.csv(Gut_Species, check.names = F) %>% mutate(across(c(r1, r2), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2) %>% arrange(BH_p)

make_supplementary_tables(name = "Supplementary_Table8",
                          sheet_name = c("A.Oral Species","B.Oral GO","C.Oral Metabolome-NEG","D.Oral Metabolome-POS","E.Plasma Metabolome-NEG","F.Plasma Metabolome-POS","G.Gut Species"),
                          data = list(df1,df2,df3,df4,df5,df6, df7),output_dir)

# Supplemental Table 9, significant differential correlation for oral species
Oral_Species <- "Outputs/Correlations_trend_analysis/oral_species_cor_diff_sig_with_trend_test.csv"
Oral_GO <- "Outputs/Correlations_trend_analysis/oral_go_cor_diff_sig_with_trend_test.csv"
Oral_MTB_NEG <- "Outputs/Correlations_trend_analysis/oral_meta.n_cor_diff_sig_with_trend_test.csv"
Oral_MTB_POS <- "Outputs/Correlations_trend_analysis/oral_meta.p_cor_diff_sig_with_trend_test.csv"
Plasma_MTB_NEG <- "Outputs/Correlations_trend_analysis/oral_plasma.meta.n_cor_diff_sig_with_trend_test.csv"
Plasma_MTB_POS <- "Outputs/Correlations_trend_analysis/oral_plasma.meta.p_cor_diff_sig_with_trend_test.csv"
Gut_Species <- "Outputs/Correlations_trend_analysis/oral_gut_species_cor_diff_sig_with_trend_test.csv"

#df1 <- read.csv(Oral_Species, check.names = F) %>% mutate(across(c(r1, r2,`r1.g1`,`r2.g2`,`r3.g3`,`r4.g4`), ~ round(.x, 2))) %>% arrange(`adj.p.trend`)
df2 <- read.csv(Oral_GO, check.names = F) %>% select(-c("num_sig_corr_t1","num_sig_corr_t2"))  %>% mutate(across(c(r1, r2,`r1.g1`,`r2.g2`,`r3.g3`,`r4.g4`), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2,`ρG1`=`r1.g1`,`ρG2`=`r2.g2`,`ρG3`=`r3.g3`,`ρG4`=`r4.g4`) %>% arrange(`adj.p.trend`)
df3 <- read.csv(Oral_MTB_NEG, check.names = F) %>% select(-c("num_sig_corr_t1","num_sig_corr_t2"))  %>% mutate(across(c(r1, r2,`r1.g1`,`r2.g2`,`r3.g3`,`r4.g4`), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2,`ρG1`=`r1.g1`,`ρG2`=`r2.g2`,`ρG3`=`r3.g3`,`ρG4`=`r4.g4`) %>% arrange(`adj.p.trend`)
df4 <- read.csv(Oral_MTB_POS, check.names = F) %>% select(-c("num_sig_corr_t1","num_sig_corr_t2"))  %>% mutate(across(c(r1, r2,`r1.g1`,`r2.g2`,`r3.g3`,`r4.g4`), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2,`ρG1`=`r1.g1`,`ρG2`=`r2.g2`,`ρG3`=`r3.g3`,`ρG4`=`r4.g4`) %>% arrange(`adj.p.trend`)
df5 <- read.csv(Plasma_MTB_NEG, check.names = F) %>% select(-c("num_sig_corr_t1","num_sig_corr_t2"))  %>% mutate(across(c(r1, r2,`r1.g1`,`r2.g2`,`r3.g3`,`r4.g4`), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2,`ρG1`=`r1.g1`,`ρG2`=`r2.g2`,`ρG3`=`r3.g3`,`ρG4`=`r4.g4`) %>% arrange(`adj.p.trend`)
df6 <- read.csv(Plasma_MTB_POS, check.names = F) %>% select(-c("num_sig_corr_t1","num_sig_corr_t2"))  %>% mutate(across(c(r1, r2,`r1.g1`,`r2.g2`,`r3.g3`,`r4.g4`), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2,`ρG1`=`r1.g1`,`ρG2`=`r2.g2`,`ρG3`=`r3.g3`,`ρG4`=`r4.g4`) %>% arrange(`adj.p.trend`)
df7 <- read.csv(Gut_Species, check.names = F) %>% select(-c("num_sig_corr_t1","num_sig_corr_t2"))  %>% mutate(across(c(r1, r2,`r1.g1`,`r2.g2`,`r3.g3`,`r4.g4`), ~ round(.x, 2))) %>% 
  dplyr::rename(ρ1=r1, ρ2=r2,`ρG1`=`r1.g1`,`ρG2`=`r2.g2`,`ρG3`=`r3.g3`,`ρG4`=`r4.g4`) %>% arrange(`adj.p.trend`)

make_supplementary_tables(name = "Supplementary_Table9",
                          sheet_name = c("A.Oral GO","B.Oral Metabolome-NEG","C.Oral Metabolome-POS","D.Plasma Metabolome-NEG","E.Plasma Metabolome-POS","F.Gut Species"),
                          data = list(df2,df3,df4,df5,df6, df7),output_dir)

# Supplemental Table 10
datasets <- c("gut_species","gut_go","gut_meta.n","gut_meta.p","gut_plasma.meta.n","gut_plasma.meta.p",
              "oral_species","oral_go","oral_meta.n","oral_meta.p","oral_plasma.meta.n","oral_plasma.meta.p","oral_gut_species")
list_res <- lapply(datasets, function(d){
  read.csv(file.path(file.path("Outputs", "Differential_Correlations","Tables"), paste0(d,"_dysco.csv"))) %>% 
    arrange(desc(t)) %>% dplyr::rename(`number of correlations` = num_corr)
})
datasets <- c("Gut Species","Gut GO","Gut MTB-NEG","Gut MTB-POS","Plasma MTB-NEG","Plasma MTB-POS",
                     "Oral Species","Oral GO","Oral MTB-NEG","Oral MTB-POS","Plasma MTB-NEG","Plasma MTB-POS", "Gut Species")


make_supplementary_tables(name = "Supplementary_Table10",
                          sheet_name = paste0(LETTERS[1:length(datasets)],".",datasets),
                          data = list_res,output_dir)



# Supplemental table 11
Fulcher_neg_species <- read.csv("Outputs/External_Dataset/Fulcher_2022/Negative/Fulcher_2022_Negative_gut_species_dysco.csv", check.names = F) %>% 
  arrange(desc(t)) %>% dplyr::rename(`number of correlations` = num_corr,taxa = species)
Fulcher_neg_go <- read.csv("Outputs/External_Dataset/Fulcher_2022/Negative/Fulcher_2022_Negative_gut_go_dysco.csv", check.names = F) %>% 
  arrange(desc(t)) %>% dplyr::rename(`number of correlations` = num_corr,taxa = species)
Fulcher_pos_species <- read.csv("Outputs/External_Dataset/Fulcher_2022/Positive/Fulcher_2022_Positive_gut_species_dysco.csv", check.names = F) %>% 
  arrange(desc(t)) %>% dplyr::rename(`number of correlations` = num_corr,taxa = species)
Fulcher_pos_go <- read.csv("Outputs/External_Dataset/Fulcher_2022/Positive/Fulcher_2022_Positive_gut_go_dysco.csv", check.names = F) %>% 
  arrange(desc(t)) %>% dplyr::rename(`number of correlations` = num_corr,taxa = species)
Garcia_species <- read.csv("Outputs/External_Dataset/Garcia_2024/Garcia_2024_gut_species_dysco.csv", check.names = F) %>% 
  arrange(desc(t)) %>% dplyr::rename(`number of correlations` = num_corr,taxa = species)
Garcia_go <- read.csv("Outputs/External_Dataset/Garcia_2024/Garcia_2024_gut_go_dysco.csv", check.names = F) %>% 
  arrange(desc(t)) %>% dplyr::rename(`number of correlations` = num_corr,taxa = species)
boston_msm <- read.csv(file.path("Outputs/External_Dataset/Rocafort_2024/boston/MSM", paste0("Rocafort_boston_MSM", "_gut_species_dysco.csv")), check.names = F) %>% 
  arrange(desc(t)) %>% dplyr::rename(`number of correlations` = num_corr,taxa = species)
boston_hetr <- read.csv(file.path("Outputs/External_Dataset/Rocafort_2024/boston/heterosexual", paste0("Rocafort_boston_heterosexual", "_gut_species_dysco.csv")), check.names = F) %>% 
  arrange(desc(t)) %>% dplyr::rename(`number of correlations` = num_corr,taxa = species)
bostwana_hetr <- read.csv(file.path("Outputs/External_Dataset/Rocafort_2024/botswana/heterosexual", paste0("Rocafort_botswana_heterosexual", "_gut_species_dysco.csv")), check.names = F) %>% 
  arrange(desc(t)) %>% dplyr::rename(`number of correlations` = num_corr,taxa = species)
uganda_hetr <- read.csv(file.path("Outputs/External_Dataset/Rocafort_2024/uganda_2/heterosexual", paste0("Rocafort_uganda_2_heterosexual", "_gut_species_dysco.csv")), check.names = F) %>%
  arrange(desc(t)) %>% dplyr::rename(`number of correlations` = num_corr,taxa = species)
no_ART <- read.csv(file.path("Outputs/External_Dataset/Armstrong_2018/MSM_HIV+ART-", paste0("Armstrong_MSM_HIV+ART-", "_gut_species_dysco.csv")), check.names = F) %>% 
  arrange(desc(t)) %>% dplyr::rename(`number of correlations` = num_corr,taxa = species)
with_ART <- read.csv(file.path("Outputs/External_Dataset/Armstrong_2018/MSM_HIV+ART", paste0("Armstrong_MSM_HIV+ART", "_gut_species_dysco.csv")), check.names = F) %>% 
  arrange(desc(t)) %>% dplyr::rename(`number of correlations` = num_corr,taxa = species)

datasets <- c(
  "Fulcher_Pre-HIV_species", "Fulcher_Pre-HIV_GO", "Fulcher_Post-HIV_species", "Fulcher_Post-HIV_GO",
  "Garcia_species", "Garcia_GO",
  "Rocafort_Boston_MSM", "Rocafort_Boston_non-MSM", "Rocafort_Bostwana_non-MSM", "Rocafort_Uganda_non-MSM",
  "Armstrong_no_ART", "Armstrong_with_ART"
)
make_supplementary_tables(name = "Supplementary_Table11",
                          sheet_name = paste0(LETTERS[1:length(datasets)],".",datasets),
                          data = list(Fulcher_neg_species,Fulcher_neg_go,Fulcher_pos_species,Fulcher_pos_go,
                                      Garcia_species,Garcia_go,
                                      boston_msm,boston_hetr,bostwana_hetr,uganda_hetr,
                                      no_ART,with_ART),output_dir)

