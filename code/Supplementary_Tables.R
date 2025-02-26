rm(list = ls())
output_dir <- file.path("Supplementary_Tables")
dir.create(output_dir)

# Supplemental Table 1
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
  filter(!grepl("oral",SCFA))


xlsx::write.xlsx(df1, file=file.path(output_dir,"Supplementary_Table1.xlsx"),
                 sheetName="A.Gut Species", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df2, file=file.path(output_dir,"Supplementary_Table1.xlsx"),
                 sheetName="B.Gut GO", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df3, file=file.path(output_dir,"Supplementary_Table1.xlsx"),
                 sheetName="C.Gut Metabolites", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df4, file=file.path(output_dir,"Supplementary_Table1.xlsx"),
                 sheetName="D.Plasma Metabolites", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df5, file=file.path(output_dir,"Supplementary_Table1.xlsx"),
                 sheetName="E.SCFA", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)


# Supplemental Table 3
Gut_Species <- "Outputs/Dysbiosis/Gut/Gut_Species/species_differential_correlations_sig.csv"
Gut_GO <- "Outputs/Dysbiosis/Gut/Gut_GO/go_differential_correlations_sig.txt"
Gut_MTB_NEG <- "Outputs/Dysbiosis/Gut/Gut_MTB_NEG/mtb_neg_differential_correlations_sig.txt"
Gut_MTB_POS <- "Outputs/Dysbiosis/Gut/Gut_MTB_POS/mtb_pos_differential_correlations_sig.txt"
Plasma_MTB_NEG <- "Outputs/Dysbiosis/Gut/Plasma_MTB_NEG/plasma_mtb_neg_differential_correlations_sig.txt"
Plasma_MTB_POS <- "Outputs/Dysbiosis/Gut/Plasma_MTB_POS/plasma_mtb_pos_differential_correlations_sig.txt"


df1 <- read.table(Gut_Species, sep = ",",header = T,check.names = F,comment.char = "",quote = "")
df2 <- read.table(Gut_GO, sep = "\t",header = T,check.names = F,comment.char = "",quote = "")
df3 <- read.table(Gut_MTB_NEG, sep = "\t",header = T,check.names = F,comment.char = "",quote = "")
df4 <- read.table(Gut_MTB_POS, sep = "\t",header = T,check.names = F,comment.char = "",quote = "")
df5 <- read.table(Plasma_MTB_NEG, sep = "\t",header = T,check.names = F,comment.char = "",quote = "")
df6 <- read.table(Plasma_MTB_POS, sep = "\t",header = T,check.names = F,comment.char = "",quote = "")


xlsx::write.xlsx(df1, file=file.path(output_dir,"Supplementary_Table3.xlsx"),
                 sheetName="A.Gut Species", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df2, file=file.path(output_dir,"Supplementary_Table3.xlsx"),
                 sheetName="B.Gut GO", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df3, file=file.path(output_dir,"Supplementary_Table3.xlsx"),
                 sheetName="C.Gut Metabolome-NEG", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df4, file=file.path(output_dir,"Supplementary_Table3.xlsx"),
                 sheetName="D.Gut Metabolome-POS", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df5, file=file.path(output_dir,"Supplementary_Table3.xlsx"),
                 sheetName="E.Plasma Metabolome-NEG", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df6, file=file.path(output_dir,"Supplementary_Table3.xlsx"),
                 sheetName="F.Plasma Metabolome-POS", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)


# Supplemental Table 4
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
  


xlsx::write.xlsx(df1, file=file.path(output_dir,"Supplementary_Table4.xlsx"),
                 sheetName="A.Oral Species", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df2, file=file.path(output_dir,"Supplementary_Table4.xlsx"),
                 sheetName="B.Oral GO", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df3, file=file.path(output_dir,"Supplementary_Table4.xlsx"),
                 sheetName="C.Oral SCFA", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)


# Supplemental Table 5
Oral_Species <- "Outputs/Dysbiosis/Oral/Oral_Species/species_differential_correlations_sig.csv"
Oral_GO <- "Outputs/Dysbiosis/Oral/Oral_GO/go_differential_correlations_sig.txt"
Oral_MTB_NEG <- "Outputs/Dysbiosis/Oral/Oral_MTB_NEG/mtb_neg_differential_correlations_sig.txt"
Oral_MTB_POS <- "Outputs/Dysbiosis/Oral/Oral_MTB_POS/mtb_pos_differential_correlations_sig.txt"
Plasma_MTB_NEG <- "Outputs/Dysbiosis/Oral/Plasma_MTB_NEG/plasma_mtb_neg_differential_correlations_sig.txt"
Plasma_MTB_POS <- "Outputs/Dysbiosis/Oral/Plasma_MTB_POS/plasma_mtb_pos_differential_correlations_sig.txt"
Gut_Species <- "Outputs/Dysbiosis/Oral/Gut_Species/species2_differential_correlations_sig.csv"

df1 <- read.table(Oral_Species, sep = ",",header = T,check.names = F,comment.char = "",quote = "")
df2 <- read.table(Oral_GO, sep = "\t",header = T,check.names = F,comment.char = "",quote = "")
df3 <- read.table(Oral_MTB_NEG, sep = "\t",header = T,check.names = F,comment.char = "",quote = "")
df4 <- read.table(Oral_MTB_POS, sep = "\t",header = T,check.names = F,comment.char = "",quote = "")
df5 <- read.table(Plasma_MTB_NEG, sep = "\t",header = T,check.names = F,comment.char = "",quote = "")
df6 <- read.table(Plasma_MTB_POS, sep = "\t",header = T,check.names = F,comment.char = "",quote = "")
df7 <- read.table(Gut_Species, sep = ",",header = T,check.names = F,comment.char = "",quote = "")


xlsx::write.xlsx(df1, file=file.path(output_dir,"Supplementary_Table5.xlsx"),
                 sheetName="A.Oral Species", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df2, file=file.path(output_dir,"Supplementary_Table5.xlsx"),
                 sheetName="B.Oral GO", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df3, file=file.path(output_dir,"Supplementary_Table5.xlsx"),
                 sheetName="C.Oral Metabolome-NEG", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df4, file=file.path(output_dir,"Supplementary_Table5.xlsx"),
                 sheetName="D.Oral Metabolome-POS", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df5, file=file.path(output_dir,"Supplementary_Table5.xlsx"),
                 sheetName="E.Plasma Metabolome-NEG", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df6, file=file.path(output_dir,"Supplementary_Table5.xlsx"),
                 sheetName="F.Plasma Metabolome-POS", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)
xlsx::write.xlsx(df7, file=file.path(output_dir,"Supplementary_Table5.xlsx"),
                 sheetName="G.Gut Species", row.names=FALSE,append = TRUE,col.names = TRUE,showNA = FALSE)

