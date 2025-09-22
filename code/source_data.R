# Creating Source Data for the figures
# Library
library(tidyverse)
library(openxlsx)

# Figure 1
fig1a <- read.table("Outputs/Differential_Abundance_Microbiome/Gut/Species/Species_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="") 
fig1b <- read.table("Outputs/Differential_Abundance_Microbiome/Gut/GO/GO_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="") 

fig1c <- read.table("Outputs/Differential_Abundance_Metabolomics/Gut/metabolites_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="") 
fig1d <- read.table("Outputs/Differential_Abundance_Metabolomics/Plasma/metabolites_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="") 

fig1e <- read.table("Outputs/Differential_Abundance_SCFA/Fecal_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="") 
fig1f <- read.table("Outputs/Differential_Abundance_SCFA/Plasma_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="") 

# Figure 2:
fig2 <- xlsx::read.xlsx(file.path( "Outputs/Diablo/Gut", "Supplementary_Table3.xlsx"),sheetIndex = 1)
fig2 <- fig2 %>% filter(abs(coefficient) >=0.8)

# Figure 3:
fig3a <- read.table("Outputs/Differential_Correlations/Networks/Gut/Holdemanella_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="") 
fig3b <- read.table("Outputs/Differential_Correlations/Networks/Gut/go_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="") 

fig3c <- read.table("Outputs/Differential_Correlations/Networks/Gut/gut_mtb.p_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="") 

fig3d <- read.table("Outputs/Differential_Correlations/Networks/Gut/plasma_mtb.n_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="") 

# Figure 4
fig4a <- read.table("Outputs/Differential_Correlations/Heatmaps/Gut_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Species")
  
fig4b <- read.table("Outputs/Differential_Correlations/Heatmaps/Oral_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Species")

# Supplementary Figure 1
supp.fig1a <- read.table("Outputs/Differential_Abundance_Microbiome_SA/Gut/Species/Species_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") 
supp.fig1b <- read.table("Outputs/Differential_Abundance_Microbiome_SA/Gut/GO/GO_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") 
supp.fig1c <- read.table("Outputs/Differential_Abundance_Metabolomics_SA/Gut/metabolites_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") 
supp.fig1d <- read.table("Outputs/Differential_Abundance_Metabolomics_SA/Plasma/metabolites_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") 


# Supplementary Figure 2
supp.fig2a <- read.table("Outputs/Differential_Abundance_Comparisons/gut_species_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Species")
supp.fig2b <- read.table("Outputs/Differential_Abundance_Comparisons/gut_go_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("GO")
supp.fig2c <- read.table("Outputs/Differential_Abundance_Comparisons/gut_metabolomics_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Metabolites")
supp.fig2d <- read.table("Outputs/Differential_Abundance_Comparisons/plasma_metabolomics_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Metabolites")

# Supplementary Figure 3
supp.fig3a <- read.table("Outputs/Differential_Correlations/Networks/Gut/Holdemanella_trend_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Feature")
supp.fig3b <- read.table("Outputs/Differential_Correlations/Networks/Gut/go_trend_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Feature")
supp.fig3c <- read.table("Outputs/Differential_Correlations/Networks/Gut/gut_mtb.p_trend_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Feature")
supp.fig3d <- read.table("Outputs/Differential_Correlations/Networks/Gut/gut_plasma.meta.n_trend_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Feature")

# Supplementary Figure 4
supp.fig4a <- read.table("Outputs/Differential_Abundance_Microbiome/Oral/Species/Species_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="",na.strings = "NA") 
supp.fig4b <- read.table("Outputs/Differential_Abundance_Microbiome/Oral/GO/GO_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="",na.strings = "NA") 
supp.fig4c <- read.table("Outputs/Differential_Abundance_SCFA/Oral_source_data.txt",
                    sep="\t",header = T, check.names = F, quote = "",comment.char ="") 


# Supplementary Figure 5
supp.fig5a <- read.table("Outputs/Differential_Abundance_Microbiome_SA/Oral/Species/Species_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") 
supp.fig5b <- read.table("Outputs/Differential_Abundance_Microbiome_SA/Oral/GO/GO_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") 

# Supplementary Figure 6
supp.fig6a <- read.table("Outputs/Differential_Abundance_Comparisons/oral_species_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Species")
supp.fig6b <- read.table("Outputs/Differential_Abundance_Comparisons/oral_go_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("GO")

# Supplementary Figure 7
supp.fig7a <- read.table("Outputs/Differential_Correlations/Networks/Oral/oral_go_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") 
supp.fig7b <- read.table("Outputs/Differential_Correlations/Networks/Oral/oral_mtb.n_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") 
supp.fig7c <- read.table("Outputs/Differential_Correlations/Networks/Oral/oral_plasma_mtb.n_source_data.txt",
                          sep="\t",header = T, check.names = F, quote = "",comment.char ="") 
supp.fig7d <- read.table("Outputs/Differential_Correlations/Networks/Oral/oral_gut_species_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") 

# Supplementary Figure 8
supp.fig8a <- read.table("Outputs/Differential_Correlations/Networks/Oral/oral_go_trend_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Feature")
supp.fig8b <- read.table("Outputs/Differential_Correlations/Networks/Oral/oral_meta.n_trend_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Feature")
supp.fig8c <- read.table("Outputs/Differential_Correlations/Networks/Oral/oral_plasma.meta.n_trend_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Feature")
supp.fig8d <- read.table("Outputs/Differential_Correlations/Networks/Oral/oral_gut_species_trend_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Feature")

# Supplementary Figure 9
supp.fig9a <- read.table("Outputs/External_Dataset/Fulcher_2022/Gut_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Taxa")
supp.fig9b <- read.table("Outputs/External_Dataset/Garcia_2024/Gut_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Taxa")
supp.fig9c <- read.table("Outputs/External_Dataset/Rocafort_2024/Gut_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Taxa")
supp.fig9d <- read.table("Outputs/External_Dataset/Armstrong_2018/Gut_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Taxa")


# Supplementary Figure 10
supp.fig10 <- read.table("Outputs/External_Dataset/Genus_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") %>%
  rownames_to_column("Genera")
# Supplementary Figure 11
supp.fig11 <- read.table("Outputs/Differential_Correlations_Validations/null_vs_real_source_data.txt",
                         sep="\t",header = T, check.names = F, quote = "",comment.char ="") 



make_supplementary_tables <- function(name,sheet_name,data,output_dir){
  wb <- createWorkbook(name)
  
  for (i in 1:length(data)){
    print(i)
    addWorksheet(wb, sheet_name[i])
    writeData(wb,sheet_name[i] ,data[[i]])
  }
  saveWorkbook(wb, file.path(output_dir,paste0(name,".xlsx")), overwrite = TRUE)
  
}


all_table_names <- c("fig1a","fig1b","fig1c","fig1d","fig1e" ,"fig1f",
                "fig2","fig3a","fig3b","fig3c","fig3d",
                "fig4a","fig4b",
                "supp.fig1a","supp.fig1b","supp.fig1c","supp.fig1d",
                "supp.fig2a","supp.fig2b","supp.fig2c","supp.fig2d",
                "supp.fig3a","supp.fig3b","supp.fig3c","supp.fig3d",
                "supp.fig4a","supp.fig4b","supp.fig4c",
                "supp.fig5a","supp.fig5b",
                "supp.fig6a","supp.fig6b",
                "supp.fig7a","supp.fig7b","supp.fig7c","supp.fig7d",
                "supp.fig8a","supp.fig8b","supp.fig8c","supp.fig8d",
                "supp.fig9a","supp.fig9b","supp.fig9c","supp.fig9d",
                "supp.fig10","supp.fig11")


all_tables <- lapply(all_table_names, function(x) get(x))
output_dir <- "Source_Data"
dir.create(output_dir,showWarnings = F)
make_supplementary_tables(name = "Source_Data",sheet_name = all_table_names, data =  all_tables,output_dir = "Source_Data")


