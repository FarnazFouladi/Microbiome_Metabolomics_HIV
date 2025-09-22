source("code/load.R")
source("code/create_heatmap.R")

input_dir <- file.path("Outputs", "Differential_Correlations","Tables")
output_dir <- file.path("Outputs", "Differential_Correlations","Heatmaps")
dir.create(output_dir,recursive = T,showWarnings = F)

########### Heatmap for disco score
for (sample_type in c("Gut", "Oral")) {
  if (sample_type == "Gut") {
    
    datasets <- c("gut_species","gut_go","gut_meta.p","gut_meta.n","gut_plasma.meta.p","gut_plasma.meta.n")
    list_res <- lapply(datasets, function(d){
      read.csv(file.path(input_dir, paste0(d,"_disco.csv")))
    })
    names(list_res) <- c("Gut Species","Gut GO","Gut MTB-POS","Gut MTB-NEG","Plasma MTB-POS","Plasma MTB-NEG")
  } else {
    
    datasets <- c("oral_species","oral_go","oral_meta.p","oral_meta.n","oral_plasma.meta.p","oral_plasma.meta.n","oral_gut_species")
    list_res <- lapply(datasets, function(d){
      read.csv(file.path(input_dir, paste0(d,"_disco.csv")))
    })
    names(list_res) <- c("Oral Species","Oral GO","Oral MTB-POS","Oral MTB-NEG","Plasma MTB-POS","Plasma MTB-NEG", "Gut Species")
  }
  
  
  create_heatmap (list_res, p_alpha = 0.01,bh_alpha = 0.1,col_names = names(list_res),top_n = 10,
                  sample_type = sample_type,output_dir = output_dir)
  
}
