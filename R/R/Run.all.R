### Run all the results

project_folder <- getwd() 
DataPath <- paste(project_folder, "/R/data/", sep ="")                 # where are the pre-processed RNA-Seq data
ResultPath <-  paste(project_folder, "/Results/", sep ="")             # where are the results from running GLMM
CodePath <-  paste(project_folder, "/R/R/", sep ="")                   # where are the code for producing the results
FigurePath <-  paste(project_folder, "/Manuscript/Figures/", sep ="")  # where the figures will be stored
SuppPath <-  paste(project_folder, "/Supplementary/", sep ="")         # where the supplementary files will be stored

source(paste(CodePath, "Run.GLMM.R", sep = "")) 
source(paste(CodePath, "Run.geNorm.R", sep =""))
source(paste(CodePath, "Run.Results.R", sep =""))
