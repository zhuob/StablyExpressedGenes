### Run all the results

# note: the following packages are required to reproduce the results

library(ggplot2)
library(reshape2)
library(scales)
library(NBPSeq)
library(GGally)
library(gridExtra)
library(grid)
library(cowplot)
library(VennDiagram)
library(lme4)



project_folder <- getwd() 
RawDataPath <- paste(project_folder, "/R/data-raw", sep ="")           # where all the raw counts data are stored
DataPath <- paste(project_folder, "/R/data/", sep ="")                 # where are the RNA-Seq data chosen after normalization filter
ResultPath <-  paste(project_folder, "/Results/", sep ="")             # where are the results from running GLMM
CodePath <-  paste(project_folder, "/R/R/", sep ="")                   # where are the code for producing the results
FigurePath <-  paste(project_folder, "/Results/Figures/", sep = "")    # where the figures will be stored
SuppPath <-  paste(project_folder, "/Supplementary/", sep ="")         # where the supplementary files will be stored

source(paste(CodePath, "SelectingData.R", sep = ""))                   # creating the three data sets
source(paste(CodePath, "Run.GLMM.R", sep = ""))                        # fit glmm to the data sets
source(paste(CodePath, "Run.geNorm.R", sep =""))                       # run geNorm (without iteration to all genes in multi-tissue group)
source(paste(CodePath, "Run.Results.R", sep =""))                      # all the figures and results needed in the paper.
