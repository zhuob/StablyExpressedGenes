### Run all the results

FigurePath <- "/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/Manuscript/Figures/"
TablePath <- "/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/Manuscript/Tables/"
CodePath <- "/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/R/R/"
DataPath <- "/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/R/data/"
ResultPath <- "/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/Results/"



source(paste(CodePath, "Run.GLMM.R", sep = "")) 
source(paste(CodePath, "Run.geNorm.R", sep =""))
source(paste(CodePath, "Run.Results.R", sep =""))
