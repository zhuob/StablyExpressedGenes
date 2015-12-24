

## Run geNorm V values for each gene
## need the CodePath and ResultPath

CodeNeeded <- paste(CodePath, "Compare.gene.R", sep="")
source(CodeNeeded)


data <- paste(ResultPath, "tissue.columbia.use.reference.iter.1.rds", sep ="")
obj <- readRDS(data)
geNorm <- rankMvalue(obj)
save.path <- paste(ResultPath, "geNorm_rank_tissue.rds", sep="")
saveRDS(geNorm2, save.path)

