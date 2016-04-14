###  run variance component analysis for multiple- tissue, leaf and seedling data


CodeNeeded <- paste(CodePath, "estimate.var.component.R", sep="")
source(CodeNeeded)

# dataname <- c("tissue", "leaf", "seedling")

dataname <- "seedling"
for ( i in 1:length(dataname)){
  
  dataobj <- paste(DataPath, dataname[i], ".rds", sep="")
  obj <- readRDS(dataobj)
  niters = 1

  ######### RUN GLMM   
var.obj <- iteration.GLMM(
obj,                                              # the data object
niter = niters,                                   # number of iteration for DESeq normalization, default is 1.
filter.factor= 3,                                 # average row means below which the genes are removed, default is 3.
topgene = 1000,                                   # number of top stably expressed genes to be used for normalization, default is 1000
save.iter0 = F
)
save.target <- paste(ResultPath, dataname[i], ".columbia.use.reference.iter.", niters, ".rds", sep="")
saveRDS(var.obj, save.target)
}


