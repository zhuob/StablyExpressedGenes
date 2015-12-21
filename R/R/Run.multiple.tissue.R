###  run variance component analysis for multiple- tissue data

library(NormqPCR)

source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/R/R/estimate.var.component.R")
source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/R/R/plot.genes.R")
source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/R/R/Compare.gene.R")
setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/R/data/")




cat("Read the RNA-Seq data ...")
tissue <- readRDS("tissue.rds")
obj <- tissue

####### calculate the variance components ----------------------------
cat("Begin iteration 0")
var.obj <- estimate.var.component(data = obj$count[1:3000, ], reference = NULL, obj$lab, obj$trt, filter.factor = 3)
refe <- var.obj$var.comp$Gene[var.obj$var.comp$Rank <=  1000]
cat("Use top 1000 stably expressed genes identified from GLMM \n to do normalization, and re-run the GLMM")
var.obj2 <- estimate.var.component(data = obj$count, reference = refe, obj$lab, obj$trt, filter.factor = 3)

var.obj2 <- var.obj
####### rank the genes by geNorm and NormFinder ---------------------
cat("Rank the genes by V values (geNorm procedure)")
obj <- var.obj2 
geNorm_Normfinder <- rankNorm(obj)


#######  Read all the data for results ----------------------\
var_tissue <- var.obj2
cze_100 <- read.table("Czechowski100.txt", header=T)
dek_50 <- read.table("Dekkers50.txt", header=T)

#  Table 3.1: basic statistics of the data
#  cat( c("# genes", "# samples", "# experiment", "# treatments"))
#  cat(c(dim(var_tissue$count), length(unique(var_tissue$lab)), length(unique(paste(var_tissue$lab, var_tissue$trt, sep="_")))))

print( match.gene(dek_50,  var_tissue, 100) )



