# stably expressed gene within seedlings

source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Code/catchToList.R")
setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Data/arab/Rsubread")
source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Code/create.geneColumn.R")

GSE36626 <- read.table("GSE36626.Rsubread.txt")
GSE36626 <- create.geneColumn(GSE36626)
trt.GSE36626 <- c(rep(1:2, 2))  

GSE39463 <- read.table("GSE39463.Rsubread.txt")
GSE39463 <- GSE39463[, -c(3, 7, 9)]  # remove duplicated runs
GSE39463 <- create.geneColumn(GSE39463)
trt.GSE39463 <- c(rep(1:4, each= 3))

GSE48235 <- read.table("GSE48235.Rsubread.txt", header = T) 
GSE48235 <- create.geneColumn(GSE48235)
trt.GSE48235 <- c(rep(1:3, each=2))

GSE51304 <- read.table("GSE51304.Rsubread.txt", header=T) # seedling
GSE51304 <- create.geneColumn(GSE51304)
trt.GSE51304 <- c(rep(1:9, each=2))

GSE54677 <- read.table("GSE54677.Rsubread.txt", header = T)
GSE54677 <- create.geneColumn(GSE54677)
trt.GSE54677 <- c(rep(1:10, 2))  ### be careful about the design structure


### leaf   ---------------------------------------------------

library(plyr)
ls <- list(GSE36626, GSE39463, GSE48235, GSE51304, GSE54677)

stable <- join_all(ls, "Gene")
stable <- stable[complete.cases(stable),]


#as.numeric(colSums(stable[, -1]))
row.names(stable) <- stable[, 1]
stable <- stable[, -1]
stable2 <- stable[rowSums(stable) >0, ]

# -------------------------------------------
# filtering by 3 for each sample
# ------------------------------------------

library(NBPSeq)
stable <- as.matrix(stable[rowSums(stable)>=3*dim(stable)[2], ])
dim(stable)
norm.factors <- estimate.norm.factors(stable)
# norm.factors <- estimate.norm.factors(stable, method=NULL)
nb.data <- prepare.nb.data(stable, norm.factors=norm.factors)
off.set <- as.numeric(nb.data$eff.lib.sizes)

# counts per million
# y <- sweep(stable, 2, off.set, "/")*1e6
# 
# id <- which(y[, ]< 1)
# head(id)
# y2 <- y[-id, ]


trt <- as.factor(c(trt.GSE36626, trt.GSE39463, trt.GSE48235,# trt.GSE42968,
                   trt.GSE51304, trt.GSE54677))
n.obs <- c(length(trt.GSE36626), length(trt.GSE39463), length(trt.GSE48235),
           length(trt.GSE51304), length(trt.GSE54677))

group <- as.factor(rep(1:length(ls), n.obs))

library(lme4)


#### run regression models for all genes 
warn <- rep(0, dim(stable)[1])
var1 <- var2 <- var3 <- c()
system.time(for ( i in 1:dim(stable)[1]){
  y <- as.numeric(stable[i, ])
  id <- 1:length(y)
  # print(i)
  mod1 <-  catchToList( glmer(y ~ 1  + (1|group) + (1|group:trt) + (1|id),
                              offset=(log(off.set)), 
                              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)),
                              family = poisson)
  ) # catch warnings or errors
  var3[i] <- as.numeric(mod1$value$varcor[1]) # variance for individual random effect
  var2[i] <- as.numeric(mod1$value$varcor[2]) # variance for trt effect
  var1[i] <- as.numeric(mod1$value$varcor[3]) # variance for group effect
  
  if (mod1$warnings[length(mod1$warnings)] != "NO"  ) ## there is warning
  {
    warn[i] = 1 # if the relative gradients are still large, then give a warning
    
  } 
})

sum <- var1 + var2 + var3
rank.sum <- rank(sum)
var.comp<- data.frame(gene= rownames(stable),individual =var3,
                      trt=var2, lab=var1, warnS=warn, sum=sum, rank.sum=rank.sum)

write.table(var.comp, "/Users/Bin/Dropbox/Zhuo/Research/Project2014/Result/leaf_5lab_Rsubread.txt")
