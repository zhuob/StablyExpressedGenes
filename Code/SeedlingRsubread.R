# stably expressed gene within seedlings

source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Code/catchToList.R")
setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Data/arab/Rsubread")
source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Code/create.geneColumn.R")

GSE32216 <- read.table("GSE32216.Rsubread.txt")
GSE32216 <- create.geneColumn(GSE32216)
trt.GSE32216 <- c(rep(1, 2), rep(2, 2))

GSE37159 <- read.table("GSE37159.Rsubread.txt", header = T) 
GSE37159 <- create.geneColumn(GSE37159)
trt.GSE37159 <- c(rep(1:4, each=2))

GSE38400 <- read.table("GSE38400.Rsubread.txt")
GSE38400 <- create.geneColumn(GSE38400)
trt.GSE38400 <- c(rep(1:4, each=3))


GSE38879 <- read.table("GSE38879.Rsubread.txt", header = T)
GSE38879 <- GSE38879[, seq(1, 23, by =2)] # 2 ILLUMINA RUNS
GSE38879 <- create.geneColumn(GSE38879)
trt.GSE38879 <- c(rep(1:4, each =3))

GSE39214 <- read.table("GSE39214.Rsubread.txt", header=T) # seedling
GSE39214 <- create.geneColumn(GSE39214)
trt.GSE39214 <- c(rep(1:4, each=3))

GSE43865 <- read.table("GSE43865.Rsubread.txt", header=T) # seedling
GSE43865 <- create.geneColumn(GSE43865)
trt.GSE43865 <- c(rep(1:2, each=3))

# I proceeded this data myself
GSE48767 <- read.table("GSE48767.Rsubread.txt", header = T)
# GSE48767 <- GSE48767[, -(1:6)]
GSE48767 <- create.geneColumn(GSE48767)
trt.GSE48767 <- rep(1:4, each=3)


GSE51119 <- read.table("GSE51119.Rsubread.txt", header = T) 
GSE51119 <- create.geneColumn(GSE51119)
trt.GSE51119 <- rep(1:5, each=2)

GSE51772 <- read.table("GSE51772.Rsubread.txt", header = T) 
GSE51772 <- create.geneColumn(GSE51772)
trt.GSE51772 <- rep(1:4, each=2)

GSE53078 <- read.table("GSE53078.Rsubread.txt", header = T) # included
GSE53078 <- create.geneColumn(GSE53078)
trt.GSE53078 <- rep(1:2, each=2)


GSE58082 <- read.table("GSE58082.Rsubread.txt")
GSE58082 <- create.geneColumn(GSE58082)
trt.GSE58082 <- rep(1:2, each =3)


GSE57806 <- read.table("GSE57806.Rsubread.txt")
GSE57806 <- create.geneColumn(GSE57806)
trt.GSE57806 <- rep(1:2, each =3)

GSE43703 <- read.table("GSE43703.Rsubread.txt")
GSE43703 <- create.geneColumn(GSE43703)
trt.GSE43703 <- rep(1:4, 2)  # this data behaves weild in normalization factor



### seedling ---------------------------------------------------

# get rid of GSE39214, because it looks strange

library(plyr)
ls <- list(#GSE32216, 
           GSE37159, GSE38400,# number 10 (this is strain, not seedling)
           GSE38879, #GSE39214,
           GSE43865, GSE48767,# add 6 more
           GSE51119, GSE51772, GSE53078, GSE58082,
           GSE57806) # number 9
          # , GSE43703)

stable <- join_all(ls, "Gene")
stable <- stable[complete.cases(stable),]
dim(stable)
#as.numeric(colSums(stable[, -1]))
row.names(stable) <- stable[, 1]
stable <- stable[, -1]


# -------------------------------------------
# filtering by 3 for each sample
# ------------------------------------------

library(NBPSeq)
stable <- as.matrix(stable[rowSums(stable)>=3*dim(stable)[2], ])
 norm.factors <- estimate.norm.factors(stable)
# norm.factors <- estimate.norm.factors(stable, method=NULL)
nb.data <- prepare.nb.data(stable, norm.factors=norm.factors)
off.set <- as.numeric(nb.data$eff.lib.sizes)

trt <- as.factor(c(trt.GSE37159, #trt.GSE38400, 
                   trt.GSE38879,# trt.GSE42968,
                   trt.GSE43865,
                   trt.GSE48767, trt.GSE51119, trt.GSE51772, trt.GSE53078, trt.GSE58082, 
                   trt.GSE57806))
n.obs <- c(length(trt.GSE37159), #length(trt.GSE38400), 
           length(trt.GSE38879),
           length(trt.GSE43865), length(trt.GSE48767), length(trt.GSE51119),
           length(trt.GSE51772), length(trt.GSE53078), length(trt.GSE58082),
           length(trt.GSE57806))

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

write.table(var.comp, "/Users/Bin/Dropbox/Zhuo/Research/Project2014/Result/seedling_10lab_Rsubread.txt")
