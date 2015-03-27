setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Data/arab/Rsubread")
source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Code/create.geneColumn.R")

#flower
GSE35288 <- read.table("GSE35288.Rsubread.txt")
GSE35288 <- create.geneColumn(GSE35288)

GSE40256 <- read.table("GSE40256.Rsubread.txt")
GSE40256 <- create.geneColumn(GSE40256)

#leaves

GSE36626 <- read.table("GSE36626.Rsubread.txt")
GSE36626 <- create.geneColumn(GSE36626)

GSE48235 <- read.table("GSE48235.Rsubread.txt")
GSE48235 <- create.geneColumn(GSE48235)

GSE54677 <- read.table("GSE54677.Rsubread.txt")
GSE54677 <- create.geneColumn(GSE54677)

# seed
GSE53952 <- read.table("GSE53952.Rsubread.txt")
GSE53952 <- create.geneColumn(GSE53952)

# root
GSE44062 <- read.table("GSE44062.Rsubread.txt")
GSE44062 <- create.geneColumn(GSE44062)

# carpel
GSE56326 <- read.table("GSE56326.Rsubread.txt")
GSE56326 <- create.geneColumn(GSE56326)

# hypocoltyl
GSE35408 <- read.table("GSE35408.Rsubread.txt")
GSE35408 <- create.geneColumn(GSE35408)


trt.GSE35288 <- c(rep(1, 3), rep(2, 3))
trt.GSE40256 <- rep(1:4, each=2)
trt.GSE36626 <- rep(1:2, 2)
trt.GSE48235 <- rep(1:3, each=2)
trt.GSE54677 <- rep(1:10, 2)
trt.GSE53952 <- rep(1:3, each=3)
trt.GSE44062 <- rep(1:4, each=2)
trt.GSE56326 <- c(rep(1, 2), rep(2:3, each=3))
trt.GSE35408 <- rep(1:5, 2)

n.obs <- c(length(trt.GSE35288),  length(trt.GSE40256),length(trt.GSE36626),
           length(trt.GSE48235), length(trt.GSE54677), length(trt.GSE53952),
           length(trt.GSE44062), length(trt.GSE56326), length(trt.GSE35408))



### seedling ---------------------------------------------------
# special filtering
GSE53952 <- GSE53952[, c(1, 2:4,11:13,20:22 )]  # the first stage was selected
trt.GSE53952 <- rep(1:3, each=3)

GSE35288 <- GSE35288[, c(1, seq(2,dim(GSE35288)[2], by=3))]

# get rid of GSE39214, because it looks strange

library(plyr)
ls <- list(#GSE32216, 
  GSE35288, #flower
  # GSE40256, #GSE39214,
  GSE48235, #leaves
  GSE53952, # seed
  #GSE44062,
  GSE56326, # carpel
  GSE35408) # hypocoltyl

stable <- join_all(ls, "Gene")
stable <- stable[complete.cases(stable),]
dim(stable)
as.numeric(colSums(stable[, -1]))
row.names(stable) <- stable[, 1]
stable <- stable[, -1]


library(NBPSeq)
stable <- as.matrix(stable[rowSums(stable)>= dim(stable)[2], ])
norm.factors <- estimate.norm.factors(stable)
nb.data <- prepare.nb.data(stable, norm.factors=norm.factors)
off.set <- as.numeric(nb.data$eff.lib.sizes)

trt <- as.factor(c(trt.GSE35288, trt.GSE48235,
                   trt.GSE53952,
                   trt.GSE56326, trt.GSE35408))
n.obs <- c(length(trt.GSE35288), length(trt.GSE48235), 
           length(trt.GSE53952), length(trt.GSE56326), length(trt.GSE35408))
          # length(trt.GSE51772), length(trt.GSE53078), length(trt.GSE58082))

group <- as.factor(rep(1:length(ls), n.obs))

library(lme4)

source("/Users/Bin/Disk1/Stat/Research/Project2014/Project1/R/catchToList.R")

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

write.table(var.comp, "/Users/Bin/Dropbox/Zhuo/Research/Project2014/Result/Tissue_5lab_Rsubread.txt")

