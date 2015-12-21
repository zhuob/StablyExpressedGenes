# library(ggplot2)
# library(reshape2)
# library(Rsubread)
# library(lme4)
# library(NBPSeq)
# library(plyr)

library(devtools)

load_all("/Users/Bin/Dropbox/Zhuo/Research/Project2014/R/") 

setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/R/data")


Czechowski <- read.table("Czechowski100.txt")
ref.100 <- data.frame(Czechowski[, 1])
colnames(ref.100) <- "Gene"

Dekkers <- read.table("Dekkers50.txt")

###  stable genes of reference using RNA-Seq data

figA <- c("At3g18780", "At5g12250","At5g60390","At4g05320","At1g13440")
figA <- toupper(figA) # traditional
figB <- c("At4g34270","At1g13320","At1g59830","At4g33380","At2g28390")
figB <- toupper(figB) # novel



xdata <- ArrangeData(type="seedling")
names(xdata)

newdata <- xdata$data[1:100, ]
saved <- estimate.var.component(data=newdata, group=xdata$group, xdata$trt, filter.factor=3)
 
mm <- match.gene(genelist=ref.100, var.comp=saved$var.comp, top=1000)

figa <- data.frame(Gene=figA)
match.gene(figa, saved$var.comp, top=1000)

norm.100 <- norm.factor(var.comp=saved$var.comp, data=saved$count, n=100)

library(stablegene)
genelist <- c("A", "B","C")
data <- rpois(90, 500)
data1 <- matrix(data, nrow=3, ncol=30)
rownames(data1) <- genelist
lower <- 0
upper <- 100

plot.gene(genelist, count=data1, figure.num="AB")


