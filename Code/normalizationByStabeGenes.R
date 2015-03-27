setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Result")
library(NBPSeq)

seedling <- read.table("seedling_9lab_Rsubread.txt")
stable <- read.table("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Data/arab/Rsubread/seedling9lab.txt")


norm.factor <- function(data1, stable, n)
  {
  genelist <- data1$gene[data1$rank.sum <= n]
  dat1 <- stable[row.names(stable) %in% genelist, ]
  estimate.norm.factors(dat1, lib.sizes= colSums(stable) )
}

norm.10 <- norm.factor(seedling, stable, 10)
norm.1e2 <- norm.factor(seedling, stable, 1e2)
norm.1e3 <- norm.factor(seedling, stable, 1e3)
norm.1e4 <- norm.factor(seedling, stable, 1e4)


x <- data.frame(x1=as.numeric(norm.10), x2 = as.numeric(norm.1e2),
                x3=as.numeric(norm.1e3), x4 = as.numeric(norm.1e4))

pairs(x, pch=20, col="red")


##### tissue data --
tissue <- read.table("Tissue_5lab_Rsubread.txt")
stable <- read.table("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Data/arab/Rsubread/tissue5lab.txt")


norm.10 <- norm.factor(tissue, stable, 10)
norm.1e2 <- norm.factor(tissue, stable, 1e2)
norm.1e3 <- norm.factor(tissue, stable, 1e3)
norm.1e4 <- norm.factor(tissue, stable, 1e4)


y <- data.frame(y1=as.numeric(norm.10), y2 = as.numeric(norm.1e2),
                y3=as.numeric(norm.1e3), y4 = as.numeric(norm.1e4))

pairs(y, col="blue", pch=20)



##### leaf data --
leaf <- read.table("leaf_5lab_Rsubread.txt")
stable <- read.table("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Data/arab/Rsubread/leaf5lab.txt")


norm.10 <- norm.factor(leaf, stable, 10)
norm.1e2 <- norm.factor(leaf, stable, 1e2)
norm.1e3 <- norm.factor(leaf, stable, 1e3)
norm.1e4 <- norm.factor(leaf, stable, 1e4)


z <- data.frame(z1=as.numeric(norm.10), z2 = as.numeric(norm.1e2),
                z3=as.numeric(norm.1e3), z4 = as.numeric(norm.1e4))

pairs(z, col="magenta", pch=20)






