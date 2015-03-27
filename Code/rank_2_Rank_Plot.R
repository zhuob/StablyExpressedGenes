library(stablegene)

setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/data")
Czech <- read.csv("Czechowski100.csv", header=T)

develop.series <- Czech[, 3:4]
colnames(develop.series) <- c("Gene", "Rank.Czech")

seedling.rank <- read.table("var.comp.seedling.new.txt", header=T)
id <- which(seedling.rank$Gene %in% develop.series[, 1])

gene.rank <- merge(develop.series, seedling.rank[id, c(1, 7)])
gene.rank2 <- gene.rank[order(gene.rank$Rank.Czech),]

plot(gene.rank$Rank, gene.rank$Rank.Czech,  pch=20)
