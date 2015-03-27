setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Data/arab/Rsubread")
source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Code/create.geneColumn.R")

GSE32216 <- read.table("GSE32216.Rsubread.txt")
GSE32216 <- create.geneColumn(GSE32216)

GSE30795 <- read.table("GSE30795.Rsubread.txt")
GSE30795 <- create.geneColumn(GSE30795)

GSE37159 <- read.table("GSE37159.Rsubread.txt", header = T) 
GSE37159 <- create.geneColumn(GSE37159)

GSE38879 <- read.table("GSE38879.Rsubread.txt", header = T)
GSE38879 <- GSE38879[, seq(1, 23, by =2)]
GSE38879 <- create.geneColumn(GSE38879)

GSE39214 <- read.table("GSE39214.Rsubread.txt", header=T) # seedling
GSE39214 <- create.geneColumn(GSE39214)

GSE40256 <- read.table("GSE40256.Rsubread.txt")
GSE40256 <- create.geneColumn(GSE40256)

GSE42968 <- read.table("GSE42968.Rsubread.txt", header = T) 
GSE42968 <- create.geneColumn(GSE42968)

GSE43865 <- read.table("GSE43865.Rsubread.txt", header=T) # seedling
GSE43865 <- create.geneColumn(GSE43865)

# I proceeded this data myself
GSE48767 <- read.table("GSE48767.Rsubread.txt", header = T)
GSE48767 <- GSE48767[, -(1:6)]
GSE48767 <- create.geneColumn(GSE48767)

GSE51119 <- read.table("GSE51119.Rsubread.txt", header = T) 
GSE51119 <- create.geneColumn(GSE51119)
# head(GSE51772)
# colSums(GSE51772[,-1]) # reasonable

GSE51772 <- read.table("GSE51772.Rsubread.txt", header = T) 
GSE51772 <- create.geneColumn(GSE51772)

GSE53078 <- read.table("GSE53078.Rsubread.txt", header = T) # included
GSE53078 <- create.geneColumn(GSE53078)

GSE58082 <- read.table("GSE58082.Rsubread.txt")
GSE58082 <- create.geneColumn(GSE58082)

trt.GSE32216 <- c(rep(1, 2), rep(2, 2))

trt.GSE37159 <- c(rep(1:4, each=2))
trt.GSE38879 <- c(rep(1:4, each =3))
trt.GSE39214 <- c(rep(1:4, each=3))
trt.GSE42957 <- c(rep(1:6, each=2)) # removed two columns
trt.GSE42968 <- c(rep(1:2, each=3))
trt.GSE43865 <- c(rep(1:2, each=3))
trt.GSE48767 <- c(rep(1:2, each=3))
trt.GSE51119 <- rep(1:5, each=2)
trt.GSE51772 <- rep(1:4, each=2)
trt.GSE53078 <- rep(1:2, each=2)

trt.GSE56922_10 <- rep(1:6, each=2)
trt.GSE56922_20 <- rep(7:10, each=2)

trt.GSE58082 <- rep(1:2, each =3)
