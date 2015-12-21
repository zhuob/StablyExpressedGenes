library(stablegene)
library(NBPSeq)
setwd("/Users/Bin/Dropbox/Zhuo/Research/RsubreadData")



GSE37159 <- read.table("GSE37159.Rsubread.txt", header = T) 
GSE37159 <- create.geneColumn(GSE37159)
trt.GSE37159 <- c(rep(1:4, each=2))

GSE38879 <- read.table("GSE38879.Rsubread.txt", header = T)
GSE38879 <- GSE38879[, seq(1, 23, by =2)] # 2 ILLUMINA RUNS
GSE38879 <- create.geneColumn(GSE38879)
trt.GSE38879 <- c(rep(11:14, each =3))

GSE43865 <- read.table("GSE43865.Rsubread.txt", header=T) 
GSE43865 <- create.geneColumn(GSE43865)
trt.GSE43865 <- c(rep(21:22, each=3))

GSE48767 <- read.table("GSE48767.Rsubread.txt", header = T)
GSE48767 <- create.geneColumn(GSE48767)
trt.GSE48767 <- rep(31:34, each=3)


GSE51119 <- read.table("GSE51119.Rsubread.txt", header = T) 
GSE51119 <- create.geneColumn(GSE51119)
trt.GSE51119 <- rep(41:45, each=2)

GSE51772 <- read.table("GSE51772.Rsubread.txt", header = T) 
GSE51772 <- create.geneColumn(GSE51772)
trt.GSE51772 <- rep(51:54, each=2)

GSE53078 <- read.table("GSE53078.Rsubread.txt", header = T)
GSE53078 <- create.geneColumn(GSE53078)
trt.GSE53078 <- rep(61:62, each=2)


GSE58082 <- read.table("GSE58082.Rsubread.txt")
GSE58082 <- create.geneColumn(GSE58082)
trt.GSE58082 <- rep(71:72, each =3)


GSE57806 <- read.table("GSE57806.Rsubread.txt")
GSE57806 <- create.geneColumn(GSE57806)
trt.GSE57806 <- rep(81:82, each =3)

ls <- list(GSE37159, GSE38879, GSE43865, GSE48767,
           GSE51119, GSE51772, GSE53078, GSE58082, GSE57806)

stable <- join_all(ls, "Gene")
stable <- stable[complete.cases(stable),]

# compare for seedling data

## Normalization first
stable2 <- stable[, -1]
norm.factors2 <- estimate.norm.factors(stable2)


# filter first
getdat <- ArrangeData("seedling")
norm.factors <- estimate.norm.factors(getdat$data)

plot(norm.factors, norm.factors2, pch=3, cex=0.5)
abline(0, 1)






