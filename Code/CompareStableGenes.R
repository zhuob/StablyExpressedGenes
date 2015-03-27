library(NBPSeq)
library(ggplot2)
library(reshape2)

# setwd("C:/Users/zhuob/Dropbox/Zhuo/Research/Project2014/Result")
setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Result")

ref.100 <- read.table("ref.100.gene.paper.txt", header=T)
# all <- read.table("var.comp.all.txt", header=T)
# tissue <- read.table("var.comp.4tissue.txt", header=T)

ref.50 <- read.table("ref.50.txt", header=T)
# seed <- read.table("var.comp.seed.txt", header=T)
# leaf <- read.table("var.comp.leaf.txt", header=T)

###  stable genes of reference using RNA-Seq data

figA <- c("At3g18780", "At5g12250","At5g60390","At4g05320","At1g13440")
figA <- toupper(figA) # traditional
figB <- c("At4g34270","At1g13320","At1g59830","At4g33380","At2g28390")
figB <- toupper(figB) # novel


### seedling ----------------------------------------

seedling9 <- read.table("seedling_9lab_Rsubread.txt")
colnames(seedling9)[7] <- "Rank"
seedling8_2 <- read.table("seedling_8lab_Rsubread.txt", header=T)
colnames(seedling8_2)[7] <- "Rank"
seedling10 <- read.table("seedling_10lab_Rsubread.txt")
colnames(seedling10)[7] <- "Rank"

nobs <- 100 # top 200
# seedling <- intersect(seedling8_2$gene[seedling8_2$Rank <= nobs],
#                          seedling8$gene[seedling8$Rank <= nobs])

se8.100 <- intersect(seedling8_2$gene[seedling8_2$Rank <= nobs],ref.100[, 1])
se8.50 <- intersect(seedling8_2$gene[seedling8_2$Rank <= nobs],ref.50[, 1])

seedling89 <- intersect(seedling9$gene[seedling9$Rank <= nobs],
                        seedling8_2$gene[seedling8_2$Rank <= nobs])

seedling10.100 <- intersect(seedling10$gene[seedling10$Rank <= nobs],ref.100[, 1])

seedling9.10 <- intersect(seedling9$gene[seedling9$Rank <= nobs],
                          seedling10$gene[seedling10$Rank <= nobs])

seedling8.10 <- intersect(seedling8_2$gene[seedling8_2$Rank <= nobs],
                          seedling10$gene[seedling10$Rank <= nobs])

c(length(se8.100), length(se9.100), length(seedling10.100))
c(length(seedling89), length(seedling9.10), length(seedling8.10))



### plot the stably expressed genes

stable <- read.table("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Data/arab/Rsubread/seedling9lab.txt")
norm.factors <- estimate.norm.factors(stable)
nb.data <- prepare.nb.data(stable, norm.factors=norm.factors)
off.set <- as.numeric(nb.data$eff.lib.sizes)


plot.gene(figA, stable, 1, 1e4, "A")
plot.gene(figB, stable, 1, 1e4, "B")
genelist <- seedling9$gene[seedling9$Rank <= 5]
plot.gene(genelist, stable, 1, 1e4, "C")



### see the rank of ref.100 or ref.50 genes


ref.100.g <- data.frame(gene= ref.100[, 1])
seedling.rank.100 <- merge(seedling9, ref.100.g, by = "gene")
seedling.rank.100 <- data.frame(seedling.rank.100)

ggplot(seedling.rank.100, aes(gene, Rank)) +
  geom_point(color="blue")+ 
  scale_y_continuous(limits=c(0, 25000)) +
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
  labs(title="A")


seedling.rank.50 <- merge(seedling9, ref.50, by= "gene")
seedling.rank.50 <- data.frame(seedling.rank.50)
colnames(seedling.rank.50)[7] <- "Rank"

ggplot(seedling.rank.50, aes(gene, Rank)) +
  geom_point(color="red")+
  scale_y_continuous(limits=c(0, 25000)) +
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
  labs(title="D")
#+
#  labs(title="Stable Genes (Derkkers) ranking in seedling data")
sum(seedling.rank.100$Rank< 2000)



# tissue -----------------------------------------------


tissue <- read.table("Tissue_5lab_Rsubread.txt")
colnames(tissue)[7] <- "Rank"
# the spearman correlation between rank of seedling and diff_tissues
genelist <- merge(tissue, seedling9, by="gene")

cor(genelist[, 7], genelist[, 13],method="spearman")


stable <- read.table("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Data/arab/Rsubread/tissue5lab.txt")
norm.factors <- estimate.norm.factors(stable)
nb.data <- prepare.nb.data(stable, norm.factors=norm.factors)
off.set <- as.numeric(nb.data$eff.lib.sizes)


plot.gene(figA, stable, 1, 1e4, "D")
plot.gene(figB, stable, 1, 1e4, "E")
genelist <- tissue$gene[tissue$Rank %in% c(10:14)]
length(genelist)
plot.gene(genelist, stable, 1, 1e4, "F")


nobs <- 100 # top 100
# seedling <- intersect(seedling8_2$gene[seedling8_2$Rank <= nobs],
#                          seedling8$gene[seedling8$Rank <= nobs])


ref.100.g <- data.frame(gene= ref.100[, 1])
tissue.rank.100 <- merge(tissue, ref.100.g, by = "gene")
tissue.rank.100 <- data.frame(tissue.rank.100)

ggplot(tissue.rank.100, aes(gene, Rank)) +
  geom_point(color="blue")+ 
  scale_y_continuous(limits=c(0, 25000)) +
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
  labs(title="B")

sum(tissue.rank.100$Rank< 1000)


tissue.rank.50 <- merge(tissue, ref.50, by= "gene")
tissue.rank.50 <- data.frame(tissue.rank.50)

ggplot(tissue.rank.50, aes(gene, Rank)) +
  geom_point(color="red")+
  scale_y_continuous(limits=c(0, 25000)) +
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
 labs(title="E")


## stable genes of reference using RNA-Seq data


id1 <- which(seedling8_2$gene %in% figB)
seedling8_2[id1, ]

id2 <- which(seedling9$gene %in% figB)
seedling9[id2, ]

id3 <- which(tissue$gene %in% figB)
tissue[id3,]


### leaf ----------------------------------------------------





leaf <- read.table("leaf_5lab_Rsubread.txt")
colnames(leaf)[7] <- "Rank"
# the spearman correlation between rank of seedling and diff_tissues
genelist <- merge(leaf, seedling9, by="gene")

cor(genelist[, 7], genelist[, 13],method="spearman")



stable <- read.table("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Data/arab/Rsubread/leaf5lab.txt")
norm.factors <- estimate.norm.factors(stable)
nb.data <- prepare.nb.data(stable, norm.factors=norm.factors)
off.set <- as.numeric(nb.data$eff.lib.sizes)

plot.gene(figA, stable, 1, 1e4, "G")
plot.gene(figB, stable, 1, 1e4,  "H")
genelist <- leaf$gene[leaf$Rank %in% c(1:5) ]
plot.gene(genelist, stable,1, 1e4,  "I")



nobs <- 100 # top 100
# seedling <- intersect(seedling8_2$gene[seedling8_2$Rank <= nobs],
#                          seedling8$gene[seedling8$Rank <= nobs])


ref.100.g <- data.frame(gene= ref.100[, 1])
leaf.rank.100 <- merge(leaf, ref.100.g, by = "gene")
leaf.rank.100 <- data.frame(leaf.rank.100)

ggplot(leaf.rank.100, aes(gene, Rank)) +
  geom_point(color="blue")+ 
  scale_y_continuous(limits=c(0, 25000)) +
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
  labs(title="C")

sum(leaf.rank.100$Rank< 1000)


leaf.rank.50 <- merge(leaf, ref.50, by= "gene")
leaf.rank.50 <- data.frame(leaf.rank.50)

ggplot(leaf.rank.50, aes(gene, Rank)) +
  geom_point(color="red")+
  scale_y_continuous(limits=c(0, 25000)) +
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
  labs(title="F")


id1 <- which(leaf$gene %in% figB)
leaf[id1, ]

id2 <- which(tissue$gene %in% figB)
tissue[id2, ]

id3 <- which(seedling9$gene %in% figB)
seedling9[id3,]



# 
# data1 <- melt(seedling9, measure.vars=c(2:4), id.vars="gene", variable.name = "sourceOfVar",
#               value.name= "var")
# p <- ggplot(data1, aes(gene, log(var)))
# p + geom_point(aes(colour = factor(sourceOfVar)))+
#   labs(title="A")
# 
# data2 <- melt(tissue, measure.vars=c(2:4), id.vars="gene", variable.name = "sourceOfVar",
#               value.name= "var")
# p <- ggplot(data2, aes(gene, log(var)))
# p + geom_point(aes(colour = factor(sourceOfVar)))+
#   labs(title="B")
# 
# data3 <- melt(leaf, measure.vars=c(2:4), id.vars="gene", variable.name = "sourceOfVar",
#               value.name= "var")
# p <- ggplot(data3, aes(gene, log(var)))
# p + geom_point(aes(colour = factor(sourceOfVar)))+
#   labs(title="C")
# 
# 

nobs <- 1000

se9.100 <- intersect(seedling9$gene[seedling9$Rank <= nobs],ref.100[, 1])
tissue.100<- intersect(tissue$gene[tissue$Rank <= nobs], ref.100[, 1] )
leaf.100<- intersect(leaf$gene[leaf$Rank <= nobs], ref.100[, 1] )

c(length(se9.100), length(tissue.100), length(leaf.100))

tissue.seedling <- intersect(tissue$gene[tissue$Rank<=nobs],
                             seedling9$gene[seedling9$Rank <=nobs])



leaf.seedling <- intersect(leaf$gene[leaf$Rank<=nobs],
                           seedling9$gene[seedling9$Rank <=nobs])


se9.50 <- intersect(seedling9$gene[seedling9$Rank <= nobs],ref.50[, 1])
tissue.50 <-  intersect(tissue$gene[tissue$Rank <= nobs], ref.50[, 1] )
leaf.50 <-  intersect(leaf$gene[leaf$Rank <= nobs], ref.50[, 1] )
c(length(se9.50), length(tissue.50), length(leaf.50))




