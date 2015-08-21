## ### create the plots for ICSA talk ### ## 
library(NBPSeq)
library(stablegene)
library(reshape2)
library(ggplot2)
setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/data")


source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/R/R/plot.gene.R")
arab <- readRDS("/Users/Bin/Dropbox/Zhuo/Research/Project2014/data/seedling.columbia.rds")
seedling <- arab$count
var.s <- arab$var.comp

#

# FIVE GENES FROM CZECHOWSKI
figA <- c("At3g18780", "At5g12250","At5g60390","At4g05320","At1g13440")
figA <- toupper(figA) # traditional
figB <- c("At4g34270","At1g13320","At1g59830","At4g33380","At2g28390")
figB <- toupper(figB) # novel 


## derkkers seed-specific 
B.dekker <- c("AT1G16790", "AT2G43370", "AT3G59990", "AT4G12590", "AT4G02080")
A.dekker <- c("AT5G60390", "AT5G44340", "AT3G18780", "AT1G49240", "AT4G25760")

norm.factors <- estimate.norm.factors(seedling)
nb.data <- prepare.nb.data(seedling, norm.factors=norm.factors)
off.set <- as.numeric(nb.data$eff.lib.sizes)


plot.gene(figA, seedling, 1, 1e4, "Traditional HKGs")
plot.gene(figB, seedling, 1, 1e4, "Stably expressed genes identified by Czechowski et. al.")

random.gene <- sample(1:100, 5)
genelist <- var.s$Gene[var.s$Rank %in% random.gene]
plot.gene(genelist, seedling, 1, 1e4, "Stable Genes identified from RNA-Seq data")








############# Scatter plot for normalization factors

norm.factor <- function(data1, stable, n)
{
  genelist <- data1$Gene[data1$Rank <= n]
  dat1 <- stable[row.names(stable) %in% genelist, ]
  estimate.norm.factors(dat1, lib.sizes= colSums(stable) )
}

norm.10 <- norm.factor(var.s, seedling, 10)
norm.1e2 <- norm.factor(var.s, seedling, 1e2)
norm.1e3 <- norm.factor(var.s, seedling, 1e3)
norm.1e4 <- norm.factor(var.s, seedling, 1e4)

norm.all <- norm.factor(var.s, seedling, dim(var.s)[1])

norm.factor1 <- estimate.norm.factors(seedling, lib.sizes=colSums(seedling))
plot(norm.1e3, norm.factor1)
plot(norm.1e3, norm.1e2)

x <- data.frame(top10=as.numeric(norm.10), top1e2 = as.numeric(norm.1e2),
                top1e3=as.numeric(norm.1e3), top1e4 = as.numeric(norm.1e4), 
                all= norm.all)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(x, pch=20, col="red", lower.panel = panel.smooth, upper.panel = panel.cor)
pairs(x, pch=20, col="red")





#############################################################
#                 RANK TO RANK PLOT
#############################################################

Mvalue <- data.frame(readRDS("Mvalue.rds"))
Rhovalue <- data.frame(readRDS("Rhovalue.rds"))
var.l <- read.table("var.comp.leaf.new.txt", head=T)
ref.100 <- read.table("Czechowski100.txt", header=T)
ref.50 <- read.table("Dekkers50.txt", header=T)


Rhovalue$rank <- rank(Rhovalue)
colnames(Rhovalue) <- c("Rho", "rank")
gene.geNORM <- rownames(Mvalue[which(Mvalue$rank <= 100),])
gene.NormFinder <- rownames(Mvalue[which(Rhovalue$rank <= 100),])

gene.l <- var.l$Gene[var.l$Rank <= 100]  ## leaf rankings

gene.cze <- ref.100[, 1]
gene.dek <- ref.50[, 1]

rank.seq <- c(100, 200, 500, 1000, 2000, 3500, 5000)
l.num <- cze.num <- dek.num <- genorm <- Normfinder<-  c()
for ( i in 1: length(rank.seq))
{
  gene.s <- var.s$Gene[var.s$Rank <= rank.seq[i]]
 # t.num[i] <- round(length(intersect(gene.s, gene.t))/length(gene.t), 3)
  l.num[i] <- round(length(intersect(gene.s, gene.l))/length(gene.l), 3)
  cze.num[i] <- round(length(intersect(gene.s, gene.cze))/length(gene.cze), 3)
  dek.num[i] <- round(length(intersect(gene.s, gene.dek))/length(gene.dek), 3)
  genorm[i] <- round(length(intersect(gene.s, gene.geNORM))/length(gene.geNORM), 3)
  Normfinder[i] <- round(length(intersect(gene.s, gene.NormFinder))/length(gene.NormFinder), 3)
  
  }

pct.data <- data.frame(rank=rank.seq, l.num, cze.num, dek.num, genorm, Normfinder)
colnames(pct.data) <- c("Rank", "Leaf(prelimilary)", 
                        "Czechowski et. al.", "Dekkers et. al.", "geNORM", "NormFinder")
pct.data <- melt(pct.data, id="Rank")


ggplot(pct.data, aes(x= Rank, y = value, colour = variable)) +
  geom_path(alpha = 0.5) + 
  geom_line(aes(variable="Czechowski et. al.")) + 
  geom_point(aes(variable="Dekkers et. al.")) + 
 # geom_line(aes(variable="Tissue"),linetype="dotted") + 
  geom_line(aes(variable="Leaf"),linetype="dotted") +
  geom_line(aes(variable="geNORM"),linetype="dotted") +
  geom_line(aes(variable="NormFinder"),linetype="dotted") +
  labs(x="top ranks in seedling", y="percentage of presence") +
  theme(legend.position=c(0.9, 0.2))


#########################################################
##  DENSITY PLOT FOR VARIANCE COMPONENT               ###
#########################################################
library(sm)


mdata <- melt(var.s[, c(1:4)], id=c("Gene"))
id <- which(mdata$value > 1e-5)
mdata <- mdata[id, ]
# create value labels 
var.source <- factor(mdata$variable, 
                     labels = c("between-sample", "between-treatment", "between-experiment")) 

# plot densities 
sm.density.compare(log(mdata$value), var.source, xlab="log-Variance", 
                   col=c("red", "blue", "black"), lwd=c(2,2,2), lty=c(1, 1, 1))

# add legend via mouse click
colfill<-col=c("red", "blue", "black")

legend("topleft", levels(var.source),fill=c("red", "blue", "black"))

head(var.s)
var.s.percent <- data.frame(matrix(nrow=dim(var.s)[1], ncol=0))
var.s.percent$Gene <- var.s$Gene
var.s.percent$sample <- var.s$individual/var.s$sum
var.s.percent$treatment <- var.s$trt/var.s$sum
var.s.percent$experiment <- var.s$lab/var.s$sum


#### percentage Plot
mdata <- melt(var.s.percent[, c(1:4)], id=c("Gene"))

id <- which(mdata$value > 1e-6)
mdata <- mdata[id, ]
# create value labels 
var.source <- factor(mdata$variable, 
                     labels = c("between-sample", "between-treatment", "between-experiment")) 

# plot densities 
sm.density.compare(mdata$value, var.source, xlab="percentage", 
                   col=c("red", "blue", "black"), lwd=c(2,2,2), lty=c(1, 1, 1))

# add legend via mouse click
colfill<- c("red", "blue", "black")

legend("topright", levels(var.source),fill=c("red", "blue", "black"))





################## USE GGPLOT2
colnames(mdata) <- c("Gene", "Source", "percentage")
mgg <- ggplot(mdata, aes(x=percentage, colour= Source)) 
mgg + geom_density(size = 1.5) +
  theme(legend.position="top")
  


#######################################################
# STACKED BAR PLOT
#######################################################
library(plyr)
gene.ids <- sample(1:dim(var.s)[1], 20)
colnames(var.s.percent)
mdata <- melt(var.s.percent[var.s$Rank %in% gene.ids, c(1:4)], id=c("Gene"))

id <- which(mdata$value > 1e-6)
mdata <- mdata[id, ]


library(scales)
ggplot(mdata,aes(x = Gene, y = value,fill = variable)) + 
  geom_bar(position = "fill",stat = "identity") + 
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
  labs( y="percentage of total variance") +
   scale_y_continuous(labels = percent_format())



mdata <- melt(var.s[var.s$Rank %in% gene.ids, c(1:4)], id=c("Gene"))

id <- which(mdata$value > 1e-6)
mdata <- mdata[id, ]


library(scales)
ggplot(mdata,aes(x = Gene, y = value,fill = variable)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
  labs( y="log-variance") +
  scale_y_continuous()



findrow <- (1:nrow(Count))[rowSums(Count==0)==0]

gene.ids <- sample(findrow, 20)



mdata <- melt(var.s[gene.ids, c(1:4)], id=c("Gene"))

Count[gene.ids,]

id <- which(mdata$value > 1e-6)
mdata <- mdata[id, ]


library(scales)
ggplot(mdata,aes(x = Gene, y = value,fill = variable)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
  labs( y="log-variance") +
  scale_y_continuous()
Count["AT4G17970",]



######################
## 
getwd()
rho <- readRDS("Rhovalue.rds")
gene.ids <- which(rho$rank <= 1000)
range(apply(seedling[gene.ids, ], 1, mean))
cor(rho$rank, var.s$Rank)









