library(stablegene)
library(ggplot2)
library(reshape2)
library(plyr)
library(NBPSeq)
library(sm)
library(scales)

## Creating all the figures and tables

# 1  (Table)

setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/data")

var.t <- read.table("var.comp.tissue.new.txt", head=T)
var.s <- read.table("var.comp.seedling.new.txt", head=T)
var.l <- read.table("var.comp.leaf.new.txt", head=T)
seedling <- read.table("SeedlingData.txt",header=T)
tissue <- read.table("TissueData.txt",header=T)
leaf <- read.table("LeafData.txt", header=T)

tissue1 <- matrix(tissue)
mean.exp <- apply(tissue1, 1, mean)

# how much variation each component explains the total variation
# on an average basis
calculate.percent <- function(data)
{
  id <- which(data$sum > 0)
  data2 <- data[id, ]
  p.lab <- data2$lab/data2$sum
  p.trt <- data2$trt/data2$sum
  p.ind <- data2$individual/data2$sum
  c(mean(p.lab), mean(p.trt), mean(p.ind))
}


pct.s <- percent(round(calculate.percent(var.s), 3))
pct.t <- percent(round(calculate.percent(var.t), 3))
pct.l <- percent(round(calculate.percent(var.l), 3))







k <- 1000
genelist <- var.s[, 1][var.s$Rank<=k]
genelist <- data.frame(Gene=genelist)

s1 <- match.gene(genelist, var.comp=var.t, top=k)
s2 <- match.gene(genelist, var.comp=var.l, top=k)
c(s1$s, s2$s)

plot.rank(genelist,var.comp=var.s,col="green",figure.num= "B")

dat.t <- var.t[, c(1, 7)]
dat.s <- var.s[, c(1, 7)]
dat.l <- var.l[, c(1, 7)]

ls <- list(dat.t, dat.s, dat.l)
rank.all <- join_all(ls, "Gene")
rank.all <- rank.all[complete.cases(rank.all),]

#  RANK TO RANK PLOT in SUPLEMENTARY 
plot(x=rank.all[,4], y = rank.all[, 3], 
     xlab="tissue", ylab="seedling", pch=20)
smart.plot(x=rank.all[,3], y = rank.all[, 4],
           xlab="seedling",ylab="leaf",resolution=32)

smart.plot(x=rank.all[,3], y = rank.all[, 2],
           xlab="seedling", ylab="tissue",resolution=32)

##############################################################################
##      VARIANCE COMPONENT FIGURE (CAUTION: TAKES TIME)                    ###
##############################################################################

data1 <- melt(var.s, measure.vars=c(2:4), id.vars="Gene", variable.name = "sourceOfVar",
              value.name= "var")
p <- ggplot(data1, aes(Gene, log(var)))
p + geom_point(aes(colour = factor(sourceOfVar)))+
  labs(title="A")

data2 <- melt(var.t, measure.vars=c(2:4), id.vars="Gene", variable.name = "sourceOfVar",
              value.name= "var")
p <- ggplot(data2, aes(Gene, log(var)))
p + geom_point(aes(colour = factor(sourceOfVar)))+
  labs(title="B")

data3 <- melt(var.l, measure.vars=c(2:4), id.vars="Gene", variable.name = "sourceOfVar",
              value.name= "var")
p <- ggplot(data3, aes(Gene, log(var)))
p + geom_point(aes(colour = factor(sourceOfVar)))+
  labs(title="C")


#########################################################
##  DENSITY PLOT FOR VARIANCE COMPONENT               ###
#########################################################


density.plot <- function(var.s, title= title)
{
  mdata <- melt(var.s[, c(1:4)], id=c("Gene"))
  id <- which(mdata$value > 1e-5)
  mdata <- mdata[id, ]
  # create value labels 
  var.source <- factor(mdata$variable, 
                       labels = c("individual", "treatment", "lab")) 
  
  # plot densities 
  sm.density.compare(log(mdata$value), var.source, xlab="Variance (log scale)")
  title(main=title)
  
  # add legend via mouse click
  colfill<-c(2:(2+length(levels(var.source)))) 
  # legend(locator(1), levels(var.source), fill=colfill)
  legend("topright", levels(var.source), fill=colfill)  
}

density.plot(var.s, title= "seedling")
density.plot(var.t, title= "tissue")
density.plot(var.l, title= "leaf")




## rank plot to see if they are consistent
var.s1 <- var.s[, c(1,7)]
var.l1 <- var.l[, c(1,7)]
var.t1 <- var.t[, c(1,7)]

ls <- list(var.s1, var.l1, var.t1)
list <- join_all(ls, "Gene")
colnames(list) <- c("Gene", "Rank.s", "Rank.l", "Rank.t")
smart.plot(x=list$Rank.s, y = list$Rank.l,pch=20, clip=32,
           xlab="seedling", ylab="leaf")
smart.plot(x=list$Rank.s, y = list$Rank.t,pch=20, clip=32,
           xlab="seedling", ylab="tissue")


##############################################################################
##                 EXPRESSION LEVEL                                        ###
##############################################################################

# FIVE GENES FROM CZECHOWSKI
figA <- c("At3g18780", "At5g12250","At5g60390","At4g05320","At1g13440")
figA <- toupper(figA) # traditional
figB <- c("At4g34270","At1g13320","At1g59830","At4g33380","At2g28390")
figB <- toupper(figB) # novel 

### plot the stably expressed genes  # SEEDLING
norm.factors <- estimate.norm.factors(seedling)
nb.data <- prepare.nb.data(seedling, norm.factors=norm.factors)
off.set <- as.numeric(nb.data$eff.lib.sizes)


plot.gene(figA, stable, 1, 1e4, "A")
plot.gene(figB, stable, 1, 1e4, "B")
genelist <- var.s$Gene[var.s$Rank %in% c(10:14)]
plot.gene(genelist, seedling, 1, 1e4, "C")


## TISSUE
norm.factors <- estimate.norm.factors(tissue)
nb.data <- prepare.nb.data(tissue, norm.factors=norm.factors)
off.set <- as.numeric(nb.data$eff.lib.sizes)


plot.gene(figA, tissue, 1, 1e4, "D")
plot.gene(figB, tissue, 1, 1e4, "E")
genelist <- var.t$Gene[var.t$Rank %in% c(10:14)]
length(genelist)
plot.gene(genelist, tissue, 1, 1e4, "F")


## LEAF
norm.factors <- estimate.norm.factors(leaf)
nb.data <- prepare.nb.data(leaf, norm.factors=norm.factors)
off.set <- as.numeric(nb.data$eff.lib.sizes)

plot.gene(figA, tissue, 1, 1e4, "G")
plot.gene(figB, tissue, 1, 1e4, "H")
genelist <- var.l$Gene[var.l$Rank %in% c(1:5)]
length(genelist)
plot.gene(genelist, leaf, 1, 1e4, "I")


##############################################################################
##                 PLOT RANK OF CZE.                                        ###
##############################################################################

## SEEDLING
ref.100 <- read.csv("Czechowski100.csv", header=T)
ref.50 <- read.table("Dekkers50.txt", header=T)
ref.100.g <- data.frame(Gene= ref.100[, 1])


seedling.rank.100 <- merge(var.s, ref.100.g, by = "Gene")
seedling.rank.100 <- data.frame(seedling.rank.100)

#### HOW MANY ARE RANKED ABOVE 5000
find.rank <- function(data, n)
{
  c(percent(sum(data$Rank<n)/dim(data)[1]), 
    sum(data$Rank<n)) 
}
find.rank(seedling.rank.100, 100)
find.rank(seedling.rank.100, 500)
find.rank(seedling.rank.100, 1000)
find.rank(seedling.rank.100, 5000)



ggplot(seedling.rank.100, aes(Gene, Rank)) +
  geom_point(color="blue")+ 
  scale_y_continuous(limits=c(0, 25000)) +
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
  labs(title="A")

seedling.rank.50 <- merge(var.t, ref.50, by= "Gene")
seedling.rank.50 <- data.frame(seedling.rank.50)

find.rank(seedling.rank.50, 100)
find.rank(seedling.rank.50, 500)
find.rank(seedling.rank.50, 1000)
find.rank(seedling.rank.50, 5000)


ggplot(seedling.rank.50, aes(Gene, Rank)) +
  geom_point(color="red")+
  scale_y_continuous(limits=c(0, 25000)) +
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
  labs(title="D")



### TISSUE
tissue.rank.100 <- merge(var.t, ref.100.g, by = "Gene")
tissue.rank.100 <- data.frame(tissue.rank.100)
#### HOW MANY ARE RANKED ABOVE 5000
find.rank(tissue.rank.100, 100)
find.rank(tissue.rank.100, 500)
find.rank(tissue.rank.100, 1000)
find.rank(tissue.rank.100, 5000)


ggplot(tissue.rank.100, aes(Gene, Rank)) +
  geom_point(color="blue")+ 
  scale_y_continuous(limits=c(0, 25000)) +
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
  labs(title="B")

tissue.rank.50 <- merge(var.t, ref.50, by= "Gene")
tissue.rank.50 <- data.frame(tissue.rank.50)
find.rank(tissue.rank.50, 100)
find.rank(tissue.rank.50, 500)
find.rank(tissue.rank.50, 1000)
find.rank(tissue.rank.50, 5000)



ggplot(tissue.rank.50, aes(Gene, Rank)) +
  geom_point(color="red")+
  scale_y_continuous(limits=c(0, 25000)) +
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
  labs(title="E")


## LEAF
leaf.rank.100 <- merge(var.l, ref.100.g, by = "Gene")
leaf.rank.100 <- data.frame(leaf.rank.100)
find.rank(leaf.rank.100, 100)
find.rank(leaf.rank.100, 500)
find.rank(leaf.rank.100, 1000)
find.rank(leaf.rank.100, 5000)





#### HOW MANY ARE RANKED ABOVE 5000

ggplot(leaf.rank.100, aes(Gene, Rank)) +
  geom_point(color="blue")+ 
  scale_y_continuous(limits=c(0, 25000)) +
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
  labs(title="C")

leaf.rank.50 <- merge(var.l, ref.50, by= "Gene")
leaf.rank.50 <- data.frame(leaf.rank.50)
find.rank(leaf.rank.50, 100)
find.rank(leaf.rank.50, 500)
find.rank(leaf.rank.50, 1000)
find.rank(leaf.rank.50, 5000)

ggplot(leaf.rank.50, aes(Gene, Rank)) +
  geom_point(color="red")+
  scale_y_continuous(limits=c(0, 25000)) +
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
  labs(title="F")


####  SIMULATE P VALUES for HYPERGEOMETRIC DISTRIBUTION

sim.p <- function(data, list, rank, nsim)
{
  list <- list[, 1]
  ref.gene <- data$Gene[data$Rank <= rank]
  obs.count <- sum(list %in% ref.gene)
  
  sim.count <- rep(0, nsim)
  for (i in 1:nsim)
  {
   samp <- sample(1:dim(data)[1], rank) 
   genes <- data$Gene[samp]
   sim.count[i] <- sum(list %in% genes)
  }
  sum(sim.count > obs.count)/nsim
}


data <- var.s
list <- ref.100.g
rank <- 1000
nsim <- 1000
sim.p(var.s, ref.100.g, rank=1000, nsim=10000)
sim.p(var.t, ref.100.g, rank=1000, nsim=10000)
sim.p(var.l, ref.100.g, rank=1000, nsim=10000)


N <- dim(var.s)[1]
N <- dim(seedling.rank.100)[1]
K <- dim(var.s)[1] - N
n <- 1000
nn <- 10000
hypergeo <- rhyper(nn, N, K, n)
sum(hypergeo>= 29)


##############################################################################
##              HOW CONSISTENT ARE THE THREE DATA SETS                     ###
##############################################################################

consistency <- function(data1, data2, top1=100, top2 = 100)
{
  ref.gene1 <- data1$Gene[data1$Rank <= top1]
  ref.gene2 <- data2$Gene[data2$Rank <= top2]
  ## how many top1 are in top 2
  top <- sum(ref.gene1 %in% ref.gene2)
  c(top, top/top1)
}

consistency(var.l, var.s, 100, 100)
consistency(var.l, var.s, 100, 500)
consistency(var.l, var.s, 100, 1000)
consistency(var.l, var.s, 100, 5000)

consistency(var.t, var.s, 100, 100)
consistency(var.t, var.s, 100, 500)
consistency(var.t, var.s, 100, 1000)
consistency(var.t, var.s, 100, 5000)

consistency(var.t, var.l, 100, 100)
consistency(var.t, var.l, 100, 500)
consistency(var.t, var.l, 100, 1000)
consistency(var.t, var.l, 100, 5000)

consistency.plot <- function(data1, data2, top1=100)
{
  n <- dim(data2)[1]
  pct <- rep(0, n)
  ref.gene1 <- data1$Gene[data1$Rank <= top1]
  for ( i in 1:n ){
    ref.gene2 <- data2$Gene[data2$Rank <= i]
    ## how many top1 are in top 2
    top <- sum(ref.gene1 %in% ref.gene2)
    pct[i] <- top/top1   
  }
  pct
}

pct.t.s <- consistency.plot(var.t, var.s)
pct.l.s <- consistency.plot(var.l, var.s)

plot(1:length(pct.t.s), pct.t.s, type="l", col="blue", ylim=c(0, 1),
     xlab="top ranked stable genes", ylab="pct of overlap")
lines(1:length(pct.l.s), pct.l.s, lty=2, col="black")

pct.s.t <- consistency.plot(var.s, var.t)
pct.l.t <- consistency.plot(var.l, var.t)
lines(1:length(pct.s.t), pct.s.t, lty=3, col="magenta")
lines(1:length(pct.l.t), pct.l.t, lty=4, col="green")

pct.s.l <- consistency.plot(var.s, var.l)
pct.t.l <- consistency.plot(var.t, var.l)
lines(1:length(pct.s.l), pct.s.l, lty=5, col="purple")
lines(1:length(pct.t.l), pct.t.l, lty=6, col="red")

legend <- c("tissue VS seedling", "leaf VS seedling", "seedling VS tissue", 
            "leaf VS tissue", "seedling VS leaf", "tissue VS leaf")
color <- c("blue", "black", "magenta", "green", "purple", "red")

legend("bottomright", legend=legend,
       lty=1:6, col=color)



##############################################################################
##                 NORMALIZATION PLOT                                      ###
##############################################################################


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

x <- data.frame(top10=as.numeric(norm.10), top1e2 = as.numeric(norm.1e2),
                top1e3=as.numeric(norm.1e3), top1e4 = as.numeric(norm.1e4))

pairs(x, pch=20, col="red")


norm.10 <- norm.factor(var.t, tissue, 10)
norm.1e2 <- norm.factor(var.t, tissue, 1e2)
norm.1e3 <- norm.factor(var.t, tissue, 1e3)
norm.1e4 <- norm.factor(var.t, tissue, 1e4)


y <- data.frame(top10=as.numeric(norm.10), top1e2 = as.numeric(norm.1e2),
                top1e3=as.numeric(norm.1e3), top1e4 = as.numeric(norm.1e4))

pairs(y, col="blue", pch=20)


norm.10 <- norm.factor(var.l, leaf, 10)
norm.1e2 <- norm.factor(var.l, leaf, 1e2)
norm.1e3 <- norm.factor(var.l, leaf, 1e3)
norm.1e4 <- norm.factor(var.l, leaf, 1e4)

z <- data.frame(top10=as.numeric(norm.10), top1e2 = as.numeric(norm.1e2),
                top1e3=as.numeric(norm.1e3), top1e4 = as.numeric(norm.1e4))

pairs(z, col="magenta", pch=20)


## apply stable genes in tissue to seedling data


norm.10 <- norm.factor(var.t, seedling, 10)
norm.1e2 <- norm.factor(var.t, seedling, 1e2)
norm.1e3 <- norm.factor(var.t, seedling, 1e3)
norm.1e4 <- norm.factor(var.t, seedling, 1e4)


z.pred <- data.frame(top10=as.numeric(norm.10), top1e2 = as.numeric(norm.1e2),
                top1e3=as.numeric(norm.1e3), top1e4 = as.numeric(norm.1e4))

pairs(z.pred, col="cyan", pch=20)


