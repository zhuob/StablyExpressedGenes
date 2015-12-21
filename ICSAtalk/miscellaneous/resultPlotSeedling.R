###  create seedling plot

library(ggplot2)
library(reshape2)
library(plyr)
library(NBPSeq)
library(sm)
setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/data")

seedling.col <- readRDS("seedling.columbia.rds")
seedling.count <- seedling.col$count

ids <- which(rownames(seedling.count) %in% seedling.col$var.comp$Gene)
seedling<- seedling.count[ids, ]
var.seedling <- seedling.col$var.comp


#########################################################
##  DENSITY PLOT FOR VARIANCE COMPONENT               ###
#########################################################
mdata <- melt(var.seedling[, c(1:4)], id=c("Gene"))
id <- which(mdata$value > 1e-5)
mdata <- mdata[id, ]
# create value labels 
var.source <- factor(mdata$variable, 
                     labels = c("between-sample", "between-treatment", "between-experiment")) 

# plot densities 
sm.density.compare(log(mdata$value), var.source, xlab="log-Variance", col=c("red", "blue", "black"))


# add legend via mouse click
colfill<- c("red", "blue", "black")

legend("topleft", levels(var.source), fill=colfill)





##############################################################################
##                 EXPRESSION LEVEL                                        ###
##############################################################################

# FIVE GENES FROM CZECHOWSKI
figA <- c("At3g18780", "At5g12250","At5g60390","At4g05320","At1g13440")
figA <- toupper(figA) # traditional
figB <- c("At4g34270","At1g13320","At1g59830","At4g33380","At2g28390")
figB <- toupper(figB) # novel 


## derkkers seed-specific 
B.dekker <- c("AT1G16790", "AT2G43370", "AT3G59990", "AT4G12590", "AT4G02080")
A.dekker <- c("AT5G60390", "AT5G44340", "AT3G18780", "AT1G49240", "AT4G25760")

source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/R/R/plot.gene.R")


### plot the stably expressed genes  # SEEDLING
norm.factors <- estimate.norm.factors(seedling)
nb.data <- prepare.nb.data(seedling, norm.factors=norm.factors)
off.set <- as.numeric(nb.data$eff.lib.sizes)


plot.gene(figA, seedling, 1, 1e4, "Traditional HKGs")
plot.gene(figB, seedling, 1, 1e4, "Reference Gene (Czechowski et. al.)")

genelist <- var.seedling$Gene[var.seedling$Rank %in% c(10:14)]
plot.gene(genelist, seedling, 1, 1e4, "Stable Genes for Seedling")



