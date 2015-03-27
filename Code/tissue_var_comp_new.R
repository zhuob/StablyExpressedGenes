library(stablegene)
setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Result/")



bar.plot <- function(data, from, nobs)
{
  if( nobs + from -1 > dim(data)[1])
  {stop("gene id out of bound")}
  data1 <- melt(data, measure.vars=c(2:4), id.vars="Gene", variable.name = "SourceOfVar",
                value.name= "var")
  gap <- dim(data)[1]
  index <-  c(1:nobs, (gap+1):(gap + nobs), (2*gap +1):(2*gap + nobs))
  
  ggplot(data1[index,], aes(x=Gene, y=var, fill=SourceOfVar)) +
    geom_bar(stat="identity") + theme_bw()+
    theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) 
  
}
  
  bar.plot(leafvar, 2100, 100)

## must have header = T
tissuevar <- read.table("var.comp.tissue.new.txt",header=T)
leafvar <- read.table("var.comp.leaf.new.txt", header=T)
seedlingvar <- read.table("var.comp.seedling.new.txt", header=T)
crez <- read.table("Czechowski100.txt", header=T)
dekker <- read.table("Dekkers50.txt", header=T)

setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/data/")
x1 <- ArrangeData("tissue")
tissue <- estimate.var.component(data=x1$data, group=x1$group, trt=x1$trt,filter.factor=0)
write.table(tissue$var.comp, "var.comp.tissue.new.txt")

x2 <- ArrangeData("seedling")
seedling <- estimate.var.component(data=x2$data, group=x2$group, trt=x2$trt,filter.factor=0)
write.table(seedling$var.comp, "var.comp.seedling.new.txt")

x3 <- ArrangeData("leaf")

figA <- c("At3g18780", "At5g12250","At5g60390","At4g05320","At1g13440")
figA <- toupper(figA) # traditional
figB <- c("At4g34270","At1g13320","At1g59830","At4g33380","At2g28390")
figB <- toupper(figB) # novel


plot.gene(figA, count=x2$data, figure.num="A")
plot.gene(figB, count=x2$data, figure.num="B")
genelist <- seedlingvar$Gene[seedlingvar$Rank <= 5]
# genelist <- tissuevar$Gene[tissuevar$Rank %in% c(10:14)]
plot.gene(as.vector(genelist), count=x2$data, figure.num="C")


plot.gene(figA, count=x1$data, figure.num="D")
plot.gene(figB, count=x1$data, figure.num="E")
genelist <- tissuevar$Gene[tissuevar$Rank%in% c(10:14)]
plot.gene(as.vector(genelist), count=x1$data, figure.num="F")



plot.gene(figA, count=x3$data, figure.num="G")
plot.gene(figB, count=x3$data, figure.num="H")
genelist <- leafvar$Gene[leafvar$Rank<=5]
plot.gene(as.vector(genelist), count=x3$data, figure.num="I")





plot.rank <- function(genelist, var.comp, col="blue", figure.num)
{
  
  merged <- merge(genelist, var.comp, by = "Gene")
  merged1 <- data.frame(merged)
  
  ggplot(merged1, aes(Gene, Rank)) +
    geom_point(color=col)+
    scale_y_continuous(limits=c(0, 25000)) +
    theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
    labs(title=figure.num)
  
}

plot.rank(genelist = data.frame(Gene=crez[, 1]), var.comp=seedlingvar,figure.num="A" )
plot.rank(genelist = data.frame(Gene=crez[, 1]), var.comp=tissuevar,figure.num="B" )
plot.rank(genelist = data.frame(Gene=crez[, 1]), var.comp=leafvar,figure.num="C" )

plot.rank(genelist = dekker, col="red",var.comp=seedlingvar,figure.num="D" )
plot.rank(genelist = dekker, col="red",var.comp=tissuevar,figure.num="E" )
plot.rank(genelist = dekker, col="red",var.comp=leafvar,figure.num="F" )



library(reshape2)
library(sm)

mdata <- melt(tissuevar[, c(1:4)], id=c("Gene"))

id <- which(mdata$value > 1e-5)
mdata <- mdata[id, ]

# Compare MPG distributions for cars with 
# 4,6, or 8 cylinders


# create value labels 
var.source <- factor(mdata$variable, 
                     labels = c("individual", "trt", "lab")) 

# plot densities 
sm.density.compare(log(mdata$value), var.source, xlab="variance (log scale)")
title(main="variance Distribution by Source (Tissue)")

# add legend via mouse click
colfill<-c(2:(2+length(levels(var.source)))) 
legend("topright", levels(var.source), fill=colfill)



