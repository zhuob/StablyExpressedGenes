## run the figures, tables and results shown in paper

## need the following directory
#   CodePath: where the code is located
# ResultPath: where the results from iteration.GLMM()  and rankVvalue() are stored
#  TablePath: where should the table be stored
# FigurePath: where should the figures be stored


CodeNeeded1 <- paste(CodePath, "Compare.gene.R", sep = "")
CodeNeeded2 <- paste(CodePath, "plot.genes.R", sep="")
source(CodeNeeded1)
source(CodeNeeded2)


######## prepare the data
var_seedling <- readRDS(paste(ResultPath, "seedling.columbia.use.reference.iter.1.rds", sep=""))
var_leaf <- readRDS(paste(ResultPath,"leaf.columbia.use.reference.iter.1.rds", sep=""))
var_tissue <- readRDS(paste(ResultPath,"tissue.columbia.use.reference.iter.1.rds", sep=""))
cze_100 <- read.table(paste(DataPath, "Czechowski100.txt", sep=""), header=T)
dek_50 <- read.table(paste(DataPath, "Dekkers50.txt", sep=""), header=T)
geNorm <- readRDS(paste(ResultPath,"geNorm_rank_tissue.rds", sep=""))

#########  Methods SECTION  ---------------------------------------

#  Table 2.1: basic statistics of the data
print( c("Group", "# genes", "# samples", "# experiment","# treatments"))
print(c("Multiple-Tissue",dim(var_tissue$count), length(unique(var_tissue$lab)), length(unique(paste(var_tissue$lab, var_tissue$trt, sep="_")))))
print(c("Leaf", dim(var_leaf$count), length(unique(var_leaf$lab)),  length(unique(paste(var_leaf$lab, var_leaf$trt, sep="_"))) ))
print(c("Seedling", dim(var_seedling$count), length(unique(var_seedling$lab)),length(unique(paste(var_seedling$lab, var_seedling$trt, sep="_")))))


#########  SECTION 1 ---------------------------------------
cat("producing results in section 1...\n")


# Figure 1: histogram of CPM for top 1000 stably-expressed genes

setEPS()
postscript(paste(FigurePath, "cpm_seedling.eps", sep=""))
plot.cpm(var_seedling,  1e3, "seedling", textsize =c(10, 30, 20, 20))
dev.off()
setEPS()
postscript(paste(FigurePath, "cpm_leaves.eps", sep=""))
plot.cpm(var_leaf,  1e3, "leaves", textsize =c(10, 30, 20, 20))
dev.off()
setEPS()
postscript(paste(FigurePath, "cpm_tissue.eps", sep=""))
plot.cpm(var_tissue,  1e3, "multiple tissues", textsize =c(10, 30, 20, 20))
dev.off()

## produce 1000 stably-expressed genes for the three groups
VarCompData <- c("tissue", "leaf", "seedling")
for ( i in 1:length(VarCompData)){
  temp <- paste(ResultPath, VarCompData[i], ".columbia.use.reference.iter.1.rds", sep ="") 
  Varcomp <- readRDS(temp) 
  writedata <- Varcomp$var.comp[, c(1, 7)]
  writedata2 <- writedata[writedata$Rank <= 1000, ]
  writedata2 <- writedata2[order(writedata2$Rank), ]
  write.csv(writedata2, paste(TablePath, VarCompData[i], "Top1000.csv", sep=""),row.names = F)
}
paste(TablePath, VarCompData[i])



#########  SECTION 2 ---------------------------------------

cat("producing results in section 2...\n")


# Figure 2: plot expression profile for the 15 genes
# five traditional house-keeping genes
figA <- c("AT3G18780", "AT5G12250","AT5G60390","AT4G05320","AT1G13440")  # HKG
# FIVE GENES FROM CZECHOWSKI
figB <- c("AT4G34270","AT1G13320","AT1G59830","AT4G33380","AT2G28390")   # NOVEL

textsize <- c(15, 1, 12, 15)  # legened, title, axis, axis.title

setEPS() 
postscript(paste(FigurePath, "A1.eps", sep=""), width = 13, height = 5)
plot.gene(figA, var_tissue, 1, 1e4,figure.num = NULL,textsize = textsize)
dev.off()

setEPS() 
postscript(paste(FigurePath, "A2.eps", sep=""), width = 13, height = 5)
plot.gene(figB, var_tissue, 1, 1e4, figure.num = NULL, textsize)
dev.off()


setEPS() 
postscript(paste(FigurePath, "A3.eps", sep=""), width = 13, height = 5)
set.seed(102)  ## 102 or 107
random.gene <- sample(1:100, 5)
genelist <- var_tissue$var.comp$Gene[var_tissue$var.comp$Rank %in% random.gene]
plot.gene(genelist,var_tissue, 1, 1e4, figure.num = NULL, textsize)
dev.off()

## table of variance components for the 15 genes
gene_q <- show_plot_gene(figA, figB, genelist, var_tissue)
colnames(gene_q)[2:4] <- c("between-sample", "between-treatment", "between-experiment")
tableVarComp <- paste(TablePath, "VarComp15Genes.csv", sep="")
write.csv(gene_q, tableVarComp, row.names = F)



#########  SECTION 3 ---------------------------------------

cat("producing results in section 3...\n")
# PLOT The figure
setEPS() 
postscript(paste(FigurePath, "rankVSrank_RNA2.eps", sep=""), width = 10, height = 5)
xtext <- c( "number of most stably expressed Genes (Multi-tissue)", "L3", "L2", "recall percentage", "L1", "L4", "L5")
TopGene(var_tissue, var_leaf, var_seedling, cze_100, dek_50, geNorm, xtext)
dev.off()

## overlap number
colnames(geNorm)[3] <- "Rank"
print(match.gene(geNorm$Gene[geNorm$Rank  <= 100], var_tissue, top = 100)$s)
print(match.gene(geNorm$Gene[geNorm$Rank  <= 100], var_tissue, top = 1000)$s)

## rank correlation
geNorm_GLMM <- merge(geNorm, var_tissue$var.comp, by="Gene")
print(cor(geNorm_GLMM$Rank.x, geNorm_GLMM$Rank.y, method="spearman"))

#########  SECTION 4 ---------------------------------------

cat("producing results in section 4...\n")
 
############# stacked bar plot 
set.seed(102)
gene_ids1 <- var_tissue$var.comp$Gene[var_tissue$var.comp$Rank %in% sample(1:1000, 20) ] 
set.seed(110)
gene_ids2 <- var_tissue$var.comp$Gene[var_tissue$var.comp$Rank %in% sample(1:20000, 20) ] 

setEPS() 
postscript(paste(FigurePath, "top1000.eps", sep=""), width = 10, height = 5)
plot.stackedBar(gene_ids1, var_tissue, percent=F, figure.num = NULL, textsize=c(20, 20, 15, 20))
dev.off()

setEPS() 
postscript(paste(FigurePath, "all.eps", sep=""), width = 10, height = 5)
plot.stackedBar(gene_ids2, var_tissue, percent=F, figure.num = NULL, textsize=c(20, 20, 15, 20))
dev.off()

 ############# density plot
setEPS() 
postscript(paste(FigurePath, "var_dens1.eps", sep=""), width = 8, height = 5)
plot.density(var_seedling, "Seedling")
dev.off()

setEPS() 
postscript(paste(FigurePath, "var_dens2.eps", sep=""), width = 8, height = 5)
plot.density(var_leaf, "Leaf")
dev.off()

setEPS() 
postscript(paste(FigurePath, "var_dens3.eps", sep=""), width = 8, height = 5)
plot.density(var_tissue, "Tissue")
dev.off()

############# variance percentage
# produce the table 
print(round(VariancePercent(var_seedling), 3))
print(round(VariancePercent(var_leaf), 3))
print(round(VariancePercent(var_tissue), 3))




#########  SECTION 5 ---------------------------------------

cat("producing results in section 5...\n")

#######  Scatter plot for normalization factors 
# legened, title, axis, axis.title

text.size <- c(20, 20, 10, 20)

setEPS() 
postscript(paste(FigurePath, "norm1.eps", sep=""), width = 8, height = 8)
plot.pair.normfactor(var_seedling, textsize = text.size)
dev.off()

setEPS() 
postscript(paste(FigurePath, "norm2.eps", sep=""), width = 8, height = 8)
plot.pair.normfactor(var_leaf, textsize = text.size)
dev.off()

setEPS() 
postscript(paste(FigurePath, "norm3.eps", sep=""), width = 8, height = 8)
plot.pair.normfactor(var_tissue, textsize = text.size)

dev.off()
## plot norm.factors of a new data set GSE66666

GSE66666 <- read.table(paste(DataPath, "GSE66666.Rsubread.txt", sep=""))
var_new <- var_seedling
var_new$count <- GSE66666

setEPS() 
postscript(paste(FigurePath, "norm4.eps", sep=""), width = 8, height = 8)
plot.pair.normfactor(var_new, text.size)
 dev.off()



