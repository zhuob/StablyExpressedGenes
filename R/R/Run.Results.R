## run the figures, tables and results shown in paper

## need the following directory
#   CodePath: where the code is located
# ResultPath: where the results from iteration.GLMM()  and rankVvalue() are stored
#  SuppPath: where should the table be stored
# FigurePath: where should the figures be stored

## textsize: 4 elements, for legend, plot.title, axis, and axis title respectively.
library(xtable)

CodeNeeded1 <- paste(CodePath, "Compare.gene.R", sep = "")
CodeNeeded2 <- paste(CodePath, "plot.genes.R", sep="")
source(CodeNeeded1)
source(CodeNeeded2)


######## prepare the data
var_seedling <- readRDS(paste(ResultPath, "seedling.columbia.use.reference.iter.1.rds", sep=""))
#var_seedling <- readRDS(paste(ResultPath, "seedling.58.iteration.1.rds", sep=""))
var_leaf <- readRDS(paste(ResultPath,"leaf.columbia.use.reference.iter.1.rds", sep=""))
var_tissue <- readRDS(paste(ResultPath,"tissue.columbia.use.reference.iter.1.rds", sep=""))
cze_100 <- read.table(paste(DataPath, "Czechowski100.txt", sep=""), header=T)
dek_50 <- read.table(paste(DataPath, "Dekkers50.txt", sep=""), header=T)
geNorm <- readRDS(paste(ResultPath,"geNorm_rank_tissue.rds", sep=""))

#########  Methods SECTION  ---------------------------------------

#  Table 2.1: basic statistics of the data
group <- c("seedling", "leaf", "multi-tissue")
n_experiment <- c(length(unique(var_seedling$lab)), length(unique(var_leaf$lab)),length(unique(var_tissue$lab)))
n_treatment <- c(length(unique(paste(var_seedling$lab, var_seedling$trt, sep="_"))), length(unique(paste(var_leaf$lab, var_leaf$trt, sep="_"))),  length(unique(paste(var_tissue$lab, var_tissue$trt, sep="_"))) )
n_samples <- c(ncol(var_seedling$count), ncol(var_leaf$count), ncol(var_tissue$count))
n_genes <- c(nrow(var_seedling$count), nrow(var_leaf$count), nrow(var_tissue$count))
summary_table <- data.frame(group, n_experiment, n_treatment, n_samples, n_genes)
print(xtable(summary_table), include.rownames=F)


#########  SECTION 1 ---------------------------------------
cat("producing results in section 1...\n")


# Figure 1: histogram of CPM for top 1000 stably-expressed genes
ts <- c(10, 25, 25, 25) # legened, title, axis, axis.title
wd <- 6; ht <- 5
A1 <- plot.cpm(var_seedling,  1e3, "seedling", textsize =ts)
ggsave(paste(FigurePath, "cpm_seedling.eps", sep =""), A1, width= wd, height= ht)
A2 <- plot.cpm(var_leaf,  1e3, "leaf", textsize =ts)
ggsave(paste(FigurePath, "cpm_leaves.eps", sep =""), A2, width= wd, height= ht)
A3 <- plot.cpm(var_tissue,  1e3, "multi-tissue", textsize =ts)
ggsave(paste(FigurePath, "cpm_tissue.eps", sep =""), A3, width= wd, height= ht)

## produce 1000 stably-expressed genes for the three groups
VarCompData <- c("tissue", "leaf", "seedling")
for ( i in 1:length(VarCompData)){
  temp <- paste(ResultPath, VarCompData[i], ".columbia.use.reference.iter.1.rds", sep ="") 
  Varcomp <- readRDS(temp) 
  writedata <- Varcomp$var.comp[, c(1:4, 7)]
  writedata2 <- writedata[writedata$Rank <= 1000, ]
  writedata2 <- writedata2[order(writedata2$Rank), ]
  colnames(writedata2)[2:4] <- c("between_sample", "bewtween_treatment", "between_experiment")
  write.csv(writedata2, paste(SuppPath, VarCompData[i], "Top1000.csv", sep=""),row.names = F)
}


## how does stably expressed genes across different tissue types
cat("Overlapped genes among the multi-tissue, leaf and seedling data")
overlap_gene <- compare.3Set(var_seedling, var_leaf, var_tissue, top =1000)
print( nrow(overlap_gene) )# 106


write.csv(overlap_gene, paste(SuppPath, "OverlappedGeneFromThreeSetsTop1000.csv", sep=""),row.names = F)
# how many genes of the 109 are overlapped with czechowski
ocz <- intersect(overlap_gene[, 1], cze_100[, 1])  #11
cat(sort(ocz),"\n")

# s1 <- sort(c("AT5G46630", "AT4G24550", "AT1G13320", "AT5G26760", "AT1G10430", 
#        "AT4G27120", "AT3G01150", "AT3G10330", "AT4G32560", "AT2G20790"))

# calculate the probability
cat(pbinom(length(ocz), size = 100, nrow(overlap_gene)/14000, lower.tail=F),"\n")
# pbinom(10, 100, 106/14000, lower.tail =F)

# what is the ranking of AT1G13320
interest_gene <- "AT1G13320"
get_rank <- function(data, interest_gene){
  
  ranks <- data$var.comp$Rank[data$var.comp$Gene == interest_gene]
  top_percent <- ranks/nrow(data$var.comp)
  c(ranks, top_percent)
}
cat(get_rank(var_seedling, interest_gene ), "\n")
cat(get_rank(var_leaf, interest_gene ), "\n")
cat(get_rank(var_tissue, interest_gene ), "\n")



# Venn diagram

require(VennDiagram)
top1000seedling <- var_seedling$var.comp$Gene[var_seedling$var.comp$Rank <=1000]

top1000tissue <- var_tissue$var.comp$Gene[var_tissue$var.comp$Rank <=1000]
top1000leaf <- var_leaf$var.comp$Gene[var_leaf$var.comp$Rank <=1000]
venn1 <- length(intersect(top1000seedling, top1000tissue))
venn2 <- length(intersect(top1000seedling, top1000leaf))
venn3 <- length(intersect(top1000leaf, top1000tissue))

# require(gplots)
# venn(list(Leaf = top1000leaf, Seedling = top1000seedling, Multi_tissue = top1000tissue))
ls <- list(Leaf = top1000leaf, Seedling = top1000seedling, Multi_tissue = top1000tissue)

venn.plot <- venn.diagram(ls,height = 300, width = 300, resolution = 100,
                          fill = c('yellow', 'purple', 'green'), 
                          paste(FigurePath, "/venn.png", sep =""), imagetype = "png")
#########  SECTION 2 ---------------------------------------

cat("producing results in section 2...\n")


# Figure 2: plot expression profile for the 15 genes
textsize <- c(20, 1, 20, 20)  # legened, title, axis, axis.title
wd <- 13; ht <- 5
# five traditional house-keeping genes
figA <- c("AT3G18780", "AT5G12250","AT5G60390","AT4G05320","AT1G13440")  # HKG
A4 <- plot.gene(figA, var_tissue, 1, 1e4,figure.num = "HKG",textsize = textsize)
ggsave(paste(FigurePath, "hkg.eps", sep =""), A4, width= wd, height= ht)

# FIVE GENES FROM CZECHOWSKI
figB <- c("AT4G34270","AT1G13320","AT1G59830","AT4G33380","AT2G28390")   # NOVEL
A5 <- plot.gene(figB, var_tissue, 1, 1e4, figure.num = "Czechowski", textsize)
ggsave(paste(FigurePath, "czechowski.eps", sep =""), A5, width= wd, height= ht)

set.seed(117)  ## set.seed
random.gene <- sample(1:100, 5)
genelist <- var_tissue$var.comp$Gene[var_tissue$var.comp$Rank %in% random.gene]
A6 <- plot.gene(genelist,var_tissue, 1, 1e4, figure.num = "GLMM", textsize)
print(A6)
ggsave(paste(FigurePath, "glmm.eps", sep=""), A6, width = wd, height=ht)
cat(as.character(genelist), "\n")

gene_15 <- data.frame(source = rep(c("HKG", "Czechowski", "GLMM"), each=5),
                      Gene = c(figA, figB, as.character(genelist)))

table_var_comp <- var_tissue$var.comp[var_tissue$var.comp$Gene %in% gene_15$Gene, c(1:4,7)]
table_var_report <- merge(gene_15, table_var_comp, by = "Gene")
table_var_report <- table_var_report[with(table_var_report, order(table_var_report$source, table_var_report$Rank)), ]
table_var_report <- data.frame(table_var_report$source, table_var_report[, c(1, 3:6)])
colnames(table_var_report)[c(1, 3:5)] <- c("source", "between-sample", "between-treatment", "between-experiment")
print(xtable(table_var_report, digits = c(0, 0, 4, 4, 4, 4, 0)), include.rownames = F)

## table of variance components for the 15 genes
gene_q <- show_plot_gene(figA, figB, genelist, var_tissue)
colnames(gene_q)[2:4] <- c("between-sample", "between-treatment", "between-experiment")
tableVarComp <- paste(SuppPath, "VarComp15Genes.csv", sep="")
write.csv(gene_q, tableVarComp, row.names = F)



#########  SECTION 3 ---------------------------------------

cat("producing results in section 3...\n")
# PLOT The figure
xtext <- c( "Number of most stably expressed Genes (Multi-tissue)", "L2", "L3", "Recall percentage", "L1", "L4", "L5")
A7 <- TopGene(var_tissue, var_seedling, var_leaf, cze_100, dek_50, geNorm, xtext)
ggsave(paste(FigurePath, "rankVSrank_RNA2.eps", sep = ""), A7, width = 10, height = 5)


## overlap number
colnames(geNorm)[3] <- "Rank"
print(match.gene(geNorm$Gene[geNorm$Rank  <= 100], var_tissue, top = 100)$s)
print(match.gene(geNorm$Gene[geNorm$Rank  <= 100], var_tissue, top = 1000)$s)

## rank correlation
geNorm_GLMM <- merge(geNorm, var_tissue$var.comp, by="Gene")
print(cor(geNorm_GLMM$Rank.x, geNorm_GLMM$Rank.y, method="spearman"))


#### toy exmaple
sample1 <- rep(1, 7)
sample2 <- c(1, 1, 1, 2, 2, 3, 4)

toy <- data.frame(sample1, sample2)
rownames(toy) <- paste("Gene", 1:7, sep="")
rank(stabMvalue(t(toy)), ties.method = "min")
rankReferenceSet(toy, print.level = 1)

#########  SECTION 4 ---------------------------------------

cat("producing results in section 4...\n")
 
############# stacked bar plot 
ts <- c(20, 20, 14, 17) # legened, title, axis, axis.title
set.seed(102)
gene_ids1 <- var_tissue$var.comp$Gene[var_tissue$var.comp$Rank %in% sample(1:1000, 20) ] 
set.seed(110)
gene_ids2 <- var_tissue$var.comp$Gene[var_tissue$var.comp$Rank %in% sample(1:20000, 20) ] 

A8 <- plot.stackedBar(gene_ids1, var_tissue, figure.num = NULL, textsize=ts)
ggsave(paste(FigurePath, "top1000.eps", sep = ""), A8, width = 10, height = 5)

A9 <- plot.stackedBar(gene_ids2, var_tissue, figure.num = NULL, textsize=ts)
ggsave(paste(FigurePath, "all.eps", sep = ""), A9, width = 10, height = 5)

#legend <-  plot.stackedBar(gene_ids2, var_tissue, figure.num = NULL, textsize=c(20, 20, 15, 20), legend = T)
#ggsave(paste(FigurePath, "stackbar_legend.eps", sep = ""), legend, width = 10, height = 1)


 ############# density plot
wd <- 7; ht <- 5
ts <- c(20, 30, 30, 30) # legened, title, axis, axis.title
A10 <- plot.density(var_seedling, "sedling", textsize = ts, legend=F)
ggsave(paste(FigurePath, "var_dens1.eps", sep = ""), A10, width = wd, height = ht)

A11 <- plot.density(var_leaf, "leaf", textsize = ts)
ggsave(paste(FigurePath, "var_dens2.eps", sep = ""), A11, width = wd, height = ht)

A12 <- plot.density(var_tissue, "multi-tissue", textsize = ts)
ggsave(paste(FigurePath, "var_dens3.eps", sep = ""), A12, width = wd, height = ht)

legend <- plot.density(var_tissue, "legend", legend=T, textsize = ts)
ggsave(paste(FigurePath, "var_dens_legend.eps", sep = ""), legend, width = 10, height = 1)

############# variance percentage
# produce the table 
pct1 <- (round(VariancePercent(var_seedling), 3))
pct2 <- (round(VariancePercent(var_leaf), 3))
pct3 <- (round(VariancePercent(var_tissue), 3))
source <- (c("between-sample", "between-treatment", "between-experiment"))

pct_table <- data.frame(source, pct1, pct2, pct3)
print(xtable(pct_table, digits = c(0, 0, 3,3, 3)), include.rownames = F)

#########  SECTION 5 ---------------------------------------

cat("producing results in section 5...\n")

#######  Scatter plot for normalization factors 
# legened, title, axis, axis.title

text.size <- c(20, 20, 80, 20)


setEPS() 
postscript(paste(FigurePath, "norm_seedling.eps", sep=""), width = 8, height = 8)
A13 <- plot.pair.normfactor(var_seedling, textsize = text.size)
print(A13)
dev.off()

setEPS() 
postscript(paste(FigurePath, "norm_leaf.eps", sep=""), width = 8, height = 8)
A14 <- plot.pair.normfactor(var_leaf, textsize = text.size)
print(A14)
dev.off()

setEPS() 
postscript(paste(FigurePath, "norm_tissue.eps", sep=""), width = 8, height = 8)
A15 <- plot.pair.normfactor(var_tissue, textsize = text.size)
print(A15)
dev.off()

## plot norm.factors of a new data set GSE66666
newdata <- readRDS(paste(DataPath, "GSE64410.rds", sep=""))
var_new <- var_tissue
var_new$count <- newdata

setEPS() 
postscript(paste(FigurePath, "norm_GSE64410.eps", sep=""), width = 8, height = 8)
A16 <- plot.pair.normfactor(var_new, text.size)
print(A16) 
dev.off()



