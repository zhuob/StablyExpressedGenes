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


#  Table 3.1: basic statistics of the data
print( c("Group", "# genes", "# samples", "# experiment","# treatments"))
print(c("Multiple-Tissue",dim(var_tissue$count), length(unique(var_tissue$lab)), length(unique(paste(var_tissue$lab, var_tissue$trt, sep="_")))))
print(c("Leaf", dim(var_leaf$count), length(unique(var_leaf$lab)),  length(unique(paste(var_leaf$lab, var_leaf$trt, sep="_"))) ))
print(c("Seedling", dim(var_seedling$count), length(unique(var_seedling$lab)),length(unique(paste(var_seedling$lab, var_seedling$trt, sep="_")))))


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





