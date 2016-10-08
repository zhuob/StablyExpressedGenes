# library(NormqPCR)
library(ggplot2)
library(reshape2)

############ Results part ----------------------------------------------


######### Section 1 -------------------------------

match.gene <- function(genelist, set, top)
{
  var.comp <- set$var.comp
  match <- var.comp[var.comp$Rank <= top, ]
  genelist <- as.vector(genelist)
  genelist <- data.frame(genelist)
  colnames(genelist) <- "Gene"
  result <- merge(genelist, match, by="Gene")
  
  lst <- list(s = dim(result)[1], 
              l = data.frame(Gene=result$Gene, Rank=result$Rank))
  return(lst)
}


cpm_result <- function(set, top  = 1000){
## convert counts to cpm
 
 count <- set$count
 ids <- set$var.comp$Rank <= top
 
 norm.factors <- estimate.norm.factors(count)
 nb.data <- prepare.nb.data(count, norm.factors=norm.factors)
 off.set <- as.numeric(nb.data$eff.lib.sizes)
 
 y <- sweep(count[ids,], 2, off.set, "/")*1e6 # CPM

}

mean_cpm <- function(set, top = 1000){
# mean cpm for each gene over all samples 
  
   y <- cpm_result(set, top= top) 
   mean_cpm = apply(y, 1, mean)
   return(mean_cpm)
 }




compare.3Set <- function(var_seedling, var_leaf, var_tissue, top = 1000){
# compare the stably expressed genes from Seedling, Leaf, and Tissue
## seedling leaf tissue are object returned from the estimate.var.comp() function
  
            var_comp1 <- var_seedling$var.comp
            var_comp2 <- var_leaf$var.comp
            var_comp3 <- var_tissue$var.comp

            seedling_list <- var_comp1$Gene[var_comp1$Rank <= top]
            leaf_list <- var_comp2$Gene[var_comp2$Rank <= top]
            tissue_list <- var_comp3$Gene[var_comp3$Rank <= top]

            common_gene <- intersect(seedling_list, intersect(leaf_list, tissue_list))
            
            d1 <- var_comp1[var_comp1[, 1] %in% common_gene, c(1, 7)];
            d2 <- var_comp2[var_comp2[, 1] %in% common_gene, c(1, 7)]
            d3 <- var_comp3[var_comp3[, 1] %in% common_gene, c(1, 7)]
            dat <- merge(merge(d1, d2, by = "Gene"), d3, by = "Gene")
            colnames(dat) <- c("Gene", "Rank_in_seedling", "Rank_in_leaf", "Rank_in_tissue")
            return(dat)
            
    }




rank_cor <- function(set1, set2, set3, genorm_normfinder, xtext){
  ### rank correlation
  
    var_comp <- set1$var.comp
    df <- merge(var_comp, genorm_normfinder, by = "Gene")
    rank_cor <- cor(data.frame(df$Rank, df$Rank_M, df$Rank_rho), method = "spearman")
    
    data2 <- join_all(list(df, set2$var.comp, set3$var.comp), by= "Gene")
    data3 <- data2[complete.cases(data2), c(7, 9, 11, 17, 23)]
    colnames(data3) <- c(xtext[1], "geNorm", "NormFinder", xtext[2], xtext[3])
    return(round(cor(data3, method="spearman"), 3))
  
    
}


VariancePercent <- function(set){
##### Table in Section 4.3 --------------------------------------
  
    var_comp <- set$var.comp
    var2 <- data.frame(sample = var_comp$individual/var_comp$sum,
                       treatment = var_comp$trt/var_comp$sum,
                       lab = var_comp$lab/var_comp$sum)
  
    per <- apply(var2, 2, mean)
    return(per)
}


show_plot_gene <- function(figA, figB, genelist, set){
## show the variance components for the 15 genes of section Stably expressed genes
#  
    
  ids <- c(which(set$var.comp$Gene %in% c(figA, figB)), genelist)
  var_plot_gene <- set$var.comp[ids, ]
  var_plot_gene$type <- NULL
  var_plot_gene$type[var_plot_gene$Gene %in% figA] <- "House-Keeping"
  var_plot_gene$type[var_plot_gene$Gene %in% figB] <- "Czechowski"
  var_plot_gene$type[is.na(var_plot_gene$type)] <- "RNA-Seq"
  #var_plot_gene[, 2:4] <- round(var_plot_gene[, 2:4]/var_plot_gene[, 6], 4)
  var_plot_gene[, 2:4] <- round(var_plot_gene[, 2:4], 4)
  var_plot_gene <- var_plot_gene[order(var_plot_gene$Rank), c(1:4, 7, 8)]
  
  return(var_plot_gene)
  
}



stabMvalue <- function (x, log = TRUE, na.rm = TRUE) {
####### calculate the V-values of the gene, algorithm of geNorm 
## input:
#       x: the expression data, with column being genes and rows being samples
#     log: whether the data need to be log transformed, default is TRUE
#   na.rm: whether na values should be removed, default is TRUE
#
# output:
#       V: the V values for each gene    

  if (!is.data.frame(x) & !is.matrix(x)) 
    stop("'x' has to of class matrix or data.frame")
  if (is.data.frame(x)) 
    x <- data.matrix(x)
  n <- ncol(x)
  if (n == 1) 
    stop("you need at least two variables (i.e., columns) for this computation")
  M <- numeric(n)
  for (j in 1:n) {
    if (log == T) 
      A <- x[, j] - x[, -j]
    else A <- log2(x[, j]/x[, -j])
    if (n > 2) {
      N <- colSums(!is.na(A))
      N[N < 1] <- NA
      Mean <- colMeans(A, na.rm = na.rm)
      M[j] <- mean(sqrt(rowSums((t(A) - Mean)^2, na.rm = na.rm)/(N - 1)))
    }
    else M[j] <- sd(A, na.rm = na.rm)
  }
  names(M) <- colnames(x)
  M
  
}


######### Section 2 -------------------------------
## rank the genes by geNorm and NormFinder
# set is an object returned from the estimate.var.comp() function

rankMvalue <- function(set){
  obj <- set
  calcuM <- stabMvalue(t(obj$count + 1), log=F, na.rm= T)
  calcuM2 <- data.frame(Gene = names(calcuM), Mvalue = calcuM, Rank = rank(calcuM))
  rownames(calcuM2) <- NULL
  return(calcuM2)
}


rankNormNotUsed <- function(set){
  
  obj <- set
  new_trt <- as.factor(noquote(paste(obj$lab, obj$trt, sep="_")))
  calcuRho <- stabMeasureRho(x = t(obj$count +1),log=F, group = new_trt)
  calcuRho2 <- data.frame( Rho = calcuRho, rank = rank(calcuRho))
  
  calcuM <- stabMeasureM(t(obj$count + 1), log=F, na.rm= T)
  calcuM2 <- data.frame(Mvalue = calcuM, rank = rank(calcuM))
  
  geNorm_NormFinder <- merge(calcuM2, calcuRho2, by = "row.names")
  colnames(geNorm_NormFinder)[c(1, 3, 5)] <- c("Gene", "Rank_M","Rank_rho")
  
  return(geNorm_NormFinder)
  
}


rankReferenceSet <- function(exprData, log = F, print.level = 0){
  
##  rank the reference set (geNorm iterative procedure)
# input:
#         exprData: the expression data with row for genes and column for samples
#              log: whether the data need to be log transformed 
#      print.level: 0 for no print, 1 for showing every elimination step
#
# output: 
#       final.rank: the ranks for genorm with iterative elimination

  n_gene <-  nrow(exprData)
  rankGene <- rep(NA, n_gene)
  stepwiseData <- exprData
  count <- 0
  
  while ( nrow(stepwiseData) > 2) {
    sortMvalue <- sort(stabMvalue(t(stepwiseData), log=log))
    if (print.level == 1){
      print( noquote( paste("step ", count + 1, ":", sep = "")))
      print(sortMvalue)}  
    rankGene[n_gene - count] <- names(tail(sortMvalue, 1))
    id <- which(rownames(stepwiseData) == names(tail(sortMvalue, 1))) 
    stepwiseData <- stepwiseData[-id, ]
    count <- count + 1
    # cat("\r", count)    # show how many iterations have been done
  }
  
  rankGene[1:2] <- rownames(stepwiseData)
  finalRank <- 1:n_gene
  finalRank[2] <- 1            # the last two genes are not ranked, so both ranked 1
  names(finalRank) <- rankGene
  return(finalRank) 
}




