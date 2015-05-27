#' regression for all genes within a dataset.
#' 
#'
#' @title estimate variance component for all genes piecewisely.
#'
#' @param data  an n-by-p data frame that combines experimental data to be analyzed, and row names are gene IDs
#' @param group a p-vector of ids,  indicating to which experiment the column should belong 
#' @param trt  a p-vector, treatment ids nested in group ids.
#' 
#' @return a list of two data frames. One is the read count matrix, and the other 
#' is a data frame described below. 
#' an n-by-7 data frame containing the columns
#' \item{gene} gene names
#' \item{individual} variance component for biological sample
#' \item{trt} variance component for treatment
#' \item{lab} variance component for experiment
#' \item{warnS} a 1-0 indicator, 1 if there is a warning and 0 if not. 
#' \item{sum}  the sum of individual, trt and lab
#' \item{Rank} the rank of sum (column 6) in ascending order 
#' 


estimate.var.component <- function(data, group, trt, filter.factor=3)
  {
  
  groups <- as.factor(group)
  treat <- as.factor(trt)
  stable <- as.matrix(data[rowSums(data)>=filter.factor*dim(data)[2], ])
  
  library(NBPSeq)
  library(lme4)
  
  norm.factors <- estimate.norm.factors(stable)
  # norm.factors <- estimate.norm.factors(stable, method=NULL)
  nb.data <- prepare.nb.data(stable, norm.factors=norm.factors)
  off.set <- as.numeric(nb.data$eff.lib.sizes)
  
  warn <- rep(0, dim(stable)[1])
  var1 <- var2 <- var3 <- c()
  
  for ( i in 1:dim(stable)[1]){ 
    y <- as.numeric(stable[i, ])  # regression for each gene
    id <- 1:length(y)
    # print(i)
    mod1 <-  catchToList( glmer(y ~ 1  + (1|group) + (1|group:trt) + (1|id),
                                offset=(log(off.set)), 
                                control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)),
                                family = poisson)) # catch warnings or errors
    
    var3[i] <- as.numeric(mod1$value$varcor[1]) # variance for individual random effect
    var2[i] <- as.numeric(mod1$value$varcor[2]) # variance for trt effect
    var1[i] <- as.numeric(mod1$value$varcor[3]) # variance for group effect
    
    if (mod1$warnings[length(mod1$warnings)] != "NO"  ) ## there is warning
    {
      warn[i] = 1 # if the relative gradients are still large, then give a warning
      
    } 
  }
  
  sum <- var1 + var2 + var3
  rank.sum <- rank(sum)
  var.comp<- data.frame(Gene= rownames(stable),individual =var3,
                        trt=var2, lab=var1, warnS=warn, sum=sum, Rank=rank.sum)
  
  return(list(count=data, var.comp=var.comp))
}

