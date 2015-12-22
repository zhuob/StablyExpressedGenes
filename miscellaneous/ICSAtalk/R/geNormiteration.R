
###########################
# SMALL SIMULATION TO SHOW EFFECT OF ITERATION IN GENORM PROCEDURE
###########################


source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/R/R/estimate.var.component.R")


one_gene <- function(mu, lab, vsd){
  n <- length(lab)
  reps <- as.vector(table(lab))
  random_f <- rnorm(length(reps), 0, sqrt(vsd))
  alpha <- rep(random_f, reps)
  new_mu <- mu + alpha
  new_mu[new_mu<0] <- 0
  expression_data <- rpois(n, exp(new_mu))
  return(expression_data)
}

generateRNASeq <- function(mu_vec, lab, vsd_vec){
  n <- length(lab)
  m <- length(mu_vec)
  RNAdat <- matrix(NA, m, n)
  
  for(i in 1: m)
  {
    RNAdat[i, ] <- one_gene(mu_vec[i], lab, vsd_vec[i])
  }
  colnames(RNAdat) <- paste("lab", lab, sep="")
  rownames(RNAdat) <- paste("Gene", 1:m, sep="")
  
  return(RNAdat)
}



rankGLMM <-  function(RNA_sim, trt, lab){
  var1 <- var2 <- var3 <-  c()
  m <- nrow(RNA_sim)
  n <- ncol(RNA_sim)
  id <- 1:n
  for ( i in 1:m){
    
    mod1 <-  catchToList( glmer(RNA_sim[i, ] ~ 1  + (1|lab) + (1|lab:trt) + (1|id),
                            control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)),
                                family = poisson)) # catch warnings or errors
    
    var1[i] <- as.numeric(mod1$value$varcor[1]) # variance for individual random effect
    var2[i] <- as.numeric(mod1$value$varcor[2]) # variance for trt effect
    var3[i] <- as.numeric(mod1$value$varcor[3]) # variance for group effect
    
  }
  
  sum <- var1 + var2 + var3
  Rank <- rank(sum)
  var.comp<- data.frame(Gene= rownames(RNA_sim),sample = var1,
                        treatment=var2, sum=sum, Rank=Rank)
  return(var.comp) 
}

stabMValue <- function(x, log = TRUE, na.rm = TRUE) {
  if (!is.data.frame(x) & !is.matrix(x)) 
    stop("'x' has to of class matrix or data.frame")
  if (is.data.frame(x)) 
    x <- data.matrix(x)
  n <- ncol(x)
  if (n == 1) 
    stop("you need at least two variables (i.e., columns) for this computation")
  M <- numeric(n)
  for (j in 1:n) {
    if (log) 
      A <- x[, j] - x[, -j]
    else A <- log2(x[, j]/x[, -j])
    if (n > 2) {
      N <- colSums(!is.na(A))
      N[N < 1] <- NA
      Mean <- colMeans(A, na.rm = na.rm)
      M[j] <- mean(sqrt(rowSums((t(A) - Mean)^2, na.rm = na.rm)/(N - 
                                                                   1)))
    }
    else M[j] <- sd(A, na.rm = na.rm)
  }
  names(M) <- colnames(x)
  M
}

##  rank the reference set
rankReferenceSet <- function(exprData, log = F, print.level = 0){
  n_gene <-  nrow(exprData)
  
  rankGene <- rep(NA, n_gene)
  stepwiseData <- exprData
  count <- 0
  while ( nrow(stepwiseData) > 2) {
    sortMvalue <- sort(stabMValue(t(stepwiseData), log=log))
    if (print.level == 1){
      print( noquote( paste("step ", count + 1, ":", sep = "")))
      print(sortMvalue)}  
    rankGene[n_gene - count] <- names(tail(sortMvalue, 1))
    id <- which(rownames(stepwiseData) == names(tail(sortMvalue, 1))) 
    stepwiseData <- stepwiseData[-id, ]
    count <- count + 1
    cat("\r", count)    # show how many iterations have been done
  }
  rankGene[1:2] <- rownames(stepwiseData)
  finalRank <- 1:n_gene
  finalRank[2] <- 1            # the last two genes are not ranked, so both ranked 1
  names(finalRank) <- rankGene
  return(finalRank) 
}

stabMeasureRho <- function(x, group, log = TRUE, na.rm = TRUE, returnAll = FALSE){
  if(class(x) == "qPCRBatch") {
    x <- t(exprs(x))
  }
  
  if(!is.data.frame(x) & !is.matrix(x))
    stop("'x' has to be of class qPCRBatch, matrix or data.frame")
  
  if(is.data.frame(x)) x <- data.matrix(x)
  k <- ncol(x) # number of variables
  n <- nrow(x) # number of samples
  if(k == 1) 
    stop("you need at least two variables (i.e., columns) for this computation")
  
  if(!log) x <- log2(x)
  if(!is.factor(group)){
    group <- factor(group)
    warning("Argument 'group' is transformed to a factor vector.")
  }
  m <- nlevels(group) # number of groups
  Levels <- levels(group)
  ngr <- as.integer(table(group))
  
  mej <- rowMeans(x, na.rm = na.rm)
  
  if(m == 1){
    warning("There is only one group. Only the variance estimates are returned!")
    mei <- colMeans(x, na.rm = na.rm)
    me <- mean(mej, na.rm = na.rm)
    N <- colSums(!is.na(x))
    N[N < 1] <- NA
    a <- rowSums((t(x-mej)-mei+me)^2, na.rm = na.rm)/(N-1)
    b <- sum(a, na.rm = na.rm)
    var.no.group <- (a - b/(k*(k-1)))/(1-2/k)
    return(var.no.group)
  }else{
    meigr <- as.matrix(aggregate(x, by = list(group), mean, na.rm = na.rm)[,-1])
    megr <- rowMeans(meigr, na.rm = na.rm)
    x.split <- split(x, f = group)
    var.group.all <- matrix(0, ncol = k, nrow = m)
    for(i in 1:m){
      x.temp <- matrix(x.split[[i]], nrow = ngr[i])
      N <- colSums(!is.na(x.temp))
      N[N < 1] <- NA
      a <- colSums((t(t(x.temp)-meigr[i,])-mej[group == Levels[i]]+megr[i])^2, na.rm = na.rm)/(N-1)
      b <- sum(a, na.rm = na.rm)
      var.group.all[i,] <- (a - b/(k*(k-1)))/(1-2/k)
    }
    if(any(var.group.all < 0)) var.group.all <- pmax(var.group.all, 0)
    
    m1i <- colMeans(meigr, na.rm = na.rm)
    m1j <- rowMeans(meigr, na.rm = na.rm)
    m1 <- mean(m1j, na.rm = na.rm)
    dif <- t(meigr - m1j) - m1i + m1
    va <- var.group.all/ngr
    tau <- max(sum(dif*dif)/((m-1)*(k-1))-mean(va), 0)
    dnew <- dif*tau/(tau+t(va))
    vanew <- t(va+tau*va/(tau+va))
    rownames(vanew) <- rownames(dnew)
    qm <- abs(dnew)+sqrt(vanew)
    qmaal <- rowMeans(qm)
    if(returnAll)
      return(list(rho = qmaal, d = dnew, v = vanew))
    else
      return(qmaal)
  }
}




x1 <- rep(1, 7)
x2 <- c(1, 1, 1, 2, 2, 3, 4)
dat <- cbind(x1, x2)
rownames(dat) <- paste("Gene", 1:7, sep="")
colnames(dat) <- paste("sample", 1:2, sep="")
rankReferenceSet(dat, print.level = 1) 
EffectDESeq <- prepare.nb.data(dat)
dataToPresent <- data.frame(cbind(dat, EffectDESeq$rel.frequencies))
colnames(dataToPresent)[c(3, 4)] <- c("normalized sample1", "normalized sample2")
dataToPresent$GLMM <- rank(as.numeric(apply(dataToPresent[, c(3, 4)], 1, var)), ties.method="min")
dataToPresent$Vvalue <- rank(as.numeric(stabMValue(t(dataToPresent[, c(1, 2)]))), ties.method = "min")
dataToPresent$geNorm <- rank(as.numeric(rankReferenceSet(dat)), ties.method = "min")
xtable(dataToPresent, digits = c(0, 0, 0, 3, 3, 0, 0, 0))
xtable(dataToPresent[, -c(3,4, 5)], digits = c(0, 0, 0, 0, 0), caption = "Toy Example", label = "table:toyexample")



var_tissue <- readRDS("/Users/Bin/Dropbox/Zhuo/Research/Project2014/DATA_FINALVERSION/tissue.columbia.use.reference.iter.1.rds")
sampRank <- sample(1:1000, ngene)

ngene <- 100
lab <- rep(1:8, each=6)
trt <- rep(1:16, each= 3)
mu_vec <- seq(2, 8, length=ngene)
mu_vec <- rowMeans(log(var_tissue$count[var_tissue$var.comp$Rank %in% sampRank, ]+1))
 vsd_vec <- rep(0, ngene)       # case 1
#vsd_vec <- seq(1.5, 0.1, length= ngene)  # case 2
# vsd_vec <- seq(0.1, 1.5, length= ngene)  # case 3
# vsd_vec <- rep(0.01,  ngene)  # case 4
# vsd_vec <- seq(1.5, 1.5, length= ngene)  # case 4
 
# vsd_vec <- var_tissue$var.comp$sum[var_tissue$var.comp$Rank %in% sampRank]  # case 5
 
names(vsd_vec)<- names(mu_vec) <- paste("Gene", 1:ngene, sep="")
RNA_sim <- generateRNASeq(mu_vec, lab, vsd_vec)
RNA_sim[RNA_sim < 1] <- 1

h1 <- stabMeasureRho(t(log(RNA_sim+1)), as.factor(trt))
normFinder <- data.frame(Gene = names(h1), NormFinderRank = rank(h1))

library(lme4)
GLM <- rankGLMM(RNA_sim, trt, lab)
#GLM[order(GLM[, 5]), c(1, 5)]
gN1 <- stabMValue(t(log(RNA_sim + 1)))
geNormRank1 <- data.frame(genormRank=rank(gN1), Gene=names(gN1))

gN2 <- rankReferenceSet(log(RNA_sim + 1),  print.level =0)
geNormRank2 <- data.frame(genormRankIter=gN2, Gene=names(gN2))
ssdata0 <- merge(GLM, geNormRank1, by = "Gene")
ssdata <- merge(ssdata0, geNormRank2, by = "Gene")
vsd <- data.frame(trueVar =vsd_vec, trueMu = mu_vec, Gene = names(mu_vec))
s2data <- merge(ssdata, vsd, by= "Gene")
colnames(s2data)[5] <- "GLMRank"

s3data <- merge(s2data, normFinder, by ="Gene")

library(dplyr)
library(reshape2)
library(ggplot2)
copa <- melt(s3data, id.vars = c(1:4, 8,9), measure.vars = c(5:7, 10), variable.name = "RankMethod",  value.name = "Rank")

ggplot(data = copa, aes(x = trueMu, y = Rank)) + 
  geom_point(aes(shape = factor(RankMethod), color = factor(RankMethod)))

cor(s3data[, c(5:10)], method='spearman')




### top 1000 from tissue
# tData <- var_tissue$count[var_tissue$var.comp$Rank <= 1000, ]
# gp <- paste(var_tissue$trt, var_tissue$lab, sep="_")
# norm.fa<- NBPSeq::estimate.norm.factors(tData)
# ax <- stabMeasureRho(t(log(tData + 1)),gp )
# ay <- stabMValue(t(log(tData + 1)))
# quantile(var_tissue$var.comp$sum[var_tissue$var.comp$Rank <= 1000])
# RhoRank <- data.frame(Gene = names(ax), NormFinderRank = rank(ax), geNom = rank(ay))
# see <- merge(RhoRank, var_tissue$var.comp, by ="Gene")
# cor(see[, c(2, 3, 9)])

