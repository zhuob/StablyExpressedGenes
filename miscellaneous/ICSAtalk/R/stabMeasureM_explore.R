
##  The stabMeasureM() function

##
library(NormqPCR)

source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/R/R/plot.genes.R")
arab <- readRDS("/Users/Bin/Dropbox/Zhuo/Research/Project2014/data/seedling.columbia.rds")
seedling <- arab$count

x <- t(seedling[1:3,])



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

stabMValue(t(data), log = F)
sort(stabMeasureM(t(data[,1:9]), log=F))


## Selecting HK genes
####  I write my own function to rank the reference set, which gives the same results
####  as the selectHKs() function in NormqPCR package.
########################################################
data(geNorm)
tissue <- as.factor(c(rep("BM", 9),  rep("FIB", 20), rep("LEU", 13),
                      rep("NB", 34), rep("POOL", 9)))
res.BM <- selectHKs(geNorm.qPCRBatch[,tissue == "BM"], method = "geNorm", 
                    Symbols = featureNames(geNorm.qPCRBatch), minNrHK = 2, log = FALSE)

data <- exprs(geNorm.qPCRBatch)
data1 <- data[, 1:9]       # this is the data used in res.BM example

rankReferenceSet(data1, print.level =1)





############################################################



library(NBPSeq)
y1 <- c(1, 1, 2, 2)
y <- rbind(y1, y1, y1)
estimate.norm.factors(y)

counts <- y
lib.sizes = colSums(counts)
r = t(t(counts)/lib.sizes)        # relative read counts
ref = exp(rowMeans(log(counts)))  # the reference gene
r.ref = ref/sum(ref)              # relative frequency
norm.factors = apply(r, 2, function(y) median((y/r.ref)[r.ref > 
                                                          0]))
norm.factors/exp(mean(log(norm.factors)))




####################  use top 1000 stably expressed genes 
#######  to get normalization and then re-run GLMM 

# would the two ranks by GLMM be different ?

library(stablegene)

estimate.var.component <- function(data, reference, group, trt, filter.factor=3)
{
  
  groups <- as.factor(group)
  treat <- as.factor(trt)
  
  ##  the reference set 
  ref <- data[row.names(data) %in% reference, ]
  
  stable <- as.matrix(data[rowSums(data)>=filter.factor*dim(data)[2], ])
  
  library(NBPSeq)
  library(lme4)
  
  # use the reference genes to do normalizations
  
  norm_factor <- estimate.norm.factors(ref, lib.sizes = colSums(stable))
  # norm2.factors <- estimate.norm.factors(stable)
  
  # norm.factors <- estimate.norm.factors(stable, method=NULL)
  nb.data <- prepare.nb.data(stable, norm.factors=norm_factor)
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
  
  return(list(count=data, var.comp=var.comp, norm.factor = norm_factor))
}


source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/R/R/plot.gene.R")
arab <- readRDS("/Users/Bin/Dropbox/Zhuo/Research/Project2014/data/seedling.columbia.rds")
seedling <- arab$count
var.s <- arab$var.comp
lab <- arab$lab
trt <- arab$trt
n_gene <- 1000

reference1 <- var.s$Gene[var.s$Rank <= n_genegene]





for (i in 1:5)
{
  data <- seedling
  var.seedling <- estimate.var.component(data, reference, group = lab, trt)
  alldata <- list(count = var.seedling$count, var.s = var.seedling$var.comp, 
                  norm.factor = var.seedling$norm.factor, trt= trt, lab = lab)
  
  reference <- alldata$var.s$Gene[alldata$var.s$Rank <= n_gene]
  
  path <-  "/Users/Bin/Dropbox/Zhuo/Research/Project2014/data/"
  data.name <- paste("seedling.columbia.use.reference.iter.", i, ".rds", sep="" )
  dest = paste(path, data.name, sep = "")
  saveRDS(alldata, dest)

}

alldat0 <- readRDS("/Users/Bin/Dropbox/Zhuo/Research/Project2014/data/seedling.columbia.rds")
alldat0 <- obj3
alldat0$var.s <- obj3$var.comp
alldat1 <- readRDS("/Users/Bin/Dropbox/Zhuo/Research/Project2014/data/seedling.columbia.use.reference.iter.1.rds")
alldat2 <- readRDS("/Users/Bin/Dropbox/Zhuo/Research/Project2014/data/seedling.columbia.use.reference.iter.2.rds")
alldat3 <- readRDS("/Users/Bin/Dropbox/Zhuo/Research/Project2014/data/seedling.columbia.use.reference.iter.3.rds")
alldat4 <- readRDS("/Users/Bin/Dropbox/Zhuo/Research/Project2014/data/seedling.columbia.use.reference.iter.4.rds")
alldat5 <- readRDS("/Users/Bin/Dropbox/Zhuo/Research/Project2014/data/seedling.columbia.use.reference.iter.5.rds")

norm_f <- cbind(alldat1$norm.factor, alldat2$norm.factor, alldat3$norm.factor, alldat4$norm.factor, alldat5$norm.factor)

length(intersect(alldat0$var.comp$Gene[alldat0$var.comp$Rank<= 1000],alldat1$var.s$Gene[alldat1$var.s$Rank<= 1000]))
length(intersect(alldat1$var.s$Gene[alldat1$var.s$Rank<= 1000],alldat2$var.s$Gene[alldat2$var.s$Rank<= 1000]))
length(intersect(alldat2$var.s$Gene[alldat2$var.s$Rank<= 1000],alldat3$var.s$Gene[alldat3$var.s$Rank<= 1000]))
length(intersect(alldat3$var.s$Gene[alldat3$var.s$Rank<= 1000],alldat4$var.s$Gene[alldat4$var.s$Rank<= 1000]))
length(intersect(alldat4$var.s$Gene[alldat4$var.s$Rank<= 1000],alldat5$var.s$Gene[alldat5$var.s$Rank<= 1000]))


compare.iter.rank <- function(set1, set2){
  
  gene1 <- set1$var.s[set1$var.s$Rank <= 1000, c(1, 7)]
  rank2 <- set2$var.s[set2$var.s$Gene %in% gene1$Gene, c(1, 7)]  
  rank2$Rank = rank(rank2$Rank)
  comp <- merge(gene1, rank2, by = "Gene")
}

c1 <- compare.iter.rank(alldat1, alldat5)
plot(c1$Rank.x, c1$Rank.y)





######## rank correlation by different method
## rank the top 1000 genes by geNorm, with stepwise elimination
# prepareCPM <- cpm_result(var_tissue, top = 1000)        # normalization used here
## about 30 minutes to compute
# geNorm1000 <- rankReferenceSet(prepareCPM, log= F)#
# genorm0 <- data.frame(geNorm1000)
# genorm0$Gene <- rownames(genorm0)
# write.table(genorm0, "geNorm1000Stepwise.txt")
# genorm1000 <- read.table("geNorm1000StepwiseCPM.txt", header=T)
# compareTop1000 <- merge(genorm1000, var_tissue$var.comp, by= "Gene")
# cor(compareTop1000$geNorm1000, compareTop1000$Rank)  # 0.195

#top1000genes <- var_tissue$var.comp$Gene[var_tissue$var.comp$Rank <= 1000]
#dat <- var_tissue$count[rownames(var_tissue$count) %in% top1000genes, ]
#genorm1000 <- rankReferenceSet(log(dat + 1), log = F)
#genorm0 <- data.frame(genorm1000)
#genorm0$Gene <- rownames(genorm0)
#write.table(genorm0, "geNorm1000StepwiseRawCount.txt")
#genorm1000 <- read.table("geNorm1000StepwiseRawCount.txt", header=T)
#compareTop1000 <- merge(genorm1000, var_tissue$var.comp, by= "Gene")
#plot(compareTop1000$genorm1000, compareTop1000$Rank)  

#################### ABOVE RESULT: NOT VERY INTERESTING #######################


# cannot see any trend
genorm <- read.table("/Users/Bin/Dropbox/Zhuo/Research/Project2014/DATA_FINALVERSION/geNorm2000StepwiseRawCount.txt", header=T)
compareTop <- merge(genorm, var_tissue$var.comp, by= "Gene")
plot(compareTop$Rank, compareTop$genorm)
gene1 <- cze_100[, 1]
gene2 <- compareTop$Gene[compareTop$genorm <=100]

ngene <- 5

id1 <- which(rownames(var_tissue$count) %in% genorm$Gene[genorm$genorm <= ngene])
TopgeNormdata <- t(var_tissue$count[id1, ])
genorm$Gene[genorm$genorm <= ngene]
log(TopgeNormdata+1)
colMeans(log(TopgeNormdata+1))
plot.gene(genorm$Gene[genorm$genorm <= ngene],var_tissue, 1, 1e4, figure.num = NULL, textsize)
a1 <- rankReferenceSet(t(log(TopgeNormdata+1)), print.level =0)
a2 <- sort(stabMValue((log(TopgeNormdata+1))))
sum(names(a1) == names(a2)) # 

ourtopgene <- var_tissue$var.comp$Gene[var_tissue$var.comp$Rank <= ngene]
id2 <- which(rownames(var_tissue$count) %in% ourtopgene) 
ourtopgene
log(t(var_tissue$count[id2, ])+1)
colMeans(log(t(var_tissue$count[id2, ])+1))
plot.gene(ourtopgene,var_tissue, 1, 1e4, figure.num =NULL, textsize)

new_trt <- as.factor(noquote(paste(obj$lab, obj$trt, sep="_")))

calcuM <- stabMeasureRho((log(TopgeNormdata+1)),group =new_trt,  log=F)



# when iterated 2000 genes, the correlation between final ranks and the original ranks is small.
data1 <- merge(geNorm_Normfinder, genorm, by ="Gene")
cor(data1[, c(3, 6)])[1, 2]

# y <- matrix(c(10, 20, 10, 100, 199, 100, 201, 200, 200, 
#              20, 20, 21, 100, 100, 100), 5, 3, byrow=T)

y <- matrix(c(9, 20, 30, 
              10, 20, 28, 
              10, 18, 31, 
              9, 21, 30, 
              30, 20, 10), 5, 3, byrow=T)
y <- data.frame(y)
rownames(y) <- c("G1", "G2", "G3", "G4", "G5")
colnames(y) <- c("Rep1", "Rep2", "Rep3")
y
ranks <- rankReferenceSet(y, print.level =1)
ranks



one_gene <- function(mu, trt, vsd){
  n <- length(trt)
  reps <- as.vector(table(trt))
  random_f <- rnorm(length(reps), 0, sqrt(vsd))
  alpha <- rep(random_f, reps)
  new_mu <- mu + alpha
  expression_data <- rpois(n, exp(new_mu))
  return(expression_data)
}

generateRNASeq <- function(mu_vec, trt, vsd_vec){
  n <- length(trt)
  m <- length(mu_vec)
  RNAdat <- matrix(NA, m, n)
  
  for(i in 1: m)
  {
    RNAdat[i, ] <- one_gene(mu_vec[i], trt, vsd_vec[i])
   }
  colnames(RNAdat) <- paste("trt", trt, sep="")
  rownames(RNAdat) <- paste("Gene", 1:m, sep="")
  
  return(RNAdat)
}

rankGLMM <-  function(RNA_sim, trt){
  var1 <- var2 <- c()
  m <- nrow(RNA_sim)
  n <- ncol(RNA_sim)
  id <- 1:n
  for ( i in 1:m){
    
    mod1 <- catchToList( glmer(RNA_sim[i, ] ~ 1  + (1|trt) + (1|id),
                               control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)),
                               family = poisson))
    var1[i] <- as.numeric(mod1$value$varcor[1]) # variance for individual random effect
    var2[i] <- as.numeric(mod1$value$varcor[2]) # variance for trt effect
    
  }
  
  sum <- var1 + var2 
  Rank <- rank(sum)
  var.comp<- data.frame(Gene= rownames(RNA_sim),sample = var1,
                        treatment=var2, sum=sum, GLMRank=Rank)
  
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




RankSimu <- function(mu_vec, trt, vsd_vec){
  
  ngene <- length(mu_vec)
  GLM <- rankGLMM(RNA_sim, trt)
  a1 <- GLM[order(GLM[, 5]), c(1, 4, 5)]
  sort(mu_vec, decreasing = T)
  #vsd_vec
  #sort(rank(vsd_vec))
  sort(rowMeans(log(RNA_sim +1)), decreasing=T)
  a2 <- rankReferenceSet(log(RNA_sim + 1), print.level =0)
  a2 <- data.frame(geNormRank = a2)
  a3 <- data.frame(trueMu =mu_vec, meanRank = rank(mu_vec))
  d1 <- merge(a2, a3, by = "row.names")
  
  colnames(d1)[1] <- "Gene"
  comb <- merge(d1, a1, by = "Gene")
  
}


ngene <- 100
trt <- rep(1:4, each=3)
#mu_vec <- runif(ngene, 1, 8)
mu_vec <- rep(c(2, 3, 5, 8), each=ngene/4)
#vsd_vec <- rep(0.2, ngene)  ## variance for random component
vsd_vec <- seq(1.5, 0.2,  length=ngene)
names(mu_vec) <- names(vsd_vec) <- paste("Gene", 1:ngene, sep="")
RNA_sim <- generateRNASeq(mu_vec, trt, vsd_vec)



ssd <- RankSimu(mu_vec, trt, vsd_vec)
vsd <- data.frame(Gene = names(vsd_vec), sd =vsd_vec)
ssd_new <- merge(ssd, vsd, by = "Gene")
ssd_new[order(ssd_new$geNormRank),]
# sort(rowMeans(log(RNA_sim+1)), decreasing = T)
# rankReferenceSet(log(RNA_sim+1), print.level = 0)
plot(ssd_new$sd, ssd_new$sum)
plot(ssd_new$trueMu, ssd_new$GLMRank, pch=3, col="red", xlab= "true log mean", ylab= "rank")
points(ssd_new$trueMu, ssd_new$geNormRank, pch=20, col="blue")
#points(ssd_new$sd, rank(ssd_new$sd), pch=5, col="blue")
corr <- cor(data.frame(trueRank = rank(ssd_new$sd), GLMRank = ssd_new$GLMRank, geNormRank=ssd_new$geNormRank))
corr

### conclusion: there's a pattern of ranking genes by genorm algorithm
### It's ranking is highly dependent on the mean values

mod1 <- catchToList( glmer(RNA_sim[1, ] ~ 1  + (1|trt) + (1|id),
                           control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)),
                           family = poisson))


mod2 <- catchToList( glmer(RNA_sim[1, ] ~ 1  + (1|trt) + (1|id),
                           control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)),
                           family = poisson))



