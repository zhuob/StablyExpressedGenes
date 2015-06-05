## simulate genes
library(NormqPCR)

generateReadCount <- function(size, N, n.stable, mu1, mu2, mu3)
{
  readCount <- data.frame(matrix(nrow= N, ncol = size))
  
  for ( i in 1:N)
  {
    if (i <= n.stable)
    { 
      readCount[i, 1:size] <- rpois(size, mu3) # the stably expressed genes
    }
    else 
    {
      readCount[i, 1:n1] <- rpois(n1, mu1[i-n.stable]) # treatment 1
      readCount[i, (n1+1):size] <- rpois(n1, mu2[i-n.stable]) # treatment 2
    }
    
  }
  
  rownames(readCount) <-  paste("Gene", 1:N, sep="")
  return(readCount)
}


## case 1 : a lot of upregulating genes 

size <- 40 # number of samples
N <- 100 # number of genes
n1 <- size/2  # number of treatment and controls
n.stable <- 20  # number of stably expressed genes

mu1 <- rep(20, N-n.stable) + 2*rnorm(N-n.stable)
mu2 <- rep(40, N-n.stable) + 2*rnorm(N-n.stable)
mu3 <- 30


readCount <- generateReadCount(size, N, n.stable, mu1, mu2, mu3)

geneNorm <- stabMeasureM(t(readCount), log=F, na.rm=T)
rank.gene <- rank(geneNorm)
stable.rank <- data.frame(Mvalue = geneNorm, rank = rank.gene)
stable.rank <- stable.rank[order(stable.rank$rank),] # sort by rank
stable.rank

## case 2:  70% stably expressed genes



size <- 40 # number of samples
N <- 100 # number of genes
n1 <- size/2  # number of treatment and controls
n.stable <- 70  # number of stably expressed genes

mu1 <- rep(20, N-n.stable) + 2*rnorm(N-n.stable)
mu2 <- rep(40, N-n.stable) + 2*rnorm(N-n.stable)
mu3 <- 30


readCount <- generateReadCount(size, N, n.stable, mu1, mu2, mu3)

geneNorm <- stabMeasureM(t(readCount), log=F, na.rm=T)
rank.gene <- rank(geneNorm)
stable.rank <- data.frame(Mvalue = geneNorm, rank = rank.gene)
stable.rank <- stable.rank[order(stable.rank$rank),]

stable.rank

## case 3, 40 % stable genes, but others are expressed in a random way

size <- 40 # number of samples
N <- 100 # number of genes
n1 <- size/2  # number of treatment and controls
n.stable <- 20  # number of stably expressed genes

mu <- rexp(1000, 1/30) + 10  # mean is 40
mu1 <- sample(mu, N-n.stable) + 2*rnorm(N-n.stable)
mu2 <- sample(mu, N-n.stable) + 2*rnorm(N-n.stable)
mu3 <- 30

readCount <- generateReadCount(size, N, n.stable, mu1, mu2, mu3)

geneNorm <- stabMeasureM(t(readCount), log=F, na.rm=T)
rank.gene <- rank(geneNorm)
stable.rank <- data.frame(Mvalue = geneNorm, rank = rank.gene, gene=1:N)
stable.rank <- stable.rank[order(stable.rank$rank),]

stable.rank
# how many are ranked top for true stably expressed genes
sum(stable.rank$gene[1:n.stable] < n.stable)


