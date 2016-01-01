
############  whether there're a lot of DE genes in Leaf data??

library(NBPSeq)

DEtest <- function(dat, trt, grp1, grp2){
  
  file <- paste("/Users/Bin/Google Drive/Study/Thesis/NCBI/Data/", dat, ".Rsubread.txt", sep ="")
  dat1 <- read.table(file)
  dat1 <- dat1[rownames(dat1) %in% genes, ]
  
  if  (dat == "GSE39463") {
      rm.39463 <- c(6, 8, 10, 13, 15:44,  51:72)-4
      dat1 <- dat1[, -rm.39463]}
  else if (dat == "GSE43865"){
      col.id.GSE43865 <- c("GSM1072464", "GSM1072465", "GSM1072466",
                         "GSM1072467", "GSM1072468", "GSM1072469") 
      dat1 <- dat1[, colnames(dat1) %in% col.id.GSE43865]
      print(dim(dat1))
  }
  else if (dat == "GSE35288")
  {
      dup <- c("GSM865212.1", "GSM865212.2", "GSM865213.1", "GSM865213.2", 
             "GSM865214.1", "GSM865214.2", "GSM865215.1", "GSM865215.2",
             "GSM865216.1", "GSM865216.2", "GSM865217.1", "GSM865217.2")
      rm.dup <- which(colnames(dat1) %in% dup)
      dat1 <- dat1[, -rm.dup]
  }
  
  nb.dat <- prepare.nb.data(as.matrix(dat1))
  res <- nbp.test(nb.dat$counts, grp.ids = trt, grp1=grp1, grp2=grp2, norm.factors = nb.dat$norm.factors)
  alpha = 0.05;
  sig.res = res$q.values < alpha;
  table(sig.res);
  return (list(a = table(sig.res), b = res))
}


DEglm <- function(dat, trt, beta){
  file <- paste("/Users/Bin/Google Drive/Study/Thesis/NCBI/Data/", dat, ".Rsubread.txt", sep ="")
  dat1 <- read.table(file)
  dat1 <- dat1[rownames(dat1) %in% genes, ]
  nb.dat <- prepare.nb.data(as.matrix(dat1))
  x <- model.matrix(~as.factor(trt))
  res <- nb.glm.test(nb.dat$counts, x, beta0 = beta, tests = "HOA")
  alpha = 0.05;
  sig.res = res$test.results$HOA$q.values < alpha
  table(sig.res);
  return (list(a = table(sig.res), b = res))
  
}

# DEglm("GSE53078",  rep(c(1:2), each=2), c(NA, 0))$a # 20430  3949 

genes <- rownames(readRDS("/Users/Bin/Dropbox/Zhuo/Research/Project2014/DATA_FINALVERSION/seedling.columbia.use.reference.iter.00.rds")$count)  
# seedling

DEtest("GSE32202",  rep(1:3, each=2), 1, 2)$a         #  
DEtest("GSE37159",  rep(1:4, each=2), 1, 2)$a         #  
DEtest("GSE38400",  rep(c(1:4), each=3), 1, 2)$a      # 
DEtest("GSE41766",  rep(1:3, each=2), 1, 2)$a         #  
DEtest("GSE43703", rep(1:2, 2), 1, 2)$a               #   
DEtest("GSE43865", rep(1:2, each=3), 1, 2)$a          #  
DEtest("GSE51119", rep(1:5, each=2), 1, 2)$a          #  
DEtest("GSE51772", rep(1:4, each=2), 1, 2)$a          #  
DEtest("GSE53078",  rep(c(1:2), each=2), 1, 2)$a      #  
DEtest("GSE60835", rep(c(1:2), each=3), 1, 2)$a       #  



## leaf 
DEtest("GSE36626", rep(1:2, 2), 1, 2)$a       #  4466 18053 
mod <- DEtest("GSE36626", rep(1:2, 2), 1, 2)
plot(mod$b)

DEtest("GSE39463", rep(1:2, each=3), 1, 2)$a  # 21410    28 
mod2 <- DEtest("GSE39463", rep(1:2, each=3), 1, 2)
plot(mod2$b)
DEtest("GSE48235", rep(1:3, each=2), 1, 2)$a  #  7444 14157 
DEtest("GSE51304", rep(1:9, each=2), 1, 2)$a  # 20684   589  
DEtest("GSE54677", rep(1:10, each=2), 1, 2)$a # 19791   447  

# tissue
DEtest("GSE35288", rep(1:2, each=3), 1, 2)$a        
DEtest("GSE35408", rep(1:5, 2), 1, 2)$a        
DEtest("GSE52966", rep(1:9, 2), 1, 2)$a        
DEtest("GSE56326", c(1,1, 2, 2, 2, 3, 3, 3), 2, 3)$a        
DEtest("GSE59167", c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4), 1, 2)$a       
DEtest("GSE59637", rep(1:2, each=2), 1, 2)$a       
DEtest("GSE62799", rep(1:2, each=3), 1, 2)$a        
DEtest("GSE63355", rep(1:8, each=2), 1, 2)$a       



resu <- DEtest("GSE36626", rep(1:2, 2), 1, 2)       #  4466 18053 
par(mfrow=c(3, 2))
plot(resu$b)

