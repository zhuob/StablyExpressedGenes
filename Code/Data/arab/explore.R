
arab.3 <- read.csv("/home/zhuob/Dropbox/Zhuo/Project_2014/arab/arabCtrl.csv", header=T)


arab.3 <- read.csv("C:/Users/zhuob/Dropbox/Zhuo/Project_2014/arab/arabCtrl.csv", header=T)

s1 <- apply(arab.3[,3:5],1, mean)
s2 <- apply(arab.3[, 6:7], 1, mean)
s3 <- apply(arab.3[, 8:10], 1, mean)
s4 <- apply(arab.3[, 11:13], 1, mean)


pdf("C:/Users/zhuob/Dropbox/Zhuo/FindData/arab/log_mean.pdf",width=7,height=5)
par(mfrow=c(2, 3))
plot(log(s1), log(s2))
plot(log(s1), log(s3))
plot(log(s1), log(s4))
plot(log(s2), log(s3))
plot(log(s2), log(s4))
plot(log(s3), log(s4))
dev.off()


library(NBPSeq)

arab.dat <- arab.3[, 3:dim(arab.3)[2]] 
# keep the count read column only

head(arab.dat)
arab.dat <- as.matrix(arab.dat)

nb.data <- prepare.nb.data(arab.dat)
names(nb.data)





arab.nozero <- arab.dat[, which(arab!=0)]

arab.dat[arab.dat==0]=NA
rho <- cor(arab.dat, use="complete.obs")
round(rho, 2)
c(min(rho), max(rho))
rho.log <- cor(log(arab.dat), use="complete.obs")
round(rho.log, 2)

round(rho.log, 2)


plot.nb.data(nb.data)
head(nb.data$counts)
grp.ids <- c(1, 1, 1, 2,2, 3,3,3, 4, 4, 4)
obj <- prepare.nbp(as.matrix(arab.dat), grp.ids)
obj2 <- estimate.disp(obj)


### estimate.dispersion is the function that gives dispersion estimate
# similaryly, we can estimate dispersion here
# specify the design matrix
x <- matrix(c(
  rep(1, 11), 
  rep(1, 3), rep(0, 8),
  rep(0, 3), 1, 1, rep(0, 6),
  rep(0, 5), 1, 1, 1, rep(0, 3),
  rep(0, 8), rep(1, 3)
  ), 11, 5
  )


# estimate dispersion by NBQ model

ptm <- proc.time()
disp <- estimate.dispersion(nb.data, x)
proc.time() - ptm
# may cause 10,000 seconds to run 
# therefore I store it in a txt file
write.table(disp$estimates, "/home/zhuob/Dropbox/Zhuo/Project_2014/arab/dispersion.txt")
dd <- read.table("/home/zhuob/Dropbox/Zhuo/Project_2014/arab/dispersion.txt")


head(disp$estimates) # pull out the first few entries
head(nb.data$rel.fre)
rel <- log( as.vector(nb.data$rel.fre))
disper <- log(disp$estimates)
length(disper)

plot(rel, disper, pch=20, cex=0.5)

## create MA plot
names(obj)


s1 <- rowSums(arab.dat[, 1:3])
s2 <- rowSums(arab.dat[, 4:5])
s3 <- rowSums(arab.dat[, 6:8])
s4 <- rowSums(arab.dat[, 9:11])
m <- obj$pseudo.lib.s[1]
log.fc.12 <- log2( (s2/2/m)/(s1/3/m)  )
mu <- rowMeans(arab.dat[, 1:5])
pi.12 <- mu/m

pi.13 <- (s1+s3)/6/m
log.fc.13 <- log2( (s3/3/m)/(s1/3/m))

pi.14 <- (s1+s4)/6/m
log.fc.14 <- log2( (s4/3/m)/(s1/3/m))

pi.23 <- (s2+s3)/5/m
log.fc.23 <- log2( (s3/3/m)/(s2/2/m))

pi.24 <- (s2+s4)/6/m
log.fc.24 <- log2( (s4/3/m)/(s2/2/m))

pi.34 <- (s3+s4)/6/m
log.fc.34 <- log2( (s4/3/m)/(s3/3/m))


pdf("/home/zhuob/Dropbox/Zhuo/Project_2014/arab/maplot.pdf"
    ,width=7,height=5)

plot(pi.12, log.fc.12, cex=0.5, pch=3, log="x", col="grey", 
     main="MA plot of group 1 vs. group2")

plot(pi.13, log.fc.13, cex=0.5, pch=3, log="x", col="cyan", 
     main="MA plot of group 1 vs. group 3")

plot(pi.14, log.fc.14, cex=0.5, pch=3, log="x", col="red", 
     main="MA plot of group 1 vs. group 4")

plot(pi.23, log.fc.23, cex=0.5, pch=3, log="x", col="blue", 
     main="MA plot of group 2 vs. group 3")

plot(pi.24, log.fc.24, cex=0.5, pch=3, log="x", col="pink", 
     main="MA plot of group 2 vs. group 4")

plot(pi.34, log.fc.34, cex=0.5, pch=3, log="x", col="black", 
     main="MA plot of group 3 vs. group 4")

dev.off()



## package glmmADMB deals with count data via glmm
library(glmmADMB)

### get a random subset of arab.dat to run this glmm model
#   treat gene id as a fixed effect, and group as random effect.

group <- as.factor(c(1, 1,1, 2,2,3,3,3, 4, 4, 4))
arab.nozero <- as.matrix(arab.dat[rowSums(arab.dat)!=0, ])

n <- 200
set.seed(21)
id <- sample (1:dim(arab.nozero)[1], n)
vari <- c()
for (i in 1:n){
  y <- arab.nozero[id[i], ]
  a <- glmmadmb(y ~1, random=  ~1 | group, zeroInflation=F, 
            link="log",  family = "nbinom")
  vari[i] <- as.numeric(a$S)  #### output the variance
  }

hist(log(vari), main="number of genes = 200")

range(vari)












