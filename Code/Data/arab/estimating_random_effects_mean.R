
arab.3 <- read.csv("/home/zhuob/Dropbox/Zhuo/Project_2014/arab/arabCtrl.csv", header=T)


arab.3 <- read.csv("C:/Users/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arabCtrl.csv", header=T)

arab.dat <- arab.3[, 3:dim(arab.3)[2]] 

## package glmmADMB deals with count data via glmm
install.packages("glmmADMB", 
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),type="source")
library(glmmADMB)

### get a random subset of arab.dat to run this glmm model
#   treat gene id as a fixed effect, and group as random effect.

group <- as.factor(c(1, 1,1, 2,2,3,3,3, 4, 4, 4))
arab.nozero <- as.matrix(arab.dat[rowSums(arab.dat)!=0, ])
arab.pre <- as.matrix(arab.dat[rowSums(arab.dat)>= dim(arab.dat)[2], ])

n <- 1000
set.seed(128)
  
vari <- c()
  # id <- sample (1:dim(arab.nozero)[1], n)
  id <- sample (1:dim(arab.pre)[1], n)
for (i in 1:n){
 # y <- arab.nozero[id[i], ]
  y <- arab.pre[id[i], ]
  a <- glmmadmb(y ~1, random=  ~1 | group, zeroInflation=F, 
                link="log",  family = "nbinom")
  vari[i] <- as.numeric(a$S)  #### output the variance
}

hist(log(vari), main=paste("number of genes =", n))

arab.pre[id[14],]

id2 <- which(log(vari)< -10)
id2
## to see which genes have estimated variance of 0
arab.pre[id[id2],]
length(id2)


# we can have glmmPQL estimate the variance, but it cannot be output.
y <- arab.pre[id[25],]
b <- glmmPQL(y ~1, random=  ~1 | group, family = quasipoisson)

a <- glmmadmb(y~1, random=  ~1 | group, zeroInflation=F, 
              link="log",  family = "nbinom")
summary(a)
## in general, a gives smaller variance than b

head(arab.dat)
y1 <- as.numeric(arab.dat[6, ])
y2 <- as.numeric(arab.dat[7502, ])
y3 <- as.numeric(arab.dat[19895, ])
y4 <- as.numeric(arab.pre[id[25],])
y5 <- as.numeric(arab.dat[23,])
y6 <- as.numeric(arab.dat[id[12],])
y7 <- as.numeric(arab.dat[7782,])
y8 <- as.numeric(arab.dat[12475,])
cbind(y1, y2, y3, y4, y5, y6, y7, y8)

a <- glmmadmb(y1~1, random=  ~1 | group, zeroInflation=F, 
              link="log",  family = "nbinom")

#  installing R-INLA package
source("http://www.math.ntnu.no/inla/givemeINLA-testing.R") 

file.show(system.file("tpl","glmmadmb.tpl",package="glmmADMB"))
y1
