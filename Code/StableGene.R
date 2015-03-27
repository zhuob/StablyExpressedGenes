
## load the data
source("/Users/Bin/Disk1/Stat/Research/Project2014/Project1/R/prepare.data.R")
source("/Users/Bin/Disk1/Stat/Research/Project2014/Project1/R/catchToList.R")

##
library(glmmADMB)

y0 <- as.numeric(stable2[13320, ])  # glmer converge only
y1 <- as.numeric(stable2[36, ])  # glmmadmb converge only
y2 <- as.numeric(stable2[15070,]) # glmer.nb converge only
y3 <- as.numeric(stable2[6124,])  # none would converge
y4 <- as.numeric(stable2[12021,]) # wield example
y <- as.numeric(stable2[4004, ]) # everyone is happy
id <- 1:dim(stable2)[2]


library(gamlss.mx)


library(lme4)
library(optimx)
##  random-effects Poisson

### mod1 and mod2 has similar results

## trt nested in group 
id <- 1:length(y)
mod1 <- glmer(y ~ 1  + (1|group/trt) + (1|id),
              offset=(log(off.set)), family = poisson, 
              control=glmerControl(optimizer="bobyqa",boundary.tol=1e-7) # advised by Ben Bolker
            # http://stackoverflow.com/questions/21344555/convergence-error-for-development-version-of-lme4
           )



relgrad <- with(mod1@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# different optimizer???
# g0.nlminb <- update(mod1,control=glmerControl(optimizer="optimx",
#                                              optCtrl=list(method="nlminb")))
# g0.LBFGSB <- update(mod1,control=glmerControl(optimizer="optimx",
#                                              optCtrl=list(method="L-BFGS-B")))



# however, it is suggested by Douglas Bates to write the nested effect this way
# (1|group/trt) = (1|group) + (1|group:trt)
# and (1|id) is used to take care of overdispersion
mod.prefer <-  glmer(y ~ 1  + (1|group) + (1|group:trt) + (1|id),
                     offset=(log(off.set)), 
                     control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)),
                     family = poisson)

## negative binomial regression in lme4 package
mod2 <- glmer.nb(y~ 1  + (1|group/trt)+offset(log(off.set)) )

dim(stable2)
phi.mom <- (apply(stable2, 1, var)-apply(stable2, 1, mean))/apply(stable2, 1, mean)^2
length(phi.mom)
id0 <- which(phi.mom < 0.5)
length(id0)

stable3 <- stable2 

n <- 100
set.seed(14285)
samp <- sample(1:dim(stable3)[1], n)

v <- v2 <- phi <- c()

for ( i in 1: 1){
  y <- as.numeric(stable3[samp[i], ])
  print(i)
  mod <- glmmadmb(y~1 + offset(log(off.set)), random = ~(1|group/trt), 
                  zeroInflation =F, link="log", family="nbinom")
 
  v[i] <- as.numeric(mod$S[1]) # the variance of random effect for group
  v2[i] <- as.numeric(mod$S[2]) # the variance of random effect for trt

  phi[i] <- 1/as.numeric(mod$alpha)
}


#### run regression models for all genes 
warn <- rep(0, dim(stable2)[1])
var1 <- var2 <- var3 <- c()
system.time(for ( i in 1:dim(stable2)[1]){
  y <- as.numeric(stable2[i, ])
  id <- 1:length(y)
  # print(i)
  mod1 <-  catchToList( glmer(y ~ 1  + (1|group) + (1|group:trt) + (1|id),
                              offset=(log(off.set)), 
                              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)),
                              family = poisson)
  ) # catch warnings or errors
  var1[i] <- sqrt(as.numeric(mod1$value$varcor[1])) # variance for individual random effect
  var2[i] <- sqrt(as.numeric(mod1$value$varcor[2])) # variance for trt effect
  var3[i] <- sqrt(as.numeric(mod1$value$varcor[3])) # variance for group effect
  
  if (mod1$warnings[length(mod1$warnings)] != "NO"  ) ## there is warning
  {
    warn[i] = 1 # if the relative gradients are still large, then give a warning
    
  } 
})

## save the estimated standard errors for all genes
var.data <- data.frame(gene = rownames(stable2), mu= mu,
                       individual = var1, trt=var2, lab=var3, warn=warn)
var.data$gene <- 
warn.index <- var.data[which(var.data[, 5]==1),][, 1] # the genes that give warnings

var.data$sum <- var.data$individual + var.data$trt + var.data$lab
var.data$rank <- rank(var.data$sum)
# var.data$sum_ind_trt <- var.data$individual + var.data$trt
# var.data$rank_ind_trt <- rank(var.data$sum_ind_trt)
# 
# expect <- exp(var.data$mu + mean(log(off.set)) + var.data$sum/2)
# std.err <- sqrt(expect + (exp(var.data$sum)-1)*expect^2 )
# var.data$CV <- std.err/expect
# var.data$CV.rank <- rank(var.data$CV)
write.table(var.data, "temp_var-8-exp-1.1.7.txt")


c("a", "B")
ref <- c("At4g27960", "At5g53300", "At1g13320","At1g59830","At1g58050", "At5g55840","At5g08290",
        "At1g13320", "At3g53090", "At4g38070","At5g60390","At1g07920", "At1g07930","At1g07940",
        "At5g46630", "At1g13320", "At5g15710","At5g12240","At5g60390", "At3g25800","At4g24550", 
        "At3g12590", "At5g65050", "At4g27960","At5g53300","At2g32170", "At1g62930","At5g55840",
        "At5g27630", "At3g26520", "At5g49720","At1g13320","At2g28390", "At1g62930","At4g26410",
        "At3g01150", "At3g28320", "At1g62930")

gene.ref <- data.frame(gene = unique(toupper(ref)))

compare <- join(var.data, gene.ref, by ="gene", type= "right")
compare[, -(2:6)]

matplot(var.data[which(var.data$lab<10), 3:5])
var.data[which(var.data$gene=="AT4G05760")]

gene.ref <- as.matrix(read.table("/Users/Bin/Dropbox/Zhuo/stably-expressed-genes.txt", header=F,
                                  stringsAsFactors=FALSE))

ids = match(gene.ref[6:24], rownames(nb.data$counts))
nb.data$ref.freq[ids,] * 1e6;

colnames(gene.ref) <- "gene"

nbdata <- nb.data$counts

mean(var.data$sum<1)

ids <- which(var.data$rank <= 100)
write.table(var.data$gene[ids], "stably-expressed-genes-100.txt", quote=F)
var.data[ids, -(2:5)]


var.data <- var.data[sort(var.data$rank),]
write.table(var.data,"/Users/Bin/Dropbox/Zhuo/rank-gene-expression-by-var.txt" )

as.numeric(stable2[15332,])

stable.9lab <- read.table("stably-expressed-genes-100.txt", header =T)
dim(stable.9lab)

overlap <- intersect(stable.9lab[,1], var.data$gene[ids])
length(overlap)


stable.expressed <- as.matrix(read.table("/Users/Bin/Dropbox/Zhuo/identified-by-jeff.txt", header=F,
                                 stringsAsFactors=FALSE))

overlap2 <- intersect(stable.expressed[, 1], var.data$gene[ids])
overlap2
range(var.data$lab)




#  compare 8 labs and 9 labs
x1 <- read.table("temp_var-8-exp-1.1.7.txt", header=T)
x2 <- read.table("temp_var1.1.7.txt", header=T)

range(x2$lab)

x2$sum <- apply(x2[, 3:5], 1, sum)
x2$rank <- rank(x2$sum)

id1 <- which(x1$rank <= 2000)
id2 <- which(x2$rank <= 200)

overlap.8.9 <- intersect(x1$gene[id1], x2$gene[id2])
length(overlap.8.9)
"AT5G44340" %in% x2$gene[id2]
as.numeric(stable2[rownames(stable2)=="AT4G17410",])

ref.paper <- read.table("/Users/Bin/Dropbox/Zhuo/Research/Project2014/manualscript/Reference/ref.100.gene.paper.txt", header=T)
head(ref.paper)

overlap <- intersect(x2$gene[id2],ref.paper[,1])
length(overlap)
quantile(x2$sum)
length(intersect(overlap, ref.paper[, 3]))


