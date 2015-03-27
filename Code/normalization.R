

setwd("/Users/Bin/Dropbox/Zhuo/Research/Project2014/Result")

ref.100 <- read.table("ref.100.gene.paper.txt", header=T)
all <- read.table("var.comp.all.txt", header=T)
tissue <- read.table("var.comp.4tissue.txt", header=T)
seedling8 <- read.table("var.comp.seedling.8lab.txt", header=T)

ref.50 <- read.table("ref.50.txt", header=T)

seedling7 <- read.table("var.comp.seedling.7lab.txt", header=T)


nobs <- 200 # top 200

lab7.8 <- intersect(seedling7$gene[seedling7$rank.sum <= nobs],
                   seedling8$gene[seedling8$rank.sum <= nobs])
length(lab7.8)

lab7.tissue <- intersect(seedling7$gene[seedling7$rank.sum <= nobs],
                    tissue$gene[tissue$rank.sum <= nobs])
length(lab7.tissue)

lab7.100 <- intersect(seedling7$gene[seedling7$rank.sum <= nobs],
                    ref.100[, 1])
length(lab7.100)



library(ggplot2)
library(reshape2)
id.plot <- c(1:5)
ids <- which(seedling7$rank.sum %in% id.plot)
# if You want to plot five genes identified by paper
# ids <- which(row.names(stable) %in% ref.100[1:5, 1])
# ids <- which(row.names(stable) %in% ref.50[1:5, 1])


y <- t(sweep(stable[ids,], 2, off.set, "/"))
x <- 1:dim(stable)[2]


data.plot <- data.frame(y, sampleID=x)
data.plot <- melt(data.plot, measure.vars=c(1:length(id.plot)), id.vars=c("sampleID"),
                  variable.name = "Gene", value.name = "ExprVal")

ggplot(data.plot, aes(as.factor(sampleID), ExprVal, group=Gene, colour=Gene)) +
  geom_line() +
# labs(title="Top 5 stable genes by Poisson Regression (seedling)") +
# labs(title="Top 5 stable genes by Czechowski (2005)") + # ref.100
 labs(title="Top 5 stable genes by Dekkers (2012)") + # ref.100
  scale_y_log10(limits=c(1e-7, 1e-3)) 


data.seedling7 <- melt(seedling7, measure.vars=c(2:4), id.vars="gene", variable.name = "sourceOfVar",
              value.name= "var")
p <- ggplot(data.seedling7, aes(gene, var))
p + geom_point(aes(colour = factor(sourceOfVar)))+
 labs(title="Variance components of Set 2")


# 
# 
# plot(1:50, y[1, ], log="y", type="l", ylim=c(1e-7, 1e-3), 
#      xlab="sample", ylab="RPM", main= "top 5 stable genes by Poisson Regression")
# lines(1:50, y[2, ], col="magenta")
# lines(1:50, y[3, ], col="blue")
# lines(1:50, y[4, ], col="green")
# lines(1:50, y[5, ], col="red")


library(lme4)
y <- as.numeric(stable[ids[3],])
id <- 1:length(y)
mod1 <- glmer(y ~ 1  + (1|group) + (1|group:trt) + (1|id), offset=(log(off.set)), 
      control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 100000)),
      family = poisson)
summary(mod1)

mu <- as.numeric(fitted(mod1))
xx <- rpois(50, mu)
plot(1:50, xx/off.set, log="y", type="l", ylim=c(1e-7, 1e-3), main="simulated data (params estimated from real data)")
xx2 <- rpois(50, mean(mu))
lines(1:50, xx2/off.set, col="red", lty=2)
lines(1:50, y[, 1], col="blue", lty=4)


## MOM  CV 
mean.MOM <- apply(sweep(stable, 2, off.set, "/"), 1, mean)
sd.MOM <- apply(sweep(stable, 2, off.set, "/"), 1, sd)

# mean.MOM <- apply(stable, 1, mean)
# sd.MOM <- apply(stable, 1, sd)
CV <- sd.MOM/mean.MOM
rank.CV <- rank(CV)
gene.MOM <- row.names(stable)


lab7.cv <- intersect(seedling7$gene[seedling7$rank.sum <= 100], row.names(stable)[rank.CV<=100])
length(lab7.cv)
cv.100 <- intersect(ref.100[, 1], row.names(stable)[rank.CV<200])
length(cv.100)
lab7.100 <- intersect(seedling7$gene[seedling7$rank.sum <= 200], ref.100[, 1])
length(lab7.100)


## 
id3 <- which(rank.CV <= 5)

y <- t(sweep(stable[id3,], 2, off.set, "/"))

data.plot <- data.frame(y, sampleID=x)
data.plot <- melt(data.plot, measure.vars=c(1:length(id.plot)), id.vars=c("sampleID"),
                  variable.name = "Gene", value.name = "ExprVal")

ggplot(data.plot, aes(as.factor(sampleID), ExprVal, group=Gene, colour=Gene)) +
  geom_line() +
  labs(title="Top 5 stable genes ranked by CV (ascending)") +
  scale_y_log10(limits=c(1e-7, 1e-4)) 
# 


id1 <- which(seedling7$gene %in% ref.100[1:5, 1])
seedling7[id1, c(1, 7)]

id2 <- which(seedling7$gene %in% ref.50[1:5, 1])
seedling7[id2, c(1, 7)]


science <- read.table("science.txt", header=T)
science$gene <- toupper(science$Gene)
write.table(science,"science1.txt")

data12 <- merge(science, seedling7, by="gene", all.x=T)
plot(1:dim(data12)[1], data12$rank.sum)





setwd("/Users/Bin/Dropbox/Zhuo")

dif <- read.table("diff_exp_genes1.txt", header=T, fill = T)
dim(dif)
dif[dif$Gene=="AT1G10660",]

x <- grep("AT" , dif[, 1])
dif2 <- dif[x, ]
length(unique(dif2[, 1]))

dif_exp.100 <- dif[dif$Gene %in% ref.100[,1], 1]
seedling.gene <- seedling7[seedling7$rank.sum <= 100, 1]

dif_exp.seedling <-  dif[dif$Gene %in% seedling.gene, 1]
dup <- duplicated(dif_exp.seedling)
dif_exp.seedling[dup]
# AT5G18500

dif_exp.100[duplicated(dif_exp.100)]
tables(dif_exp.100)
##  AT2G16860


