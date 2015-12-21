## GuestLecture

library(NBPSeq)



# load the data
data(arab);
nb.data <- prepare.nb.data(arab)
print(nb.data);
## Specify treatment structure
grp.ids = as.factor(c(1, 1, 1, 0, 0, 0));
x = model.matrix(~grp.ids);
## Specify the null hypothesis
## The null hypothesis is beta[1]=0 (beta[1] is the log fold change).
beta0 = c(NA, 0);
## Fit NB regression model and perform large sample tests.
## The step can take long if the number of genes is large
fit = nb.glm.test(arab, x, beta0, subset=1:50);
  ## The result contains the data, the dispersion estimates and the test results
print(str(fit));
## Show HOA test results for top ten genes
subset = order(fit$test.results$HOA$p.values)[1:10];
cbind(fit$data$counts[subset,], fit$test.results$HOA[subset,]);
## Show LR test results
subset = order(fit$test.results$LR$p.values)[1:10];
cbind(fit$data$counts[subset,], fit$test.results$LR[subset,]);





## generate normal 
n <- 100
x1 <- rnorm(n, mean=20, sd= sqrt(21))
# x1 <- rpois(n, 20)
# var(x2) <- mu + mu^2/size
x2 <- rnbinom(n, size=20, mu=20)

df <- data.frame(obs= 1:n, Normal=x1, NB= x2)
library(reshape2)
df <- melt(df, id="obs")
library(ggplot2)
ggplot(data=df, aes(x= value, fill=variable)) + 
    geom_histogram(binwidth=2, alpha=0.5, position="identity")









## Load Arabidopsis data
data(arab);
## Specify treatment groups
## grp.ids = c(1, 1, 1, 2, 2, 2); # Numbers or strings are both OK
grp.ids = rep(c("mock", "hrcc"), each=3);
## Estimate normalization factors
norm.factors = estimate.norm.factors(arab);
print(norm.factors);
## Prepare an NBP object, adjust the library sizes by thinning the
## counts. For demonstration purpose, only use the first 100 rows of
## the arab data.
set.seed(999);
obj = prepare.nbp(arab[1:100,], grp.ids, lib.size=colSums(arab), norm.factors=norm.factors);
print(obj);
## Fit a dispersion model (NBQ by default)
obj = estimate.disp(obj);
plot(obj);
## Perform exact NB test
## grp1 = 1;
## grp2 = 2;
grp1 = "mock";
grp2 = "hrcc";
obj = exact.nb.test(obj, grp1, grp2);
## Print and plot results
print(obj);
par(mfrow=c(3,2));
plot(obj);
