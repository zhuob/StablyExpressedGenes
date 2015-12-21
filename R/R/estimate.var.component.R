#'   run poisson regression with random effect , and keep track of warnings
#'   or errors of regression for a single gene
#'
#'
#' @title Poisson regression with random effect
#'
#' @param expr  a regression expression
#' @return a list
#' \item{value} a list returned by the regression model
#' \item{warnings} warning message. If none, "NO" is returned.
#' \item{error} error message. If none, "NO" is returned.
#'


## this function is used to save warnings and errors as output
# http://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function

catchToList <- function(expr) {
    val <- NULL
    myWarnings <- "NO"
    wHandler <- function(w) {
        myWarnings <<- c(myWarnings, w$message)
        invokeRestart("muffleWarning")
    }
    myError <- "NO"
    eHandler <- function(e) {
        myError <<- e$message
        NULL
    }
    val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
    list(value = summary(val), warnings = myWarnings, error=myError)
} 









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

estimate.var.component <- function(data, reference = NULL, group, trt, filter.factor=3)
{
    
    groups <- as.factor(group)
    treat <- as.factor(trt)
    
    stable <- as.matrix(data[rowSums(data)>=filter.factor*dim(data)[2], ])
    
    
    ##  the reference set
    if (is.null(reference)) {
        ref <- stable
    }
    else
    {
       ref <- data[row.names(data) %in% reference, ]
    }
    
    
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
    
    cat(paste("total number of genes to be evaluated: ", nrow(stable)))
    cat("\n runing the GLMM for genes...  \n")
    for ( i in 1:dim(stable)[1]){
        
        cat('\r', i)  # show how many iterations have been executed
        
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
    
    return(list(count=stable, var.comp=var.comp, norm.factor = norm_factor, trt= trt, lab = group))
}
