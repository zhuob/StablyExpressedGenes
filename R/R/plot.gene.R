#'   Plot the Counts Per Million (CPM) of a given gene list from RNA-Seq data 
#' 
#' 
#' @title Plot CPM of a gene list
#'
#' @param genelist a list of gene IDs
#' @param data read count data set to be plotted
#' @param lower  lower limit for the y-axis
#' @param upper  upper limit for the y-axis
#' @param figure.num  A string that is used to name the plot
#' @return a figure
#'



plot.gene <- function(genelist, count, lower=1, upper=1e4, figure.num)
{
    library(ggplot2)
    library(reshape2)
    
    norm.factors <- estimate.norm.factors(count)
    nb.data <- prepare.nb.data(count, norm.factors=norm.factors)
    off.set <- as.numeric(nb.data$eff.lib.sizes)
    
  ids <- which(row.names(count) %in% genelist)
  
  x <- 1:dim(count)[2]
  y <- t(sweep(count[ids,], 2, off.set, "/"))*1e6 # CPM
  
  data.plot <- data.frame(y, Sample=x)
  data.plot <- melt(data.plot, measure.vars=c(1:length(ids)), id.vars=c("Sample"),
                    variable.name = "Gene", value.name = "CPM")
  data.plot$Sample <- as.factor(data.plot$Sample)
  
  p <- ggplot(data.plot, aes(Sample, CPM, group=Gene, colour=Gene))
  p + geom_line() +
    theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
    labs(title=figure.num) +
    scale_y_log10(limits=c(lower, upper)) +
    ## theme on top
   theme(legend.position="top")
  #+ theme(legend.position = c(2, 0.8*upper))
  #  + geom_text( aes(1, 9000, label = "A"))
  
}
