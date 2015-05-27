#'  estimate norm factors by most stably expressed genes 
#' 
#' @title Estimate norm. factor
#'
#' @param var.comp ranked genes by sum of variance components.
#' @param data  data frame of read count
#' @param n  number of top stably expressed genes
#' @return a list
#' \item{s} the number of genes matched
#' \item{Gene} the gene list successfully matched
#' \item{Rank} the corresponding ranks 
#'

norm.factor <- function(var.comp, data, n)
{
  genelist <- var.comp$Gene[var.comp$Rank <= n]
  dat1 <- data[row.names(data) %in% genelist, ]
  estimate.norm.factors(dat1, lib.sizes= colSums(data) )
}