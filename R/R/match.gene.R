#'  match the top ranked genes 
#' 
#' 
#' @title Match the top ranked genes
#'
#' @param genelist a data frame containing the top gene IDs
#' @param var.comp ranked genes by sum of variance components.
#' @param top  how many genes are to be matched 
#' @return a list
#' \item{s} the number of genes matched
#' \item{Gene} the gene list successfully matched
#' \item{Rank} the corresponding ranks 
#'

match.gene <- function(genelist, var.comp, top)
{
  match <- var.comp[var.comp$Rank <= top, ]
  genelist <- as.vector(genelist)
  result <- merge(genelist, match, by="Gene")
  
  lst <- list(s = dim(result)[1], 
              l = data.frame(Gene=result$Gene, Rank=result$Rank))
  return(lst)
}