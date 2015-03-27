##' @title Converting a data.frame whose row.names are gene into  a new data
##' containing the row.names as a new column, facilitating to merge by gene
##' @param  data  input data
##' 
##' @return  return a new data set
##'
##' @author
##' 

create.geneColumn <- function(data)
{
  Gene <- row.names(data) # get the gene names
  d2 <- data.frame(Gene, data)
  row.names(d2) <- 1:dim(data)[1]
  return(d2)
}


