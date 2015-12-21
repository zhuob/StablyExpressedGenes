#'   Plot the rank of a given gene list 
#' 
#' 
#' @title Plot the rank of genes
#'
#' @param genelist a data frame containing the gene IDs
#' @param var.comp ranked genes by sum of variance components.
#' @param col  the color of the dot
#' @param figure.num  A string that is used to name the plot
#' @return a figure
#'


plot.rank <- function(genelist, var.comp, col="blue", figure.num)
{
    
    merged <- merge(genelist, var.comp, by = "Gene")
    merged1 <- data.frame(merged)
    
    ggplot(merged1, aes(Gene, Rank)) +
    geom_point(color=col)+
    scale_y_continuous(limits=c(0, 25000)) +
    theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1)) +
    labs(title=figure.num)
    
}

