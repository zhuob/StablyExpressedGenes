#'   Plot the Counts Per Million (CPM) of a given gene list from RNA-Seq data
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

##### Figure in Section 4.1 --------------------------------------

plot.cpm <- function(set, top=1e3, figure.num, breaks=50, textsize =rep(20, 4)){
  
  cpm0 <- as.vector(mean_cpm(set, top = top))
  cpm1 = cbind.data.frame(cpm0, obs = 1:length(cpm0))
  print(range(cpm0))
  f <- ggplot(data =cpm1, aes(x = cpm0)  )  + 
    geom_histogram(binwidth=max(cpm0)/breaks, colour="black", fill="white") + 
    theme(legend.position="top",
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]), 
          axis.text=element_text(size=textsize[3]), 
          axis.title=element_text(size=textsize[4],face="bold")) +
    labs(x = "CPM", y = "Frequency", title = figure.num)
  return(f)
}


plot.gene <- function(genelist, set, lower=1, upper=1e4, figure.num, textsize = rep(20, 4)){
    
    lab <- set$lab
    count <- set$count

    
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
    data.plot$lab <- lab

    
    p <- ggplot(data.plot) +
    geom_line(aes(x=Sample,y=CPM,group=Gene, color=Gene, linetype= Gene), size=1)+
    facet_grid(.~lab,scales="free", space = "free") +   # divide the experiments
    theme(axis.text.x = element_text(size=15,angle = 90, hjust = 1)) +
    labs(title=figure.num) +
    scale_y_log10(limits=c(lower, upper)) +    # log scale
    theme(legend.position="top", 
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]), 
          axis.text=element_text(size=textsize[3]), 
          axis.title=element_text(size=textsize[4],face="bold")) +
     scale_x_discrete(breaks = seq(1, length(x), 2) ) + 
      guides(fill = guide_legend(title = "", keywidth = 1, keyheight = 1),
             linetype=guide_legend(keywidth = 3, keyheight = 1),
             colour=guide_legend(keywidth = 3, keyheight = 1)) +
    #+
          scale_color_brewer(palette="Set1")
          # scale_fill_continuous(name = "MYTGE")
    #+ theme(legend.position = c(2, 0.8*upper))
    #  + geom_text( aes(1, 9000, label = "A"))
    p
}




##### Figure in Section 4.2 --------------------------------------

TopGene <- function(set1, set2, set3, cze, dek, geNorm, xtext, textsize = rep(20, 4)){
  
  
  var_set1 <-set1$var.comp
  var_set2 <- set2$var.comp
  var_set3 <- set3$var.comp
  
  # top stably expressed genes
  gene.cze <- cze[, 1]
  gene.dek <- dek[, 1]
  
  #gene.geNORM <- geNorm_NormFinder_CV$Gene[geNorm_NormFinder_CV$rank <=100]
  
  gene.geNORM <- geNorm$Gene[geNorm$Rank <= 100]
  #gene.NormFinder <- geNormNormFinder$Gene[geNormNormFinder$Rank_rho <= 100]
  
  gene.set2 <- var_set2$Gene[var_set2$Rank <= 100]  ##  
  gene.set3 <- var_set3$Gene[var_set3$Rank <= 100]  ##  
  
  
  #rank.seq <- c(100, 200, 500, 1000, 2000, 3500, 5000)
  rank.seq <- seq(100, 5000, by = 100)
  
  set2.num <- set3.num <- cze.num <- dek.num <- genorm <- Normfinder<-  c()
  
  for ( i in 1: length(rank.seq))
  {
    gene.set1 <- var_set1$Gene[var_set1$Rank <= rank.seq[i]]
    set2.num[i] <- round(length(intersect(gene.set1, gene.set2))/length(gene.set2), 3)
    set3.num[i] <- round(length(intersect(gene.set1, gene.set3))/length(gene.set3), 3)
    cze.num[i] <- round(length(intersect(gene.set1, gene.cze))/length(gene.cze), 3)
    dek.num[i] <- round(length(intersect(gene.set1, gene.dek))/length(gene.dek), 3)
    genorm[i] <- round(length(intersect(gene.set1, gene.geNORM))/length(gene.geNORM), 3)
  #  Normfinder[i] <- round(length(intersect(gene.set1, gene.NormFinder))/length(gene.NormFinder), 3)
    
  }
  
  
  pct.data <- data.frame(rank=rank.seq, genorm, set2.num, set3.num, cze.num, dek.num )#, Normfinder)
  colnames(pct.data) <- c("Rank",xtext[5], xtext[2], xtext[3],xtext[6], xtext[7])#, "NormFinder")
  pct.data <- melt(pct.data, id="Rank", variable.name = "List",  value.name = "Percentage")
  
  
  ggplot(pct.data,  aes(x= Rank, y = Percentage, color = List,  linetype = List)) + 
    geom_line(size = 1.5) +
    labs(x=xtext[1], y=xtext[4]) +
    theme(legend.position="top", 
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]), 
          axis.text=element_text(size=textsize[3]), 
          axis.title=element_text(size=textsize[4],face="bold")) +
       scale_y_continuous(labels = percent)  +
    guides(fill = guide_legend(keywidth = 1, keyheight = 1),
             linetype=guide_legend(keywidth = 3, keyheight = 1),
             colour=guide_legend(keywidth = 3, keyheight = 1))
    
    #  scale_fill_continuous(guide = "legend") + 
  #   guides(fill = guide_legend(keywidth = 10, keyheight = 1))
  }








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

##### Figure in Section 4.4 --------------------------------------


norm.factor <- function(target_var, target_count, n)
{
    genelist <- target_var$Gene[target_var$Rank <= n]
    dat1 <- target_count[row.names(target_count) %in% genelist, ]
    estimate.norm.factors(dat1, lib.sizes= colSums(target_count) )
}



plot.pair.normfactor <- function(set, textsize = rep(20,4)){
  
  target_var <- set$var.comp
  target_count <- set$count
  
  norm.10 <- norm.factor(target_var, target_count, 10)
  norm.1e2 <- norm.factor(target_var, target_count, 1e2)
  norm.1e3 <- norm.factor(target_var, target_count, 1e3)
  norm.1e4 <- norm.factor(target_var, target_count, 1e4)
  
  
  norm.all <- norm.factor(target_var, target_count, dim(target_var)[1])
  
  norm.factor1 <- estimate.norm.factors(target_count, lib.sizes=colSums(target_count))
  
  x <- data.frame(top10=as.numeric(norm.10), top1e2 = as.numeric(norm.1e2),
                  top1e3=as.numeric(norm.1e3), top1e4 = as.numeric(norm.1e4),
                  all= as.numeric(norm.all))
  
  x <- round(x, 2)
  
  abcd <- ggpairs(x, upper = "blank", axisLabels = "none") +
    theme(legend.text = element_text(size = textsize[1]),
          panel.grid.major = element_blank(), 
          axis.ticks=element_blank(), 
          panel.border = element_rect(linetype = "dashed", colour = "blue", fill = NA), 
          plot.title = element_text(size = textsize[2]), 
          axis.text=element_text(size=textsize[3]), 
          axis.title=element_text(size=textsize[4],face="bold")) 
  
  ## fill the uppper part manually
  for ( i in 1:(ncol(x)-1)){
    for ( j in (i + 1) :ncol(x) ){
      x1 <- x[, i]
      x2 <- x[, j]
      r <- round(cor(x1, x2), 3)
      add_plot <- ggally_text(r, size = 10)
      abcd[i, j] <- add_plot
    }
  }
  abcd 
}



#######################################################
# STACKED BAR PLOT
#######################################################
## target_var:  the variance component data set, returned by estimate.var.component()
## gene_ids: the genes to be ploted
## percent:  whether use percentage of variances


plot.stackedBar <- function(gene_ids, set, figure.num, textsize = rep(20,4), legend=F){
  
    target_var <- set$var.comp
    library(scales)

  
    colnames(target_var)[2:4] <- c("sample", "treatment", "experiment")
    mdata <- melt(target_var[target_var$Gene %in% gene_ids, c(1:4)], id=c("Gene"), 
                      variable.name = "Source",  value.name = "value")
    id <- which(mdata$value > 1e-6)
    mdata <- mdata[id, ]
      
   p1 <-  ggplot(mdata,aes(x = Gene, y = value,fill = Source)) +
        geom_bar(stat = "identity") +
        theme( legend.text = element_text(size = textsize[1]),
               legend.position="top",
               plot.title = element_text(size = textsize[2]), 
               axis.text.x = element_text(size=textsize[3],angle = 45, hjust = 1),
               axis.text.y = element_text(size = textsize[3]),
               axis.title=element_text(size=textsize[4],face="bold")) +
          labs( y="Variance", title=figure.num) +
          scale_y_continuous() +
          scale_fill_grey() + 
          guides(linetype = guide_legend(keywidth = 3, keyheight = 1), 
                 colour = guide_legend(keywidth = 3, keyheight = 1))
   
  p1
}



# this function is used to extract the legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}



plot.density <- function(set, figure.num, y_label="Density", textsize = rep(30, 4), legend = F)
{
    var.s <- set$var.comp
    var.s.percent <- data.frame(matrix(nrow=dim(var.s)[1], ncol=0))
    var.s.percent$Gene <- var.s$Gene
    var.s.percent$sample <- var.s$individual/var.s$sum
    var.s.percent$treatment <- var.s$trt/var.s$sum
    var.s.percent$experiment <- var.s$lab/var.s$sum

    mdata <- melt(var.s.percent[, c(1:4)], id=c("Gene"))
    colnames(mdata) <- c("Gene", "Source", "percentage")
    mdata$Source <- paste("between-", mdata$Source, sep= "")

    

    
    A1 <- ggplot(mdata, aes(x=percentage, color = Source,  linetype = Source)) +
         geom_density(size = 1.5, show.legend=FALSE) +
      stat_density(aes(x=percentage, colour=Source),
                   geom="line",position="identity") +
       theme(legend.position="top",
            legend.text = element_text(size = textsize[1]),
            plot.title = element_text(size = textsize[2]), 
            axis.text=element_text(size=textsize[3]), 
            axis.title=element_text(size=textsize[4],face="bold") ) +
      labs(title=figure.num, x = "Percentage", y = y_label) + 
      guides(fill = guide_legend(keywidth = 1, keyheight = 1), 
             linetype = guide_legend(keywidth = 2, keyheight = 1), 
             colour = guide_legend(keywidth = 3, keyheight = 1))

   
    
    if(legend){    # plot the legend in a separate figure
      mylegend<-g_legend(A1)
      p_fig <- grid.arrange(mylegend)
    }
    
    else{
      p_fig <- ggplot(mdata, aes(x=percentage, color = Source,  linetype = Source)) +
        geom_density(size = 1.5) +
        theme(legend.position="none",
              # legend.text = element_text(size = textsize[1]),
              plot.title = element_text(size = textsize[2]), 
              axis.text=element_text(size=textsize[3]), 
              axis.title=element_text(size=textsize[4],face="bold") ) +
        guides(guide_legend(show= F) )+ 
        labs(title=figure.num, x = "Percentage", y = y_label)
      print(p_fig)
    }
    
    return(p_fig)
}





#####  specific for addressing format required by the magzine
# figure 1
plot_cpm <- function(set, top=1e3, figure.num, y_label="",  breaks=50, textsize =rep(20, 4)){
  
  cpm0 <- as.vector(mean_cpm(set, top = top))
  cpm1 = cbind.data.frame(cpm0, obs = 1:length(cpm0))
  print(range(cpm0))
  f <- ggplot(data =cpm1, aes(x = cpm0)  )  + 
    geom_histogram(binwidth=max(cpm0)/breaks, colour="black", fill="white") + 
    theme(legend.position="top",
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]), 
          axis.text=element_text(size=textsize[3]), 
          axis.title=element_text(size=textsize[4],face="bold")) +
    labs(x = "CPM", y = y_label, title = figure.num)
  return(f)
}


## figure 2
create_stable_genes <- function(genelist, set){
  
  lab <- set$lab
  count <- set$count
  
  
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
  data.plot$lab <- lab
  
  return(data.plot)
}

five_genes_plot <- function(data, lower=1, upper=1e4, figure.num, x_label = "", textsize = rep(20, 4)){
  nx <- length(unique(data$Sample))
  p <- ggplot(data) +
    geom_line(aes(x=Sample,y=CPM,group=Gene, color=Gene, linetype= Gene), size=1)+
    facet_grid(.~lab,scales="free", space = "free") +   # divide the experiments
    # theme(axis.text.x = element_text(size=15,angle = 90, hjust = 1)) +
    labs(title=figure.num, x = x_label) +
    scale_y_log10(limits=c(lower, upper)) +    # log scale
    theme(legend.position="top", 
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]), 
          axis.text.y=element_text(size=textsize[3]), 
          axis.text.x=element_text(size=textsize[3], angle=90, hjust =1),
          axis.title=element_text(size=textsize[4],face="bold")) +
    scale_x_discrete(breaks = seq(1, nx, 2) ) + 
    guides(fill = guide_legend(title = "", keywidth = 1, keyheight = 1),
           linetype=guide_legend(keywidth = 3, keyheight = 1),
           colour=guide_legend(keywidth = 3, keyheight = 1))  +   scale_color_brewer(palette="Set1")
  
  return(p)
}




stacked_bar <- function(gene_ids, set, figure.num, x_label = "", textsize = rep(20,4)){
  
  target_var <- set$var.comp
  library(scales)
  
  
  colnames(target_var)[2:4] <- c("sample", "treatment", "experiment")
  mdata <- melt(target_var[target_var$Gene %in% gene_ids, c(1:4)], id=c("Gene"), 
                variable.name = "Source",  value.name = "value")
  id <- which(mdata$value > 1e-6)
  mdata <- mdata[id, ]
  
    
    p1 <-  ggplot(mdata,aes(x = Gene, y = value,fill = Source)) +
      geom_bar(stat = "identity") +
      theme(
          # legend.text = element_text(size = textsize[1]),
             legend.position="none",
             plot.title = element_text(size = textsize[2]), 
             axis.text.x = element_text(size=textsize[3],angle = 45, hjust = 1),
             axis.text.y = element_text(size = textsize[3]),
             axis.title=element_text(size=textsize[4],face="bold")) +
      labs( y="Variance", title=figure.num, x = x_label) +
      scale_y_continuous() +
      scale_fill_grey() + 
      guides(linetype = guide_legend(keywidth = 3, keyheight = 1), 
             colour = guide_legend(keywidth = 3, keyheight = 1))

  
  p1
}





