#' Combine the different GSE data sets to make Set 1 , Set 2 and Set 3. 
#' 
#' 
#' @title combine different data sets togather
#'
#' @param type  a string which can be "seedling", "tissue" or "leaf"
#' @return a list
#' \item{data} a n-by-p data frame of counts, with row name being the gene IDs
#' \item{trt} a p-vector indicating the treatment factor
#' \item{group} a p-vector of experimental group IDs
#'


ArrangeData <- function(type)
{
  library(plyr)
 
  if (type == "seedling")  
  {  
    GSE37159 <- read.table("GSE37159.Rsubread.txt", header = T) 
    GSE37159 <- create.geneColumn(GSE37159)
    trt.GSE37159 <- c(rep(1:4, each=2))
    
    GSE38879 <- read.table("GSE38879.Rsubread.txt", header = T)
    GSE38879 <- GSE38879[, seq(1, 23, by =2)] # 2 ILLUMINA RUNS
    GSE38879 <- create.geneColumn(GSE38879)
    trt.GSE38879 <- c(rep(11:14, each =3))
    
    GSE43865 <- read.table("GSE43865.Rsubread.txt", header=T) 
    GSE43865 <- create.geneColumn(GSE43865)
    trt.GSE43865 <- c(rep(21:22, each=3))
    
    GSE48767 <- read.table("GSE48767.Rsubread.txt", header = T)
    GSE48767 <- create.geneColumn(GSE48767)
    trt.GSE48767 <- rep(31:34, each=3)
    
    
    GSE51119 <- read.table("GSE51119.Rsubread.txt", header = T) 
    GSE51119 <- create.geneColumn(GSE51119)
    trt.GSE51119 <- rep(41:45, each=2)
    
    GSE51772 <- read.table("GSE51772.Rsubread.txt", header = T) 
    GSE51772 <- create.geneColumn(GSE51772)
    trt.GSE51772 <- rep(51:54, each=2)
    
    GSE53078 <- read.table("GSE53078.Rsubread.txt", header = T)
    GSE53078 <- create.geneColumn(GSE53078)
    trt.GSE53078 <- rep(61:62, each=2)
    
    
    GSE58082 <- read.table("GSE58082.Rsubread.txt")
    GSE58082 <- create.geneColumn(GSE58082)
    trt.GSE58082 <- rep(71:72, each =3)
    
    
    GSE57806 <- read.table("GSE57806.Rsubread.txt")
    GSE57806 <- create.geneColumn(GSE57806)
    trt.GSE57806 <- rep(81:82, each =3)
 
    ls <- list(GSE37159, GSE38879, GSE43865, GSE48767,
      GSE51119, GSE51772, GSE53078, GSE58082, GSE57806)

    stable <- join_all(ls, "Gene")
    stable <- stable[complete.cases(stable),]

    row.names(stable) <- stable[, 1]
    stable <- stable[, -1]
    
    
    trt <- as.factor(c(trt.GSE37159, trt.GSE38879, trt.GSE43865, 
             trt.GSE48767, trt.GSE51119, trt.GSE51772, 
             trt.GSE53078, trt.GSE58082, trt.GSE57806))
    
    n.obs <- c(length(trt.GSE37159), length(trt.GSE38879), length(trt.GSE43865),
               length(trt.GSE48767), length(trt.GSE51119), length(trt.GSE51772), 
               length(trt.GSE53078), length(trt.GSE58082), length(trt.GSE57806))
    
    group <- rep(1:length(ls), n.obs)
    
    
    index.id <- c()
    for ( i in 1: dim(stable)[1])
    {
      data2 <- as.matrix(stable)
      m <- tapply(data2[i, ], trt, sum)
      # get the sum for each gene by group
      # x <- c(1, 2, 3, 1, 3, 0, 0, 0,0)
      # trt <- rep(c(1, 2, 3), each=3)
      # tapply(x, trt, sum)  --->  6, 4, 0
      
      index.id[i] <- any(m==0)
    }
    

    stable <- stable[index.id==FALSE, ]
    
  }
  
  
  ### leaf   ---------------------------------------------------
  
  else if (type=="leaf")
  {
    
    GSE36626 <- read.table("GSE36626.Rsubread.txt")
    GSE36626 <- create.geneColumn(GSE36626)
    trt.GSE36626 <- c(rep(1:2, 2))  
    
    GSE39463 <- read.table("GSE39463.Rsubread.txt")
    GSE39463 <- GSE39463[, -c(3, 7, 9)]  # remove duplicated runs
    GSE39463 <- create.geneColumn(GSE39463)
    trt.GSE39463 <- c(rep(11:14, each= 3))
    
    GSE48235 <- read.table("GSE48235.Rsubread.txt", header = T) 
    GSE48235 <- create.geneColumn(GSE48235)
    trt.GSE48235 <- c(rep(21:23, each=2))
    
    GSE51304 <- read.table("GSE51304.Rsubread.txt", header=T)
    GSE51304 <- create.geneColumn(GSE51304)
    trt.GSE51304 <- c(rep(31:39, each=2))
    
    GSE54677 <- read.table("GSE54677.Rsubread.txt", header = T)
    GSE54677 <- create.geneColumn(GSE54677)
    trt.GSE54677 <- c(rep(41:50, 2))  ### be careful about the design structure
    
 
    ls <- list(GSE36626, GSE39463, GSE48235, GSE51304, GSE54677)
    
    stable <- join_all(ls, "Gene")
    stable <- stable[complete.cases(stable),]
    
    
  
    
    #as.numeric(colSums(stable[, -1]))
    row.names(stable) <- stable[, 1]
    stable <- stable[, -1]
  
    
    #  group and treatment
    trt <- as.factor(c(trt.GSE36626, trt.GSE39463, trt.GSE48235, trt.GSE51304, trt.GSE54677))
    n.obs <- c(length(trt.GSE36626), length(trt.GSE39463), length(trt.GSE48235),
               length(trt.GSE51304), length(trt.GSE54677))
    
    
    ## if any treatment have all 0 read counts, then this gene will be removed from the data.
    index.id <- c()
    for ( i in 1: dim(stable)[1])
    {
        data2 <- as.matrix(stable)
        m <- tapply(data2[i, ], trt, sum)
        index.id[i] <- any(m==0)
    }
    
    stable <- stable[index.id==FALSE, ]  
    group <- rep(1:length(ls), n.obs)
    
  }
  
  
  else if (type =="tissue")
  {
    #flower
    GSE35288 <- read.table("GSE35288.Rsubread.txt")
    GSE35288 <- create.geneColumn(GSE35288)
    
#     GSE40256 <- read.table("GSE40256.Rsubread.txt")
#     GSE40256 <- create.geneColumn(GSE40256)
    
    #leaves
    
#     GSE36626 <- read.table("GSE36626.Rsubread.txt")
#     GSE36626 <- create.geneColumn(GSE36626)
    
    GSE48235 <- read.table("GSE48235.Rsubread.txt")
    GSE48235 <- create.geneColumn(GSE48235)
    
#     GSE54677 <- read.table("GSE54677.Rsubread.txt")
#     GSE54677 <- create.geneColumn(GSE54677)
    
    # seed
    GSE53952 <- read.table("GSE53952.Rsubread.txt")
    GSE53952 <- create.geneColumn(GSE53952)
    
    # root
#     GSE44062 <- read.table("GSE44062.Rsubread.txt")
#     GSE44062 <- create.geneColumn(GSE44062)
    
    # carpel
     GSE56326 <- read.table("GSE56326.Rsubread.txt")
     GSE56326 <- create.geneColumn(GSE56326)
    
    # hypocoltyl
    GSE35408 <- read.table("GSE35408.Rsubread.txt")
    GSE35408 <- create.geneColumn(GSE35408)
    
    
    trt.GSE35288 <- c(rep(1, 3), rep(2, 3))
#    trt.GSE40256 <- rep(1:4, each=2)
#    trt.GSE36626 <- rep(1:2, 2)
    trt.GSE48235 <- rep(11:13, each=2)
#    trt.GSE54677 <- rep(1:10, 2)
#    trt.GSE53952 <- rep(21:23, each=3)
#    trt.GSE44062 <- rep(1:4, each=2)
    trt.GSE56326 <- c(rep(31, 2), rep(32:33, each=3))
    trt.GSE35408 <- rep(41:45, 2)
    
    
    ### tissue ---------------------------------------------------
    # special filtering
    GSE53952 <- GSE53952[, c(1, 2:4,11:13,20:22 )]  # the first stage was selected
    trt.GSE53952 <- rep(41:43, each=3)
    
    GSE35288 <- GSE35288[, c(1, seq(2,dim(GSE35288)[2], by=3))]
    

    ls <- list(
      GSE35288, #flower
      GSE48235, #leaves
      GSE53952, # seed
      GSE56326, # carpel
      GSE35408) # hypocoltyl
    
    stable <- join_all(ls, "Gene")
    stable <- stable[complete.cases(stable),]
    dim(stable)
    as.numeric(colSums(stable[, -1]))
    row.names(stable) <- stable[, 1]
    stable <- stable[, -1]
    

    trt <- as.factor(c(trt.GSE35288, trt.GSE48235, trt.GSE53952,
             trt.GSE56326, trt.GSE35408))
    n.obs <- c(length(trt.GSE35288), length(trt.GSE48235), 
               length(trt.GSE53952), length(trt.GSE56326), length(trt.GSE35408))


    index.id <- c()
    for ( i in 1: dim(stable)[1])
    {
        data2 <- as.matrix(stable)
        m <- tapply(data2[i, ], trt, sum)
        index.id[i] <- any(m==0)
    }
      

    stable <- stable[index.id==FALSE, ]

    group <- rep(1:length(ls), n.obs)
  }
  
  return (list(data= stable, trt=trt, group=group))
  
}






