## selecting data sets


library(NBPSeq)

setwd(RawDataPath)


decideData <- function(datalist){

  num.labs <- length(datalist)
  data.list <- list()
  size <- c()
  
  
  ## repeated over time. chose samples at time 0.
  col.id.GSE43865 <- c("GSM1072470", "GSM1072471", "GSM1072472", "GSM1072473",
    "GSM1072474", "GSM1072475", "GSM1072476", "GSM1072477", "GSM1072478",
    "GSM1072479", "GSM1072480", "GSM1072481", "GSM1072482", "GSM1072483",
    "GSM1072484", "GSM1072485", "GSM1072486", "GSM1072487", "GSM1072488",
    "GSM1072489", "GSM1072490", "GSM1072491", "GSM1072492", "GSM1072493",
   "GSM1072494", "GSM1072495", "GSM1072496", "GSM1072497", "GSM1072498",
   "GSM1072499", "GSM1072500", "GSM1072501", "GSM1072502", "GSM1072503",
   "GSM1072504", "GSM1072505")

  col.id.GSE60835 <- paste("GSM14896", 38:43, sep ="")
  
  ## remove duplicated column of GSE35288
  col.id.GSE35288 <- c("GSM865212.1", "GSM865212.2", "GSM865213.1", "GSM865213.2", 
                       "GSM865214.1", "GSM865214.2", "GSM865215.1", "GSM865215.2",
                       "GSM865216.1", "GSM865216.2", "GSM865217.1", "GSM865217.2")
  
  ## choose samples at hour 6
  col.id.GSE39463 <- c("GSM969663.1", "GSM969664.1", "GSM969665.1", "GSM969667.1", "GSM969668.1",
                       "GSM969669",   "GSM969669.1", "GSM969669.2", "GSM969670",   "GSM969670.1",
                       "GSM969671",   "GSM969671.1", "GSM969672",   "GSM969672.1", "GSM969673",  
                       "GSM969673.1", "GSM969674",   "GSM969674.1", "GSM969675",   "GSM969675.1",
                       "GSM969676",   "GSM969677",   "GSM969678",   "GSM969679",   "GSM969679.1",
                       "GSM969680",   "GSM969681",   "GSM969682",   "GSM969683",   "GSM969684",  
                       "GSM969685",   "GSM969685.1", "GSM969686",   "GSM969686.1", "GSM969693",  
                       "GSM969694",   "GSM969695",   "GSM969696",   "GSM969696.1", "GSM969697" , 
                       "GSM969698",   "GSM969698.1", "GSM969699",   "GSM969700",   "GSM969700.1",
                       "GSM969701",   "GSM969701.1", "GSM969702",   "GSM969703",   "GSM969704" , 
                       "GSM969705",   "GSM969706",   "GSM969707",   "GSM969708",   "GSM969709",  
                       "GSM969710")  
  
  select.index <- c(col.id.GSE39463, col.id.GSE35288, col.id.GSE43865, col.id.GSE60835)
  
  for ( k  in 1: num.labs)
  {
    data.name <- paste(datalist[k], ".Rsubread.txt", sep="")
    
    data.txt <- read.table(data.name, header=T) 
   

    
 
    keep <- ! names(data.txt) %in% select.index
    
    rm.index <- c()
    
    data.txt <- data.txt[, keep]
    data.list$data[[k]] <- data.txt
    print(datalist[k])
    size[k] <- ncol(data.txt)
    
    data.list$GEOname[k] <- datalist[k]
  }
  
  df <- data.list$data[[1]]
  nrep <- c()
  nrep[1] <- dim(df)[2]
  for ( i in 2: (num.labs))
  {
    df <- cbind(df, data.list$data[[i]])
    nrep[i] <- dim(data.list$data[[i]])[2]
  }
  
  
  results <- list(count= df, GEO = datalist, size = size)
  return(results)
  
}    


rm.experiment <- function(df, datalist, norm.threshold = c(0.7, 1.3)){
  
  df$lab <- rep(df$GEO, df$size)
  ## check out the norm factors
  norm.factors <- estimate.norm.factors(df$count)
  nf <- data.frame(lab = df$lab, sample = names(norm.factors), norm.factors = norm.factors)

  rm.id <- which(nf$norm.factors < norm.threshold[1] | norm.factors > norm.threshold[2]) 
  rm.exp <- unique(nf$lab[rm.id])
  new.exp <- datalist[! datalist %in% rm.exp ]
  df <- decideData(new.exp)
  nf <- nf[nf$lab %in% new.exp, ]
  
  result <- list(df= df, nf = nf)
  return(result)
}







#########-------------------------------- Seedling data-------------------------


# GEO.seedling <- c("GSE30720", "GSE32202", "GSE32216", "GSE37159", "GSE38286", 
#                   "GSE38400", "GSE38879", "GSE39214", "GSE41432", "GSE41776",
#                   "GSE43703", "GSE43865", "GSE48767", "GSE50597",
#                   "GSE51119", "GSE51772", "GSE53078", "GSE55482",
#                   "GSE57806", "GSE58082", "GSE58662", "GSE58974", "GSE60835",
#                   "GSE64870", "GSE66666")

# 0. GSE30720 	two biological replicates for col-0
# col.GSE30720 <- c("GSM762078", "GSM762079")


Col.ecotype <- sort(c( "GSE32202",# "GSE32216", 
                       "GSE37159", # "GSE38286", 
                  # "GSE38400",
                  "GSE39214", # "GSE41432",
                  "GSE41766", # "GSE43703" ,
                  "GSE43865", # "GSE50597",  
                  "GSE51119", "GSE51772", 
                  "GSE53078", #"GSE58974",
                  "GSE60835", 
                  #"GSE64870",
                  "GSE66666"))

# 1. GSE43865 has 42 samples, we only extracted the first 6
# 2.  discard GSE58974, because the mapping quantity is too low (less than 20%)
# 3. the ecotypes of GSE38879 , GSE39214, GSE57806, GSE58082, 
#    GSE58662, are not Col 0
#   
# 4. GSE38286 to be examined 
#
# 5. If GSE38879 is chosen, need to filter the duplicate runs

# 6. to be checked GSE32216, GSE40256, GSE43703, GSE50597 
# the map quantity  >=59.3%           23% ~ 68%, 31% ~ 58%

# 7. GSE60835: GSM1489632	- GSM1489637 c(1, 1, 1, 0, 0, 0)
# 8. GSE38879 not  Col 0
# 9 
# GSE38286 NOT YET processed
# GSE41432 NOT YET PROCESSED
# GSE43703 high counts  (removed)
# GSE32216  high counts (need to decide)
# GSE50597  (need to decide)
# GSE58974  low counts (removed)
# GSE64870  not Columbia
# GSE32216  only 2 columbia


cat("selecting seedling experiments\n")
df <- decideData(Col.ecotype)

threshold <- c(0.5,1.5)
rm1 <- rm.experiment(df, Col.ecotype, threshold)  
library(ggplot2)

ggplot(data= rm1$nf, aes(x=lab, y= norm.factors, group = lab, col = as.factor(lab), shape= as.factor(lab))) + 
  geom_point() + scale_shape_manual(values = 1:length(rm1$df$GEO)) + 
  theme(axis.text.x = element_text(angle = 45, size =10))

# dts <-  as.vector(unique(rm1$nf$lab))
# rm2 <- rm.experiment(rm1$df, dts, threshold)  
# 
# ggplot(data= rm2$nf, aes(x=lab, y= norm.factors, group = lab, col = as.factor(lab), shape= as.factor(lab))) + 
#   geom_point() + scale_shape_manual(values = 1:length(dts)) + 
#   theme(axis.text.x = element_text(angle = 45, size =10))


seedling <- rm1$df

# remove GSE32216                  COL: 7-10
# remove GSE39214                  COL: 31-42
# remove Col 3-4, 7-8 of GSE43703  COL:  51 52 55 56
#@ remove GSE43703                 COL: 49:56
# remove GSE50597                  COL: 63-70
# remove GSE58974                  COL: 93:102
# remove col 7-12 of GSE60835 NOT COLUMBIA 0      COL: 109:114
# remove GSE64870                  COL: 115:136
# remove GSE66666                  COL: 137:142
# 
# remove GSE38400                  COL: 19:30    (low quality)


# treatment Structure 
col.trt32202 <- rep(1:3, each=2)
col.trt37159 <- c(rep(c(1:4), each=2))
#@ col.trt38400 <- rep(c(1:4), each=3)
col.trt41766 <- c(rep(c(1:3), each=2))
# col.trt43703 <- rep(1:2, 2)
col.trt43865 <- rep(c(1:2), each=3)
col.trt51119 <- rep(c(1:5), each=2)
col.trt51772 <- rep(c(1:4), each=2)
col.trt53078 <- rep(c(1:2), each=2)
col.trt60835 <- rep(c(1:2), each=3)
col.trt66666 <- rep(c(1, 2), each= 3)
trt <- c(col.trt32202, col.trt37159, # col.trt38400,
         #@ trt <- c(col.trt32202, col.trt37159, # col.trt38400,
         col.trt41766,
         # col.trt43703, 
         col.trt43865, col.trt51119, 
         col.trt51772, col.trt53078, col.trt60835, 
         col.trt66666)

lab <- rep(1:length(seedling$size), seedling$size)

seedling$lab <- lab; seedling$trt <- trt


saveRDS(seedling, paste(DataPath, "seedling.rds", sep = ""))



#####------------------------ LEAF DATA --------------------------------------------------

leaf_Col0 <- c("GSE36626", "GSE39463", "GSE43983", #"GSE45989", 
               "GSE48235", "GSE51304", "GSE54677", "GSE55884", 
               "GSE58029")
#,"GSE67777", "GSE67956") 

# GSE36626   c(1, 0, 1, 0)
# GSE39463  c(1:4, each=3)
# GSE43983  No replicate   CONSIDER rep(1:4, 2)
# GSE48235: rep(1:3, each=2)
# GSE45989: use GSM1121350 -GSM1121365	
# GSE51304: rep(1:9, each=2)
# GSE54677: rep(1:10, each=2)
# GSE55884: use GSM1347943- GSM1347948 rep(1:2, each=3)
# GSE58029: no replicate
# GSE67777: rosettle leaves: RNA-seq data at different stages of leaf aging. (not enough samples)
# GSE67956: rosettle leaves  GSM1659609	-GSM1659644  (0h, D1, D2, D3, D6h, L3) X (hml, ivd, wt) X 2 (not enought samples)
cat("selecting leaf experiments\n")
df <- decideData(leaf_Col0)


threshold <- c(0.5,1.5)
rm1 <- rm.experiment(df, leaf_Col0, threshold)  
library(ggplot2)

ggplot(data= rm1$nf, aes(x=lab, y= norm.factors, group = lab, col = as.factor(lab), shape= as.factor(lab))) + 
  geom_point() + scale_shape_manual(values = 1:length(rm1$df$GEO)) + 
  theme(axis.text.x = element_text(angle = 45, size =10))


leaf <- rm1$df


trt.GSE36626 <- rep(1:2, 2)
trt.GSE39463 <- rep(1:4, each=3)
# trt.GSE45989 <- rep(1:4, each=4)
trt.GSE48235 <- rep(1:3, each=2)
trt.GSE51304 <- rep(1:9, each=2)
trt.GSE54677 <- rep(1:10, 2)  ### be careful about the design structure
trt.GSE55884 <- rep(1:2, each=3)
trt.GSE67777 <- rep(1:4, each=3)
trt.GSE67956 <- c(1, 1, 2, 2, 3, 3)

trt <- c(trt.GSE36626, trt.GSE39463, trt.GSE48235, trt.GSE51304, 
         trt.GSE54677)
lab <- rep(1:5, c(4, 12, 6, 18, 20))

leaf$trt <- trt
leaf$lab <- lab

#saveRDS(leaf, "/Users/Bin/Dropbox/Zhuo/Research/Project2014/StablyExpressedGenes/R/data/leaf.rds")
saveRDS(leaf, paste(DataPath, "leaf.rds", sep = ""))



#####-------------------  MULTIPLE TISSUE DATA ---------------------
# aerial tissue:  GSE62799: GSM1533521- GSM1533526 (1, 1, 1, 0, 0, 0)
#                 GSE58856: Not sure about ecotype
# carpel:  GSE56326  GSM1359146- GSM1359153 (1,1, 2, 2, 2, 3, 3, 3)
# epidermis: GSE60183 GSM1466934-GSM1466939 (1,1, 1, 0, 0, 0)
# FLOWER: GSE35288 (COLUMBIA), GSM865212-GSM865217 (1,1 1, 0, 0,0)
#         GSE40256 GSM989339-GSM989346 (1, 1, 2, 2, 3, 3, 4, 4)  USE COL 1-4
#         GSE51069(COLUMBIA) GSM1237330-GSM1237333  (1, 1, 0, 0)
# unopened flower bud:  GSE57215, col GSM1377347-GSM1377369, USE RNA-Seq GSM1377353	- GSM1377358 (1, 1, 1, 0, 0, 0)
# hypocotyl:   GSE35408  GSM867674－ GSM867678， GSM951964－GSM951968 （1，2， 3， 4， 5， 1， 2， 3， 4， 5）
# inflorescences & siliques: GSE59637 GSM1441333-GSM1441336 (1, 1, 2, 2)
# leaves:  see above USE 48235
# root:   GSE52966  GSM1279534-GSM1279551  rep(1:9, 2)
#         GSE44062 GSM1077432-GSM1077439  rep(1:4, each =2)
#         GSE64381 (Maybe low quanlity)
#         GSE64410: May choose GSM1576088-GSM1576135(0 hour: col 1,2, 3 hour: 3-4, 8 hour 5-6, 8 treatments)
# root tip: GSE59167 GSM1429679-GSM1429693  c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4)
# roots, leaves, flowers, and siliques:  GSE43983 WT X clf28 No replicate
# seed:   GSE42957  Hybrid
#         GSE53952 GSM1303953-GSM1303979  7-8Hour: 1, 2, 3, 9-10Hour: 4, 5,6, 11-12Hour: 7, 8, 9. TRT: fae1 CL37 PDAT           
#         GSE61061 GSM1496098-GSM1496103  c(1, 1, 1, 0, 0, 0)

#  shoot apical meristem    GSE63355  GSM1546538-GSM1546553 rep(1:8, each=2)

# Tissue_col0 <- sort(c("GSE62799", #"GSE58856",             #aerial tissue
#                  "GSE56326",                          # carpel
#                  "GSE60183",                          # epidermis
#                  "GSE35288", # "GSE40256", "GSE51069",# flower 
#                  # "GSE57215",                        # flower bud
#                  "GSE35408",                          # hypocotyl
#                  "GSE59637",                          # infuoresence
#                  "GSE52966", "GSE44062", "GSE64381", "GSE64410",  # root
#                  "GSE59167", #"GSE42957","GSE53952",  # root tip
#                #   "GSE51304",                          # leaf
#                  "GSE61061",                          # seed
#                  "GSE63355"))                        # shoot apical
# 
# # GES53952: maybe NOT COLUMBIA TYPE
# # GSE57215: removed low quantity  < 50%
# # GSE51069: removed low quantity
# # GSE64381: removed low quantity
# # GSE60183: low quantity

# Tissue_col0 <- sort(Tissue_col0)


Tissue_col1 <- sort(c("GSE62799",              #aerial tissue
                 "GSE56326",              # carpel
                 "GSE35288",              # flower 
                 "GSE60183",              # epidermis
                 # epidermis
                 # flower bud
                 "GSE35408",              # hypocotyl
                 "GSE59637",              # infuoresence
                 "GSE52966",              # root
                 "GSE59167",              # root tip
                 #"GSE51304",                          # leaf
                 "GSE61061",              # seed
                 "GSE63355") )             # shoot apical



cat("selecting multi-tissue experiments\n")
df <- decideData(Tissue_col1)

threshold <- c(0.5, 1.5)
rm1 <- rm.experiment(df, Tissue_col1, threshold)  
library(ggplot2)

ggplot(data= rm1$nf, aes(x=lab, y= norm.factors, group = lab, col = as.factor(lab), shape= as.factor(lab))) + 
  geom_point() + scale_shape_manual(values = 1:length(rm1$df$GEO)) + 
  theme(axis.text.x = element_text(angle = 45, size =10))

tissue <- rm1$df

trt.35288 <- rep(1:2, each=3)
trt.35408 <- rep(1:5, 2)
trt.52966 <- rep(1:9, 2)
trt.56326 <- c(1,1, 2, 2, 2, 3, 3, 3)
trt.59167 <- c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4)
trt.59637 <- c(1, 1, 2, 2)
trt.60183 <- rep(1:2, each=3)
trt.61061 <- rep(1:2, each=3)
trt.62799 <- rep(1:2, each=3)
trt.63355 <- rep(1:8, each=2)

trt <- c(trt.35288, trt.35408, trt.52966, trt.56326, trt.59167, 
         trt.59637, trt.60183, trt.61061, trt.62799, trt.63355)

lab <- rep(1:length(tissue$size), tissue$size)


tissue$trt <- trt; tissue$lab <- lab
ids <- tissue$lab == 5
tissue$count[1:4, ids]
saveRDS(tissue, paste(DataPath, "tissue.rds", sep = ""))






