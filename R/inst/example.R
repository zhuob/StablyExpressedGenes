##  info need to be prepared
#  create a folder (500G to be safe) and put all the following material under this folder
# 1 SRA Toolkit installed 
# 2 Runinfo.csv
# 3 annotation file  .gtf
# 4 compressed reference genome .fa


source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/R/R/SRAtoFASTQ.R")
source("/Users/Bin/Dropbox/Zhuo/Research/Project2014/R/R/FASTQtoCount.R")
source("/home/zhuob/Dropbox/Zhuo/Research/Project2014/R/R/SRAtoFASTQ.R")
source("/home/zhuob/Dropbox/Zhuo/Research/Project2014/R/R/FASTQtoCount.R")

## all 1-4 files must be in the current working directory.
setwd("/Users/Bin/Disk1/protocol/Rsubread/test/")
setwd("/home/zhuob/protocol/Rsubread/test/")

path <- "/Users/Bin/Disk1/protocol/Rsubread/test/sratoolkit.2.3.5-2-mac64/bin/"
path="/home/zhuob/protocol/Rsubread/test/sratoolkit.2.3.5-2-ubuntu64/bin/"

file <- "GSE32216SraRunInfo.csv"

GSE32216 <- SRAtoFASTQ(file, path)



reference.genome <- "Arabidopsis_thaliana.TAIR10.22.dna.toplevel.fa"
annotation <- "Arabidopsis_thaliana.TAIR10.22.gtf"

coun <- FASTQtoCount(GSE32216, reference.genome, annotation)



