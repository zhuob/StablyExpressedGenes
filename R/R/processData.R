source("/home/stats/zhuob/data/alignment/NCBI/SRAtoCount.R")

setwd("/home/stats/zhuob/data/alignment/FASTQ/")   # where the reference index and FASTQ files are stored

file <- "GSE67956"


## Read the RunInfo file 
## all: if True, process all the samples in an experiment
#	if False, process the given samples specified by obs 


  information <- ReadRunInfo(file, all = T)
# information <- ReadRunInfo(file, all = F, obs = c(1, 2))


for (i in 1: length(information$fs) ){

	SRA2FASTQ(information, i)
	FASTQ2BAM(information, i)
}


count <- BAM2COUNTS(information)


dest.path <- "/home/stats/zhuob/data/alignment/Count/"

dest.file <- paste(dest.path, file, ".Rsubread.txt", sep="")
write.table(count, dest.file)





##########
##########  Old Version ########################

# source("/home/stats/zhuob/data/alignment/NCBI/SRAtoCount.R")

# file <- "GSE64381"


# filedownload <- SRAtoFASTQ(file)
# FASTQtoCount(filedownload)

