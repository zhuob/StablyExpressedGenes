
# step 1: download SRA file
# step 2: convert SRA to FASTQ
# step 3: align() by Rsubread
# step 4: featureCounts() to generate counts


library(Rsubread)


# this function allows me to process one sample a time, or all simultaneously
# input:  if all = True, then I'll process all the samples
#         otherwise, I'll just process the sample specified as "obs = " AND obs could be a vector of several samples


ReadRunInfo <- function(file, all = T, obs = 1 ){


	file1 <- paste("/home/stats/zhuob/data/alignment/NCBI/", file, "SraRunInfo.csv", sep="")   # where are the RunInfo File stored
	RunInfo <- read.csv(file1, stringsAsFactors=F)
	id <- which(RunInfo$LibraryStrategy=="RNA-Seq"
              & RunInfo$LibrarySelection=="cDNA"
                &RunInfo$LibraryLayout=="SINGLE"
                 &RunInfo$LibrarySource=="TRANSCRIPTOMIC")


	if (all == T){
	RunInfo <- RunInfo[id, ]
	}

	else {
	RunInfo <- RunInfo[id[obs], ]
	}

	
 	(fs <- basename(RunInfo$download_path))
 	filepath <- RunInfo$download_path
 	fsname <- RunInfo$SampleName
 	fastq.name <- RunInfo$Run
	layout <- RunInfo$LibraryLayout

	obj <- list(fs = fs, filepath= filepath, fsname = fsname, fastq.name = fastq.name)

	return(obj)
}




# input: returned object from ReadRunInfo()
#	  i: the ith sample to be processed
# output: a FASTQ file

## Note: need to set the the directory where FASTQ files locates as the current working directory
SRA2FASTQ <- function(obj, i){
	
#	setwd("/home/stats/zhuob/data/alignment/FASTQ/")
	SRA.path ="/home/stats/zhuob/sratoolkit.2.4.2-centos_linux64/bin/"

	filepath = obj$filepath	
	fs <- obj$fs
	fastq.name <- obj$fastq.name
 	fsname <- obj$fsname
  	datafile <- obj$filename	

	setwd("/home/stats/zhuob/data/alignment/FASTQ/")
	# download SRA files
	download.file(url= filepath[i], destfile = paste("/home/stats/zhuob/data/alignment/FASTQ/", fs[i], sep=""))
	## generate the terminal command line, Converting SRA to FASTQ
	fastq <- paste(SRA.path, "fastq-dump --split-3", sep="")
	cmd <- paste(fastq, fs[i])
	cat(cmd, "\n")
	
	# call the command line code
	system(cmd)

}



## Convert FASTQ to bam files  
# input: returned object from ReadRunInfo()
# output: bam file

FASTQ2BAM <- function(obj, i){
	
	fastq.name <- obj$fastq.name
 	fsname <- obj$fsname
 	datafile <- obj$filename

	# reference genome
  	ref.genome <- "Arabidopsis_thaliana.TAIR10.22.dna.toplevel.fa"

#	 build reference genome index IF YOU HAVEN'T DONE SO.  ONLY NEED TO DO ONCE!!!!
#	 in this function, the index file should be in the same directory as the FASTQ files 
# 	 buildindex(base="arab_index", reference=ref.genome)


	file1 <- paste(fastq.name[i], ".fastq", sep="")
    	out <- paste(fastq.name[i], ".bam", sep="")

	align(index="arab_index", readfile1=file1, input_format ="FASTQ",  output_file=out)

	# remove the fastq files and sra files
	file2 <- paste(fastq.name[i], ".sra", sep="")
	cmd <- paste("rm", file1, file2)
   	system(cmd)

}

## Convert FASTQ to bam files
# input: returned object from ReadRunInfo()
# output: summarized read counts

BAM2COUNTS <- function(obj){

        annotation <- "Arabidopsis_thaliana.TAIR10.22.gtf"       # annotation file

	fastq.name <- obj$fastq.name
	fsname <- obj$fsname

	bamfile <- paste(fastq.name, ".bam", sep="")

	# a list of samples can be processed simultaneously.
 	fc <- featureCounts(files=bamfile, annot.ext=annotation,
                      isGTFAnnotationFile=T)

 	cun <- fc$counts
  	cun <- data.frame(cun[order(row.names(cun)),])
  	colnames(cun) <- fsname

	return(cun)
}

