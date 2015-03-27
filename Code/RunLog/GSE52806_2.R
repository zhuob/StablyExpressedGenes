

library(Rsubread)
library(edgeR)
library(limma)

datafile <- "GSE52806"

setwd("/home/zhuob/protocol/Rsubread")
read.file <- paste(datafile, "SraRunInfo.csv", sep="" )
RunInfo <- read.csv(read.file, stringsAsFactors=F)[(11:20),]
(fs <- basename(RunInfo$download_path))
filepath <- RunInfo$download_path
fsname <- RunInfo$SampleName
fastq.name <- RunInfo$Run
layout <- RunInfo$LibraryLayout


for (i in 1:length(fs))
{
  path="/home/zhuob/protocol/software/sratoolkit.2.3.5-2-ubuntu64/bin/"
  fastq <- paste(path, "fastq-dump --split-3", sep="")
  
  fileSRA <- paste("/home/zhuob/protocol/Rsubread/", fs[i], sep="")
  cmd <- paste(fastq, fileSRA)
  cat(cmd, "\n")
  system(cmd)
  
  cmd2 <- paste("rm", fs[i])
  system(cmd2)
}

# buildindex(base="arab_index", reference="Arabidopsis_thaliana.TAIR10.22.dna.toplevel.fa")
for (i in fastq.name)
{
    file1 <- paste(i, ".fastq", sep="")
    out <- paste(i, ".bam", sep="")
    
    align(index="arab_index", readfile1=file1, input_format ="FASTQ",  output_file=out)
   remove.file <- paste("/home/zhuob/protocol/Rsubread/", file1, sep="")
    cmd <- paste("rm", file1)
   system(cmd)
}



bamfile <- paste(fastq.name, ".bam", sep="")

## The options for featureCounts are discussed here. see #78
# http://seqanswers.com/forums/showthread.php?t=30258&page=4-

fc <- featureCounts(files=bamfile, annot.ext="Arabidopsis_thaliana.TAIR10.22.gtf",
                    isGTFAnnotationFile=T)

cun <- fc$counts
cun <- data.frame(cun[order(row.names(cun)),])
colnames(cun) <- fsname

dest.path <- "/home/zhuob/Dropbox/Zhuo/Research/Project2014/Data/arab/Rsubread/"

dest.file <- paste(dest.path, datafile, ".Rsubread.txt", sep="")
write.table(cun, dest.file)





