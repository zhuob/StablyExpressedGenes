#' Process FASTQ file to RNA-Seq count table
#' 
#' 
#'
#' @note need to install Rsubread first. In the "RunInfo" file, choose \code{LibraryStrategy=="RNA-Seq"},
#' \code{LibrarySelection=="cDNA"}, and \code{RunInfo$LibraryLayout=="SINGLE"}.
#'
#' @title RNA-Seq alignment to obtain read count table
#'
#' @param obj  An object returned by function \code{SRAtoFASTQ}
#' @param ref.genome  the reference genome to which the FASTQ file should be aligned
#' @param annotation  annotation file of this species
#' @return a data frame 
#' \item{count} an n-by-p data frame, where n is the number of genes and is the number of samples
#'

FASTQtoCount <- function(obj, ref.genome, annotation)
{
  library(Rsubread)
  fastq.name <- obj$fastq.name
  fsname <- obj$fsname
  
  # build reference genome index
  buildindex(base="build_index", reference=ref.genome)
  
  # do the alignments by converting .fastq file into .bam file
  for (i in fastq.name)
  {
    file1 <- paste(i, ".fastq", sep="")
    out <- paste(i, ".bam", sep="")
    
    align(index="build_index", readfile1=file1, input_format ="FASTQ",  output_file=out)
  
    # remove the fastq files and sra files
    file2 <- paste(i, ".sra", sep="")
    
   cmd <- paste("rm", file1, file2)
   system(cmd)
  }
  
  
  bamfile <- paste(fastq.name, ".bam", sep="")
  
  ## The options for featureCounts are discussed here. see #78
  # http://seqanswers.com/forums/showthread.php?t=30258&page=4-
  
  # a list of samples can be processed simultaneously.
  fc <- featureCounts(files=bamfile, annot.ext=annotation,
                      isGTFAnnotationFile=T)
  
  cun <- fc$counts
  cun <- data.frame(cun[order(row.names(cun)),])
  colnames(cun) <- fsname
  
  # system("rm  ../*.bam")
    
  return(cun)
}



