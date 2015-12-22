#' Download SRA file from NCBI and convert it into FASTQ
#' 
#' 
#' SRA files can be downloaded automatically from NCBI, and then are converted into "FASTQ" file 
#' using SRA Toolkit (version 2.3.5-2)
#'
#' @note need to install SRA toolkit first. In the "RunInfo" file, choose \code{LibraryStrategy=="RNA-Seq"},
#' \code{LibrarySelection=="cDNA"}, and \code{RunInfo$LibraryLayout=="SINGLE"}.
#'
#' @title downloading SRA file and unzip it to FASTQ
#'
#' @param file  A RunInfo file (.csv format) generated from NCBI, which should be in current working directory.
#' @param SRA.path  where the SRA Toolkit "bin" folder locates
#' @param      the path where FASTQ file should be stored
#' @return the returned objects will be used in \code{FASTQtoCount} function
#' 



SRAtoFASTQ <- function(file, SRA.path)
  {
  RunInfo <- read.csv(file, stringsAsFactors=F)
  id <- which(RunInfo$LibraryStrategy=="RNA-Seq" 
              & RunInfo$LibrarySelection=="cDNA"
              &RunInfo$LibraryLayout=="SINGLE")
  RunInfo <- RunInfo[id, ]
  
  
  
  (fs <- basename(RunInfo$download_path))
  filepath <- RunInfo$download_path
  fsname <- RunInfo$SampleName
  fastq.name <- RunInfo$Run
  layout <- RunInfo$LibraryLayout
  
  for (i in 1:length(fs))
  {

    ## generate the terminal command line 
    fastq <- paste(SRA.path, "fastq-dump --split-3", sep="")
    # cmd <- paste(fastq, fileSRA)
    cmd <- paste(fastq, fs[i])
    cat(cmd, "\n")
    
    # call the command line code
    system(cmd)

  }
  
  # save this for alignment functions
  return(list(fastq.name=fastq.name, fsname =fsname))
  
  }