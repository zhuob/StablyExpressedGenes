% Generated by roxygen2 (4.0.2): do not edit by hand
\name{FASTQtoCount}
\alias{FASTQtoCount}
\title{RNA-Seq alignment to obtain read count table}
\usage{
FASTQtoCount(obj, ref.genome, annotation)
}
\arguments{
\item{obj}{An object returned by function \code{SRAtoFASTQ}}

\item{ref.genome}{the reference genome to which the FASTQ file should be aligned}

\item{annotation}{annotation file of this species}
}
\value{
a data frame
\item{count} an n-by-p data frame, where n is the number of genes and is the number of samples
}
\description{
Process FASTQ file to RNA-Seq count table
}
\note{
need to install Rsubread first. In the "RunInfo" file, choose \code{LibraryStrategy=="RNA-Seq"},
\code{LibrarySelection=="cDNA"}, and \code{RunInfo$LibraryLayout=="SINGLE"}.
}

