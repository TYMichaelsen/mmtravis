% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_loadMetaT.R
\name{mt_loadMetaT}
\alias{mt_loadMetaT}
\title{Load data produced by the \href{https://github.com/TYMichaelsen/MetaT}{MetaT} pipeline into R.}
\usage{
mt_loadMetaT(counts.txt,seqstat.txt,mtmeta = NULL)
}
\arguments{
\item{counts.txt}{(\emph{required}) The count table as outputted from the \href{https://github.com/TYMichaelsen/MetaT}{MetaT} pipeline}

\item{seqstat.txt}{(\emph{required}) The seqstat table as outputted from the \href{https://github.com/TYMichaelsen/MetaT}{MetaT} pipeline}

\item{mtmeta}{(\emph{optional}) data.frame with metadata associated to samples (columns). See \code{\link{mt_load}} for details.}
}
\value{
A \code{mmt} object.
}
\description{
mt_loadMetaT
}
\details{
This function is a wrapper to transform .txt files outputted from the \href{https://github.com/TYMichaelsen/MetaT}{MetaT} pipeline into R.
}
\author{
Thomas Yssing Michaelsen \email{tym@bio.aau.dk}
}
