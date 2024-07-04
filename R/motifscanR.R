#' motifscanR: motifscan in R
#'
#' The motifscanR package is designed for analyzing many sequences and many
#' motifs to find which sequences contain which motifs.
#'
#'
#' @useDynLib motifscanR, .registration = TRUE
#' @docType package
#' @import Matrix SummarizedExperiment methods
#' @importFrom Rcpp sourceCpp
#' @importFrom parallel makeCluster stopCluster parLapply parSapply
#' @importFrom S4Vectors DataFrame
#' @importFrom Rsamtools getSeq
#' @importFrom GenomicFeatures transcripts
#' @importFrom GenomicRanges GRanges GRangesList promoters
#' @importFrom Biostrings getSeq letterFrequency DNAString DNAStringSet reverseComplement
#' @importFrom TFBSTools PWMatrixList PWMatrix name bg ID tags matrixClass
#' @importFrom IRanges IRanges IRangesList
#' @importFrom BSgenome getBSgenome as.list
#' @importFrom GenomeInfoDb genome seqlengths
#' @importMethodsFrom GenomicRanges seqnames start
#' @importMethodsFrom TFBSTools as.matrix
#' @importClassesFrom Biostrings DNAString DNAStringSet
#' @importClassesFrom GenomicRanges GenomicRanges
#' @importClassesFrom GenomicFeatures TxDb
#' @importClassesFrom TFBSTools PWMatrix PFMatrix PWMatrixList PFMatrixList
#' @author Louis Kwok
#' @name motifscanR
NULL



.onUnload <- function(libpath) {
  library.dynam.unload("motifscanR", libpath)
}

# Display startup message
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("motifscanR 0.9.3 2024-07-04")
}
