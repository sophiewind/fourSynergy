#' @title fourSynergy: Ensemble based interaction calling in 4C-seq data
"_PACKAGE"
#' @description fourSynergy is an R Bioconductor package for 4C-seq data
#' analysis.
#'
#' @details
#' fourSynergy is an ensemble algorithm for the analysis of 4C-seq data. It is
#' based on r3c-seq, fourSig, peakC and R.4C-ker.
#'
#' @section Main Functions:
#' \itemize{
#'   \item{\code{createIa()}: Reads in all relevant files from the fourSynergy
#'   pipeline and stores informaion in fourSynergy object.}
#'   \item{\code{consensusIa()}: Runs the weighted voting based
#'   ensemble algorithm to find interacting regions.}
#'   \item{\code{diffIa()}: Performs differential interaction analysis based on
#'   DEseq2.}
#' }
#'
#' @section Installation:
#' To install this package, use:
#' \code{BiocManager::install("fourSynergy")}
#'
#' @name fourSynergy
#' @aliases fourSynergy-package
NULL
