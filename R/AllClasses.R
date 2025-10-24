#' @title fourSynergy Class
#' @description S4 class storing data collected from 4C-seq analyses.
#' @slot metadata Experimental metadata from config file.
#' @slot expInteractions Base tool interactions found in the experiment.
#' @slot ctrlInteractions Base tool interactions found in the control.
#' @slot expConsensus Consensus interactions found in the experiment.
#' @slot ctrlConsensus Consensus interactions found in the control.
#' @slot vp Viewpoint position.
#' @slot vfl Virtual fragment library.
#' @slot tracks Path to the alignment files.
#' @slot differential Results of differential interaction calling.
#'
#' @rdname ia4C-class
#' @export
setClass(
    "fourSynergy",
    representation = list(
        metadata = "list",
        expInteractions = "GRangesList",
        ctrlInteractions = "GRangesList",
        expConsensus = "GRanges",
        ctrlConsensus = "GRanges",
        vp = "GRanges",
        vfl = "GRanges",
        tracks = "character",
        differential = "ANY"
    ),
    prototype = prototype(
        metadata = list(),
        expInteractions = GRangesList(),
        ctrlInteractions = GRangesList(),
        expConsensus = GRanges(),
        ctrlConsensus = GRanges(),
        vp = GRanges(),
        vfl = GRanges(),
        tracks = character(0),
        differential = NULL
    ),
    validity = function(object) {
        # Check if all slots are of correct class
        if (!(is.list(object@metadata))) {
            return("Metadata must be a list")
        }
        if (!(inherits(object@expInteractions, "GRangesList"))) {
            return("Exp interactions must be a GRangesList")
        }
        if (!(inherits(object@ctrlInteractions, "GRangesList"))) {
            return("Ctrl interactions must be a GRangesList")
        }
        if (!(inherits(object@expConsensus, "GRanges"))) {
            return("Exp consensus must be a GRanges")
        }
        if (!(inherits(object@ctrlConsensus, "GRanges"))) {
            return("Ctrl consensus must be a GRanges")
        }
        if (!(inherits(object@vp, "GRanges"))) {
            return("VP must be a GRanges")
        }
        if (!(inherits(object@vfl, "GRanges"))) {
            return("VFL must be a GRanges")
        }
        if (!(is.character(object@tracks))) {
            return("Tracks must be a character vector")
        }
        return(TRUE)
    }
)
