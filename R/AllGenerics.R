#' Set differential attribute
#'
#' @param object fourSynergy object with interactions from all base tools.
#' @param value DESeqResults-object.
#'
#' @return Modified fourSynergy-object
#' @keywords internal
setGeneric("setDifferential", function(object, value) {
    standardGeneric("setDifferential")
})

#' @rdname setDifferential
setMethod(
    "setDifferential", signature(object = "fourSynergy"),
    function(object, value) {
        if (inherits(value, "DESeqResults")) {
            object@differential <- value
        } else {
            stop("Value muss ein DESeqResults-Objekt sein")
        }
        return(object)
    }
)

#' Set dds attribute
#'
#' @param object fourSynergy object with interactions from all base tools.
#' @param value DESeqDataSet-object.
#'
#' @return Modified fourSynergy-object
#' @keywords internal
setGeneric("setDds", function(object, value) {
    standardGeneric("setDds")
})

#' @rdname setDds
setMethod(
    "setDds", signature(object = "fourSynergy"),
    function(object, value) {
        if (inherits(value, "DESeqDataSet")) {
            object@dds <- value
        } else {
            stop("Value must be DESeqDataSet object.")
        }
        return(object)
    }
)


setGeneric("plotDiffIa", function(ia) {
    standardGeneric("plotDiffIa")
})

setGeneric("consensusIa", function(ia, model = "F1") {
    standardGeneric("consensusIa")
})

setGeneric("createIa", function(res_path = character(), config = list(),
                                tracks = "") {
    standardGeneric("createIa")
})

setGeneric("differentialAnalysis", function(ia, fitType = "local") {
    standardGeneric("differentialAnalysis")
})

setGeneric("plotConsenusIa", function(ia = GRangesList(), genes_of_interest =
                                        NULL) {
    standardGeneric("plotConsenusIa")
})

setGeneric("plotIaIndiviualTools", function(ia = GRangesList(),
                                            genes_of_interest = NULL) {
    standardGeneric("plotIaIndiviualTools")
})

setGeneric("plotBaseTracks", function(sia = GRangesList(),
                                    highlight_regions = NULL) {
    standardGeneric("plotBaseTracks")
})

setGeneric("plotConsensusTracks", function(sia = GRangesList(),
                                        highlight_regions = NULL) {
    standardGeneric("plotConsensusTracks")
})

setGeneric("checkConfig", function(config = character()) {
    standardGeneric("checkConfig")
})

