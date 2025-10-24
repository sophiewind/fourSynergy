#' Set differential attribute
#'
#' @param object fourSynergy object with interactions from all base tools.
#' @param value DESeqResults-object.
#'
#' @return Modified fourSynergy-object
#' @export
setGeneric("setDifferential", function(object, value) {
    standardGeneric("setDifferential")
})

#' Set differential attribute
#'
#' @param object fourSynergy object with interactions from all base tools
#' @param value DESeqResults-object
#'
#' @return Modified fourSynergy-object
#' @export
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
