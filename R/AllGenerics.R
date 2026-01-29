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

