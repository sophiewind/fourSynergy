#' plotConsensusIa
#'
#' This function creates a karyotype plot displaying the interaction calls from
#' the consensus approach.
#'
#' @param ia fourSynergy object with interactions from all base tools
#' (peakC, r3c-seq, fourSig, r4cker) and other relevant information.
#' @param genes_of_interest Vector with genes of interest.
#' @param cex.chr character expansion of chromosome label.
#' @param cex.ideo character expansion base numbers of ideogram.
#' @param cex.y.track character expansion y axis track.
#' @param cex.vp character expansion viewpoint label.
#'
#' @return karyoplot with calling results.
#'
#' @examples
#' config <- system.file("extdata", "Datasets", "Demo", "info.yaml",
#'     package = "fourSynergy"
#' )
#' res_path <- system.file("extdata", "results", "Demo",
#'     package = "fourSynergy"
#' )
#' tracks <- system.file("extdata", "results", "Demo", "alignment",
#'     package = "fourSynergy"
#' )
#' sia <- createIa(res_path = res_path, config = config, tracks = tracks)
#' sia <- consensusIa(ia = sia, model = "AUPRC")
#' plotConsensusIa(ia = sia)
#' @export
plotConsensusIa <- function(ia = GRangesList(), genes_of_interest = NULL,
                            cex.chr = 1, cex.ideo = 0.6, cex.y.track = 0.6,
                            cex.vp = 1) {
    TxDb <- switch(ia@metadata$organism,
        "mm9" = TxDb.Mmusculus.UCSC.mm9.knownGene,
        "mm10" = TxDb.Mmusculus.UCSC.mm10.knownGene,
        "hg19" = TxDb.Hsapiens.UCSC.hg19.knownGene,
        "hg38" = TxDb.Hsapiens.UCSC.hg38.knownGene
    )

    if (is.null(genes_of_interest)) {
        kp <- createKaryoplot(ia, cex = cex.chr)
    } else {
        kp <- createKaryoplot(ia, type = "2", cex = cex.chr)
    }

    # Tracks
    bgs <- readBedGraph(ia)
    plotTracks(ia, kp, bgs,
        r0 = 0.5, r1 = 1, cex.vp = cex.vp,
        cex.y.track = cex.y.track
    )

    # Interactions
    if (length(ia@expConsensus) == 0) {
        warning(
            "fourSynergy found no interactions in condition. Did you run",
            " `consensusIa()`?"
        )
    } else if (length(ia@expConsensus[ia@expConsensus$significance > 0]) == 0) {
        warning("fourSynergy found no interactions in condition.")
    } else {
        kpPlotRegions(kp, ia@expConsensus[ia@expConsensus$significance > 0],
            col = "firebrick4", r0 = 0.45, r1 = 0.225
        )
    }

    if (length(ia@ctrlConsensus) == 0) {
        warning(
            "fourSynergy found no interactions in control. Did you run",
            "`consensusIa()`?"
        )
    } else if (length(ia@ctrlConsensus[ia@ctrlConsensus$significance > 0]) ==
        0) {
        warning("fourSynergy found no interactions in control.")
    } else {
        kpPlotRegions(kp, ia@ctrlConsensus[ia@ctrlConsensus$significance > 0],
            col = "darkblue", r0 = 0, r1 = 0.225
        )
    }

    # Genes
    if (!is.null(genes_of_interest)) {
        plot_genes(ia, kp, genes_of_interest, TxDb)
    }

    # VP
    kpAbline(kp, v = start(ia@vp), col = "black", lty = 2, data.panel = 1)
    kpAddBaseNumbers(kp,
        tick.dist = 100000, tick.len = 10, cex = cex.ideo,
        minor.tick.col = "gray"
    )

    legend(0.85, 0.8,
        legend = c("condition", "control"),
        fill = c("firebrick4", "darkblue"), border = NA,
        bty = "o", cex = 0.6
    )
}
