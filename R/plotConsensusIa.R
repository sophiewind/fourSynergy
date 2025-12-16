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
#' @param cex.y.lab character expansion for y labels.
#' @param cex.y.track character expansion y axis track.
#' @param cex.vp character expansion viewpoint label.
#' @param cex.leg character expansion for legend.
#' @param highlight_regions regions to highlight in the plot.
#' @param plot_spider plotting connections from VP to interactions.
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
                            cex.chr = 1, cex.ideo = 0.6, cex.y.lab = 0.6,
                            cex.y.track = 0.6, cex.vp = 1, cex.leg = 0.6,
                            highlight_regions = NULL, plot_spider = FALSE) {
    TxDb <- switch(ia@metadata$organism,
                "mm10" = TxDb.Mmusculus.UCSC.mm10.knownGene,
                "hg19" = TxDb.Hsapiens.UCSC.hg19.knownGene)

    if (is.null(genes_of_interest)) {
        kp <- createKaryoplot(ia, cex = cex.chr)
    } else {
        kp <- createKaryoplot(ia, type = "2", cex = cex.chr)
    }

    # Regions
    if (!is.null(highlight_regions)){
        plotRegions(ia, kp, highlight_regions)
    }

    # Tracks
    bgs <- readBedGraph(ia)
    plotTracks(ia, kp, bgs,
            r0 = 0.5, r1 = 1, cex.vp = cex.vp, cex.y.track = cex.y.track)

    if (plot_spider == TRUE){
        r_height <- c(0.4, 0.225, 0.175)
        kpAbline(kp, chr = paste0("chr", ia@metadata$VPchr), h = r_height[2])

        start.reg <- GRanges(seqnames = seqnames(ia@vp),
            ranges = start(ia@vp))
        start.regs <- rep(start.reg,
                        length(ia@expConsensus[
                            ia@expConsensus$significance > 0]))
        kpPlotLinks(kp, start.regs,
                ia@expConsensus[ia@expConsensus$significance > 0],
                col = "firebrick4",r0 = r_height[1], r1 = 0.45)
        start.regs_ctrl <-
            rep(start.reg,
                length(ia@ctrlConsensus[ia@ctrlConsensus$significance > 0]))
        kpPlotLinks(kp, start.regs_ctrl,
                    ia@ctrlConsensus[ia@ctrlConsensus$significance > 0],
                    col = "darkblue",r0 = r_height[3], r1 = r_height[2])
    } else {
        r_height <- c(0.45, 0.225, 0.225)
    }

    # Structure plot
    kpAbline(kp, chr = paste0("chr", ia@metadata$VPchr), h = r_height[1])
    kpAbline(kp, chr = paste0("chr", ia@metadata$VPchr), h = r_height[3])

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
                    col = "firebrick4", r0 = r_height[2], r1 = r_height[1])
        kpAddLabels(kp,
                    labels = "condition", r0 = r_height[2], r1 = r_height[1],
                    cex = cex.y.lab, data.panel = 1)
    }

    if (length(ia@ctrlConsensus) == 0) {
        warning(
            "fourSynergy found no interactions in control. Did you run",
            "`consensusIa()`?"
        )
    } else if (length(ia@ctrlConsensus[ia@ctrlConsensus$significance > 0]) ==
            0){
        warning("fourSynergy found no interactions in control.")
    } else {
        kpPlotRegions(kp, ia@ctrlConsensus[ia@ctrlConsensus$significance > 0],
                    col = "darkblue", r0 = 0, r1 = r_height[3])
        kpAddLabels(kp, labels = "control", r0 = 0, r1 = r_height[3],
                    cex = cex.y.lab, data.panel = 1)
    }

    # Genes
    if (!is.null(genes_of_interest)) {
        plot_genes(ia, kp, genes_of_interest, TxDb)
    }

    # VP
    kpAbline(kp, v = start(ia@vp), col = "black", lty = 2, data.panel = 1)
    kpAddBaseNumbers(kp, tick.dist = 100000, tick.len = 10, cex = cex.ideo,
                    minor.tick.col = "gray")

    # Legend
    has_control <-!is.null(ia@metadata$control)
    lg <- if (has_control) c("condition", "control") else "condition"
    col <- if (has_control) c("firebrick4", "darkblue") else "firebrick4"
    legend(0.85, 0.8, legend = lg, fill = col, border = NA, bty = "o",
        cex = cex.leg)
}
