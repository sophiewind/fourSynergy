#' plotDiffIa
#'
#' This function creates a karyoplot with the differential interactions calls.
#'
#' @param ia fourSynergy object with interactions from all base tools
#' (peakC, r3c-seq, fourSig, r4cker) and more relevant information.
#' @param genes_of_interest Vector with genes of interest. Set to `all` if you
#' want to plot all genes in this area.
#' @param cex.chr character expansion of chromosome label.
#' @param cex.ideo character expansion base numbers of ideogram.
#' @param cex.y.track character expansion y axis track.
#' @param cex.vp character expansion viewpoint label.
#' @param plot_spider plotting connections from VP to interactions
#' @param cex.leg character expansion for legend.
#' @param highlight_regions regions to highlight in the plot
#' @param cex.y.lab character expansion for y labels.
#' @param gene.name.cex character expansion for gene names.
#'
#' @return DESeq2 results of differential interaction calling.
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
#' sia <- differentialAnalysis(ia = sia, fitType = "mean")
#' plotDiffIa(ia = sia)
#' @export
plotDiffIa <- function(ia, genes_of_interest = NULL, cex.chr = 1,
                    cex.y.lab = 0.6, cex.ideo = 0.6, cex.y.track = 0.6,
                    cex.vp = 1, cex.leg = 0.6,
                    plot_spider = FALSE, highlight_regions = NULL,
                    gene.name.cex = 1) {
    # Check input
    if (is.null(ia@differential)) {
        stop("Empty @differential slot. Did you run `differentialAnalysis()`?")
    }

    # Remove NAs and transform to GRanges
    res <- ia@differential
    res.gr <- res %>%
        as.data.frame() %>%
        na.omit() %>%
        rownames_to_column("seqnames") %>%
        tidyr::separate(., seqnames,
                        sep = "[:-]",
                        c("seqnames", "start", "end")
        ) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    if (nrow(ia@differential) == 0 || length(res.gr[res.gr$padj < 0.05]) == 0) {
        warning("No significant differential interaction were found.")
    }

    # Tracks and baseplot
    bgs <- readBedGraph(ia)
    if (!is.null(genes_of_interest)) {
        panel <- "2"
    } else {
        panel <- "1"
    }

    # Plot highlighted regions
    kp <- createKaryoplot(ia, type = panel, cex = cex.chr)
    if (!is.null(highlight_regions)){
        plotRegions(ia, kp, highlight_regions)
    }

    # Plot interactions
    kpPlotRegions(kp, res.gr[res.gr$padj < .05], col = "grey", r0 = 0.4,
            r1 = 1)
    plotTracks(ia, kp, bgs, r0 = 0.4, r1 = 1, cex.vp = cex.vp,
            cex.y.track = cex.y.track)

    if (plot_spider == TRUE){
        r_height <- c(0.35, 0.275, 0.235, 0.16, 0.15, 0.075)
        kpAbline(kp, chr = paste0("chr", ia@metadata$VPchr), h = r_height[3])

        start.reg <- GRanges(
            seqnames = seqnames(ia@vp),
            ranges = start(ia@vp)
        )
        start.regs <-
            rep(start.reg,
                length(ia@expConsensus[ia@expConsensus$significance > 0]))
        kpPlotLinks(kp, start.regs,
                    ia@expConsensus[ia@expConsensus$significance > 0],
                    col = "firebrick4",r0 = 0.35, r1 = 0.39)
        start.regs_ctrl <-
            rep(start.reg,
                length(ia@ctrlConsensus[ia@ctrlConsensus$significance > 0]))
        kpPlotLinks(kp, start.regs_ctrl,
                    ia@ctrlConsensus[ia@ctrlConsensus$significance > 0],
                    col = "darkblue",r0 = 0.235, r1 = 0.275)
    } else {
        r_height <- c(0.39, 0.293, 0.293, 0.196, 0.194, 0.097)
    }
    # Structure plot
    kpAbline(kp, chr = paste0("chr", ia@metadata$VPchr), h = r_height[1])
    kpAbline(kp, chr = paste0("chr", ia@metadata$VPchr), h = r_height[2])
    kpAbline(kp, chr = paste0("chr", ia@metadata$VPchr), h = r_height[4])
    kpAbline(kp, chr = paste0("chr", ia@metadata$VPchr), h = r_height[6])

    kpPlotRegions(kp, ia@expConsensus[ia@expConsensus$significance > 0],
                col = "firebrick4", r0 = r_height[2], r1 = r_height[1])
    kpAddLabels(kp,
                labels = "condition", r0 = r_height[2], r1 = r_height[1],
                cex = cex.y.lab, data.panel = 1)
    kpPlotRegions(kp, ia@ctrlConsensus[ia@ctrlConsensus$significance > 0],
                col = "darkblue", r0 = r_height[4], r1 = r_height[3])

    kpAddLabels(kp,
                labels = "control", r0 = r_height[4], r1 = r_height[3],
                cex = cex.y.lab, data.panel = 1)
    max_value <- max(abs(res.gr$log2FoldChange))
    max_panel_range <- 0.5
    scaled_y1 <- (res.gr$log2FoldChange / max_value)* max_panel_range

    # Plot L2FC
    kpBars(kp,
        chr = paste0("chr", ia@metadata$VPchr), x0 = start(res.gr),
        x1 = end(res.gr), y0 = 0.5, y1 = 0.5 + scaled_y1,
        r0 = 0, r1 = r_height[5], col = "yellow")

    kpAddLabels(kp,
                labels = "L2FC", r0 = 0, r1 = r_height[5],
                cex = cex.y.lab, data.panel = 1)

    # Plot base L2FC
    kpAbline(kp, v = start(ia@vp), col = "black", lty = 2, lwd = 1)

    # Ideogram labels
    kpAddBaseNumbers(kp, tick.dist = 100000, tick.len = 10, cex = cex.ideo,
                    minor.tick.col = "gray")

    # Plot genes
    if (!is.null(genes_of_interest)) {
        plot_genes(ia, kp, genes_of_interest, gene.name.cex = gene.name.cex)
    }

    # Legend
    has_control <-!is.null(ia@metadata$control)
    lg <- if (has_control) c("condition", "control") else "condition"
    col <- if (has_control) c("firebrick4", "darkblue") else "firebrick4"
    legend(0.85, 0.8, legend = lg, fill = col, border = NA, bty = "o",
        cex = cex.leg)
}
