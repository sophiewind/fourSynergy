#' plotDiffIa
#'
#' This function creates a karyoplot with the differential interactions calls.
#'
#' @param ia fourSynergy object with interactions from all base tools
#' (peakC, r3c-seq, fourSig, r4cker) and more relevant information.
#' @param genes_of_interest Vector with genes of interest.
#' @param cex.chr character expansion of chromosome label.
#' @param cex.ideo character expansion base numbers of ideogram.
#' @param cex.y.track character expansion y axis track.
#' @param cex.vp character expansion viewpoint label.
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
                        cex.ideo = 0.6, cex.y.track = 0.6, cex.vp = 1) {
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
    kp <- createKaryoplot(ia, type = panel, cex = cex.chr)

    # Structure plot
    kpAbline(kp, chr = paste0("chr", ia@metadata$VPchr), h = 0.25)
    kpAbline(kp, chr = paste0("chr", ia@metadata$VPchr), h = 0.0625)
    kpAbline(kp, chr = paste0("chr", ia@metadata$VPchr), h = 0.25 - 0.0625)
    kpAbline(kp, chr = paste0("chr", ia@metadata$VPchr), h = 0.125)

    # Plot interactions
    kpPlotRegions(kp, res.gr[res.gr$padj < .05], col = "grey", r0 = 0.25,
                    r1 = 1)
    plotTracks(ia, kp, bgs,
        r0 = 0.25, r1 = 1, cex.vp = cex.vp,
        cex.y.track = cex.y.track
    )
    kpPlotRegions(kp, ia@expConsensus[ia@expConsensus$significance > 0],
        col = "firebrick4", r0 = 0.25 - 0.0625, r1 = 0.25
    )
    kpPlotRegions(kp, ia@ctrlConsensus[ia@ctrlConsensus$significance > 0],
        col = "darkblue", r0 = 0.125, r1 = 0.125 + 0.0625
    )

    max_value <- max(abs(res.gr$log2FoldChange))
    scaled_y1 <- (res.gr$log2FoldChange / max_value)

    # Plot L2FC
    kpBars(kp,
        chr = paste0("chr", ia@metadata$VPchr), x0 = start(res.gr),
        x1 = end(res.gr), y0 = 0.125 - 0.0625, y1 = 0.125 - 0.0625 + scaled_y1,
        r0 = 0.125 - 0.0625, r1 = 0.125, col = "yellow"
    )

    # Plot base L2FC
    kpAbline(kp, v = start(ia@vp), col = "black", lty = 2, lwd = 1)

    # Ideogram labels
    kpAddBaseNumbers(kp,
        tick.dist = 100000, tick.len = 10, cex = cex.ideo,
        minor.tick.col = "gray"
    )

    # Plot genes
    if (!is.null(genes_of_interest)) {
        plot_genes(ia, kp, genes_of_interest)
    }

    # Legend
    legend(0.85, 0.8,
        legend = c("condition", "control"), bty = "o", cex = 0.6,
        fill = c("firebrick4", "darkblue"), border = NA
    )
}
