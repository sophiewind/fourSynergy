#' Internal function to create karyoplot
#'
#' @param ia fourSynergy object with interactions from all base tools
#' (peakC, r3c-seq, fourSig, r4cker) and other relevant information.
#' @param type Plottype.
#' @param cex character expansion
#' @param cex.axis character expansion axis
#' @param cex.lab character expansion labels
#' @param cex.main character expansion main
#'
#' @return karyoplot base
#' @keywords internal
createKaryoplot <- function(ia, type = 1, cex = 1, cex.axis = 1,
                            cex.lab = 1, cex.main = 1) {
    # Adjust plot settings
    pp <- getDefaultPlotParams(plot.type = type)

    pp$ideogramheight <- 5
    pp$data1inmargin <- 5
    pp$data2inmargin <- 5
    pp$topmargin <- 1
    pp$bottommargin <- 1
    pp$data1height <- 200
    pp$data2height <- 100
    pp$leftmargin <- 0.15

    # Display the plot
    plotKaryotype(
        genome = ia@metadata$organism,
        chromosomes = seqnames(ia@vp),
        plot.type = type,
        zoom = makeGRangesFromDataFrame(data.frame(
            seqnames = seqnames(ia@vp),
            start = start(ia@vfl[1]),
            end = end(ia@vfl[length(ia@vfl)])
        )), plot.params = pp,
        cex = cex, cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main
    )
}

#' Internal function to read bedGraphs
#' @param ia fourSynergy object with interactions from all base tools.
#' (peakC, r3c-seq, fourSig, r4cker) and other relevant information.
#'
#' @return GrangesList of bedGraph content
#'
#' @keywords internal
readBedGraph <- function(ia) {
    bgs <- GRangesList()
    if (!is.null(ia@tracks)) {
        for (i in ia@metadata$conditionRep) {
            try(cond <- read.delim(
                paste0(
                    ia@tracks, ia@metadata$condition, "_", i,
                    "_sorted.bedGraph"
                ),
                header = FALSE
            ) %>%
                `colnames<-`(c("seqnames", "start", "end", "reads")))
            try(cond <- cond[startsWith(cond$seqnames, "chr"), ] %>%
                makeGRangesFromDataFrame(keep.extra.columns = TRUE))
            try(bgs[[paste0(ia@metadata$condition, "_", i)]] <- cond)
        }

        for (i in ia@metadata$controlRep) {
            ctrl <- read.delim(
                paste0(
                    ia@tracks, ia@metadata$control, "_", i,
                    "_sorted.bedGraph"
                ),
                header = FALSE
            ) %>%
                `colnames<-`(c("seqnames", "start", "end", "reads"))
            ctrl <- ctrl[startsWith(ctrl$seqnames, "chr"), ] %>%
                makeGRangesFromDataFrame(keep.extra.columns = TRUE)
            bgs[[paste0(ia@metadata$control, "_", i)]] <- ctrl
        }
    }
    return(bgs)
}

#' Internal function to plot tracks
#'
#' @param ia fourSynergy object with interactions from all base tools
#' (peakC, r3c-seq, fourSig, r4cker) and other relevant information.
#' @param kp Kayroplot object.
#' @param bgs GrangesList of bedGraphs.
#' @param r0 r0 defines the vertical range of the data panel.
#' @param r1 r1 defines the vertical range of the data panel.
#'
#' @return karyoplot with tracks
#'
#' @keywords internal
plotTracks <- function(ia, kp, bgs, r0 = 0, r1 = 1, cex.vp = 1,
                        cex.y.track = 0.6) {
    for (bg in names(bgs)) {
        if (grepl(ia@metadata$condition, bg)) {
            kpPlotDensity(kp,
                        data.panel = 1, bgs[[bg]],
                        col = adjustcolor("firebrick4", alpha.f = 1 /
                                            length(bgs) * 2),
                        r0 = (r0 + r1) / 2, r1 = r1, window.size = 10000)
        } else {
            kpPlotDensity(kp,
                        data.panel = 1, bgs[[bg]],
                        col = adjustcolor("darkblue", alpha.f = 1 /
                                                length(bgs) * 2),
                        r0 = (r0 + r1) / 2, r1 = r0, window.size = 10000)
        }
    }
    kpText(kp,
        chr = paste0("chr", ia@metadata$VPchr), x = start(ia@vp),
        y = (r1 - 0.05), labels = "VP", data.panel = 1, cex = cex.vp,
        srt = 0, pos = 2, offset = 0.2)

    # Axis
    max_reads_per_replica <- vapply(
        bgs, function(gr) max(gr$reads, na.rm = TRUE),
        numeric(1))
    max_count <- max(max_reads_per_replica, na.rm = TRUE)
    pos.labels <- pretty(c(0,max_count), n = 10)
    plot_max <- max(pos.labels)
    ticks_to_use <- seq(-plot_max, plot_max, by = plot_max / 2)
    labels_to_show <- abs(ticks_to_use)
    kpAxis(kp, r0 = r0, r1 = 1, ymin = -plot_max, ymax = plot_max,
        numticks = 5, labels = labels_to_show, cex = cex.y.track, side = 1)
}

#' plot_genes
#'
#' Internal function to plot_genes.
#'
#' @param ia fourSynergy object with interactions from all base tools
#' @param kp karyoplot object
#' @param genes_of_interest genes of special interest
#' @keywords internal
#' @return karyoplot layer with genes
#' @noRd
plot_genes <- function(ia, kp, genes_of_interest, TxDb, panel = "2") {
    TxDb <- switch(ia@metadata$organism,
                "mm10" = TxDb.Mmusculus.UCSC.mm10.knownGene,
                "hg19" = TxDb.Hsapiens.UCSC.hg19.knownGene)

    genes.data <- makeGenesDataFromTxDb(TxDb, karyoplot = kp)
    genes.data <- addGeneNames(genes.data)
    genes.data <- mergeTranscripts(genes.data)
    sg <- genes.data$genes[mcols(genes.data$genes)$name %in%
                            genes_of_interest]

    if (length(sg) == 0 && genes_of_interest != "all") {
        warning("Genes of interest not found in the nearbait area.")
    }

    plotg <- function(data) {
        plot_params <- list(
            karyoplot = kp,
            data = data,
            gene.name.cex = 0.6,
            gene.name.position = "right",
            r0 = 0.2,
            data.panel = "2")
        do.call(kpPlotGenes, plot_params)
    }

    if (length(genes_of_interest) == 1 && genes_of_interest == "all") {
        plotg(genes.data)
    } else {
        for (i in sg$name) {
            gene <- sg[sg$name == i]
            # zr <- GRanges(
            #     seqnames = unique(as.character(seqnames(gene))),
            #     ranges = IRanges(
            #         start = max(1, start(gene) - 1000),
            #         end = end(gene) + 1000
            #     ),
            #     strand = strand(gene)
            # )
            zr <- GRanges(
                seqnames = unique(as.character(seqnames(gene))),
                ranges = IRanges(
                    start = max(1, start(gene)) +
                        ((end(gene) - max(1, start(gene)))/2),
                    end = max(1, start(gene)) +
                        ((end(gene) - max(1, start(gene)))/2) + 1
                ),
                strand = strand(gene)
            )
            pdf(NULL)
            kp_zoom <- plotKaryotype(
                genome = ia@metadata$organism,
                chromosomes = as.character(seqnames(gene)),
                zoom = zr,
                plot.type = 1
            )
            gene_data <- makeGenesDataFromTxDb(TxDb, karyoplot = kp_zoom)
            gene_data <- addGeneNames(gene_data)
            gene_data <- mergeTranscripts(gene_data)
            dev.off()
            plotg(gene_data)
        }
    }
}

#' Internal function to highlight regions in karyoplot.
#'
#' @param ia fourSynergy object with interactions from all base tools
#' (peakC, r3c-seq, fourSig, r4cker) and other relevant information.
#' @param kp karyoplot.
#' @param highlight_regions regions to highlight in the plot.
#'
#' @return karyoplot base
#' @keywords internal
plotRegions <- function(ia, kp, highlight_regions){
    if (str_detect(highlight_regions,
                "chr[1-9XYM]+\\:\\d+\\-\\d+(,\\s*chr[1-9XYM]+:\\d+-\\d+)*")){
        reg <- highlight_regions %>%
            stringr::str_split_1(',') %>%
            trimws() %>%
            as.data.frame() %>%
            separate('.', into = c("seqnames", "start", "end")) %>%
            makeGRangesListFromDataFrame()
    } else if (endsWith(highlight_regions, '.bed')){
        b.f <- read.delim(highlight_regions, header = FALSE)
        if (b.f[1,2] == "start"){
            b.f <- read.delim(highlight_regions, header = TRUE)
        }
        b.f <- b.f [,seq(1, 3)]
        colnames(b.f) <- c('seqnames', 'start', 'end')
        reg <- b.f %>%
            makeGRangesListFromDataFrame()
    } else {
        stop("The regions have to be either provided as string with following ",
            "pattern: 'chr3:33000000-34000000' or as comma separated string ",
            "(chr3:33000000-34000000, chr3:35000000-35500000) or as valid ",
            ".bed file.")
    }

    for (i in seq(1, length(reg))){
        kpBars(kp,
            chr = paste0("chr", ia@metadata$VPchr),
            x0 = start(reg[[i]]),
            x1 = end(reg[[i]]),
            y0 = 0, y1 = 1,
            r0 = 0.5, r1 = 1, col = adjustcolor("grey", alpha.f = 0.4),
            border = NA)
    }
}

