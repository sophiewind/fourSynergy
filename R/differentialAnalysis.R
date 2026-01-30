#' differentialAnalysis
#'
#' This function performs differential analysis to identify differential
#' interacting regions using DESeq2.
#'
#' @param ia fourSynergy object with interactions from all base tools
#' (peakC, r3c-seq, fourSig, r4cker) and other relevant information.
#' @param fitType Parameter for DESeq2s estimateDispersions(). Should be either
#' "parametric", "local", "mean", or "glmGamPoi" for the type of fitting of
#' dispersions to the mean intensity.
#'
#' @references https://doi.org/10.1186/s13059-014-0550-8
#' @return sia object with GRanges of DESeq results in the diff slot.
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
#' @export
differentialAnalysis <- function(ia, fitType = "local") {
    # No consensus present
    if (length(getExpConsensus(ia)) == 0) {
        stop(
            "fourSynergy found no interactions in condition. Did you run",
            " `consensusIa()`?"
        )
    } else if (length(getCtrlConsensus(ia)) == 0) {
        stop(
            "fourSynergy found no interactions in control. Did you run",
            " `consensusIa()`?"
        )
    }

    # No interactions in consensus
    else if (length(getExpConsensus(ia)[
        getExpConsensus(ia)$significance > 0]) == 0) {
        stop(
            "fourSynergy found no interactions in condition. A differential ",
            "analysis is not possible"
        )
    } else if (length(getCtrlConsensus(ia)[
        getCtrlConsensus(ia)$significance > 0]) == 0) {
        stop(
            "fourSynergy found no interactions in control. A differential ",
            "analysis is not possible"
        )
    } else {
        # Create colData -------------------------------------------------------
        coldata <- data.frame("condition" = rep(
            c("condition", "control"),
            each = length(getMetadata(ia)$conditionRep)
        ))
        rownames(coldata) <- c(paste(getMetadata(ia)$condition,
                                    getMetadata(ia)$conditionRep, sep = "_"),
            paste(getMetadata(ia)$control,
                getMetadata(ia)$controlRep, sep = "_"))

        ## Create consensus peaks
        cond <- getExpConsensus(ia)[
            which(getExpConsensus(ia)$significance > 0), ]
        ctrl <- getCtrlConsensus(ia)[
            which(getCtrlConsensus(ia)$significance > 0), ]
        cons <- GenomicRanges::union(na.omit(cond), na.omit(ctrl))

        ## Create count matrices
        if (ia@tracks != "") {
            bams <- paste0(ia@tracks, rownames(coldata), "_sorted.bam")
        } else {
            bams <- list.files(
                paste0("./results/", getMetadata(ia)$author, "/alignment/"),
                "_sorted.bam$"
            )
        }

        counts_coll <- lapply(bams, function(exp) {
            bamCount(exp, cons, shift = getMetadata(ia)$readLength / 2)
        })

        # Transform list to df
        samples <- as.data.frame(do.call(cbind, counts_coll))
        rownames(samples) <- as.character(cons)
        colnames(samples) <- gsub("_sorted.bam", "", basename(bams))

        #try(heatmap(as.matrix(samples)))

        ## Run DESeq2
        dds <- DESeqDataSetFromMatrix(
            countData = samples,
            colData = coldata,
            design = ~condition
        )
        resultsNames(dds)
        dds <- estimateSizeFactors(dds)
        tryCatch(
            {
                dds <- estimateDispersions(dds, fitType = fitType)
            },
            error = function(e) {
                message(
                    "DESeq estimateDispersions crashed. You can try to adjust ",
                    "the `fitType` param."
                )
            }
        )

        dds <- nbinomWaldTest(dds)
        res <- results(dds, contrast = c("condition", "condition", "control"))
        norm_counts <- counts(dds, normalized = TRUE)
        norm_counts_log <- log(norm_counts + 1, 10)
        ia <- setDifferential(ia, res)
        ia <- setDds(ia, dds)

        return(ia)
    }
}
