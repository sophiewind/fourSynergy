#' Set differential attribute
#'
#' @param object fourSynergy object with interactions from all base tools.
#' @param value DESeqResults-object.
#'
#' @return Modified fourSynergy-object
#' @noRd
#' @keywords internal
setDifferential <- function(object, value) {
    if (inherits(value, "DESeqResults")) {
        object@differential <- value
    } else {
        stop("Value muste be DESeqResults object.")
    }
    return(object)
}

#' Set dds attribute
#'
#' @param object fourSynergy object with interactions from all base tools.
#' @param value DESeqDataSet-object.
#'
#' @return Modified fourSynergy-object
#' @noRd
#' @keywords internal
setDds <- function(object, value) {
    if (inherits(value, "DESeqDataSet")) {
        object@dds <- value
    } else {
        stop("Value must be DESeqDataSet object.")
    }
    return(object)
}

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
    if (length(ia@expConsensus) == 0) {
        stop(
            "fourSynergy found no interactions in condition. Did you run",
            " `consensusIa()`?"
        )
    } else if (length(ia@ctrlConsensus) == 0) {
        stop(
            "fourSynergy found no interactions in control. Did you run",
            " `consensusIa()`?"
        )
    }

    # No interactions in consensus
    else if (length(ia@expConsensus[ia@expConsensus$significance > 0]) == 0) {
        stop(
            "fourSynergy found no interactions in condition. A differential ",
            "analysis is not possible"
        )
    } else if (length(ia@ctrlConsensus[ia@ctrlConsensus$significance > 0]) ==
        0) {
        stop(
            "fourSynergy found no interactions in control. A differential ",
            "analysis is not possible"
        )
    } else {
        # Create colData -------------------------------------------------------
        coldata <- data.frame("condition" = rep(
            c("condition", "control"),
            each = length(ia@metadata$conditionRep)
        ))
        rownames(coldata) <- c(
            paste(ia@metadata$condition, ia@metadata$conditionRep, sep = "_"),
            paste(ia@metadata$control, ia@metadata$controlRep, sep = "_")
        )

        ## Create consensus peaks
        cond <- ia@expConsensus[which(ia@expConsensus$significance > 0), ]
        ctrl <- ia@ctrlConsensus[which(ia@ctrlConsensus$significance > 0), ]
        cons <- GenomicRanges::union(na.omit(cond), na.omit(ctrl))

        ## Create count matrices
        if (ia@tracks != "") {
            bams <- paste0(ia@tracks, rownames(coldata), "_sorted.bam")
        } else {
            bams <- list.files(
                paste0("./results/", ia@metadata$author, "/alignment/"),
                "_sorted.bam$"
            )
        }

        counts_coll <- lapply(bams, function(exp) {
            bamCount(exp, cons, shift = ia@metadata$readLength / 2)
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

        # plotMA(res, ylim = c(-2, 2))
        # plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")

        ia <- setDifferential(ia, res)
        ia <- setDds(ia, dds)

        return(ia)
    }
}
