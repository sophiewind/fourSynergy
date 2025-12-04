#' consensusIa
#'
#' This function performs an optimized weighted voting of 4C-seq tools.
#'
#' @param ia fourSynergy object with interactions from all base tools
#' (peakC, r3c-seq, fourSig, r4cker) and other relevant information.
#' @param model Selected optimization model. Either 'F1' or 'AUPRC'.
#'
#' @return fourSynergy object with interactions from all base tools and weighted
#' voting results.
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
#' @export
consensusIa <- function(ia, model = "F1") {
    dummy.gr <- makeGRangesFromDataFrame(
        data.frame(
            seqnames = paste0("chr1"),
            start = 1,
            end = 1,
            nReads = -1,
            significance = 0
        ),
        keep.extra.columns = TRUE,
        seqinfo = GenomeInfoDb::Seqinfo(genome = ia@metadata$organism)
    )
    vfl <- ia@vfl
    greater.than <- 0.5

    # Assign weights to vfl
    if (model %in% c("F1", "AUPRC")) {
        # Read weights
        weights <- utils::read.csv(system.file(
            "extdata",
            paste0("weighting_", tolower(model), ".csv"),
            package = "fourSynergy"
        )) %>%
            column_to_rownames("X")
        ia_tmp <- ia@expInteractions

        # Add weights and rating to intervals each tool is in one list element
        for (i in seq(1, nrow(weights))) {
            gr <- ia@expInteractions[startsWith(
                prefix = rownames(weights)[i],
                x = names(ia_tmp)
            )][[1]]
            mcols(ia_tmp[startsWith(
                prefix = rownames(weights)[i],
                x = names(ia_tmp)
            )][[1]])$weight <-
                weights[i, ]
        }

        # Calculate rating
        ia_tmp <- GRangesList(lapply(ia_tmp, \(x) {
            mcols(x)$rating <- mcols(x)$significance * mcols(x)$weight
            x
        }))

        if (length(ia@ctrlInteractions) > 0) {
            ia_tmp_ctrl <- ia@ctrlInteractions

            # Add weights and rating to intervals
            for (i in seq(1, nrow(weights))) {
                gr <- ia@ctrlInteractions[
                    startsWith(
                        prefix = rownames(weights)[i],
                        x = names(ia_tmp_ctrl)
                    )
                ][[1]]
                mcols(ia_tmp_ctrl[
                    startsWith(
                        prefix = rownames(weights)[i],
                        x = names(ia_tmp_ctrl)
                    )
                ][[1]])$weight <-
                    weights[i, ]
            }

            ia_tmp_ctrl <- GRangesList(lapply(ia_tmp_ctrl, \(x) {
                mcols(x)$rating <- mcols(x)$significance * mcols(x)$weight
                x
            }))
        }

        # Add rates to vfl
        for (i in seq_along(ia_tmp)) {
            ov <- findOverlaps(vfl, ia_tmp[[i]])
            name <- paste0("rate_", names(ia_tmp[i]))

            # Add name col
            mcols(vfl)[[name]] <- 0

            # Map ratings to vfl
            mcols(vfl[queryHits(ov)])[[name]] <-
                ia_tmp[[i]][subjectHits(ov)]$rating
        }

        # Add rates to vfl if control is available
        if (length(ia@ctrlInteractions) > 0) {
            for (i in seq_along(ia_tmp_ctrl)) {
                ov <- findOverlaps(vfl, ia_tmp_ctrl[[i]])
                name <- paste0("rate_", names(ia_tmp_ctrl[i]))

                # Add name col
                mcols(vfl)[[name]] <- 0

                # Map ratings to vfl
                mcols(vfl[queryHits(ov)])[[name]] <-
                    ia_tmp_ctrl[[i]][subjectHits(ov)]$rating
            }
        }

        # Sum up ratings and add to vfl
        rate.sum.cond <- mcols(vfl) %>%
            as.data.frame() %>%
            dplyr::select(matches("rate.*condition$")) %>%
            apply(1, function(x) {
                sum(unlist(x))
            })

        rate.sum.ctrl <- mcols(vfl) %>%
            as.data.frame() %>%
            dplyr::select(matches("rate.*control$")) %>%
            apply(1, function(x) {
                sum(unlist(x))
            })

        # condition
        maj.cond <- vfl
        mcols(maj.cond)$rate_total_condition <- rate.sum.cond

        # control
        maj.ctrl <- vfl
        mcols(maj.ctrl)$rate_total_control <- rate.sum.ctrl

        # Set maj. significance to 0 bc its old value
        maj.cond$significance <- 0
        if (!is.null(ia@ctrlInteractions)) {
            maj.ctrl$significance <- 0
        }

        # condition
        if (length(maj.cond) > 0 &&
            sum(maj.cond$rate_total_condition > greater.than) > 0) {
            mcols(maj.cond[maj.cond$rate_total_condition >
                greater.than, ])$significance <- 3
        } else if (model == "auprc") {
            maj.cond <- maj.cond
        } else {
            maj.cond <- dummy.gr
            maj.cond$rate_total_condition <- 0
        }

        # control
        if (length(ia@ctrlInteractions) > 0) {
            if (length(maj.ctrl) > 0 &&
                sum(maj.ctrl$rate_total_control > greater.than) > 0) {
                mcols(maj.ctrl[maj.ctrl$rate_total_control >
                    greater.than, ])$significance <- 3
            } else if (model == "auprc") {
                maj.ctrl <- maj.ctrl
            } else {
                maj.ctrl <- dummy.gr
                maj.ctrl$rate_total_control <- 0
            }
        }
        maj.cond$tool <- "weighted_voting"

        if (!is.null(ia@ctrlInteractions)) {
            maj.ctrl$tool <- "weighted_voting"
        }
        ia@expConsensus <- maj.cond

        if (!is.null(ia@ctrlInteractions)) {
            ia@ctrlConsensus <- maj.ctrl
        }
        return(ia)
    } else {
        stop("'model' should be either 'F1' or 'AUPRC'.")
    }
}
