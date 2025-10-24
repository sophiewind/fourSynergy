#' createIa
#'
#' This function reads the interaction bed files created by the pipeline and
#' transfers this information into an GrangesList.
#'
#' @param res_path Path to results created by the pipeline.
#' Typically stored in the `results/[dataset]/nearbait_area.bed`.
#' @param config Path of config file.
#' @param tracks Path to alignment files.
#'
#' @return fourSynergy object with interactions from all base tools.
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
#' ia <- createIa(res_path = res_path, config = config, tracks = tracks)
#'
#' @export
createIa <- function(res_path = character(), config = list(), tracks = "") {
    if (!dir.exists(res_path)) {
        stop("The path '", res_path, "' does not exist.", call. = FALSE)
    }
    if (!file.exists(config)) {
        stop("The config file '", config, "' does not exist.", call. = FALSE)
    }
    if (!dir.exists(tracks)) {
        stop("The tracks path '", tracks, "' does not exist.", call. = FALSE)
    }
    config <- yaml::read_yaml(config)
    dataset <- config$author
    ia <- GRangesList()
    vfl <- read.delim(paste0(res_path, "/nearbait_area.bed"), sep = "\t")
    area.frags <- makeGRangesFromDataFrame(vfl)
    bed_path <- paste0(res_path, "/sia/", dataset)

    if (!is.null(tracks)) tracks <- paste0(res_path, "/alignment/")
    files_to_read <- list(
        foursig = c(1, 3, 5, 11), peakcSig = c(11, 21, 31, 51),
        r3c = c(2000, 5000, 10000), r4cker = c("nearbait")
    )

    for (type in names(files_to_read)) {
        if (type == "r4cker") {
            suffixes <- c("_condition", "_control")
            for (suffix in suffixes) {
                f.p <- paste0(
                    bed_path, "_r4cker_nearbait", suffix,
                    "_nearbait.bed"
                )
                try(ia[[paste0("rep.r4cker_nearbait", gsub("_", ".", suffix))]]
                <- readAndTag(
                        file_path = f.p,
                        paste0("rep.r4cker_nearbait", suffix),
                        config$organism
                    ))
            }
        } else {
            for (i in files_to_read[[type]]) {
                suffixes <- c("_condition", "_control")
                for (suffix in suffixes) {
                    f.p <- paste0(
                        bed_path, "_", type, "_", i, suffix,
                        "_nearbait.bed"
                    )
                    try(ia[[paste0(
                        "rep.", gsub("Sig", "", type), "_", i,
                        gsub("_", ".", suffix)
                    )]] <-
                        readAndTag(
                            file_path = f.p, paste0(
                                "rep.", type,
                                "_", i, suffix
                            ),
                            config$organism
                        ))
                }
            }
        }
    }
    ia <- new("fourSynergy",
        metadata = config,
        expInteractions = ia[grepl("condition$", names(ia))],
        ctrlInteractions = ia[grepl("control$", names(ia))],
        vfl = area.frags,
        vp = GRanges(
            paste0("chr", config$VPchr),
            IRanges(config$VPpos, config$VPpos)
        ),
        tracks = tracks
    )
    return(ia)
}
