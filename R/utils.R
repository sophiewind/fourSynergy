#' readAndTag
#'
#' Internal function to read base tools solutions in FourSynergy format.
#'
#' @param file_path Path to the .bed files
#' (`results/[dataset]/nearbait_area.bed`).
#' @param tag tool name.
#' @param org organism.
#' @return GRanges with interaction calls and tool name as mcol.
#' @keywords internal
readAndTag <- function(file_path, tag, org) {
    # Create dummy file
    dummy.gr <- makeGRangesFromDataFrame(
        data.frame(
            seqnames = paste0("chr1"), start = 1, end = 1, nReads = -1,
            significance = 0
        ),
        keep.extra.columns = TRUE,
        seqinfo = GenomeInfoDb::Seqinfo(genome = org)
    )

    # Read and tag data
    if (file.exists(file_path)) {
        f <- read.delim(file_path, header = TRUE)
        f$tool <- tag
        f <- f[f$significance > 0, ]
        if (nrow(f) > 0) {
            return(makeGRangesFromDataFrame(f, keep.extra.columns = TRUE))
        } else {
            dummy.gr
        }
    } else {
        message("File not found: ", file_path)
        return(NULL)
    }
}

#' plotpreTracks
#'
#' @param ia fourSynergy object with interactions from all base tools.
#' @param highlight_regions regions to highlight in the plot.
#'
#' @returns list with read counts and Granges of bedfiles
#' @keywords internal
plotpreTracks <- function(ia, highlight_regions = NULL) {
    bgs <- readBedGraph(ia)
    tmp <- ia@vfl

    # Collect reads
    collect <- list()
    for (i in seq(1, length(bgs))){
        ov <- findOverlaps(tmp, bgs[[i]])
        tmp$reads <- 0
        tmp[queryHits(ov)]$reads <- bgs[[i]][subjectHits(ov)]$reads
        collect[[names(bgs)[i]]] <- tmp$reads
    }
    coll.df <- as.data.frame(collect)

    # Average reads
    cond.reads <- coll.df %>%
        dplyr::select(., matches(ia@metadata$condition)) %>%
        mutate(mean_cond = rowMeans(.)) %>%
        dplyr::select(mean_cond)
    tmp$cond_reads <- cond.reads
    if(!is.null(ia@metadata$control)){
        ctrl.reads <- coll.df %>%
            dplyr::select(., matches(ia@metadata$control)) %>%
            mutate(mean_ctrl = rowMeans(.)) %>%
            dplyr::select(mean_ctrl)
        tmp$ctrl_reads <- ctrl.reads
    }

    if (!is.null(highlight_regions)){
        if (grepl("chr[1-9XYM]+\\:\\d+\\-\\d+(,\\s*chr[1-9XYM]+:\\d+-\\d+)*",
                highlight_regions)){
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
            warning("The regions have to be either provided as string with ",
                    "following pattern: 'chr3:33000000-34000000' or as comma ",
                    "separated string (chr3:33000000-34000000, chr3:35000000-",
                    "35500000) or as valid .bed file.")
            reg <- NULL
        }
    } else {
        reg <- NULL
    }
    return(list(tmp = tmp, reg = reg))
}

#' checkConfig
#'
#' @param config config file with path.
#'
#' @returns TRUE if config is valid.
#' @export
#'
#' @examples
#' config <- system.file("extdata", "Datasets", "Demo", "info.yaml",
#'     package = "fourSynergy")
#' checkConfig(config)
checkConfig <- function(config) {
    if (is.null(config) || !(grepl("\\.yaml$", config, ignore.case = TRUE))) {
        stop("Info file is missing or not of type .yaml")
    }
    pattern <- "^[A-Za-z0-9_ \\+\\-]+$"
    config <- read_yaml(config, fileEncoding = "UTF-8")
    if (is.null(config$author) || !(is.character(config$author)) || !(grepl(pattern, config$author))) {
        stop("Author is missing or has an incorrect format.")
    }
    organism <- tolower(config$organism)
    pattern <- "^(hg19|hg38|mm9|mm10|[a-z]{2}\\d{1,2})$"
    if (is.null(config$organism) || !(grepl(pattern, organism))) {
        stop("Organism is missing or has an incorrect format.")
    }
    vpchr <- as.character(config$VPchr)
    pattern <- "^(\\d{1,2}|X|Y)$"
    if (is.null(config$VPchr) || !(grepl(pattern, vpchr, ignore.case = TRUE))) {
        stop("VPchr is missing or has an incorrect format.")
    }
    if (is.null(config$VPpos) || !(is.numeric(config$VPpos))) {
        stop("VPpos is missing or has an incorrect format.")
    }
    check_reps <- function(reps) {
        if (!is.list(reps) && !is.vector(reps) || length(reps) < 2) {
            return(FALSE)
        }
        if (!all(sapply(reps, is.numeric))) {
            return(FALSE)
        }
        return(TRUE)
    }
    if (is.null(config$control)) {
        warning("control is missing.")
    }
    if (is.null(config$controlRep) || !(check_reps(config$controlRep))) {
        warning("controlRep is missing or has an incorrect format.")
    }
    if (!is.null(config$condition)) {
        if (is.null(config$conditionRep) || !(check_reps(config$conditionRep))) {
            stop("condition is provided, but conditionRep is missing or has an incorrect format.")
        }
    }
    if (is.null(config$readLength) || !(is.numeric(config$readLength))) {
        stop("readLength is missing or has an incorrect format.")
    }
    primer_f <- config$PrimerF
    primer_r <- config$PrimerR
    pattern <- "^[GATCatgc]+$"
    if (is.null(primer_f) || !(grepl(pattern, primer_f)) || is.null(primer_r) || !(grepl(pattern, primer_r))) {
        stop("One or both primer are missing or have an incorrect format.")
    }
    if (is.null(config$REEnz) || !(is.list(config$REEnz)) && !(is.vector(config$REEnz))) {
        stop("REEnz is missing or not a list.")
    }
    if (!all(sapply(config$REEnz, is.character))) {
        stop("REEnz must be a list where all entries are strings.")
    }
    if (!is.null(config$RESeq) && (is.list(config$RESeq) || is.vector(config$RESeq))) {
        is_valid_seq <- all(sapply(config$RESeq, function(seq) grepl("^[GATCatgc]+$", seq)))
        if (!(is_valid_seq)) {
            stop("RESeq contains invalid symbols")
        }
    } else {
        stop("RESeq is missing or not a list.")
    }
    return(TRUE)
}




# string input
# input_a <- "chr3:1000-4000"
# input_b <- "chr3:1000-4000, chr2:1000-2000"
#
# input_b %>% stringr::str_split_1(', ') %>%
#     as.data.frame() %>%
#     separate('.', into = c("seqnames", "start", "end")) %>%
#     makeGRangesListFromDataFrame()
#
# input_a %>% stringr::str_split_1(', ') %>%
#     as.data.frame() %>%
#     separate('.', into = c("seqnames", "start", "end")) %>%
#     makeGRangesListFromDataFrame()
#
# # bed
# x <- read.delim('./test.bed') %>%
#     select(any_of(c('chr', 'seqnames', 'start', 'end')))
# colnames(x) <- c('seqnames', 'start', 'end')
#
# x %>%
#     makeGRangesListFromDataFrame()
#
#
# bed_check <- function(bed){
#     if (!file.exists(bed)){
#         stop("The bed file '", bed, "' does not exist.", call. = FALSE)
#     }
#
#     tryCatch({
#         b.f <- read.delim(bed, header = FALSE)
#         b.f <- b.f [,1:3]
#         colnames(b.f) <- c('seqnames', 'start', 'end')
#         b.f <- b.f %>%
#             makeGRangesListFromDataFrame()
#         return(bf)},
#         stop('bed not formatted as expected.')
#     )
# }
