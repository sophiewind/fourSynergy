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
