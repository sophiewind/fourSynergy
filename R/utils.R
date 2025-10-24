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
