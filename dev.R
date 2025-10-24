# library(yaml)  # read_yaml
# library(GenomicRanges)
# library(dplyr)
# library(karyoploteR)
# library(bamsignals)
# library(DESeq2)
# library(tibble)
# sapply(list.files('/home/rstudio/fourSynergy///R', full.names = T),
#        source)
aut <- "Mantis_Ldlrad4_TGFb"

# install.packages('../fourSynergy_0.0.0.9000.tar.gz', source = T)
library(fourSynergy)


config <- read_yaml(paste0("/host/Datasets/", aut, "/info.yaml"))

# into globals or utils
dummy.gr <- makeGRangesFromDataFrame(
    data.frame(
        seqnames = paste0("chr", config$VPchr),
        start = 0,
        end = 0,
        nReads = -1,
        significance = 0
    ),
    keep.extra.columns = T,
    # seqinfo = Seqinfo(genome = config$organism)
)


ia <- createIa(
    res_path = paste0("/host/results/", aut),
    config = paste0("/host/Datasets/", aut, "/info.yaml"),
    tracks = paste0("/host/results/", aut, "/alignment")
)
ia@metadata
plotIaIndiviualTools(ia = ia, show_track = TRUE)

# F1
ia_cons <- consensusIa(ia = ia, model = "AUPRC")

if (config$author %in% c("Harris_c3_c4_VP-20")) {
    try(res <- differentialAnalysis(ia = ia_consens, fitType = "mean"))
    ## mean for DESEq
} else {
    try(res <- differentialAnalysis(ia = ia_cons))
    ## mean for DESEq
}
# try(saveRDS(res, paste0('./', config$author, '_diff.rds')))
try(plotDiffIa(
    res = res,
    ia = ia_consens,
    paste0("./", config$author, "_differential.png")
))
rm(res)

# AUPRC
ia_consens <- consensusIa(
    ia = ia,
    "x",
    out = TRUE,
    metric = "auprc"
)
if (config$author %in% c()) {
    try(res <- differentialAnalysis(ia = ia_consens, fitType = "mean"))
    ## mean for DESEq
} else {
    try(res <- differentialAnalysis(ia = ia_consens))
    ## mean for DESEq
}

try(plotDiffIa(
    res = res,
    ia = ia_consens,
    paste0("./", config$author, "_differential_auprc.png")
))
# rm(res)
