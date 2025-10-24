# Setup data -------------------------------------------------------------------
# TODO load testdata
# ia hier einmal createn bzw. irgendeins im env haben
library(yaml)
library(DESeq2)

config <- system.file("extdata", "Datasets", "Demo", "info.yaml",
    package = "fourSynergy"
)
res_path <- system.file("extdata", "results", "Demo",
    package = "fourSynergy"
)
tracks <- system.file("extdata", "results", "Demo", "alignment",
    package = "fourSynergy"
)
sia <- createIa(res_path = res_path, config = config, tracks = tracks)
# createIa ---------------------------------------------------------------------
test_that("createIa inserts metadata in the right way.", {
    expect_equal(sia@metadata, read_yaml(system.file("extdata", "Datasets",
        "Demo", "info.yaml",
        package = "fourSynergy"
    )))
})

test_that("createIa reads all sia.bed files (experiment).", {
    expect_equal(
        names(sia@expInteractions),
        c(
            "rep.foursig_1.condition",
            "rep.foursig_3.condition",
            "rep.foursig_5.condition",
            "rep.foursig_11.condition",
            "rep.peakc_11.condition",
            "rep.peakc_21.condition",
            "rep.peakc_31.condition",
            "rep.peakc_51.condition",
            "rep.r3c_2000.condition",
            "rep.r3c_5000.condition",
            "rep.r3c_10000.condition",
            "rep.r4cker_nearbait.condition"
        )
    )
})

test_that("createIa reads all sia.bed files (control).", {
    expect_equal(
        names(sia@ctrlInteractions),
        c(
            "rep.foursig_1.control",
            "rep.foursig_3.control",
            "rep.foursig_5.control",
            "rep.foursig_11.control",
            "rep.peakc_11.control",
            "rep.peakc_21.control",
            "rep.peakc_31.control",
            "rep.peakc_51.control",
            "rep.r3c_2000.control",
            "rep.r3c_5000.control",
            "rep.r3c_10000.control",
            "rep.r4cker_nearbait.control"
        )
    )
})

test_that("createIa inserts significance in the right way.", {
    expect_equal(
        sum(sia@expInteractions$rep.peakc_11.condition$significance),
        45
    ) # TODO adjust path
})

test_that("createIa handles wrong input", {
    expect_error(
        createIa(res_path = "testpath", config = "config"),
        "The path 'testpath' does not exist."
    )
})

# test_that("createIa handles wrong input", {
#   expect_error(
#     createIa(res_path = "./inst/extdata/Dataset_test/", config = 'config'),
#     # TODO adjust path
#     "The config file 'config' does not exist."
#   )
# })

## consensusIa -----------------------------------------------------------------
test_that("consensusIa handles wrong input", {
    expect_error(
        consensusIa(sia, model = "X"),
        "'model' should be either 'F1' or 'AUPRC'."
    )
})

# Function tests ---------------------------------------------------------------
test_that("Differntial analysis requires consensusIa() results.", {
    expect_error(
        differentialAnalysis(sia),
        paste0(
            "fourSynergy found no interactions in condition. ",
            "Did you run `consensusIa()`?"
        )
    )
})


# Get consens with testweights
sia <- consensusIa(ia = sia, model = "AUPRC")

# right weights assigned?
# On simulated data no calls of no tool in firist row
test_that("consenusIa assigns right weights", {
    expect_equal(mcols(sia@expConsensus)$rate_total_condition[1], 0)
})

# Read in weights
weights <- read.csv(system.file("extdata",
    "weighting_auprc.csv",
    package = "fourSynergy"
)) %>%
    column_to_rownames("X")
weight.fs <- weights["rep.foursig_1", ]


# FourSig_1 has a significance of 1 in Range chr1 67333203-67333721,
# the weighting for fourSig is 0.004975252 so the signficance score
# should be 0.004975252.
test_that("consensusIa assigns right weights (condition).", {
    target_range <- GRanges("chr18", IRanges(67333203, 67333721))
    overlapping_consensus <- sia@expConsensus[sia@expConsensus %over%
        target_range]

    expect_equal(
        mcols(overlapping_consensus)$rate_rep.foursig_1.condition,
        1 * weight.fs
    )
})

# test_that("consenusIa assigns right weights (control)", {
#   target_range <- GRanges("chr1", IRanges(10000, 10500))
#   overlapping_consensus <- sia@ctrlConsensus %over% target_range
#
#   expect_equal(
#     mcols(overlapping_consensus)$rate_rep.peakc_21.control,
#     0
#   )
# })

test_that("Warnings is raised if no genes of interest are available.", {
    expect_warning(
        plotIaIndiviualTools(sia, "TP53"),
        "Genes of interest not found in the nearbait area."
    )
})

test_that("Warnings is raised if no genes of interest are available.", {
    expect_warning(
        plotConsensusIa(sia, "TP53"),
        "Genes of interest not found in the nearbait area."
    )
})

# Manipulate sia
sia_tmp <- sia
sia_tmp@expConsensus <- GRanges()
test_that("Warnings is raised if no genes of interest are available.", {
    expect_warning(
        plotConsensusIa(sia_tmp, "TP53"),
        "fourSynergy found no interactions in condition. Did you run",
        "`consensusIa()`?"
    )
})

# test_that("Warnings is raised if no genes of interest are available.", {
#   expect_warning(plotConsensusIa(sia_tmp),
#                  'fourSynergy found no interactions in condition.',
#                  'Either no consensus interactions could be found or you did not run',
#                  ' `consensusIa()`.')
# })


test_that("plotDiffIa requires differential results.", {
    expect_error(
        plotDiffIa(sia),
        paste0("Empty @differential slot. Did you run `differentialAnalysis()`?")
    )
})

sia <- suppressWarnings(differentialAnalysis(sia))

test_that("Differential analysis returns DESeqResults object.", {
    expect_true(inherits(sia@differential, "DESeqResults"))
})

kp <- createKaryoplot(sia)
test_that("Warnings is raised if no genes of interest are available.", {
    expect_warning(plot_genes(sia, kp, "X"))
})


test_that("plotDifferentialIa creates a plot.", {
    tmp <- tempfile(fileext = ".png")
    png(tmp)
    suppressWarnings(plotDiffIa(sia))
    dev.off()
    expect_true(file.exists(tmp))
})

# sia_no_diff <- sia
# sia_no_diff@differential <- DESeqResults(DataFrame = data.frame())
# test_that(paste0("plotDiffIa breaks if no significant differential ",
#           "interactions were found"), {
#             expect_warning(plotDiffIa(sia_no_diff),
#             "No significant differential interaction were found")
#           })
