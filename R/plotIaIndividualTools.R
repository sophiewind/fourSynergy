#'
#' This function creates a karyoplot with the interactions calls of the
#' individual tools.
#'
#' @param ia fourSynergy object with interactions from all base tools
#' (peakC, r3c-seq, fourSig, r4cker) and other relevant information.
#' @param genes_of_interest Vector with genes of interest.
#' @param cex.chr character expansion of chromosome label.
#' @param cex.ideo character expansion base numbers of ideogram.
#' @param cex.y.track character expansion y axis track.
#' @param cex.y.lab character expansion y lab.
#' @param cex.vp character expansion viewpoint label.
#' @param highlight_regions regions to highlight in the plot
#'
#' @return karyoplot with calling results.
#' @export
#'
#' @examples
#'
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
#' plotIaIndiviualTools(ia = sia, genes_of_interest = c("Ldlrad4", "Cep76"))
plotIaIndiviualTools <- function(ia = GRangesList(), genes_of_interest = NULL,
                                 cex.chr = 1, cex.ideo = 0.6, cex.y.track = 0.6,
                                 cex.y.lab = 0.6, cex.vp = 1,
                                 highlight_regions = NULL) {
  # Create baseplot
  if (is.null(genes_of_interest)) {
    kp <- createKaryoplot(ia, cex = cex.chr)
  } else {
    kp <- createKaryoplot(ia, type = "2", cex = cex.chr)
  }
  
  # Tracks
  bgs <- readBedGraph(ia)
  if (!is.null(highlight_regions)){
    for (i in seq(1, length(highlight_regions))){
      kpBars(kp,
             chr = paste0("chr", ia@metadata$VPchr), 
             x0 = start(highlight_regions[[i]]),
             x1 = end(highlight_regions[[i]]),
             y0 = 0, y1 = 1,
             r0 = 0.4, r1 = 1, col = adjustcolor("grey", alpha.f = 0.4),
             border = NA
      )
    }
  }
  
  plotTracks(ia, kp, bgs,
             r0 = 0.5, r1 = 1, cex.vp = cex.vp,
             cex.y.track = cex.y.track)
  
  # Interactions
  groups <- list(
    list(prefix = "rep.foursig_", nums = c(1, 3, 5, 11), col = "#3288bd"),
    list(prefix = "rep.peakc_", nums = c(11, 21, 31, 51), col = "#9c0142"),
    list(prefix = "rep.r3c_", nums = c(2000, 5000, 10000), col = "#fdae61"),
    list(prefix = "rep.r4cker_", nums = "nearbait", col = "forestgreen"))
  
  total_regions <- sum(vapply(
    groups, function(grp) length(grp$nums),
    integer(1)))
  slot_height <- 1 / total_regions
  r0_slots <- seq(0, (total_regions - 1)) * 0.5
  r1_slots <- seq(1, (total_regions)) * 0.5
  counter <- 1
  
  for (grp in groups) {
    for (i in seq_along(grp$nums)) {
      n <- grp$nums[i]
      varname <- paste0(grp$prefix, n, ".condition")
      label <- paste0(sub("^rep\\.", "", grp$prefix), n)
      r0 <- as.numeric(r0_slots[counter] * slot_height)
      r1 <- as.numeric(r1_slots[counter] * slot_height) * 0.95
      kpPlotRegions(kp, ia@expInteractions[[varname]],
                    col = grp$col, r0 = r0 + (r1 - r0) / 2, r1 = r1
      )
      kpAddLabels(kp,
                  labels = label, r0 = r0 + (r1 - r0) / 2, r1 = r1,
                  cex = cex.y.lab
      )
      counter <- counter + 1
    }
  }
  
  # Genes
  if (!is.null(genes_of_interest)) {
    plot_genes(ia, kp, genes_of_interest)
  }
  
  # VP
  kpAbline(kp, v = start(ia@vp), col = "black", lty = 2, data.panel = 1)
  kpAddBaseNumbers(kp,
                   tick.dist = 100000, tick.len = 10, cex = cex.ideo,
                   minor.tick.col = "gray"
  )
  
  # Legend
  legend(0.85, 0.8,
         legend = c("condition", "control"),
         fill = c("firebrick4", "darkblue"), border = NA,
         bty = "o", cex = 0.6
  )
}