#' plotBaseTracks
#'
#' @param ia fourSynergy object with interactions from all base tools
#' @param highlight_regions regions to highlight in the plot
#' @param max_range maximum plotting range
#'
#' @returns Track-plots for all treatments with interactions from base tools
#' @export
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
#' plotBaseTracks(sia)
plotBaseTracks <- function(ia, highlight_regions = NULL, max_range = 3000) {
    if (is.null(ia) || !(is(ia, "fourSynergy"))){
        stop("Object is not a fourSynergy object or not provided.")
    }
    res <- plotpreTracks(ia, highlight_regions)
    if(is.null(res$tmp)){
        stop("res does not containt tmp. Check plotpreTracks")
    }
    p.col <- c("rep.r3c_2000" = "#fdae61", "rep.r4cker_nearbait" = "#228B22",
            "rep.peakc_21" = "#9c0142", "rep.foursig_1" = "#3288bd")

    if (!is.null(ia@metadata$control)){
        tools <- paste0('rep.', rep(c('peakc_21', 'foursig_1', 'r3c_2000',
                             'r4cker_nearbait'), each = 2),
               c('.condition', '.control'))
    } else {
        tools <- paste0('rep.', rep(c('peakc_21', 'foursig_1', 'r3c_2000',
                                      'r4cker_nearbait')), '.condition')
    }
    for (p in tools) {
        if (endsWith(p, 'condition')) {
            inter <- ia@expInteractions[[p]]
        } else {
            try(inter <- ia@ctrlInteractions[[p]])
        }

        ov <- findOverlaps(res$tmp, inter)
        mcols(res$tmp)[[p]] <- 0
        mcols(res$tmp[queryHits(ov)])[[p]] <- inter[subjectHits(ov)]$significance
    }
    tmp.df <- as.data.frame(res$tmp)
    create_plot <- function(df, treat, mean_col) {
        set <- c(mean_col, paste0("rep.r3c_2000.", treat),
                paste0("rep.r4cker_nearbait.", treat),
                paste0("rep.peakc_21.", treat),
                paste0("rep.foursig_1.", treat))
        p <- df %>%  dplyr::select(start, end, all_of(set)) %>%
            pivot_longer(!c(all_of(mean_col), 'start', 'end'),
                    names_to = "peak") %>%
            mutate(short_peak = sub("\\.(condition|control)$", "", peak),
                value = as.factor(value)) %>% ggplot()

        if (!is.null(res$reg)) {
            p <- p + geom_rect(as.data.frame(res$reg),
                            aes(xmin = start, ymin =  0, xmax = end,
                                ymax = Inf), alpha = 0.5)
        }
        p <- p + geom_segment(aes(x = start, xend = end, y = 0,
                                yend = .data[[mean_col]]), color = 'grey40') +
            ylim(0, max_range) +
            geom_segment(aes(x = start, xend = end, y = 0,
                            yend = .data[[mean_col]],
                            alpha = value, color = short_peak)) +
            scale_alpha_manual(values = c(0, 1, 1, 1)) +
            facet_grid(. ~ peak) +
            scale_color_manual(values = p.col) +
            guides(alpha = "none", color = 'none') +
            labs(y = 'Counts', x = 'Position') +
            theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
            geom_vline(aes(xintercept = start(ia@vp)), linetype = 'dashed')
        return(p)
    }
    p.cond <- create_plot(tmp.df, "condition", "mean_cond")
    if(!is.null(ia@metadata$control)){
        p.ctrl <- create_plot(tmp.df, "control", "mean_ctrl")
    } else {
        p.ctrl <- NULL
    }
    return(list(p.cond, p.ctrl))
}
