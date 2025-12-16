#' plotConsensusTracks
#'
#' @param ia fourSynergy object with interactions from all base tools
#' @param highlight_regions regions to highlight in the plot
#' @param max_range maximum plotting range
#'
#' @returns Track-plots for all treatments with interactions from consensus tool
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
#' sia <- consensusIa(sia, model = "AUPRC")
#' plotConsensusTracks(sia)
plotConsensusTracks <- function(ia, highlight_regions = NULL,
                                max_range = 3000) {
    if(is.null(ia) || !(is(ia, "fourSynergy"))){
        stop("Object is not a fourSynergy object or not provided.")
    }
    if (length(ia@expConsensus) == 0) {
        stop("fourSynergy found no interactions in condition. Did you run",
            " `consensusIa()`?"
        )
    }
    res <- plotpreTracks(ia, highlight_regions)
    if(is.null(res$tmp)){
        stop("res does not containt tmp. Check plotpreTracks")
    }
    tmp <- res$tmp
    tmp$ens_cond <- ia@expConsensus$significance
    try(tmp$ens_ctrl <- ia@ctrlConsensus$significance)
    tmp.df <- as.data.frame(tmp)

    create_plot_ens <- function(df, mean_col, ens_col){
        p <- df %>% dplyr::select(start, end, all_of(mean_col),
                                all_of(ens_col)) %>%
            pivot_longer(!c(all_of(mean_col), 'start', 'end'),
                        names_to = "peak") %>%
            mutate(value = as.factor(value)) %>%
            ggplot()
        if (!is.null(res$reg)) {
            p <- p + geom_rect(data = as.data.frame(res$reg),
                            mapping = aes(xmin =  start, ymin =  0,
                                            xmax = end, ymax = Inf),
                            alpha = 0.5)
        }
        p <- p + geom_segment(aes(x = start, xend = end, y = 0,
                                yend = .data[[mean_col]]),
                            color = 'grey40') +
            ylim(c(0, max_range)) +
            geom_segment(aes(x = start, xend = end, y = 0,
                            yend = .data[[mean_col]],
                            alpha = value, color = peak)) +
            scale_alpha_manual(values = c(0, 1, 1, 1)) +
            facet_grid(. ~ peak) +
            scale_color_manual(values = setNames("purple", ens_col)) +
            guides(alpha = "none", color = 'none') +
            labs(y = 'Counts?', x = 'Position') +
            theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
            geom_vline(aes(xintercept = start(ia@vp)), linetype = 'dashed')
    }
    p.cond <- create_plot_ens(tmp.df, "mean_cond", "ens_cond")

    if (!is.null(ia@metadata$control)){
        p.ctrl <- create_plot_ens(tmp.df, "mean_ctrl", "ens_ctrl")
    } else {
        p.ctrl <- NULL
    }
    return(list(p.cond, p.ctrl))
}
