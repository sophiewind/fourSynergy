# Laden der erforderlichen Bibliotheken
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyBS)
library(shinyvalidate)
library(dplyr)
library(ggplot2)
library(reshape2)
library(DT)
library(stringr)
library(tidyr)
library(fourSynergy)
library(bslib)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(biomaRt)
library(org.Mm.eg.db)
library(cowplot)
library(UpSetR)
library(gVenn)
library(shinycustomloader)
library(DESeq2)  #
library(bamsignals)  #
library(karyoploteR)  #
library(tibble)  #
source('/home/rstudio/fourSynergy/utils.R')
source('/home/rstudio/fourSynergy/R/karyoplots.R')
source('/home/rstudio/fourSynergy/R/plotIaIndividualTools.R')
source('/home/rstudio/fourSynergy/R/plotConsensusIa.R')
source('/home/rstudio/fourSynergy/R/plotDifferentialIa.R')
source('/home/rstudio/fourSynergy/R/differentialAnalysis.R')

# Definition der UI
ui <- shinyUI(dashboardPage(
  skin = "green",
  title = "fourSynergy",
  dashboardHeader(
    title = "fourSynergy"
  ),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload", tabName = "upload"),
      menuItem("QC", tabName = "qc"),
      menuItem("Base Tools", tabName = "base_tools"),
      menuItem("Ensemble", tabName = "ensemble"),
      menuItem("Differential Analysis", tabName = "diff_analysis"),
      menuItem("Pipeline and Support", tabName = "pipeline"))
  ),

  # Inhalt
  dashboardBody(
    useShinyjs(),
    tags$style(HTML("


.box.box-solid.box-primary>.box-header {
  color:#FFFFFF;
  background:#00a65a
                    }

.box.box-solid.box-primary{
border-bottom-color:#666666;
border-left-color:#666666;
border-right-color:#666666;
border-top-color:#666666;
}

                                    ")),
    style = "background-color: #FFFFFF;",
    tabItems(
      tabItem(tabName = "upload",
              fluidRow(
                column(12,
                       fileInput("config", "Upload config:", accept = ".yaml"),
                       fileInput("zip", "Zip and upload results/ directory:", multiple = TRUE),
                       fileInput("sia", "Upload results/sia directory:", multiple = TRUE)),
              )
      ),

      tabItem(tabName = "qc",
              fluidRow(
                column(12,
                       h2("Metadata"),
                       fluidRow(column(6,
                                       value_box(
                                         title = "Viewpoint",
                                         value = textOutput("meta_vp"),
                                         showcase = icon("eye"),
                                         theme = "green",
                                         showcase_layout = "top right",
                                       )),
                                column(6,
                                       value_box(
                                         title = "First and 2nd restiction enzyme",
                                         value = textOutput("meta_re"),
                                         showcase = icon("scissors"),
                                         theme = "green",
                                         showcase_layout = "top right",
                                         margin = margin(t = 10, b = 10, l = 10,
                                                         r = 10),
                                       ))),
                       div(style = "height: 20px;"),
                       fluidRow(column(6,
                                       # valueBoxOutput("meta_org")),
                                       value_box(
                                         title = "Organism",
                                         value = textOutput("meta_orga"),
                                         showcase = icon("compass"),
                                         theme = "green",
                                         showcase_layout = "top right",
                                         margin = margin(t = 10, b = 10, l = 10, r = 10)
                                       )),
                                column(6,
                                       value_box(
                                         title = "Read length",
                                         value = textOutput("meta_rl"),
                                         showcase = icon("ruler"),
                                         theme = "green",
                                         showcase_layout = "top right",
                                         full_screen = TRUE
                                       ),),
                       ),
                       div(style = "height: 20px;"),
                       fluidRow(column(6,
                                       value_box(
                                         title = "Condition",
                                         subtitle = textOutput("meta_rep_co"),
                                         value = textOutput("meta_cond_1"),
                                         showcase = icon("capsules"),
                                         theme = "green",
                                         showcase_layout = "top right",
                                         full_screen = TRUE
                                       ),),
                                column(6,
                                       value_box(
                                         title = "Control",
                                         subtitle = textOutput("meta_rep_ct"),
                                         value = textOutput("meta_ctrl_2"),
                                         showcase = icon("x"),
                                         theme = "green",
                                         showcase_layout = "top right",
                                         full_screen = TRUE
                                       ))),

                       hr(),

                       #dataTableOutput("metadata"),
                       h2("FASTQC statistics"),
                       h3("Sequencing statistics"),
                       dataTableOutput("fastqc_tab"),
                       column(6,
                              h3("Sequence counts"),
                              withLoader(plotOutput("seq_plot"), loader = "dnaspin")),
                       column(6,
                              h3("Per sequence quality scores"),
                              withLoader(plotOutput("seq_qual"), loader = "dnaspin")),
                       hr(),
                       h2("Alignment statistics"),
                       dataTableOutput("align_tab"),
                       hr(),
                       h2("4C-seq statistics"),
                       dataTableOutput("basic_tab")
                )
              )
      ),

      # Base Tools ---------------------------
      tabItem(tabName = "base_tools",
              h2("Base tools karyoplot"),
              fluidRow(
                column(6, withLoader(plotOutput("karyo_base", height = 800),
                                     loader = "dnaspin")),
                column(6,
                       box(title = 'Karyoplot settings',
                           solidHeader = TRUE,
                           status = "primary",
                           width = 12,
                           h3("Select genes of interest"),
                           p("Please be aware that plotting genes can take some time..."),
                           withLoader(uiOutput("genes_select_base"),
                                       loader = "dnaspin"),
                           textInput("in_regions", "Please enter genomic regions to highlight. If you want to highlight multiple regions please separate the regions using ','.", "chr:start-end"),
                           #checkboxInput("spider_base", "Show spiderplot."),
                           actionButton("button_genes_base", "Update karyoplot", icon = icon("reload")))),
                hr(),
              ),
              h2("Base tool track plot"),
              withLoader(plotOutput("plot_bt", height = 800), loader = "dnaspin"),
              h2('Interaction calls per tool'),
              #fluidRow(
              h3('Condition'),
              dataTableOutput("tab_base_cond"),
              #column(6,
              h3('Control'),
              dataTableOutput("tab_base_ctrl"),
              h2("Compare condition and control"),
              box(title = "UpsetR",
                  solidHeader = TRUE,
                  width = 12,
                  status = 'primary',
                  radioButtons("rb_upset", "Select what you want to compare",
                               choices = c("Tools in condition" = "tico",
                                           "Tools in control"= "tict",
                                           "Condition and control"= "coco"), selected = "tico"),
                  hr(),
                  checkboxGroupInput("cb_tools", label = "Tools to compare",
                                     choices = c("foursig_1",
                                                 "foursig_3",
                                                 "foursig_5",
                                                 "foursig_11",
                                                 "r3c_2000",
                                                 "r3c_5000",
                                                 "r3c_10000",
                                                 "peakc_11",
                                                 "peakc_21",
                                                 "peakc_31",
                                                 "peakc_51",
                                                 "r4cker_nearbait"), inline = T),


              ),

              plotOutput("upset_base", height = 800),
              hr()
      ),

      # Ensemble -----------------------------
      tabItem(tabName = "ensemble",
              fluidRow(
                column(6,
                       box(title = "Ensemble interaction calling settings",
                           solidHeader = TRUE,
                           status = "primary",
                           width = 12,

                       radioButtons("rb_model",
                                    "Choose F1 or AUPRC optimized model.",
                                    choices = c("F1", "AUPRC")),
                       bsButton(inputId = "pu",
                                label = "",
                                icon = icon("info"),
                                style = "info",
                                size = "extra-small"),
                       bsPopover(id = "pu",
                                 title = "Model",
                                 content = paste0(
                                   "Weighted voting is either optimized on F1 score or Area Under Precision-Recall Curve."
                                 ),
                                 placement = "right",
                                 trigger = "hover",
                                 options = list(container = "body")
                       ),
                       actionButton("run_ens", "Start interaction calling.", icon = icon("play"))),
                ),

                column(6,
                       box(title = 'Karyoplot settings',
                           solidHeader = TRUE,
                           status = "primary",
                           h3("Select genes of interest"), width = 12,
                           p("Please be aware that plotting genes can take some time..."),
                           withLoader(uiOutput("genes_select_ens"), loader = "dnaspin"),
                           textInput("in_regions_ens", "Please enter genomic regions to highlight.", "chr:start-end"),
                           checkboxInput("spider_ens", "Show spiderplot.", FALSE),
                           actionButton("button_genes_ens", "Update karyoplot", icon = icon("reload"))),
                )
              ),


              withLoader(plotOutput("karyo_ens", height = 1200), loader = "dnaspin"),
              h2('Interaction calls fourSynergy'),
              fluidRow(
                column(6, h3('Condition'),
                       dataTableOutput("tab_ens_cond")),
                column(6,
                       h3('Control'),
                       dataTableOutput("tab_ens_ctrl"))),
              h2("Compare tools"),
              plotOutput("upset_ens"),
              selectInput("area_venn", "Select area from Venn diagram to see which regions are included:",
                          choices = c("condition",
                                      "overlap",
                                      "control")),
              dataTableOutput("tab_ens")
      ),

      # Differential Analysis ----
      tabItem(tabName = "diff_analysis",
              fluidRow(
                column(6,
                       column(6,
                              value_box(
                                title = "Condition",
                                subtitle = textOutput("meta_rep_co"),
                                value = textOutput("meta_cond_2"),
                                showcase = icon("capsules"),
                                theme = "green",
                                showcase_layout = "top right",
                                full_screen = TRUE
                              ),),
                       column(6,
                              value_box(
                                title = "Control",
                                subtitle = textOutput("meta_rep_ct"),
                                value = textOutput("meta_ctrl_1"),
                                showcase = icon("x"),
                                theme = "green",
                                showcase_layout = "top right",
                                full_screen = TRUE
                              ))),
                column(6,
                       box(h3("Select genes of interest"), width = 6,
                           p("Please be aware that plotting genes can take some time..."),
                           uiOutput('genes_select_diff'),
                           textInput("in_regions_diff", "Please enter genomic regions to highlight.", "chr:start-end"),
                           checkboxInput("spider_diff", "Show spiderplot.", FALSE),
                           actionButton("button_genes_diff", "Update karyoplot", icon = icon("reload"))))),
                withLoader(plotOutput('karyo_diff', height = 800), loader = "dnaspin"),
                fluidRow(
                  column(6,
                         h2("Heatmap"),
                         plotOutput("hm_diff")),
                  column(6,
                         h2("MA"),
                         plotOutput("plot_ma"))),
              h2("Results differential analysis (DESeq2)"),
                dataTableOutput("tab_diff"),


      ),
      tabItem("pipeline",
              #fluidRow(
                  h1("Pipeline"),
                  h1("Support"),
                  h2("Code availibilty"),
                  "The fourSynergy R package is available here: TODO\n",
                  "The fourSynergy snakemake pipeline can be found here: TODO\n",
                  h2("Documentation"),
                  "A comprehensive description of fourSynergy is in the Vignette (TODO).\n",
                  h2("Bugs"),
              "Bugs can be reported here: https://github.com/sophiewind/scafari/issues",
              h2("Help"),
              "For questions, feature requests and errors, so not hestitate to contact Sophie Wind (sophie.wind@uni-muenster.de).\n"
              )
    )
  )
))


# Definition des Servers
server <- function(input, output, session) {
    options(shiny.maxRequestSize = 500 * 1024^2)

    # Input regions control
    iv <- InputValidator$new()

    # Add rules for each region input
    iv$add_rule("in_regions", sv_regex("chr[1-9XYM]+:\\d+-\\d+|chr:start-end",
                                            "The genomic region should follow the scheme chr:start-end. If you do not want to highlight a region enter `chr:start-end`."))
    iv$add_rule("in_regions_ens", sv_regex("chr[1-9XYM]+:\\d+-\\d+|chr:start-end",
                                            "The genomic region should follow the scheme chr:start-end. If you do not want to highlight a region enter `chr:start-end`."))
    iv$add_rule("in_regions_diff", sv_regex("chr[1-9XYM]+:\\d+-\\d+|chr:start-end",
                                            "The genomic region should follow the scheme chr:start-end. If you do not want to highlight a region enter `chr:start-end`."))

    iv$enable()

    # General validation function
    is_valid_region <- function(region_input) {
        grepl("chr[1-9XYM]+:\\d+-\\d+|chr:start-end", region_input)
    }

    # Reactives for validation
    is_valid_base <- reactive({ is_valid_region(input$in_regions) })
    is_valid_ens  <- reactive({ is_valid_region(input$in_regions_ens) })
    is_valid_diff <- reactive({ is_valid_region(input$in_regions_diff) })

    # Observe and enable/disable buttons
    observeEvent(is_valid_base(), {
        if (is_valid_base()) {
            shinyjs::enable("button_genes_base")
        } else {
            shinyjs::disable("button_genes_base")
        }
    })

    observeEvent(is_valid_ens(), {
        if (is_valid_ens()) {
            shinyjs::enable("button_genes_ens")
        } else {
            shinyjs::disable("button_genes_ens")
        }
    })

    observeEvent(is_valid_diff(), {
        if (is_valid_diff()) {
            shinyjs::enable("button_genes_diff")
        } else {
            shinyjs::disable("button_genes_diff")
        }
    })


  # Create sia ----
  # sia <- createIa('/host/results/Geeven_sox/',
  #                 '/host/Datasets/Geeven_sox/info.yaml',
  #                 '/host/results/Geeven_sox/alignment/')
  # sia <- consensusIa(sia, model = 'AUPRC')
  # sia <- differentialAnalysis(sia)

  # Metadata ----
  output$meta_org <- renderValueBox({
    valueBox(
      value = sia@metadata$organism,
      subtitle = "organism",
      icon = if (grepl("^mm", sia@metadata$organism)) {
        icon("mouse")
      } else {
        icon("person")
      }
      # color = "blue",
    )
  })
  output$meta_re <- renderText({
    paste(sia@metadata$REEnz, collapse = ", ")
  })

  output$meta_vp <- renderText({
    paste0("chr:", sia@metadata$VPchr, "\n", start(sia@vp))
  })
  output$meta_rl <- renderText({
    sia@metadata$readLength
  })
  output$meta_cond_1 <- renderText({
    sia@metadata$condition
  })

  output$meta_cond_2 <- renderText({
      sia@metadata$condition
  })

  output$meta_rep_co <- renderText({
    max(sia@metadata$conditionRep)
  })

  output$meta_orga <- renderText({
    sia@metadata$organism
  })
  output$meta_ctrl_1 <- renderText({
    sia@metadata$control
  })

  output$meta_ctrl_2 <- renderText({
      sia@metadata$control
  })
  output$meta_rep_ct <- renderText({
    max(sia@metadata$controlRep)
  })


  # Upload ----
  output$upload_text <- renderText({
    paste("You selected:", input$datei)
  })

  output$upload_text <- renderText({
    paste("You selected:", input$datei)
  })

  output$metadata <- renderDataTable({
    file.out <- paste0('fourSyerngy_metadata_', sia@metadata$author)
    sia@metadata %>%
      unlist() %>%
      as.matrix(ncol = 2) %>%
      as.data.frame() %>%
      datatable(colnames = "",
                extensions = "Buttons",
                options = list(pageLength = 21,
                               dom = "Bfrtip",
                               buttons = list(
                                   list(extend = "csv", filename = file.out),
                                   list(extend = "excel", filename = file.out),
                                   list(extend = "copy", filename = file.out))))
  }, server = FALSE)

  # Seq ----
  output$fastqc_tab <- renderDataTable({
      file.out <- paste0('fourSyerngy_fastqc_', sia@metadata$author)

    # Quality table
    read.delim("/host/results/Geeven_sox/multiqc_data/multiqc_fastqc.txt") %>%
      DT::datatable(., extensions = "Buttons",
                    options = list(
                        dom = "Bfrtip",
                        buttons = list(
                            list(extend = "csv", filename = file.out),
                            list(extend = "excel", filename = file.out),
                            list(extend = "copy", filename = file.out)
                        )))
  }, server = FALSE)

  # TODO change order, etc. -> see multiqc
  output$seq_plot <- renderPlot({
    read.delim("/host/results/Geeven_sox/multiqc_data/fastqc_sequence_counts_plot.txt") %>%
      melt() %>%
      ggplot(., aes(x = value, y = Sample, fill = variable)) +
      geom_bar(stat = 'identity')
  })

  output$seq_qual <- renderPlot({
    read.delim("/host/results/Geeven_sox/multiqc_data/fastqc_per_sequence_quality_scores_plot.txt") %>%
      pivot_longer(cols = starts_with("X"), names_to = "time", values_to = "val") %>%
      mutate(
        val = gsub("\\(|\\)", "", val),
      ) %>%
      separate(val, c('x', 'y'), sep = ', ') %>%
      mutate(x = as.numeric(x),
             y = as.numeric(y)) %>%
      dplyr::select(Sample, time, x, y) %>%
      ggplot(., aes(x = x, y = y, color = Sample)) +
      geom_rect(aes(xmin = -Inf, xmax = 20, ymin = -Inf, ymax = Inf),
                fill = "pink", alpha = 0.03, color = NA) +
      geom_rect(aes(xmin = 20, xmax = 28.5, ymin = -Inf, ymax = Inf),
                fill = "yellow", alpha = 0.03, color = NA) +
      geom_rect(aes(xmin = 28.5, xmax = Inf, ymin = -Inf, ymax = Inf),
                fill = "palegreen", alpha = 0.03, color = NA) +
      geom_line() +
      xlim(c(0, 40)) +
      labs(x = "Mean sequence quality (Phred score)",
           y = "Count")
  })


  output$align_tab <- renderDataTable({
      file.out <- paste0('fourSyerngy_flagstat_', sia@metadata$author)

    read.delim('/host/results/Geeven_sox/multiqc_data_1/samtools-flagstat-pct-table.txt') %>%
      datatable(colnames = paste0(gsub('\\.', ' ', colnames(.)), " (%)"),
                rownames = FALSE,
                extensions = "Buttons",
                options = list(
                    dom = "Bfrtip",
                    buttons = list(
                        list(extend = "csv", filename = file.out),
                        list(extend = "excel", filename = file.out),
                        list(extend = "copy", filename = file.out)
                    )))
  }, server = FALSE)

  output$basic_tab <- renderDataTable({
      file.out <- paste0('fourSyerngy_basic4cseq_', sia@metadata$author)

    files <- list.files("/host/results/Geeven_sox/basic4cseq/stats/",
                        pattern = "*.txt", full.names = TRUE)
    dfs <- lapply(files, read.delim, sep = ":", header = FALSE)

    # Finde die gemeinsamen Spalten
    coll <- Reduce(rbind, dfs)
    coll$Sample <- rep(gsub("_stats.txt", "", basename(files)), each = 4)

    # coll$ratio <- coll$V2 %>%
    #     str_extract(., "\\([^)]+\\)") %>%
    #     gsub("\\(|\\)", "", .) %>%
    #     str_extract("\\d+.\\d+") %>%
    #     as.numeric()
    coll <- coll %>%
      pivot_wider(id_cols = V1, names_from = Sample, values_from = V2)

    colnames(coll)[1] <- ''
    df <- coll %>%
      t() %>%
      as.data.frame()

    colnames(df) <- df[1, ]

    df <- df[-1, ] # remove the first row which was used as headers

    datatable(df, extensions = "Buttons",
              options = list(
                  dom = "Bfrtip",
                  buttons = list(
                      list(extend = "csv", filename = file.out),
                      list(extend = "excel", filename = file.out),
                      list(extend = "copy", filename = file.out)
                  )))
  }, server = FALSE)


  ## Base Tools ----
  output$base_tools_text <- renderText({
    ""
  })

  # BAsic karyoplot
  output$karyo_base <- renderPlot({
    plotIaIndiviualTools(
      sia, cex.chr = 2, cex.ideo = 1,
      cex.y.track = 1, cex.y.lab = 1, cex.vp = 1.5)
  })

  observeEvent(input$button_genes_base, {

    if (is.null(input$genes_select_base)) {
      # Behandlung wenn keine Gene ausgewählt sind
      genes <- NULL
    } else {
      genes <- input$genes_select_base
      message(paste("Ausgewählte Gene:", paste(genes, collapse = ", ")))
    }

    if (input$in_regions == "chr:start-end") {
      output$karyo_base <- renderPlot({
        plotIaIndiviualTools(
          sia, cex.chr = 2, cex.ideo = 1,
          cex.y.track = 1, cex.y.lab = 1, cex.vp = 1.5,
          genes_of_interest = genes
        )
      })
    } else {
      output$karyo_base <- renderPlot({
        plotIaIndiviualTools(
          sia, cex.chr = 2, cex.ideo = 1,
          cex.y.track = 1, cex.y.lab = 1, cex.vp = 1.5,
          genes_of_interest = genes,
          highlight_regions = input$in_regions
        )
      })
    }
  })

  # output$karyo_base <- renderPlot({
  #   if (is.null(input$genes_select_base)) {
  #     plotIaIndiviualTools(sia, cex.chr = 2, cex.ideo = 1,
  #                          cex.y.track = 1, cex.y.lab = 1, cex.vp = 1.5)
  #   } else {
  #     genes <- input$genes_select_base
  #     print(genes)
  #     plotIaIndiviualTools(sia, cex.chr = 2, cex.ideo = 1,
  #                          cex.y.track = 1, cex.y.lab = 1, cex.vp = 1.5,
  #                          genes_of_interest = genes)
  #   }
  # })



  output$plot_bt <- renderPlot({
    p <- trackPlots(sia)
    plot_grid(p[[1]] +
                theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(),
                      text = element_text(size = 16)),
              p[[2]] +
                theme(text = element_text(size = 16)), ncol = 1)
  })

  observeEvent(input$button_genes_base, {
      if (input$in_regions == "chr:start-end") {
          output$plot_bt <- renderPlot({
          p <- trackPlots(sia)
          plot_grid(p[[1]] +
                        theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.title.x = element_blank(),
                              text = element_text(size = 16)),
                    p[[2]] +
                        theme(text = element_text(size = 16)), ncol = 1)})
    } else {
        message('regions')
      reg <- input$in_regions %>%
          stringr::str_split_1(', ') %>%
          as.data.frame() %>%
          separate('.', into = c("seqnames", "start", "end")) %>%
          makeGRangesListFromDataFrame()
      output$plot_bt <- renderPlot({
      p <- trackPlots(sia, reg)
      plot_grid(p[[1]] +
                  theme(axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        text = element_text(size = 16)),
                p[[2]] +
                  theme(text = element_text(size = 16)), ncol = 1)})
    }
  })

  output$genes_select_base <- get_gene_selection("base")
  output$genes_select_ens <- get_gene_selection("ens")
  output$genes_select_diff <- get_gene_selection("diff")

  # Condition
  output$tab_base_cond <- renderDataTable({
    file.out <- paste0('fourSyerngy_base_algorithm_condition_',
                       sia@metadata$author)

    vfl <- sia@vfl
    for (i in seq(1, length(sia@expInteractions))){
      ov <- findOverlaps(vfl, sia@expInteractions[[i]])
      mcols(vfl)[[sia@expInteractions[[i]]$tool[1]]] <- '-'
      mcols(vfl[queryHits(ov),])[sia@expInteractions[[i]]$tool[1]] <- 'X'
    }

    vfl.df <- vfl %>% as.data.frame(keep.extra.columns = TRUE)# %>%
    colnames(vfl.df) <- gsub('000', 'k', gsub('Sig', '', gsub('_condition', '', gsub('rep.', '', colnames(vfl.df)))))


    # sia@expInteractions %>%  unlist() %>%
    #   as.data.frame(keep.extra.columns = TRUE) %>%
    #   dplyr::select(seqnames, start, end, tool) %>%
    #   `colnames<-`(c('Seqnames', 'Start', 'End', 'Tool')) %>%
    vfl.df %>%
      dplyr::select(-width, -strand) %>%
      filter(if_any(everything(.), ~. == "X")) %>%
      datatable(rownames = FALSE,
                extensions = "Buttons",
                options = list(
                    dom = "Bfrtip",
                    buttons = list(
                        list(extend = "csv", filename = file.out),
                        list(extend = "excel", filename = file.out),
                        list(extend = "copy", filename = file.out)
                    )))
  }, server = FALSE)

  output$tab_base_ctrl <- renderDataTable({
    file.out <- paste0('fourSyerngy_base_algorithm_control_', sia@metadata$author)

    vfl <- sia@vfl
    for (i in seq(1, length(sia@ctrlInteractions))){
      ov <- findOverlaps(vfl, sia@ctrlInteractions[[i]])
      mcols(vfl)[[sia@ctrlInteractions[[i]]$tool[1]]] <- '-'
      mcols(vfl[queryHits(ov),])[sia@ctrlInteractions[[i]]$tool[1]] <- 'X'
    }

    vfl.df <- vfl %>% as.data.frame(keep.extra.columns = TRUE)
    colnames(vfl.df) <- gsub('000', 'k',
                             gsub('Sig', '',
                                  gsub('_control', '',
                                       gsub('rep.', '', colnames(vfl.df)))))


    vfl.df %>%
      dplyr::select(-width, -strand) %>%
      filter(if_any(everything(.), ~. == "X")) %>%
      datatable(rownames = FALSE,
                extensions = "Buttons",
                options = list(
                    dom = "Bfrtip",
                    buttons = list(
                        list(extend = "csv", filename = file.out),
                        list(extend = "excel", filename = file.out),
                        list(extend = "copy", filename = file.out)
                    )))
  }, server = FALSE)

  ## Upset ---
  # genomic_overlaps <- computeOverlaps(sia@expInteractions)
  # plotUpSet(genomic_overlaps)
  observeEvent(input$rb_upset, {
    updateCheckboxGroupInput(session, "cb_tools", selected = NULL)
  })

  output$upset_base <- renderPlot({
    if (input$rb_upset == "tico"){
      sel <- input$cb_tools
      if (length(sel) <= 2){
        ggplot(data.frame(), aes()) +
          theme_void() +
          labs(title =  "Please select at least three tools to create an upset plot.")
      }
      else{
        sia@expInteractions[names(sia@expInteractions) %in% paste0('rep.', sel, '.condition')] %>%
          as.list() %>%
          lapply(function(x) as.data.frame(x, stringsAsFactors = FALSE)) %>%
          lapply(function(x) x %>%
                   dplyr::mutate(region = paste0(seqnames, ":", start, "-", end))) %>%
          lapply(function(x) x %>%
                   dplyr::select(region) %>%
                   unlist) %>%
          fromList(.) %>%
          upset(.,order.by = "freq", text.scale = 2, nsets = length(sel))
      }
    } else if (input$rb_upset == "tict"){
      sel <- input$cb_tools
      if (length(sel) <= 2){
        ggplot(data.frame(), aes()) +
          theme_void() +
          labs(title =  "Please select at least three tools to create an upset plot.")
      }
      else{
        sia@ctrlInteractions[names(sia@ctrlInteractions) %in% paste0('rep.', sel, '.control')] %>%
          as.list() %>%
          lapply(function(x) as.data.frame(x, stringsAsFactors = FALSE)) %>%
          lapply(function(x) x %>%
                   dplyr::mutate(region = paste0(seqnames, ":", start, "-", end))) %>%
          lapply(function(x) x %>%
                   dplyr::select(region) %>%
                   unlist) %>%
          fromList(.) %>%
          upset(.,order.by = "freq", text.scale = 2, nsets = length(sel))
      }
    } else {
      sel <- input$cb_tools
      if (length(sel) <= 1){
        ggplot(data.frame(), aes()) +
          theme_void() +
          labs(title =  "Please select at least two tools to create an upset plot.")
      }
      else{
        sel <- input$cb_tools
        comb <- c(sia@expInteractions, sia@ctrlInteractions)
        comb[names(comb) %in% c(paste0("rep.", sel, ".condition"), paste0("rep.", sel, ".control"))] %>%
          as.list() %>%
          lapply(function(x) as.data.frame(x, stringsAsFactors = FALSE)) %>%
          lapply(function(x) x %>%
                   dplyr::mutate(region = paste0(seqnames, ":", start, "-", end))) %>%
          lapply(function(x) x %>%
                   dplyr::select(region) %>%
                   unlist) %>%
          fromList(.) %>%
          upset(.,order.by = "freq", nsets = length(sel)*2, text.scale = 2)
      }
    }
  })

  # Plot kary ens ----
  output$karyo_ens <- renderPlot({
      plotConsensusIa(sia, cex.chr = 2, cex.ideo = 1,
                                 cex.y.track = 1, cex.vp = 1.5)
  })

  observeEvent(input$button_genes_ens, {
    if (is.null(input$genes_select_ens)) {
      genes <- NULL
    } else {
      genes <- input$genes_select_ens
    }

    if (input$in_regions_ens == "chr:start-end") {
      message(genes)
      if (input$spider_ens == FALSE){
        output$karyo_ens <- renderPlot({
          plotConsensusIa(sia, cex.chr = 2, cex.ideo = 1,
                          cex.y.track = 1, cex.vp = 1.5,
                          genes_of_interest = genes)
        })
      } else {
          # Record plot
          k_ens <- reactive({
              plotConsensusIa(sia, cex.chr = 2, cex.ideo = 1,
                              cex.y.track = 1, cex.vp = 1.5, plot_spider = TRUE)
              recordPlot()
          })
          # Display plot
          output$karyo_ens <- renderPlot({
              replayPlot(k_ens())
          })
      }

    } else {
      if (input$spider_ens == FALSE){
      output$karyo_ens <- renderPlot({
        plotConsensusIa(sia, cex.chr = 2, cex.ideo = 1,
                        cex.y.track = 1, cex.vp = 1.5,
          genes_of_interest = genes,
          highlight_regions = input$in_regions_ens
        )
      })
      } else {
          k_ens <- reactive({
              plotConsensusIa(sia, cex.chr = 2, cex.ideo = 1,
                              cex.y.track = 1, cex.vp = 1.5,
                              genes_of_interest = genes,
                              highlight_regions = input$in_regions_ens,
                              plot_spider = TRUE
              )
              recordPlot()
          })

          output$karyo_ens <- renderPlot({replayPlot(k_ens())})

      }
    }
  })
  #---

  output$tab_ens_cond <- renderDataTable({
    file.out <- paste0('fourSyerngy_ensemble_condition_', sia@metadata$author)
    sia@expConsensus %>%
      as.data.frame(keep.extra.columns = TRUE) %>%
      dplyr::select(seqnames, start, end, tool) %>%
      `colnames<-`(c('Seqnames', 'Start', 'End')) %>%
      datatable(rownames = FALSE, extensions = "Buttons",
                options = list(
                    dom = "Bfrtip",
                    buttons = list(
                        list(extend = "csv", filename = file.out),
                        list(extend = "excel", filename = file.out),
                        list(extend = "copy", filename = file.out)
                    )))
  }, server = FALSE)

  output$tab_ens_ctrl <- renderDataTable({
    file.out <- paste0('fourSyerngy_ensemble_control_', sia@metadata$author)
    sia@ctrlConsensus %>%
      as.data.frame(keep.extra.columns = TRUE) %>%
      dplyr::select(seqnames, start, end, tool) %>%
      `colnames<-`(c('Seqnames', 'Start', 'End')) %>%
      datatable(rownames = FALSE, extensions = "Buttons",
                options = list(
                    dom = "Bfrtip",
                    buttons = list(
                        list(extend = "csv", filename = file.out),
                        list(extend = "excel", filename = file.out),
                        list(extend = "copy", filename = file.out)
                    )))
  }, server = FALSE)

  output$upset_ens <- renderPlot({
    gVenn::computeOverlaps(GRangesList(Condition = sia@expConsensus[sia@expConsensus$significance > 0], Control = sia@ctrlConsensus[sia@ctrlConsensus$significance > 0])) %>%
      plotVenn(fontsize = 11)
  })


  output$tab_ens <- renderDataTable({
    file.out <- paste0('fourSyerngy_ensemble_', input$area_venn,
                       '_', sia@metadata$author)
    ov <- gVenn::computeOverlaps(GRangesList(Condition = sia@expConsensus[sia@expConsensus$significance > 0], Control = sia@ctrlConsensus[sia@ctrlConsensus$significance > 0]))
    if (input$area_venn == "condition"){
      set <- ov$reduced_regions[ov$reduced_regions$intersect_category == '10']
    } else if (input$area_venn == 'overlap'){
      set <- ov$reduced_regions[ov$reduced_regions$intersect_category == '11']
    } else {
      set <- ov$reduced_regions[ov$reduced_regions$intersect_category == '01']
    }
    set %>%
      as.data.frame() %>%
      dplyr::select(-width, -strand, -intersect_category) %>%
      datatable(., extensions = "Buttons",
                rownames = FALSE,
                options = list(
                    dom = "Bfrtip",
                    buttons = list(
                        list(extend = "csv", filename = file.out),
                        list(extend = "excel", filename = file.out),
                        list(extend = "copy", filename = file.out)
                    )))
  }, server = FALSE)

  # Diff ----
  output$karyo_diff <- renderPlot({
    plotDiffIa(sia)  #, cex.chr = 2, cex.ideo = 1, cex.y.track = 1, cex.vp = 1.5
  })

  # observeEvent(input$button_genes_diff, {
  #
  #   if (is.null(input$genes_select_diff)) {
  #     genes <- NULL
  #   } else {
  #     genes <- input$genes_select_diff
  #   }
  #
  #   if (input$in_regions_diff == "chr:start-end") {
  #     if (input$spider_diff == FALSE){
  #
  #         k_diff <- reactive({
  #       plotDiffIa(
  #         sia,
  #         #cex.chr = 2, cex.ideo = 1,
  #         #cex.y.track = 1, cex.vp = 1.5,
  #         genes_of_interest = genes
  #       )
  #             recordPlot()
  #
  #     })
  #     } else {
  #         k_diff <- reactive({
  #             plotDiffIa(
  #                 sia,
  #                 #cex.chr = 2, cex.ideo = 1,
  #                 #cex.y.track = 1, cex.vp = 1.5,
  #                 genes_of_interest = genes, plot_spider = TRUE
  #             )
  #             recordPlot()
  #         })
  #
  #     }
  #       output$karyo_diff <- renderPlot({
  #           replayPlot(k_diff())
  #       })
  #   } else {
  #       if (input$in_regions_diff == "chr:start-end") {
  #           k_diff <- reactive({
  #         plotDiffIa(
  #           sia,
  #           #cex.chr = 2, cex.ideo = 1,
  #           #cex.y.track = 1, cex.vp = 1.5,
  #           genes_of_interest = genes,
  #           highlight_regions = input$in_regions_diff, spider = TRUE)
  #           recordPlot()
  #           })
  #
  #       } else {
  #           k_diff <- reactive({
  #               plotDiffIa(
  #                   sia,
  #                   #cex.chr = 2, cex.ideo = 1,
  #                   #cex.y.track = 1, cex.vp = 1.5,
  #                   genes_of_interest = genes, plot_spider = TRUE
  #               )
  #               recordPlot()
  #           })
  #
  #       }
  #       output$karyo_diff <- renderPlot({
  #           replayPlot(k_diff())
  #       })
  #   }
  # })
  observeEvent(input$button_genes_diff, {

      # Pick genes from the correct input
      if (is.null(input$genes_select_diff)) {
          genes <- NULL
      } else {
          genes <- input$genes_select_diff     # corrected here!
      }

      # TRUE if user specified a region for highlighting
      region_highlight <- input$in_regions_diff != "chr:start-end"

      k_diff <- reactive({
          plotDiffIa(
              sia,
              genes_of_interest = genes,
              highlight_regions = if (region_highlight) input$in_regions_diff else NULL,
              plot_spider = input$spider_diff
          )
          recordPlot()
      })

      output$karyo_diff <- renderPlot({
          replayPlot(k_diff())
      })
  })
  #---

  output$tab_diff <- renderDataTable({
    file.out <- paste0('fourSyerngy_differential_', sia@metadata$author)
    sia@differential %>%
      as.data.frame() %>%
      datatable(., extensions = "Buttons",
                options = list(
                    dom = "Bfrtip",
                    buttons = list(
                        list(extend = "csv", filename = file.out),
                        list(extend = "excel", filename = file.out),
                        list(extend = "copy", filename = file.out)
                    )))
  }, server = FALSE)

  output$plot_ma <- renderPlot({
    sia@differential %>%
      plotMA()
  })

  output$hm_diff <- renderPlot({
    sia@dds %>%
      counts() %>%
          heatmap
  })

}

# Starten der App
shinyApp(ui = ui, server = server)
