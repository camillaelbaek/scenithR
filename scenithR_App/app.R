# app.R — ScenithR - Shiny app for analysis of SCENITH assay 
# Sequential gating workflow: Singlets → Live/Dead → Signal
# Channel roles configured via preset or custom selection
# Metadata: well_code, genotype, perturbation, treatment (opt), time (opt)
required_cran <- c(
  "shiny", "dplyr", "tidyr", "ggplot2", "DT",
  "stringr", "purrr", "scales", "readr",
  "readxl", "ggridges", "ggpubr", "viridis",
  "ggbeeswarm", "sp", "magick", "imager"
)

required_bioc <- c(
  "flowCore", "flowViz", "ggcyto", "openCyto"
)

installed <- installed.packages()[, "Package"]

missing_cran <- setdiff(required_cran, installed)

if (length(missing_cran) > 0) {
  install.packages(
    missing_cran,
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

installed <- installed.packages()[, "Package"]

missing_bioc <- setdiff(required_bioc, installed)

if (length(missing_bioc) > 0) {
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

suppressPackageStartupMessages({
  library(shiny)
  library(flowCore)
  library(flowViz)
  library(ggcyto)
  library(openCyto)
  library(scales)
  library(tidyverse)
  library(ggpubr)
  library(viridis)
  library(DT)
  library(stringr)
  library(sp)
  library(grid)
  library(readxl)
  library(readr)
})

# ── Panel presets ──────────────────────────────────────────────────────────────
PANEL_PRESETS <- list(
  "Scenith – APC (puromycin) / FITC (live-dead)"  = list(scatter_x="FSC.A", scatter_y="FSC.H", live_dead="FITC.A",  signal="APC.A"),
  "Scenith – APC (puromycin) / BV421 (live-dead)" = list(scatter_x="FSC.A", scatter_y="FSC.H", live_dead="BV421.A", signal="APC.A"),
  "Scenith – PE (puromycin) / FITC (live-dead)"   = list(scatter_x="FSC.A", scatter_y="FSC.H", live_dead="FITC.A",  signal="PE.A"),
  "No live-dead stain – APC signal"                = list(scatter_x="FSC.A", scatter_y="FSC.H", live_dead=NULL,      signal="APC.A"),
  "Custom – configure below"                       = list(scatter_x=NULL,    scatter_y=NULL,    live_dead=NULL,      signal=NULL)
)

DEFAULT_VERTICES <- list(x1=5000, y1=0, x2=8000, y2=20000, x3=120000, y3=90000, x4=120000, y4=68000)

perturbation_cols_default <- c(
  "Co"="black", "DG"="#7B1FA2", "O"="#D32F2F",
  "DGO"="#00897B", "UNST"="#1976D2", "DMSO_25uL"="#9E9E9E"
)
genotype_cols_default <- c(
  "WT"="#1B9E77", "WT_HSV"="#7570B3", "KO"="#D95F02", "KO_HSV"="#E7298A"
)

# ── Helpers ────────────────────────────────────────────────────────────────────
normalize_well_code <- function(w) {
  w        <- toupper(trimws(as.character(w)))
  row_part <- substr(w, 1, 1)
  col_part <- suppressWarnings(as.integer(substr(w, 2, nchar(w))))
  ifelse(is.na(col_part), NA_character_, sprintf("%s%02d", row_part, col_part))
}

auto_palette <- function(values, known) {
  new_v <- setdiff(unique(na.omit(as.character(values))), names(known))
  if (!length(new_v)) return(known)
  c(known, setNames(scales::hue_pal(l=55, c=70)(length(new_v)), new_v))
}

theme_alba <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face="bold", size=base_size+2),
      plot.subtitle    = element_text(size=base_size),
      strip.background = element_rect(fill="grey95"),
      strip.text       = element_text(face="bold"),
      panel.grid.minor = element_blank()
    )
}

geomfi        <- function(x) { x <- x[x > 0]; exp(mean(log(x), na.rm=TRUE)) }
save_plot_png <- function(p, f, w=7, h=5) ggsave(f, p, width=w, height=h, dpi=300, units="in")
empty_gg      <- function(msg) ggplot() + annotate("text", x=.5, y=.5, label=msg, size=4.5, colour="grey50") + theme_void()

# as.numeric() cast prevents "expected double, got integer" from numericInput
make_singlet_gate <- function(x1, y1, x2, y2, x3, y3, x4, y4, xchan, ychan) {
  m <- matrix(as.numeric(c(x1,y1, x2,y2, x3,y3, x4,y4)), ncol=2, byrow=TRUE)
  colnames(m) <- c(xchan, ychan)
  polygonGate(filterId="Singlets", gate=m)
}

make_thresh_gate <- function(channel, lo=-Inf, hi=Inf, id="gate") {
  p <- list(c(lo, hi)); names(p) <- channel
  do.call(rectangleGate, c(list(filterId=id), p))
}

# ── UI ─────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .small-note { color:#555; font-size:.92em; }
    .ok   { color:#176; font-weight:600; }
    .busy { color:#a33; font-weight:600; }
    .tag-box { background:#f4f4f4; border-radius:6px; padding:8px 12px;
               font-size:.9em; margin-top:6px; }
    .gate-tbl td, .gate-tbl th { padding:3px 6px; }
    .step-note { color:#555; font-style:italic; font-size:.9em; margin-bottom:6px; }
  "))),
  titlePanel("Scenith analysis"),

  sidebarLayout(
    sidebarPanel(width = 3,
      h4("1) FCS files"),
      fileInput("fcs_files", NULL, multiple=TRUE, accept=".fcs"),
      tags$p(class="small-note", "Select all .fcs files from your experiment."),
      hr(),

      h4("2) Plate metadata"),
      fileInput("meta_file", NULL, accept=c(".csv",".xlsx",".xls")),
      downloadButton("dl_template", "Download CSV template", class="btn-sm btn-default"),
      tags$p(class="small-note",
        "Required: ", tags$b("well_code"), ", ", tags$b("genotype"), ", ",
        tags$b("perturbation"), ".", br(),
        "Optional: ", tags$b("treatment"), " (e.g. starvation, drug), ",
        tags$b("time"), " (e.g. 4h, 24h).", br(),
        "Extra columns are kept in tables."),
      uiOutput("meta_summary_box"),
      hr(),

      actionButton("run", "Run analysis", class="btn-primary btn-block"),
      tags$p(class="small-note",
             "Set up channels and adjust gates in the tabs, then click Run.")
    ),

    mainPanel(width = 9,
      tabsetPanel(id="tabs",

        # ── Overview ────────────────────────────────────────────────────────────
        tabPanel("Overview",
          h4("Status"),
          uiOutput("status_box"),
          hr(),
          h4("Workflow"),
          tags$ol(
            tags$li("Upload FCS files and metadata in the sidebar."),
            tags$li("Configure channel roles in the ", tags$b("Channels"), " tab."),
            tags$li("Review and adjust the singlet gate in ", tags$b("Gate 1: Singlets"), "."),
            tags$li("Review and adjust the live/dead gate in ", tags$b("Gate 2: Live/Dead"), "."),
            tags$li("Review and adjust the signal gate in ", tags$b("Gate 3: Signal"), "."),
            tags$li("Click ", tags$b("Run analysis"), " in the sidebar to generate results.")
          ),
          hr(),
          h4("Perturbation conditions for Scenith parameters"),
          tags$p("The following perturbation labels trigger parameter calculations — use these exact names in the metadata:"),
          tags$ul(
            tags$li(tags$b("Co"), " — control"),
            tags$li(tags$b("DG"), " — 2-deoxy-D-glucose (glycolysis inhibitor)"),
            tags$li(tags$b("O"),  " — oligomycin A (OXPHOS inhibitor)"),
            tags$li(tags$b("DGO"), " — DG + oligomycin (combined)")
          )
        ),

        # ── Metadata ────────────────────────────────────────────────────────────
        tabPanel("Metadata",
          h4("Uploaded metadata"),
          DTOutput("meta_tbl"),
          hr(),
          h4("FCS ↔ metadata coverage"),
          tags$p(class="small-note", "Green = matched, red = unmatched."),
          DTOutput("coverage_tbl"),
          hr(),
          h4("Perturbation overview"),
          tags$p(class="small-note",
                 "Well counts per perturbation, broken down by available grouping variables."),
          plotOutput("p_perturbation_dist", height="320px"),
          hr(),
          DTOutput("perturbation_count_tbl")
        ),

        # ── Channels ────────────────────────────────────────────────────────────
        tabPanel("Channels",
          h4("Panel / channel configuration"),
          tags$p(class="step-note",
                 "Select a preset for common Scenith panels, or choose 'Custom' ",
                 "to assign channels manually from your FCS file."),
          fluidRow(
            column(5,
              selectInput("preset", "Panel preset",
                          choices = names(PANEL_PRESETS), width="100%"),
              hr(),
              conditionalPanel("input.preset == 'Custom \u2013 configure below'",
                h5("Custom channel assignment"),
                uiOutput("custom_channel_ui")
              ),
              conditionalPanel("input.preset != 'Custom \u2013 configure below'",
                uiOutput("preset_channel_summary")
              ),
              hr(),
              tags$p(class="small-note",
                "You can also upload a panel CSV to document your staining.",
                br(), "Required columns: channel, role (scatter_x, scatter_y, live_dead, signal, other).",
                br(), "This file is for reference — channel roles above take precedence."),
              fileInput("panel_file", "Upload panel CSV (optional)",
                        accept=c(".csv"), width="100%"),
              uiOutput("panel_tbl_ui")
            ),
            column(7,
              h5("Available channels in uploaded FCS files"),
              uiOutput("channel_list_ui")
            )
          )
        ),

        # ── Gate 1: Singlets ────────────────────────────────────────────────────
        tabPanel("Gate 1: Singlets",
          tags$p(class="step-note",
            "Adjust the polygon gate to select single cells on the scatter plot. ",
            "The gate updates the preview immediately. Click Run to apply to the full analysis."),
          fluidRow(
            column(4,
              h5("Polygon vertices"),
              tags$p(class="small-note",
                     "Coordinates in the scatter channel units (e.g. FSC-A and FSC-H)."),
              tags$table(class="gate-tbl", style="width:100%",
                tags$tr(tags$th("Vertex"), tags$th("X"), tags$th("Y")),
                tags$tr(
                  tags$td("Bottom-left"),
                  tags$td(numericInput("g1_x1", NULL, DEFAULT_VERTICES$x1, step=1000, width="100%")),
                  tags$td(numericInput("g1_y1", NULL, DEFAULT_VERTICES$y1, step=1000, width="100%"))
                ),
                tags$tr(
                  tags$td("Top-left"),
                  tags$td(numericInput("g1_x2", NULL, DEFAULT_VERTICES$x2, step=1000, width="100%")),
                  tags$td(numericInput("g1_y2", NULL, DEFAULT_VERTICES$y2, step=1000, width="100%"))
                ),
                tags$tr(
                  tags$td("Top-right"),
                  tags$td(numericInput("g1_x3", NULL, DEFAULT_VERTICES$x3, step=1000, width="100%")),
                  tags$td(numericInput("g1_y3", NULL, DEFAULT_VERTICES$y3, step=1000, width="100%"))
                ),
                tags$tr(
                  tags$td("Bottom-right"),
                  tags$td(numericInput("g1_x4", NULL, DEFAULT_VERTICES$x4, step=1000, width="100%")),
                  tags$td(numericInput("g1_y4", NULL, DEFAULT_VERTICES$y4, step=1000, width="100%"))
                )
              ),
              actionButton("g1_reset", "Reset to defaults", class="btn-sm btn-default",
                           style="margin-top:8px;"),
              hr(),
              h5("Preview sample"),
              uiOutput("g1_sample_ui"),
              hr(),
              h5("Gate statistics"),
              uiOutput("g1_stats_ui")
            ),
            column(8,
              plotOutput("p_g1_preview", height="450px")
            )
          )
        ),

        # ── Gate 2: Live/Dead ────────────────────────────────────────────────────
        tabPanel("Gate 2: Live/Dead",
          tags$p(class="step-note",
            "Set a threshold on the live/dead channel. Cells below the threshold are kept as live. ",
            "Skipped automatically if no live/dead channel is configured."),
          uiOutput("g2_content")
        ),

        # ── Gate 3: Signal ───────────────────────────────────────────────────────
        tabPanel("Gate 3: Signal",
          tags$p(class="step-note",
            "Set the minimum signal threshold (e.g. puromycin APC-A). ",
            "Only cells at or above this value are used in the Scenith analysis."),
          fluidRow(
            column(4,
              h5("Signal threshold"),
              numericInput("g3_threshold", "Keep cells \u2265", value=80, min=0, step=10, width="100%"),
              hr(),
              h5("Preview samples"),
              uiOutput("g3_sample_ui"),
              hr(),
              h5("Gate statistics"),
              uiOutput("g3_stats_ui")
            ),
            column(8,
              plotOutput("p_g3_preview", height="450px")
            )
          )
        ),

        # ── QC: Cell counts ─────────────────────────────────────────────────────
        tabPanel("Cell counts",
          h4("Cells retained after each gating step"),
          plotOutput("p_cell_counts", height=560),
          downloadButton("dl_cell_counts", "Download PNG"),
          hr(),
          DTOutput("cell_counts_tbl")
        ),

        # ── Summary ──────────────────────────────────────────────────────────────
        tabPanel("Summary",
          h4("Mean signal per sample"),
          plotOutput("p_mean_both", height=720),
          downloadButton("dl_mean_puro", "Download PNG"),
          hr(),
          DTOutput("puro_summary_tbl")
        ),

        # ── Scenith parameters ──────────────────────────────────────────────────
        tabPanel("Scenith parameters",
          h4("Derived Scenith metabolic parameters"),
          DTOutput("scenith_tbl"),
          hr(),
          h4("Glucose dependence and FAO/AAO capacity (Co vs DG vs DGO)"),
          plotOutput("p_co_dg_dgo", height=650),
          downloadButton("dl_co_dg_dgo", "Download PNG"),
          hr(),
          h4("Mitochondrial dependence and glycolytic capacity (Co vs O vs DGO)"),
          plotOutput("p_co_o_dgo", height=650),
          downloadButton("dl_co_o_dgo", "Download PNG"),
          hr(),
          h4("Geometric mean signal per genotype \u00d7 perturbation"),
          plotOutput("p_puro_bar", height=420),
          downloadButton("dl_puro_bar", "Download PNG")
        )
      )
    )
  )
)

# ── Server ─────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Metadata ──────────────────────────────────────────────────────────────────
  meta_df <- reactive({
    req(input$meta_file)
    ext <- tolower(tools::file_ext(input$meta_file$name))
    df  <- tryCatch({
      if (ext %in% c("xlsx","xls")) readxl::read_excel(input$meta_file$datapath)
      else                          readr::read_csv(input$meta_file$datapath, show_col_types=FALSE)
    }, error = function(e) NULL)
    validate(need(!is.null(df), "Could not read metadata. Upload a valid CSV or XLSX."))
    required <- c("well_code","genotype","perturbation")
    missing  <- setdiff(required, colnames(df))
    validate(need(length(missing)==0,
      paste0("Metadata missing: ", paste(missing, collapse=", "),
             ". Download the template from the sidebar.")))
    df %>% mutate(well_code = normalize_well_code(well_code))
  })

  output$meta_summary_box <- renderUI({
    if (is.null(input$meta_file)) return(NULL)
    tryCatch({
      df     <- meta_df()
      extras <- setdiff(colnames(df), c("well_code","genotype","perturbation"))
      tags$div(class="tag-box",
        tags$span(class="ok", paste0("\u2713 ", nrow(df), " wells loaded")), br(),
        tags$span("Genotypes: ",    tags$b(paste(sort(unique(na.omit(df$genotype))),     collapse=", "))), br(),
        tags$span("Perturbations: ",tags$b(paste(sort(unique(na.omit(df$perturbation))), collapse=", "))),
        if ("treatment" %in% extras) tagList(br(), tags$span("Treatments: ", tags$b(paste(sort(unique(na.omit(df$treatment))), collapse=", ")))),
        if ("time"      %in% extras) tagList(br(), tags$span("Timepoints: ", tags$b(paste(sort(unique(na.omit(df$time))),      collapse=", "))))
      )
    }, error = function(e) tags$p(class="busy", "\u26a0 Metadata error \u2014 check column names."))
  })

  output$dl_template <- downloadHandler(
    filename = "scenith_metadata_template.csv",
    content  = function(file) {
      expand.grid(well_row=LETTERS[2:8], well_col=1:12, stringsAsFactors=FALSE) %>%
        arrange(well_row, well_col) %>%
        mutate(well_code    = sprintf("%s%02d", well_row, well_col),
               genotype     = "",
               perturbation = "",
               treatment    = "",
               time         = "") %>%
        select(well_code, genotype, perturbation, treatment, time) %>%
        readr::write_csv(file)
    }
  )

  # ── FCS loading ───────────────────────────────────────────────────────────────
  fcs_raw <- reactive({
    req(input$fcs_files)
    withProgress(message="Loading FCS files\u2026", {
      fs <- read.flowSet(files=input$fcs_files$datapath,
                         alter.names=TRUE, truncate_max_range=FALSE)
      sampleNames(fs) <- input$fcs_files$name
      fs
    })
  })

  avail_channels <- reactive({
    req(fcs_raw())
    colnames(fcs_raw())
  })

  # ── Channel configuration ─────────────────────────────────────────────────────
  channel_cfg <- reactive({
    preset_name <- req(input$preset)
    preset      <- PANEL_PRESETS[[preset_name]]
    if (preset_name != "Custom \u2013 configure below") return(preset)
    list(
      scatter_x = if (!is.null(input$ch_scatter_x)) input$ch_scatter_x else "FSC.A",
      scatter_y = if (!is.null(input$ch_scatter_y)) input$ch_scatter_y else "FSC.H",
      live_dead = if (!is.null(input$ch_live_dead) && input$ch_live_dead != "None") input$ch_live_dead else NULL,
      signal    = if (!is.null(input$ch_signal))    input$ch_signal    else "APC.A"
    )
  })

  output$custom_channel_ui <- renderUI({
    chs <- c("None", avail_channels())
    tagList(
      selectInput("ch_scatter_x", "Scatter X (singlet gate X-axis)", choices=chs, selected=chs[chs=="FSC.A"][1]),
      selectInput("ch_scatter_y", "Scatter Y (singlet gate Y-axis)", choices=chs, selected=chs[chs=="FSC.H"][1]),
      selectInput("ch_live_dead", "Live/Dead channel (or None)",     choices=chs, selected="None"),
      selectInput("ch_signal",    "Signal channel (puromycin etc.)", choices=chs[chs!="None"])
    )
  })

  output$preset_channel_summary <- renderUI({
    cfg <- PANEL_PRESETS[[input$preset]]
    if (is.null(cfg)) return(NULL)
    tags$div(class="tag-box",
      tags$span(class="ok", "\u2713 Preset channels:"), br(),
      tags$span("Scatter X: ",  tags$b(cfg$scatter_x)), br(),
      tags$span("Scatter Y: ",  tags$b(cfg$scatter_y)), br(),
      tags$span("Live/Dead: ",  tags$b(if(is.null(cfg$live_dead)) "None" else cfg$live_dead)), br(),
      tags$span("Signal: ",     tags$b(cfg$signal))
    )
  })

  output$channel_list_ui <- renderUI({
    if (is.null(input$fcs_files))
      return(tags$p(class="small-note", "Upload FCS files to see available channels."))
    chs <- avail_channels()
    tags$div(class="tag-box",
      tags$b(paste0(length(chs), " channels detected:")), br(),
      tags$code(paste(chs, collapse="   "))
    )
  })

  output$panel_tbl_ui <- renderUI({
    if (is.null(input$panel_file)) return(NULL)
    df <- tryCatch(readr::read_csv(input$panel_file$datapath, show_col_types=FALSE), error=function(e) NULL)
    if (is.null(df)) return(tags$p(class="busy", "Could not read panel file."))
    tagList(hr(), h5("Uploaded panel"), renderDT(datatable(df, options=list(pageLength=15, dom="t"))))
  })

  # ── Gating helpers ─────────────────────────────────────────────────────────────
  cur_singlet_gate <- reactive({
    cfg <- channel_cfg()
    req(cfg$scatter_x, cfg$scatter_y)
    make_singlet_gate(
      input$g1_x1, input$g1_y1, input$g1_x2, input$g1_y2,
      input$g1_x3, input$g1_y3, input$g1_x4, input$g1_y4,
      cfg$scatter_x, cfg$scatter_y
    )
  })

  observeEvent(input$g1_reset, {
    updateNumericInput(session, "g1_x1", value=DEFAULT_VERTICES$x1)
    updateNumericInput(session, "g1_y1", value=DEFAULT_VERTICES$y1)
    updateNumericInput(session, "g1_x2", value=DEFAULT_VERTICES$x2)
    updateNumericInput(session, "g1_y2", value=DEFAULT_VERTICES$y2)
    updateNumericInput(session, "g1_x3", value=DEFAULT_VERTICES$x3)
    updateNumericInput(session, "g1_y3", value=DEFAULT_VERTICES$y3)
    updateNumericInput(session, "g1_x4", value=DEFAULT_VERTICES$x4)
    updateNumericInput(session, "g1_y4", value=DEFAULT_VERTICES$y4)
  })

  # ── Sample selectors ──────────────────────────────────────────────────────────
  output$g1_sample_ui <- renderUI({
    req(fcs_raw())
    selectInput("g1_sample", NULL, choices=sampleNames(fcs_raw()), width="100%")
  })
  output$g2_sample_ui <- renderUI({
    req(fcs_raw())
    selectInput("g2_sample", NULL, choices=sampleNames(fcs_raw()), width="100%")
  })
  output$g3_sample_ui <- renderUI({
    req(fcs_raw())
    selectInput("g3_samples", NULL, choices=sampleNames(fcs_raw()),
                multiple=TRUE, selected=sampleNames(fcs_raw())[1:min(4,length(fcs_raw()))],
                width="100%")
  })

  # ── Gate 1 preview ────────────────────────────────────────────────────────────
  output$p_g1_preview <- renderPlot({
    req(fcs_raw(), input$g1_sample)
    cfg <- channel_cfg(); req(cfg$scatter_x, cfg$scatter_y)
    pg  <- cur_singlet_gate()
    idx <- which(sampleNames(fcs_raw()) == input$g1_sample)
    if (!length(idx)) return(empty_gg("Sample not found."))
    ggcyto(fcs_raw()[idx], aes_string(x=cfg$scatter_x, y=cfg$scatter_y)) +
      geom_hex(bins=60) +
      geom_gate(pg, colour="red", size=0.7) +
      geom_stats() +
      scale_fill_viridis_c(option="magma") +
      theme_alba() +
      labs(title    = "Singlet gate preview",
           subtitle = paste0("Sample: ", input$g1_sample),
           x=cfg$scatter_x, y=cfg$scatter_y, fill="Count")
  })

  output$g1_stats_ui <- renderUI({
    req(fcs_raw(), input$g1_sample)
    tryCatch({
      cfg <- channel_cfg(); req(cfg$scatter_x, cfg$scatter_y)
      pg  <- cur_singlet_gate()
      idx <- which(sampleNames(fcs_raw()) == input$g1_sample)
      fs_sub <- Subset(fcs_raw()[idx], pg)
      n_raw  <- nrow(exprs(fcs_raw()[[idx]]))
      n_keep <- nrow(exprs(fs_sub[[1]]))
      tags$div(class="tag-box",
        tags$span(paste0("Raw events: ",  n_raw)),  br(),
        tags$span(paste0("After gate: ",  n_keep)), br(),
        tags$span(paste0("Retained: ",    round(100*n_keep/n_raw, 1), "%"))
      )
    }, error = function(e) tags$p(class="small-note", "Run FCS upload to see stats."))
  })

  # ── Gate 2: Live/Dead ─────────────────────────────────────────────────────────
  output$g2_content <- renderUI({
    cfg <- channel_cfg()
    if (is.null(cfg$live_dead)) {
      return(tags$div(style="margin-top:20px;",
        tags$p(class="ok", "\u2713 No live/dead channel configured \u2014 this step is skipped."),
        tags$p(class="small-note", "To enable, select a channel in the Channels tab.")))
    }
    fluidRow(
      column(4,
        h5("Live/Dead threshold"),
        tags$p(class="small-note",
               paste0("Channel: ", cfg$live_dead,
                      ". Cells BELOW the threshold are kept as live.")),
        numericInput("g2_threshold",
                     paste0("Keep cells < (", cfg$live_dead, ")"),
                     value=4000, min=0, step=500, width="100%"),
        hr(),
        h5("Preview sample"),
        uiOutput("g2_sample_ui"),
        hr(),
        h5("Gate statistics"),
        uiOutput("g2_stats_ui")
      ),
      column(8,
        plotOutput("p_g2_preview", height="450px")
      )
    )
  })

  output$p_g2_preview <- renderPlot({
    cfg <- channel_cfg()
    if (is.null(cfg$live_dead)) return(invisible(NULL))
    req(fcs_raw(), input$g2_sample, input$g2_threshold)
    pg  <- cur_singlet_gate()
    idx <- which(sampleNames(fcs_raw()) == input$g2_sample)
    if (!length(idx)) return(empty_gg("Sample not found."))
    fs_sing <- tryCatch(Subset(fcs_raw()[idx], pg), error=function(e) NULL)
    if (is.null(fs_sing)) return(empty_gg("Could not apply singlet gate."))
    thresh   <- input$g2_threshold
    fitc_max <- max(exprs(fs_sing[[1]])[, cfg$live_dead], na.rm=TRUE)
    ggcyto(fs_sing, aes_string(x=cfg$live_dead)) +
      annotate("rect", xmin=-Inf,   xmax=thresh,   ymin=-Inf, ymax=Inf, fill="#1976D2", alpha=.15) +
      annotate("rect", xmin=thresh, xmax=fitc_max, ymin=-Inf, ymax=Inf, fill="red",     alpha=.10) +
      geom_density(fill="#1976D2", alpha=.6) +
      geom_vline(xintercept=thresh, colour="red", linetype="dashed", linewidth=.8) +
      scale_x_log10() + theme_alba() +
      labs(title    = "Live/Dead gate preview",
           subtitle = paste0("Blue = live (kept)   |   Red = dead (discarded)   |   threshold = ", thresh),
           x=paste0(cfg$live_dead, " (log10)"), y="Density")
  })

  output$g2_stats_ui <- renderUI({
    cfg <- channel_cfg()
    if (is.null(cfg$live_dead) || is.null(input$g2_sample) || is.null(input$g2_threshold)) return(NULL)
    tryCatch({
      pg      <- cur_singlet_gate()
      idx     <- which(sampleNames(fcs_raw()) == input$g2_sample)
      fs_sing <- Subset(fcs_raw()[idx], pg)
      lg      <- make_thresh_gate(cfg$live_dead, hi=input$g2_threshold, id="Live")
      fs_live <- Subset(fs_sing, lg)
      n_sing  <- nrow(exprs(fs_sing[[1]]))
      n_live  <- nrow(exprs(fs_live[[1]]))
      tags$div(class="tag-box",
        tags$span(paste0("After singlet gate: ", n_sing)), br(),
        tags$span(paste0("After live gate: ",    n_live)), br(),
        tags$span(paste0("Retained: ", round(100*n_live/n_sing, 1), "%"))
      )
    }, error=function(e) tags$p(class="small-note", "Preview not yet available."))
  })

  # ── Gate 3: Signal ────────────────────────────────────────────────────────────
  output$p_g3_preview <- renderPlot({
    cfg <- channel_cfg(); req(cfg$scatter_x, cfg$signal)
    req(fcs_raw(), input$g3_samples, input$g3_threshold)
    pg  <- cur_singlet_gate()
    idx <- which(sampleNames(fcs_raw()) %in% input$g3_samples)
    if (!length(idx)) return(empty_gg("No samples selected."))
    fs_sing <- tryCatch(Subset(fcs_raw()[idx], pg), error=function(e) NULL)
    if (is.null(fs_sing)) return(empty_gg("Could not apply singlet gate."))
    fs_live <- if (!is.null(cfg$live_dead) && !is.null(input$g2_threshold)) {
      lg <- make_thresh_gate(cfg$live_dead, hi=input$g2_threshold, id="Live")
      tryCatch(Subset(fs_sing, lg), error=function(e) fs_sing)
    } else fs_sing
    thresh <- input$g3_threshold
    ggcyto(fs_live, aes_string(x=cfg$scatter_x, y=cfg$signal)) +
      geom_hex(bins=60) +
      scale_y_log10(limits=c(1, NA)) +
      geom_hline(yintercept=thresh, colour="red", linetype="dashed", linewidth=.8) +
      scale_fill_viridis_c(option="magma") +
      facet_wrap(~ name) +
      theme_alba() +
      labs(title    = "Signal gate preview",
           subtitle = paste0("Dashed line = threshold (", cfg$signal, " \u2265 ", thresh, ")"),
           x=cfg$scatter_x, y=paste0(cfg$signal, " (log10)"), fill="Count")
  })

  output$g3_stats_ui <- renderUI({
    cfg <- channel_cfg()
    if (is.null(cfg$signal) || is.null(input$g3_samples) || !length(input$g3_samples)) return(NULL)
    tryCatch({
      pg      <- cur_singlet_gate()
      idx     <- which(sampleNames(fcs_raw()) == input$g3_samples[1])
      fs_sing <- Subset(fcs_raw()[idx], pg)
      fs_live <- if (!is.null(cfg$live_dead) && !is.null(input$g2_threshold)) {
        lg <- make_thresh_gate(cfg$live_dead, hi=input$g2_threshold, id="Live")
        tryCatch(Subset(fs_sing, lg), error=function(e) fs_sing)
      } else fs_sing
      sg     <- make_thresh_gate(cfg$signal, lo=input$g3_threshold, id="Signal")
      fs_sig <- tryCatch(Subset(fs_live, sg), error=function(e) NULL)
      n_live <- nrow(exprs(fs_live[[1]]))
      n_sig  <- if (!is.null(fs_sig)) nrow(exprs(fs_sig[[1]])) else NA
      tags$div(class="tag-box",
        tags$span(paste0("After live gate (sample 1): ", n_live)), br(),
        tags$span(paste0("After signal gate: ",          n_sig)),  br(),
        tags$span(paste0("Retained: ", round(100*n_sig/n_live, 1), "%"))
      )
    }, error=function(e) tags$p(class="small-note", "Preview not yet available."))
  })

  # ── Main analysis ─────────────────────────────────────────────────────────────
  analysis <- eventReactive(input$run, {
    req(fcs_raw(), meta_df())
    meta <- meta_df()
    cfg  <- channel_cfg()
    validate(need(!is.null(cfg$scatter_x), "Configure channels before running."))
    validate(need(!is.null(cfg$signal),    "Signal channel not configured."))

    withProgress(message="Running analysis\u2026", value=0, {

      setProgress(.1, detail="Applying gates\u2026")
      fs      <- fcs_raw()
      pg      <- cur_singlet_gate()
      fs_sing <- Subset(fs, pg)

      fs_live <- if (!is.null(cfg$live_dead) && !is.null(input$g2_threshold)) {
        lg <- make_thresh_gate(cfg$live_dead, hi=input$g2_threshold, id="Live")
        Subset(fs_sing, lg)
      } else fs_sing

      sg     <- make_thresh_gate(cfg$signal, lo=input$g3_threshold, id="Signal")
      fs_sig <- Subset(fs_live, sg)

      sample_map <- tibble(sample=sampleNames(fs)) %>%
        mutate(well_code = normalize_well_code(str_extract(sample, "[A-H]\\d{1,2}"))) %>%
        left_join(meta, by="well_code")

      p_cols <- auto_palette(meta$perturbation, perturbation_cols_default)
      g_cols <- auto_palette(meta$genotype,     genotype_cols_default)

      setProgress(.3, detail="Building cell counts\u2026")
      cell_counts <- tibble(
        sample         = sampleNames(fs),
        n_raw          = map_int(seq_along(fs),      ~ nrow(exprs(fs[[.x]]))),
        n_singlets     = map_int(seq_along(fs_sing), ~ nrow(exprs(fs_sing[[.x]]))),
        n_live_singlet = map_int(seq_along(fs_live), ~ nrow(exprs(fs_live[[.x]])))
      ) %>% left_join(sample_map, by="sample")

      setProgress(.5, detail="Summarising signal\u2026")
      puro_summary <- map_df(seq_along(sampleNames(fs_live)), function(i) {
        ff_live <- fs_live[[i]]
        ff_sig  <- fs_sig[[i]]
        tibble(
          sample          = sampleNames(fs_live)[i],
          n_live          = nrow(exprs(ff_live)),
          n_signal        = nrow(exprs(ff_sig)),
          pct_signal      = ifelse(nrow(exprs(ff_live))==0, NA_real_,
                                   100*nrow(exprs(ff_sig))/nrow(exprs(ff_live))),
          mean_signal_all = mean(exprs(ff_live)[, cfg$signal], na.rm=TRUE),
          mean_signal_pos = mean(exprs(ff_sig)[,  cfg$signal], na.rm=TRUE)
        )
      }) %>% left_join(sample_map, by="sample")

      setProgress(.65, detail="Building cell-level frame\u2026")
      sig_thresh   <- input$g3_threshold
      genotype_lvl <- sort(unique(na.omit(meta$genotype)))

      cell_level <- map_df(seq_along(sampleNames(fs_live)), function(i) {
        as_tibble(exprs(fs_live[[i]])) %>% mutate(sample=sampleNames(fs_live)[i])
      }) %>% left_join(sample_map, by="sample")

      cell_filtered <- cell_level %>%
        mutate(genotype = factor(genotype, levels=genotype_lvl)) %>%
        filter(.data[[cfg$signal]] >= sig_thresh)

      grp_vars <- "genotype"
      if ("treatment" %in% colnames(meta) &&
          n_distinct(na.omit(meta$treatment)) > 1) grp_vars <- c(grp_vars, "treatment")
      if ("time"      %in% colnames(meta) &&
          n_distinct(na.omit(meta$time))      > 1) grp_vars <- c(grp_vars, "time")

      setProgress(.8, detail="Computing Scenith parameters\u2026")
      geo_means <- cell_level %>% #was cell_filtered
        mutate(genotype = factor(genotype, levels=genotype_lvl)) %>%
        filter(perturbation %in% c("Co","DG","O","DGO")) %>%
        group_by(across(all_of(c(grp_vars, "perturbation")))) %>%
        summarise(geo_mean = geomfi(.data[[cfg$signal]]), .groups="drop")

      has_dg <- all(c("Co","DG","DGO") %in% unique(geo_means$perturbation))
      has_o  <- all(c("Co","O","DGO")  %in% unique(geo_means$perturbation))

      scenith_dg <- if (has_dg) {
        geo_means %>% filter(perturbation %in% c("Co","DG","DGO")) %>%
          pivot_wider(names_from=perturbation, values_from=geo_mean) %>%
          mutate(glucose_dependence = 100*((Co-DG)/(Co-DGO)),
                 fao_aao_capacity   = 100-glucose_dependence)
      } else NULL

      scenith_o <- if (has_o) {
        geo_means %>% filter(perturbation %in% c("Co","O","DGO")) %>%
          pivot_wider(names_from=perturbation, values_from=geo_mean) %>%
          mutate(mito_dependence     = 100*((Co-O)/(Co-DGO)),
                 glycolytic_capacity = 100-mito_dependence)
      } else NULL

      scenith_summary <- if (!is.null(scenith_dg) && !is.null(scenith_o)) {
        full_join(scenith_dg, scenith_o, by=grp_vars, suffix=c("_dg","_o")) %>%
          select(all_of(grp_vars), Co_dg, DG, DGO_dg, glucose_dependence, fao_aao_capacity,
                 Co_o, O, DGO_o, mito_dependence, glycolytic_capacity)
      } else tibble(note="Co, DG, O and DGO perturbations required for Scenith parameters.")

      setProgress(.9, detail="Building plots\u2026")

      p_cell_counts <- cell_counts %>%
        ggplot(aes(x=reorder(sample, n_live_singlet), y=n_live_singlet, fill=genotype)) +
        geom_col() + coord_flip() +
        scale_fill_manual(values=g_cols, na.value="grey80") + theme_alba() +
        labs(title="Cells retained after singlet + live gating",
             x="Sample", y="# Live singlet cells", fill="Genotype")

      p_mean_all <- ggplot(puro_summary,
                           aes(x=reorder(sample, mean_signal_all),
                               y=mean_signal_all, fill=perturbation)) +
        geom_col() + coord_flip() +
        scale_fill_manual(values=p_cols, na.value="grey80") + theme_alba() +
        labs(title="Mean signal \u2013 all live cells",
             x="Sample", y=paste0("Mean ", cfg$signal), fill="Perturbation")

      p_mean_pos <- ggplot(puro_summary,
                           aes(x=reorder(sample, mean_signal_pos),
                               y=mean_signal_pos, fill=perturbation)) +
        geom_col() + coord_flip() +
        scale_fill_manual(values=p_cols, na.value="grey80") + theme_alba() +
        labs(title="Mean signal \u2013 signal-positive cells",
             x="Sample", y=paste0("Mean ", cfg$signal, " (signal+)"), fill="Perturbation")

      p_mean_both <- ggpubr::ggarrange(p_mean_all, p_mean_pos, ncol=1)

      # ── Scenith density plots ──────────────────────────────────────────────────
      p_co_dg_dgo <- if (!is.null(scenith_dg)) {
        tryCatch({
          cell_filtered %>% filter(perturbation %in% c("Co","DG","DGO")) %>%
            ggplot(aes(x=.data[[cfg$signal]], fill=perturbation)) +
            geom_density(alpha=.6, aes(color=perturbation)) +
            scale_x_log10(limits=c(1,10e6)) +
            facet_wrap(~ genotype, ncol=2) +
            scale_fill_manual(values=p_cols) +
            scale_colour_manual(values=p_cols, guide="none") +
            geom_vline(
              data        = geo_means %>% filter(perturbation %in% c("Co","DG","DGO")),
              aes(xintercept=geo_mean, colour=perturbation),
              linetype="dashed", linewidth=.6, inherit.aes=FALSE, show.legend=FALSE
            ) +
            geom_segment(
              data        = scenith_dg %>% filter(is.finite(glucose_dependence)),
              aes(x=DG, xend=Co, y=.5, yend=.5),
              arrow=arrow(length=unit(.15,"cm")), inherit.aes=FALSE
            ) +
            geom_label(
              data        = scenith_dg %>% filter(is.finite(glucose_dependence)),
              aes(x=exp((log(DG)+log(Co))/2), y=.2,
                  label=paste0("1. Glc dep = ",round(glucose_dependence,1),"%")),
              size=3.5, inherit.aes=FALSE
            ) +
            geom_segment(
              data        = scenith_dg %>% filter(is.finite(fao_aao_capacity)),
              aes(x=DGO, xend=DG, y=1, yend=1),
              arrow=arrow(length=unit(.15,"cm")),
              inherit.aes=FALSE, colour="purple4"
            ) +
            geom_label(
              data        = scenith_dg %>% filter(is.finite(fao_aao_capacity)),
              aes(x=exp((log(DGO)+log(DG))/2), y=.8,
                  label=paste0("4. FAO/AAO cap = ",round(fao_aao_capacity,1),"%")),
              size=3.5, inherit.aes=FALSE, colour="purple4"
            ) +
            theme_alba(14) +
            labs(title    = "Glucose dependence and FAO/AAO capacity",
                 subtitle = "Co vs DG vs DGO; arrows = distances between geometric means",
                 x=paste0(cfg$signal," (log10)"), y="Density", fill="Perturbation")
        }, error=function(e) empty_gg(paste("Plot error:", conditionMessage(e))))
      } else empty_gg("Co, DG and DGO perturbations required.")

      p_co_o_dgo <- if (!is.null(scenith_o)) {
        tryCatch({
          cell_filtered %>% filter(perturbation %in% c("Co","O","DGO")) %>%
            ggplot(aes(x=.data[[cfg$signal]], fill=perturbation)) +
            geom_density(alpha=.6, aes(color=perturbation)) +
            scale_x_log10(limits=c(1,10e6)) +
            facet_wrap(~ genotype, ncol=2) +
            scale_fill_manual(values=p_cols) +
            scale_colour_manual(values=p_cols, guide="none") +
            geom_vline(
              data        = geo_means %>% filter(perturbation %in% c("Co","O","DGO")),
              aes(xintercept=geo_mean, colour=perturbation),
              linetype="dashed", linewidth=.6, inherit.aes=FALSE, show.legend=FALSE
            ) +
            geom_segment(
              data        = scenith_o %>% filter(is.finite(mito_dependence)),
              aes(x=O, xend=Co, y=.5, yend=.5),
              arrow=arrow(length=unit(.15,"cm")),
              inherit.aes=FALSE, colour="blue"
            ) +
            geom_label(
              data        = scenith_o %>% filter(is.finite(mito_dependence)),
              aes(x=exp((log(O)+log(Co))/2), y=.2,
                  label=paste0("2. Mito dep = ",round(mito_dependence,1),"%")),
              size=3.5, inherit.aes=FALSE, colour="blue"
            ) +
            geom_segment(
              data        = scenith_o %>% filter(is.finite(glycolytic_capacity)),
              aes(x=DGO, xend=O, y=1, yend=1),
              arrow=arrow(length=unit(.15,"cm")), inherit.aes=FALSE
            ) +
            geom_label(
              data        = scenith_o %>% filter(is.finite(glycolytic_capacity)),
              aes(x=exp((log(DGO)+log(O))/2), y=.8,
                  label=paste0("3. Glyc cap = ",round(glycolytic_capacity,1),"%")),
              size=3.5, inherit.aes=FALSE
            ) +
            theme_alba(14) +
            labs(title    = "Mitochondrial dependence and glycolytic capacity",
                 subtitle = "Co vs O vs DGO; arrows = distances between geometric means",
                 x=paste0(cfg$signal," (log10)"), y="Density", fill="Perturbation")
        }, error=function(e) empty_gg(paste("Plot error:", conditionMessage(e))))
      } else empty_gg("Co, O and DGO perturbations required.")

      bar_data <- geo_means %>%
        filter(perturbation %in% c("Co","DG","O","DGO")) %>%
        mutate(perturbation=factor(perturbation, levels=c("Co","DG","O","DGO")))

      facet_vars <- intersect(c("treatment","time"), grp_vars)
      facet_f    <- if (length(facet_vars))
        as.formula(paste("~", paste(c("genotype", facet_vars), collapse=" + ")))
      else ~ genotype

      p_puro_bar <- ggplot(bar_data, aes(x=perturbation, y=geo_mean, fill=perturbation)) +
        geom_col(width=.7) +
        facet_wrap(facet_f, ncol=4) +
        scale_fill_manual(values=p_cols) +
        theme_alba(14) +
        labs(title    = paste0("Signal per condition (", cfg$signal, ")"),
             subtitle = "Geometric mean in live singlet cells",
             x="Perturbation", y=paste0("Geometric mean ", cfg$signal), fill="Perturbation")

      setProgress(1)

      list(
        sample_map=sample_map, cell_counts=cell_counts,
        puro_summary=puro_summary, scenith_summary=scenith_summary,
        g_cols=g_cols, p_cols=p_cols,
        plots=list(
          p_cell_counts=p_cell_counts, p_mean_both=p_mean_both,
          p_co_dg_dgo=p_co_dg_dgo, p_co_o_dgo=p_co_o_dgo, p_puro_bar=p_puro_bar
        )
      )
    })
  }, ignoreInit=TRUE)

  # ── Status box ────────────────────────────────────────────────────────────────
  output$status_box <- renderUI({
    n_fcs    <- if (!is.null(input$fcs_files)) nrow(input$fcs_files) else 0L
    has_meta <- !is.null(input$meta_file)
    tags$ul(
      tags$li(if (n_fcs>0) tags$span(class="ok",   paste0("\u2713 ", n_fcs, " FCS files loaded"))
              else         tags$span(class="busy",  "\u2b06 Upload FCS files")),
      tags$li(if (has_meta) tags$span(class="ok",  "\u2713 Metadata loaded")
              else          tags$span(class="busy", "\u2b06 Upload metadata")),
      tags$li(tags$span(class="small-note", paste0("Panel: ", input$preset))),
      tags$li(tags$span(class="small-note",
                paste0("Singlet vertices: (",
                       paste(c(input$g1_x1,input$g1_y1,input$g1_x2,input$g1_y2,
                               input$g1_x3,input$g1_y3,input$g1_x4,input$g1_y4),
                             collapse=", "), ")"))),
      tags$li(tags$span(class="small-note",
                paste0("Live/Dead threshold: ",
                       if (!is.null(input$g2_threshold)) input$g2_threshold else "n/a"))),
      tags$li(tags$span(class="small-note", paste0("Signal threshold: ", input$g3_threshold)))
    )
  })

  # ── Metadata tab outputs ──────────────────────────────────────────────────────
  output$meta_tbl <- renderDT({
    req(meta_df())
    datatable(meta_df(), options=list(pageLength=15), caption="Uploaded plate metadata")
  })

  output$coverage_tbl <- renderDT({
    req(fcs_raw(), meta_df())
    sm <- tibble(sample=sampleNames(fcs_raw())) %>%
      mutate(well_code=normalize_well_code(str_extract(sample,"[A-H]\\d{1,2}"))) %>%
      left_join(meta_df(), by="well_code") %>%
      mutate(matched=!is.na(genotype) & !is.na(perturbation)) %>%
      select(sample, well_code, genotype, perturbation, any_of(c("treatment","time")), matched)
    datatable(sm, options=list(pageLength=15)) %>%
      formatStyle("matched", backgroundColor=styleEqual(c(TRUE,FALSE), c("#d4edda","#f8d7da")))
  })

  output$p_perturbation_dist <- renderPlot({
    req(meta_df())
    df  <- meta_df()
    pal <- auto_palette(df$perturbation, perturbation_cols_default)

    # facet by grouping columns that are present and have >1 unique value
    facet_vars <- intersect(c("genotype","treatment","time"), colnames(df))
    facet_vars <- facet_vars[
      sapply(facet_vars, function(v) n_distinct(na.omit(df[[v]])) > 1)
    ]

    p <- df %>%
      count(perturbation, across(any_of(facet_vars))) %>%
      ggplot(aes(x=perturbation, y=n, fill=perturbation)) +
      geom_col(width=.7) +
      geom_text(aes(label=n), vjust=-.4, size=3.5, colour="grey30") +
      scale_fill_manual(values=pal, na.value="grey80") +
      theme_alba() +
      theme(legend.position="none") +
      labs(x="Perturbation", y="Number of wells", title=NULL)

    if (length(facet_vars) > 0)
      p <- p + facet_wrap(as.formula(paste("~", paste(facet_vars, collapse="+"))),
                          scales="free_y")
    p
  })

  output$perturbation_count_tbl <- renderDT({
    req(meta_df())
    df       <- meta_df()
    grp_cols <- intersect(c("genotype","perturbation","treatment","time"), colnames(df))
    df %>%
      group_by(across(all_of(grp_cols))) %>%
      summarise(n_wells=n(), .groups="drop") %>%
      arrange(across(all_of(grp_cols))) %>%
      datatable(options=list(pageLength=20, dom="ft"),
                caption="Wells per group")
  })

  # ── Analysis result outputs ───────────────────────────────────────────────────
  output$p_cell_counts <- renderPlot({ req(analysis()); analysis()$plots$p_cell_counts })
  output$p_mean_both   <- renderPlot({ req(analysis()); analysis()$plots$p_mean_both })
  output$p_co_dg_dgo   <- renderPlot({ req(analysis()); analysis()$plots$p_co_dg_dgo })
  output$p_co_o_dgo    <- renderPlot({ req(analysis()); analysis()$plots$p_co_o_dgo })
  output$p_puro_bar    <- renderPlot({ req(analysis()); analysis()$plots$p_puro_bar })

  output$cell_counts_tbl <- renderDT({
    req(analysis())
    datatable(analysis()$cell_counts %>% arrange(genotype, perturbation, sample),
              options=list(pageLength=10), caption="Cell counts per sample")
  })
  output$puro_summary_tbl <- renderDT({
    req(analysis())
    datatable(analysis()$puro_summary %>% arrange(genotype, perturbation, sample),
              options=list(pageLength=10), caption="Signal summary per sample")
  })
  output$scenith_tbl <- renderDT({
    req(analysis())
    datatable(analysis()$scenith_summary, options=list(pageLength=10),
              caption="Scenith-derived parameters")
  })

  # ── Downloads ─────────────────────────────────────────────────────────────────
  output$dl_cell_counts <- downloadHandler("cell_counts.png",
    function(f) save_plot_png(analysis()$plots$p_cell_counts, f, 10, 7))
  output$dl_mean_puro   <- downloadHandler("mean_signal.png",
    function(f) save_plot_png(analysis()$plots$p_mean_both, f, 9, 8))
  output$dl_co_dg_dgo   <- downloadHandler("dist_co_dg_dgo.png",
    function(f) save_plot_png(analysis()$plots$p_co_dg_dgo, f, 12, 8))
  output$dl_co_o_dgo    <- downloadHandler("dist_co_o_dgo.png",
    function(f) save_plot_png(analysis()$plots$p_co_o_dgo, f, 12, 8))
  output$dl_puro_bar    <- downloadHandler("bar_signal.png",
    function(f) save_plot_png(analysis()$plots$p_puro_bar, f, 12, 6))
}

shinyApp(ui, server)
