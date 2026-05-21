# app.R
# Scenith analysis – Shiny app version
# Author: celbaek (converted to Shiny)
#
# Sample-to-group mapping is driven by an uploaded CSV/XLSX metadata file
# instead of hardcoded plate-layout rules. Upload FCS files + metadata, click Run.

suppressPackageStartupMessages({
  library(shiny)
  library(flowCore)
  library(flowViz)
  library(ggcyto)
  library(ggridges)
  library(openCyto)
  library(scales)
  library(tidyverse)
  library(ggpubr)
  library(viridis)
  library(ggbeeswarm)
  library(DT)
  library(stringr)
  library(sp)      # point.in.polygon
  library(grid)    # unit()
  library(readxl)  # XLSX metadata support
  library(readr)   # CSV metadata support
})

# -------------------------
# Default color schemes
# Auto-extended for any new values found in uploaded metadata
# -------------------------
treatment_cols_default <- c(
  "Co"        = "#4D4D4D",
  "DG"        = "#7B1FA2",
  "O"         = "#D32F2F",
  "DGO"       = "#00897B",
  "UNST"      = "#1976D2",
  "DMSO_25uL" = "#9E9E9E"
)

genotype_cols_default <- c(
  "WT"     = "#1B9E77",
  "WT_HSV" = "#7570B3",
  "KO"     = "#D95F02",
  "KO_HSV" = "#E7298A"
)

# -------------------------
# Helper functions
# -------------------------

# Normalize well codes to "A01" format.
# Accepts "B1", "b01", "B01" and returns "B01".
normalize_well_code <- function(w) {
  w        <- toupper(trimws(as.character(w)))
  row_part <- substr(w, 1, 1)
  col_part <- suppressWarnings(as.integer(substr(w, 2, nchar(w))))
  ifelse(is.na(col_part), NA_character_, sprintf("%s%02d", row_part, col_part))
}

# Extend a named color vector with auto-generated colors for unknown values.
auto_palette <- function(values, known_cols) {
  vals     <- unique(na.omit(as.character(values)))
  new_vals <- setdiff(vals, names(known_cols))
  if (length(new_vals) == 0) return(known_cols)
  new_pal  <- setNames(scales::hue_pal(l = 55, c = 70)(length(new_vals)), new_vals)
  c(known_cols, new_pal)
}

theme_fifi <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2),
      plot.subtitle    = element_text(size = base_size),
      plot.caption     = element_text(size = base_size - 2, colour = "grey40"),
      strip.background = element_rect(fill = "grey95"),
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

geomfi <- function(x) {
  x <- x[x > 0]
  exp(mean(log(x), na.rm = TRUE))
}

save_plot_png <- function(plot_obj, file, width = 7, height = 5, dpi = 300) {
  ggsave(filename = file, plot = plot_obj,
         width = width, height = height, dpi = dpi, units = "in")
}

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .small-note  { color: #555; font-size: 0.93em; }
      .busy        { color: #a33; font-weight: 600; }
      .ok          { color: #176; font-weight: 600; }
      .tag-box     { background: #f4f4f4; border-radius: 6px;
                     padding: 8px 12px; font-size: 0.9em; margin-top: 6px; }
    "))
  ),
  titlePanel("Scenith analysis – Shiny app"),

  sidebarLayout(
    sidebarPanel(

      # ---- 1) FCS files -------------------------------------------------------
      h4("1) Load FCS files"),
      fileInput("fcs_files", "Upload FCS files (multiple)",
                multiple = TRUE, accept = ".fcs"),
      tags$p(class = "small-note",
             "Select all .fcs files from your experiment folder."),
      hr(),

      # ---- 2) Metadata --------------------------------------------------------
      h4("2) Load plate metadata"),
      fileInput("meta_file", "Upload metadata (CSV or XLSX)",
                accept = c(".csv", ".xlsx", ".xls")),
      downloadButton("dl_template", "Download CSV template", class = "btn-sm btn-default"),
      tags$p(class = "small-note",
             "Required columns: ", tags$b("well_code"), " (e.g. B01), ",
             tags$b("genotype"), ", ", tags$b("treatment"), ".", br(),
             "Extra columns are carried through and appear in the tables."),
      uiOutput("meta_summary_box"),
      hr(),

      # ---- 3) Gating ----------------------------------------------------------
      h4("3) Gating parameters"),
      numericInput("fitc_threshold",
                   "Live/Dead FITC-A threshold (keep cells BELOW this)",
                   value = 4000, min = 1),
      numericInput("puro_threshold",
                   "Puromycin APC-A threshold (keep cells AT OR ABOVE this)",
                   value = 80, min = 0),
      hr(),

      # ---- 4) QC plot options -------------------------------------------------
      h4("4) QC plot options"),
      numericInput("rep_singlet_idx",
                   "Representative sample index (singlet QC)", value = 1, min = 1),
      numericInput("rep_live_idx",
                   "Representative sample index (live-gate QC)", value = 1, min = 1),
      textInput("rep_puro_idx",
                "Sample indices for puromycin QC (comma-separated)",
                value = "1,2,3,4,5,6"),
      checkboxInput("make_all_fsc_plot",
                    "Render FSC-A vs FSC-H for ALL samples (can be slow)",
                    value = FALSE),
      hr(),

      actionButton("run", "Run / Recompute analysis", class = "btn-primary"),
      tags$p(class = "small-note",
             "Click Run after uploading files or changing any parameter.")
    ),

    mainPanel(
      tabsetPanel(

        # ---- Overview tab -----------------------------------------------------
        tabPanel("Overview",
          h4("Experimental overview"),
          tags$ul(
            tags$li(tags$b("DG"), " – 2-deoxy-D-glucose (glycolysis inhibitor)"),
            tags$li(tags$b("O"), " – oligomycin A (ATP synthase / OXPHOS inhibitor)"),
            tags$li(tags$b("DGO"), " – DG + oligomycin (combined)"),
            tags$li(tags$b("Puromycin"), " – APC-A; translation / ATP-availability readout"),
            tags$li(tags$b("Live/Dead dye"), " – FITC-A")
          ),
          tags$p("Upload FCS files and a metadata file, then click ", tags$b("Run.")),
          hr(),
          uiOutput("status_box")
        ),

        # ---- Metadata tab -----------------------------------------------------
        tabPanel("Metadata",
          h4("Uploaded plate metadata"),
          tags$p(class = "small-note",
                 "Check that all wells you expect are present and correctly spelled."),
          DTOutput("meta_tbl"),
          hr(),
          h4("FCS files × metadata coverage"),
          tags$p(class = "small-note",
                 "Green rows matched successfully; red rows could not be matched."),
          DTOutput("coverage_tbl")
        ),

        # ---- Sample mapping tab -----------------------------------------------
        tabPanel("Sample mapping",
          h4("Sample-to-well mapping (filename → metadata)"),
          DTOutput("sample_map_tbl"),
          uiOutput("unmatched_warn")
        ),

        # ---- Gating QC tab ---------------------------------------------------
        tabPanel("Gating QC",
          h4("Step 1: Singlets gate (FSC-A vs FSC-H)"),
          plotOutput("p_qc_singlets", height = 420),
          downloadButton("dl_qc_singlets", "Download PNG"),

          hr(),
          h4("Step 1 (optional): FSC-A vs FSC-H for ALL samples"),
          uiOutput("all_fsc_notice"),
          plotOutput("p_fsc_all", height = 700),
          downloadButton("dl_fsc_all", "Download PNG"),

          hr(),
          h4("Step 2: Live-cell gate (FITC-A)"),
          plotOutput("p_qc_live", height = 420),
          downloadButton("dl_qc_live", "Download PNG"),

          hr(),
          h4("Step 3: Puromycin gate (APC-A vs FSC-A)"),
          plotOutput("p_qc_puro", height = 520),
          downloadButton("dl_qc_puro", "Download PNG")
        ),

        # ---- Cell counts tab -------------------------------------------------
        tabPanel("Cell counts",
          h4("QC: Cells remaining after singlet + live gating"),
          plotOutput("p_cell_counts", height = 560),
          downloadButton("dl_cell_counts", "Download PNG"),
          hr(),
          DTOutput("cell_counts_tbl")
        ),

        # ---- Puromycin summary tab -------------------------------------------
        tabPanel("Puromycin summary",
          h4("Mean puromycin signal per sample"),
          plotOutput("p_mean_both", height = 720),
          downloadButton("dl_mean_puro", "Download PNG"),
          hr(),
          DTOutput("puro_summary_tbl")
        ),

        # ---- Scenith parameters tab ------------------------------------------
        tabPanel("Scenith parameters",
          h4("Scenith-derived parameters per genotype"),
          DTOutput("scenith_tbl"),
          hr(),
          h4("Glucose dependence and FAO/AAO capacity (Co vs DG vs DGO)"),
          plotOutput("p_co_dg_dgo_params", height = 650),
          downloadButton("dl_co_dg_dgo", "Download PNG"),
          hr(),
          h4("Mitochondrial dependence and glycolytic capacity (Co vs O vs DGO)"),
          plotOutput("p_co_o_dgo_params", height = 650),
          downloadButton("dl_co_o_dgo", "Download PNG"),
          hr(),
          h4("Bar plots: geometric mean APC-A per genotype × treatment"),
          plotOutput("p_puro_bar", height = 420),
          downloadButton("dl_puro_bar", "Download PNG")
        )
      )
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output, session) {

  # --------------------------------------------------------------------------
  # Metadata reactive: re-parses whenever the uploaded file changes.
  # --------------------------------------------------------------------------
  meta_df <- reactive({
    req(input$meta_file)
    ext <- tolower(tools::file_ext(input$meta_file$name))

    df <- tryCatch({
      if (ext %in% c("xlsx", "xls")) {
        readxl::read_excel(input$meta_file$datapath)
      } else {
        readr::read_csv(input$meta_file$datapath, show_col_types = FALSE)
      }
    }, error = function(e) NULL)

    validate(need(!is.null(df),
      "Could not read the metadata file. Please upload a valid CSV or XLSX."))

    validate(need("well_code" %in% colnames(df),
      paste0("Metadata must contain a 'well_code' column (e.g. B01). ",
             "Download the template from the sidebar.")))

    df <- df %>% mutate(well_code = normalize_well_code(well_code))

    missing_cols <- setdiff(c("genotype", "treatment"), colnames(df))
    validate(need(length(missing_cols) == 0,
      paste0("Metadata is missing required column(s): ",
             paste(missing_cols, collapse = ", "),
             ". Download the template from the sidebar.")))
    df
  })

  # --------------------------------------------------------------------------
  # Template download: full rows B-H × cols 1-12 with empty group columns.
  # --------------------------------------------------------------------------
  output$dl_template <- downloadHandler(
    filename = "scenith_metadata_template.csv",
    content  = function(file) {
      template <- expand.grid(
        well_row = LETTERS[2:8],
        well_col = 1:12,
        stringsAsFactors = FALSE
      ) %>%
        arrange(well_row, well_col) %>%
        mutate(
          well_code = sprintf("%s%02d", well_row, well_col),
          genotype  = "",
          treatment = ""
        ) %>%
        select(well_code, genotype, treatment)
      readr::write_csv(template, file)
    }
  )

  # --------------------------------------------------------------------------
  # Sidebar metadata summary box.
  # --------------------------------------------------------------------------
  output$meta_summary_box <- renderUI({
    if (is.null(input$meta_file)) return(NULL)
    tryCatch({
      df         <- meta_df()
      n_wells    <- nrow(df)
      genotypes  <- paste(sort(unique(na.omit(df$genotype))),  collapse = ", ")
      treatments <- paste(sort(unique(na.omit(df$treatment))), collapse = ", ")
      tags$div(class = "tag-box",
        tags$span(class = "ok", paste0("✓ ", n_wells, " wells loaded")), br(),
        tags$span("Genotypes: ",  tags$b(genotypes)),  br(),
        tags$span("Treatments: ", tags$b(treatments))
      )
    }, error = function(e) {
      tags$p(class = "busy", "⚠ Metadata error — check column names.")
    })
  })

  # --------------------------------------------------------------------------
  # Main analysis: triggered only by the Run button.
  # --------------------------------------------------------------------------
  results <- eventReactive(input$run, {

    req(input$fcs_files)
    validate(need(nrow(input$fcs_files) > 0, "Please upload FCS files first."))
    validate(need(!is.null(input$meta_file),
                  "Please upload a metadata file (CSV or XLSX)."))

    meta <- meta_df()

    # -- Load flowSet (from browser-temp paths; original names restored)
    df <- read.flowSet(
      files              = input$fcs_files$datapath,
      alter.names        = TRUE,
      truncate_max_range = FALSE
    )
    sampleNames(df) <- input$fcs_files$name

    # -- Sample map: extract well code from filename, join to metadata
    sample_map <- tibble(sample = sampleNames(df)) %>%
      mutate(
        well_code = normalize_well_code(
          str_extract(sample, "[A-H]\\d{1,2}")
        )
      ) %>%
      left_join(meta, by = "well_code")

    # -- Dynamic color palettes (extend defaults with any new values)
    t_cols <- auto_palette(meta$treatment, treatment_cols_default)
    g_cols <- auto_palette(meta$genotype,  genotype_cols_default)

    # -- Gates
    sqrcut <- matrix(
      c(8000, 20000, 5000, 0, 120000, 68000, 120000, 90000),
      ncol = 2, byrow = TRUE
    )
    colnames(sqrcut) <- c("FSC.A", "FSC.H")

    pg              <- polygonGate(filterId = "Singlets", gate = sqrcut)
    singleCell      <- Subset(df, pg)
    fitc_threshold  <- input$fitc_threshold
    rg_fitc         <- rectangleGate("FITC.A" = c(0, fitc_threshold))
    singleCell_live <- Subset(singleCell, rg_fitc)
    puro_threshold  <- input$puro_threshold
    gate_puro       <- rectangleGate("APC.A" = c(puro_threshold, Inf))
    puro_cells      <- Subset(singleCell_live, gate_puro)

    # -- Cell counts
    samples        <- sampleNames(df)
    raw_counts     <- purrr::map_int(seq_along(df),              ~ nrow(exprs(df[[.x]])))
    singlet_counts <- purrr::map_int(seq_along(singleCell),      ~ nrow(exprs(singleCell[[.x]])))
    live_counts    <- purrr::map_int(seq_along(singleCell_live), ~ nrow(exprs(singleCell_live[[.x]])))

    cell_counts <- tibble(
      sample         = samples,
      n_raw          = raw_counts,
      n_singlets     = singlet_counts,
      n_live_singlet = live_counts
    ) %>% left_join(sample_map, by = "sample")

    # -- Puromycin summary per sample
    sample_ids   <- sampleNames(singleCell_live)
    puro_summary <- purrr::map_df(seq_along(sample_ids), function(i) {
      ff_live <- singleCell_live[[i]]
      ff_puro <- puro_cells[[i]]
      tibble(
        sample         = sample_ids[i],
        n_live         = nrow(exprs(ff_live)),
        n_puro         = nrow(exprs(ff_puro)),
        pct_puro       = ifelse(nrow(exprs(ff_live)) == 0, NA_real_,
                                100 * nrow(exprs(ff_puro)) / nrow(exprs(ff_live))),
        mean_puro_live = mean(exprs(ff_live)[, "APC.A"], na.rm = TRUE),
        mean_puro_pos  = mean(exprs(ff_puro)[, "APC.A"], na.rm = TRUE)
      )
    }) %>% left_join(sample_map, by = "sample")

    # -- Cell-level data frame (live singlets)
    sample_ids_live <- sampleNames(singleCell_live)
    genotype_levels <- sort(unique(na.omit(meta$genotype)))

    cell_level <- purrr::map_df(seq_along(sample_ids_live), function(i) {
      ff <- singleCell_live[[i]]
      as_tibble(exprs(ff)) %>% mutate(sample = sample_ids_live[i])
    }) %>% left_join(sample_map, by = "sample")

    cell_level_filtered <- cell_level %>%
      mutate(genotype = factor(genotype, levels = genotype_levels)) %>%
      filter(APC.A >= puro_threshold)

    geo_means <- cell_level_filtered %>%
      mutate(genotype = factor(genotype, levels = genotype_levels)) %>%
      filter(treatment %in% c("Co", "DG", "DGO", "O")) %>%
      group_by(genotype, treatment) %>%
      summarise(geo_mean = geomfi(APC.A), .groups = "drop")

    # -- Scenith parameters (require Co/DG/O/DGO; graceful fallback otherwise)
    has_dg_set <- all(c("Co","DG","DGO") %in% unique(geo_means$treatment))
    has_o_set  <- all(c("Co","O","DGO")  %in% unique(geo_means$treatment))

    scenith_dg <- if (has_dg_set) {
      geo_means %>%
        filter(treatment %in% c("Co","DG","DGO")) %>%
        pivot_wider(names_from = treatment, values_from = geo_mean) %>%
        mutate(
          glucose_dependence = 100 * ((Co - DG) / (Co - DGO)),
          fao_aao_capacity   = 100 - glucose_dependence
        )
    } else NULL

    scenith_o <- if (has_o_set) {
      geo_means %>%
        filter(treatment %in% c("Co","O","DGO")) %>%
        pivot_wider(names_from = treatment, values_from = geo_mean) %>%
        mutate(
          mito_dependence     = 100 * ((Co - O) / (Co - DGO)),
          glycolytic_capacity = 100 - mito_dependence
        )
    } else NULL

    scenith_summary <- if (!is.null(scenith_dg) && !is.null(scenith_o)) {
      full_join(scenith_dg, scenith_o, by = "genotype", suffix = c("_dg","_o")) %>%
        select(genotype,
               Co_dg, DG, DGO_dg, glucose_dependence, fao_aao_capacity,
               Co_o,  O,  DGO_o,  mito_dependence, glycolytic_capacity)
    } else {
      tibble(note = "Scenith parameters require treatments: Co, DG, O and DGO.")
    }

    # ---- QC plots ------------------------------------------------------------
    n_samp      <- length(df)
    rep_singlet <- max(1L, min(as.integer(input$rep_singlet_idx), n_samp))
    rep_live    <- max(1L, min(as.integer(input$rep_live_idx),    length(singleCell)))

    p_qc_singlets <- ggcyto(df[rep_singlet], aes(x = FSC.A, y = FSC.H)) +
      geom_hex(bins = 60) +
      geom_gate(pg, colour = "red", size = 0.6) +
      geom_stats() +
      scale_fill_viridis_c(option = "magma") +
      theme_fifi() +
      labs(title    = "Step 1: Singlet gate on FSC-A vs FSC-H",
           subtitle = paste0("Representative sample: ", sampleNames(df)[rep_singlet]),
           x = "FSC-A", y = "FSC-H", fill = "Cell density")

    # Optional all-samples FSC plot
    p_fsc_all <- NULL
    if (isTRUE(input$make_all_fsc_plot)) {
      all_fsc <- purrr::map_df(seq_along(df), function(i) {
        tibble(FSC.A  = exprs(df[[i]])[,"FSC.A"],
               FSC.H  = exprs(df[[i]])[,"FSC.H"],
               sample = sampleNames(df)[i])
      }) %>% left_join(sample_map, by = "sample")

      pct_df <- purrr::map_df(seq_along(df), function(i) {
        dat    <- as_tibble(exprs(df[[i]])[, c("FSC.A","FSC.H")])
        inside <- sp::point.in.polygon(dat$FSC.A, dat$FSC.H,
                                       sqrcut[,1], sqrcut[,2]) > 0
        tibble(sample = sampleNames(df)[i],
               pct_singlets = 100 * mean(inside), n_raw = nrow(dat))
      })

      gate_df <- rbind(as.data.frame(sqrcut), as.data.frame(sqrcut)[1,])
      colnames(gate_df) <- c("FSC.A","FSC.H")

      all_fsc2 <- all_fsc %>%
        left_join(pct_df, by = "sample") %>%
        mutate(sample_label = paste0(sample, " | ", genotype,
                                     " | ", treatment, " | n=", n_raw))

      p_fsc_all <- ggplot(all_fsc2, aes(x = FSC.A, y = FSC.H)) +
        geom_hex(bins = 35) +
        geom_polygon(data = gate_df, aes(x = FSC.A, y = FSC.H),
                     inherit.aes = FALSE, color = "red", linewidth = 0.4, fill = NA) +
        scale_fill_viridis_c(option = "magma") +
        facet_wrap(~ sample_label, scales = "free", ncol = 8) +
        geom_text(
          data        = dplyr::distinct(all_fsc2, sample_label, pct_singlets),
          aes(x = -Inf, y = Inf,
              label = paste0(round(pct_singlets,1),"% singlets")),
          inherit.aes = FALSE,
          hjust = -0.1, vjust = 1.2, size = 2.5, color = "red", fontface = "bold"
        ) +
        theme_fifi(8) +
        labs(title    = "FSC-A vs FSC-H for all samples (raw data)",
             subtitle = "Hex-binned; red polygon = singlet gate; text = % singlets per sample",
             x = "FSC-A", y = "FSC-H", fill = "Cell density")
    }

    fitc_max  <- max(exprs(singleCell[[rep_live]])[,"FITC.A"], na.rm = TRUE)

    p_qc_live <- ggcyto(singleCell[rep_live], aes(x = FITC.A)) +
      annotate("rect", xmin = 1, xmax = fitc_threshold,
               ymin = -Inf, ymax = Inf, fill = "#1976D2", alpha = 0.15) +
      annotate("rect", xmin = fitc_threshold, xmax = fitc_max,
               ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.10) +
      annotate("segment",
               x    = seq(fitc_threshold, fitc_max, length.out = 20),
               xend = seq(fitc_threshold, fitc_max, length.out = 20) + 500,
               y = -Inf, yend = Inf, colour = "red", alpha = 0.15, linewidth = 0.5) +
      geom_density(alpha = 0.6, fill = "#1976D2") +
      geom_vline(xintercept = fitc_threshold, colour = "red",
                 linetype = "dashed", linewidth = 0.7) +
      scale_x_log10() +
      theme_fifi() +
      labs(title    = "Step 2: Live-cell gate on FITC-A",
           subtitle = paste0("Representative sample: ", sampleNames(singleCell)[rep_live]),
           x = "FITC-A (log10)", y = "Density")

    # Puromycin QC
    parse_idx <- function(txt, max_n) {
      x <- suppressWarnings(as.integer(str_split(txt, ",", simplify = TRUE)))
      x <- x[!is.na(x) & x >= 1 & x <= max_n]
      if (length(x) == 0) x <- 1L
      unique(x)
    }
    plot_list <- parse_idx(input$rep_puro_idx, length(singleCell_live))

    p_qc_puro <- ggcyto(singleCell_live[plot_list], aes(x = FSC.A, y = APC.A)) +
      geom_hex(bins = 60) +
      scale_y_log10(limits = c(1, NA)) +
      geom_hline(yintercept = puro_threshold, color = "red", linetype = "dashed") +
      scale_fill_viridis_c(option = "magma") +
      theme_fifi() +
      facet_wrap(~ name) +
      labs(title    = "Step 3: Puromycin gate on APC-A",
           subtitle = paste0("Dashed line = APC-A threshold (", puro_threshold, ")"),
           x = "FSC-A", y = "APC-A (log10)", fill = "Cell density")

    # Cell counts plot
    p_cell_counts <- cell_counts %>%
      ggplot(aes(x = reorder(sample, n_live_singlet), y = n_live_singlet,
                 fill = genotype)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = g_cols, na.value = "grey80") +
      theme_fifi() +
      labs(title    = "QC: Cells remaining after singlet + live gating",
           subtitle = "Samples ordered by retained live singlets; colored by genotype",
           x = "Sample", y = "# Live singlet cells", fill = "Genotype")

    # Puromycin mean plots
    p_mean_live <- ggplot(puro_summary,
                          aes(x = reorder(sample, mean_puro_live),
                              y = mean_puro_live, fill = treatment)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = t_cols, na.value = "grey80") +
      theme_fifi() +
      labs(title = "Mean puromycin signal (all live cells)",
           x = "Sample", y = "Mean APC-A", fill = "Treatment")

    p_mean_pos <- ggplot(puro_summary,
                         aes(x = reorder(sample, mean_puro_pos),
                             y = mean_puro_pos, fill = treatment)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = t_cols, na.value = "grey80") +
      theme_fifi() +
      labs(title = "Mean puromycin signal (puro+ cells only)",
           x = "Sample", y = "Mean APC-A (puro+)", fill = "Treatment")

    p_mean_both <- ggpubr::ggarrange(p_mean_live, p_mean_pos, ncol = 1)

    # Scenith density plots – graceful fallback if treatments are absent
    empty_plot <- function(msg) {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = msg, size = 5, color = "grey40") +
        theme_void()
    }

    p_co_dg_dgo_params <- if (!is.null(scenith_dg)) {
      cell_level_filtered %>%
        filter(treatment %in% c("Co","DG","DGO")) %>%
        ggplot(aes(x = APC.A, fill = treatment)) +
        geom_density(alpha = 0.6, aes(color = treatment)) +
        scale_x_log10(limits = c(1, 10e6)) +
        facet_wrap(~ genotype, ncol = 2) +
        scale_fill_manual(values = t_cols) +
        scale_colour_manual(values = t_cols, guide = "none") +
        geom_vline(
          data        = geo_means %>% filter(treatment %in% c("Co","DG","DGO")),
          aes(xintercept = geo_mean, colour = treatment),
          linetype = "dashed", linewidth = 0.6,
          inherit.aes = FALSE, show.legend = FALSE
        ) +
        geom_segment(
          data = scenith_dg,
          aes(x = DG, xend = Co, y = 0.5, yend = 0.5),
          arrow = arrow(length = unit(0.15,"cm")), inherit.aes = FALSE
        ) +
        geom_label(
          data = scenith_dg,
          aes(x = (DG + Co) / 0.8, y = 0.2,
              label = paste0("1. Glc dep = ", round(glucose_dependence, 1), "%")),
          size = 4, inherit.aes = FALSE
        ) +
        geom_segment(
          data = scenith_dg,
          aes(x = DGO, xend = DG, y = 1, yend = 1),
          arrow = arrow(length = unit(0.15,"cm")),
          inherit.aes = FALSE, color = "purple4"
        ) +
        geom_label(
          data = scenith_dg,
          aes(x = (DGO + DG) / 8, y = 0.8,
              label = paste0("4. FAO/AAO cap = ", round(fao_aao_capacity, 1), "%")),
          size = 4, inherit.aes = FALSE, color = "purple4"
        ) +
        theme_fifi(14) +
        labs(title    = "Scenith parameters – Glucose dependence and FAO/AAO capacity",
             subtitle = "Densities pooled across replicates; arrows = distances between geometric means",
             x = "APC-A (log10)", y = "Density", fill = "Treatment")
    } else {
      empty_plot("Co, DG and DGO treatments required for this plot.")
    }

    p_co_o_dgo_params <- if (!is.null(scenith_o)) {
      cell_level_filtered %>%
        filter(treatment %in% c("Co","O","DGO")) %>%
        ggplot(aes(x = APC.A, fill = treatment)) +
        geom_density(alpha = 0.6, aes(color = treatment)) +
        scale_x_log10(limits = c(1, 10e6)) +
        facet_wrap(~ genotype, ncol = 2) +
        scale_fill_manual(values = t_cols) +
        scale_colour_manual(values = t_cols, guide = "none") +
        geom_vline(
          data        = geo_means %>% filter(treatment %in% c("Co","O","DGO")),
          aes(xintercept = geo_mean, colour = treatment),
          linetype = "dashed", linewidth = 0.6,
          inherit.aes = FALSE, show.legend = FALSE
        ) +
        geom_segment(
          data = scenith_o,
          aes(x = O, xend = Co, y = 0.5, yend = 0.5),
          arrow = arrow(length = unit(0.15,"cm")),
          inherit.aes = FALSE, color = "blue"
        ) +
        geom_label(
          data = scenith_o,
          aes(x = (O + Co) / 0.8, y = 0.2,
              label = paste0("2. Mito dep = ", round(mito_dependence, 1), "%")),
          size = 4, inherit.aes = FALSE, color = "blue"
        ) +
        geom_segment(
          data = scenith_o,
          aes(x = DGO, xend = O, y = 1, yend = 1),
          arrow = arrow(length = unit(0.15,"cm")), inherit.aes = FALSE
        ) +
        geom_label(
          data = scenith_o,
          aes(x = (DGO + O) / 5, y = 0.8,
              label = paste0("3. Glyc cap = ", round(glycolytic_capacity, 1), "%")),
          size = 4, inherit.aes = FALSE
        ) +
        theme_fifi(14) +
        labs(title    = "Scenith parameters – Mitochondrial dependence and glycolytic capacity",
             subtitle = "Densities pooled across replicates; arrows = distances between geometric means",
             x = "APC-A (log10)", y = "Density", fill = "Treatment")
    } else {
      empty_plot("Co, O and DGO treatments required for this plot.")
    }

    # Bar plot
    bar_puro <- geo_means %>%
      filter(treatment %in% c("Co","DG","O","DGO")) %>%
      mutate(treatment = factor(treatment, levels = c("Co","DG","O","DGO")))

    p_puro_bar <- ggplot(bar_puro,
                         aes(x = treatment, y = geo_mean, fill = treatment)) +
      geom_col(width = 0.7) +
      facet_wrap(~ genotype, ncol = 4) +
      scale_fill_manual(values = t_cols) +
      theme_fifi(14) +
      labs(title    = "Translation per condition and genotype",
           subtitle = "Geometric mean APC-A (puromycin MFI) in live singlet cells",
           x = "Condition", y = "Geometric mean APC-A", fill = "Treatment")

    list(
      df = df, sample_map = sample_map, meta = meta,
      sqrcut = sqrcut, pg = pg,
      singleCell = singleCell, singleCell_live = singleCell_live,
      puro_cells = puro_cells,
      cell_counts = cell_counts, puro_summary = puro_summary,
      geo_means = geo_means, scenith_summary = scenith_summary,
      t_cols = t_cols, g_cols = g_cols,
      plots = list(
        p_qc_singlets      = p_qc_singlets,
        p_fsc_all          = p_fsc_all,
        p_qc_live          = p_qc_live,
        p_qc_puro          = p_qc_puro,
        p_cell_counts      = p_cell_counts,
        p_mean_both        = p_mean_both,
        p_co_dg_dgo_params = p_co_dg_dgo_params,
        p_co_o_dgo_params  = p_co_o_dgo_params,
        p_puro_bar         = p_puro_bar
      )
    )
  }, ignoreInit = TRUE)

  # --------------------------------------------------------------------------
  # Status box (Overview tab)
  # --------------------------------------------------------------------------
  output$status_box <- renderUI({
    n_fcs    <- if (!is.null(input$fcs_files)) nrow(input$fcs_files) else 0L
    has_meta <- !is.null(input$meta_file)
    fcs_msg  <- if (n_fcs == 0) "⬆ Upload FCS files." else
                  paste0("✓ ", n_fcs, " FCS file(s) loaded.")
    meta_msg <- if (!has_meta) "⬆ Upload metadata CSV/XLSX." else
                  "✓ Metadata loaded."
    if (n_fcs == 0 || !has_meta) {
      tags$p(class = "busy", fcs_msg, " ", meta_msg)
    } else {
      tags$p(class = "ok",
             paste0(fcs_msg, "  ", meta_msg,
                    "  → Click 'Run / Recompute analysis'."))
    }
  })

  # --------------------------------------------------------------------------
  # Metadata tab outputs
  # --------------------------------------------------------------------------
  output$meta_tbl <- renderDT({
    req(meta_df())
    datatable(meta_df(), options = list(pageLength = 15),
              caption = "Uploaded plate metadata")
  })

  output$coverage_tbl <- renderDT({
    req(results())
    sm <- results()$sample_map %>%
      mutate(matched = !is.na(genotype) & !is.na(treatment)) %>%
      select(sample, well_code, genotype, treatment, matched) %>%
      arrange(matched, sample)
    datatable(sm, options = list(pageLength = 15),
              caption = "FCS files and their metadata match status") %>%
      formatStyle("matched",
                  backgroundColor = styleEqual(c(TRUE, FALSE),
                                               c("#d4edda", "#f8d7da")))
  })

  # --------------------------------------------------------------------------
  # Sample mapping tab
  # --------------------------------------------------------------------------
  output$sample_map_tbl <- renderDT({
    req(results())
    datatable(results()$sample_map, options = list(pageLength = 10),
              caption = "FCS files mapped to metadata")
  })

  output$unmatched_warn <- renderUI({
    req(results())
    n_unmatched <- sum(is.na(results()$sample_map$genotype) |
                       is.na(results()$sample_map$treatment))
    if (n_unmatched > 0) {
      tags$p(class = "busy",
             paste0("⚠ ", n_unmatched, " sample(s) could not be matched. ",
                    "Check that well codes in filenames match those in your metadata."))
    }
  })

  # --------------------------------------------------------------------------
  # Gating QC plots
  # --------------------------------------------------------------------------
  output$p_qc_singlets <- renderPlot({ req(results()); results()$plots$p_qc_singlets })
  output$p_qc_live     <- renderPlot({ req(results()); results()$plots$p_qc_live })
  output$p_qc_puro     <- renderPlot({ req(results()); results()$plots$p_qc_puro })

  output$all_fsc_notice <- renderUI({
    req(results())
    if (!isTRUE(input$make_all_fsc_plot))
      tags$p(class = "small-note",
             "Enable via the sidebar checkbox to render this plot.")
  })
  output$p_fsc_all <- renderPlot({
    req(results())
    if (!isTRUE(input$make_all_fsc_plot)) return(invisible(NULL))
    results()$plots$p_fsc_all
  })

  # --------------------------------------------------------------------------
  # Cell counts
  # --------------------------------------------------------------------------
  output$p_cell_counts <- renderPlot({ req(results()); results()$plots$p_cell_counts })
  output$cell_counts_tbl <- renderDT({
    req(results())
    datatable(results()$cell_counts %>% arrange(genotype, treatment, sample),
              options = list(pageLength = 10),
              caption = "Cell counts per sample at each gating step")
  })

  # --------------------------------------------------------------------------
  # Puromycin summary
  # --------------------------------------------------------------------------
  output$p_mean_both <- renderPlot({ req(results()); results()$plots$p_mean_both })
  output$puro_summary_tbl <- renderDT({
    req(results())
    datatable(results()$puro_summary %>% arrange(genotype, treatment, sample),
              options = list(pageLength = 10),
              caption = "Puromycin summary per sample")
  })

  # --------------------------------------------------------------------------
  # Scenith parameters
  # --------------------------------------------------------------------------
  output$p_co_dg_dgo_params <- renderPlot({ req(results()); results()$plots$p_co_dg_dgo_params })
  output$p_co_o_dgo_params  <- renderPlot({ req(results()); results()$plots$p_co_o_dgo_params })
  output$p_puro_bar         <- renderPlot({ req(results()); results()$plots$p_puro_bar })
  output$scenith_tbl <- renderDT({
    req(results())
    datatable(results()$scenith_summary, options = list(pageLength = 10),
              caption = "Scenith-derived parameters per genotype")
  })

  # --------------------------------------------------------------------------
  # Downloads
  # --------------------------------------------------------------------------
  output$dl_qc_singlets <- downloadHandler(
    filename = function() "qc_singlets_plot.png",
    content  = function(f) save_plot_png(results()$plots$p_qc_singlets, f, 7, 5))
  output$dl_fsc_all <- downloadHandler(
    filename = function() "fsc_all_samples.png",
    content  = function(f) {
      req(isTRUE(input$make_all_fsc_plot))
      save_plot_png(results()$plots$p_fsc_all, f, 18, 12)
    })
  output$dl_qc_live <- downloadHandler(
    filename = function() "qc_live_plot.png",
    content  = function(f) save_plot_png(results()$plots$p_qc_live, f, 7, 5))
  output$dl_qc_puro <- downloadHandler(
    filename = function() "qc_puro_plot.png",
    content  = function(f) save_plot_png(results()$plots$p_qc_puro, f, 10, 6))
  output$dl_cell_counts <- downloadHandler(
    filename = function() "cell_counts_plot.png",
    content  = function(f) save_plot_png(results()$plots$p_cell_counts, f, 10, 7))
  output$dl_mean_puro <- downloadHandler(
    filename = function() "mean_puro_plots.png",
    content  = function(f) save_plot_png(results()$plots$p_mean_both, f, 9, 8))
  output$dl_co_dg_dgo <- downloadHandler(
    filename = function() "dist_co_dg_dgo_params.png",
    content  = function(f) save_plot_png(results()$plots$p_co_dg_dgo_params, f, 12, 8))
  output$dl_co_o_dgo <- downloadHandler(
    filename = function() "dist_co_o_dgo_params.png",
    content  = function(f) save_plot_png(results()$plots$p_co_o_dgo_params, f, 12, 8))
  output$dl_puro_bar <- downloadHandler(
    filename = function() "puro_bar_genotype.png",
    content  = function(f) save_plot_png(results()$plots$p_puro_bar, f, 12, 6))
}

shinyApp(ui, server)
