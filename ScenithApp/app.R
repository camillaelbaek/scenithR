# app.R
# Scenith analysis – Shiny app version
# Author: celbaek (converted to Shiny)

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
  library(sp)   # point.in.polygon
  library(grid) # unit()
})

# -------------------------
# Consistent color scheme
# -------------------------
treatment_cols <- c(
  "Co"         = "#4D4D4D",
  "DG"         = "#7B1FA2",
  "O"          = "#D32F2F",
  "DGO"        = "#00897B",
  "UNST"       = "#1976D2",
  "DMSO_25uL"  = "#9E9E9E"
)

genotype_cols <- c(
  "WT"      = "#1B9E77",
  "WT_HSV"  = "#7570B3",
  "KO"      = "#D95F02",
  "KO_HSV"  = "#E7298A"
)

theme_alba <- function(base_size = 12){
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

# Helper to save plots via downloadHandler
save_plot_png <- function(plot_obj, file, width = 7, height = 5, dpi = 300) {
  ggsave(filename = file, plot = plot_obj, width = width, height = height, dpi = dpi, units = "in")
}

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .small-note { color: #555; font-size: 0.95em; }
      .busy { color: #a33; font-weight: 600; }
    "))
  ),
  titlePanel("Scenith analysis – Shiny app (Alba Perez Arribas)"),

  sidebarLayout(
    sidebarPanel(
      h4("1) Load data"),
      fileInput(
        "fcs_files",
        "Upload FCS files (multiple)",
        multiple = TRUE,
        accept = c(".fcs")
      ),
      tags$p(class = "small-note",
             "Tip: select all FCS files from your experiment folder and upload them together."
      ),
      hr(),

      h4("2) Gating parameters"),
      numericInput("fitc_threshold", "Live/Dead FITC-A threshold (keep < threshold)", value = 4000, min = 1),
      numericInput("puro_threshold", "Puromycin APC-A threshold (keep ≥ threshold)", value = 80, min = 0),
      hr(),

      h4("3) QC plot options"),
      numericInput("rep_singlet_idx", "Representative sample index (singlet QC)", value = 51, min = 1),
      numericInput("rep_live_idx",   "Representative sample index (live QC)",   value = 39, min = 1),
      textInput("rep_puro_idx", "Representative sample indices (puromycin QC; comma-separated)", value = "15,27,39,51,1,4"),
      checkboxInput("make_all_fsc_plot", "Render the ALL-samples FSC-A vs FSC-H plot (can be slow)", value = FALSE),

      hr(),
      actionButton("run", "Run / Recompute analysis", class = "btn-primary"),
      tags$p(class = "small-note",
             "Any change above requires clicking Run to recompute."
      )
    ),

    mainPanel(
      tabsetPanel(
        tabPanel(
          "Overview",
          h4("Experimental overview"),
          tags$ul(
            tags$li(tags$b("DG"), " (2-deoxy-D-glucose; glycolysis inhibitor)"),
            tags$li(tags$b("Oligomycin A"), " (ATP synthase inhibitor; blocks OXPHOS)"),
            tags$li(tags$b("DGO"), " (DG + oligomycin; combined inhibition)"),
            tags$li(tags$b("Puromycin"), " (APC-A; translation activity)"),
            tags$li(tags$b("Live/Dead dye"), " (FITC-A)")
          ),
          tags$p(
            "This app reproduces your report workflow: import → sample mapping → singlet gate → live gate → puromycin gate → QC + summary + Scenith parameters."
          ),
          tags$hr(),
          uiOutput("status_box")
        ),

        tabPanel(
          "Sample mapping",
          h4("Sample-to-well mapping (from filename)"),
          DTOutput("sample_map_tbl"),
          tags$p(class = "small-note",
                 "Mapping rules match your Rmd: columns 1–3 WT, 4–6 WT_HSV, 7–9 KO, 10–12 KO_HSV; rows B/C/D/E/G/A map to Co/DG/O/DGO/UNST/DMSO_25uL."
          )
        ),

        tabPanel(
          "Gating QC",
          h4("Step 1: Singlets gate (FSC-A vs FSC-H)"),
          plotOutput("p_qc_singlets", height = 420),
          downloadButton("dl_qc_singlets", "Download PNG"),

          hr(),
          h4("Step 1 (optional): FSC-A vs FSC-H for ALL samples (raw)"),
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

        tabPanel(
          "Cell counts",
          h4("QC: Number of cells remaining after singlet + live gating"),
          plotOutput("p_cell_counts", height = 560),
          downloadButton("dl_cell_counts", "Download PNG"),
          hr(),
          h4("Interactive table"),
          DTOutput("cell_counts_tbl")
        ),

        tabPanel(
          "Puromycin summary",
          h4("Mean puromycin signal per sample"),
          plotOutput("p_mean_both", height = 720),
          downloadButton("dl_mean_puro", "Download PNG"),
          hr(),
          h4("Interactive table"),
          DTOutput("puro_summary_tbl")
        ),

        tabPanel(
          "Scenith parameters",
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

  # A single eventReactive for everything (keeps the app responsive)
  results <- eventReactive(input$run, {

    req(input$fcs_files)
    validate(need(nrow(input$fcs_files) > 0, "Please upload .fcs files first."))

    # ---- Load flowSet from uploaded files
    # input$fcs_files$datapath are temp file paths
    df <- read.flowSet(
      files              = input$fcs_files$datapath,
      alter.names        = TRUE,
      truncate_max_range = FALSE
    )

    # Attempt to keep original filenames for sample names
    # flowSet sampleNames are derived from file basenames; use uploaded original names
    sampleNames(df) <- input$fcs_files$name

    # ---- Sample map (same rules as your Rmd)
    sample_map <- tibble(sample = sampleNames(df)) %>%
      mutate(
        well_code = str_extract(sample, "[A-H]\\d{2}"),
        well_row  = str_sub(well_code, 1, 1),
        well_col  = suppressWarnings(as.integer(str_sub(well_code, 2, 3)))
      ) %>%
      mutate(
        genotype = case_when(
          well_col %in% 1:3   ~ "WT",
          well_col %in% 4:6   ~ "WT_HSV",
          well_col %in% 7:9   ~ "KO",
          well_col %in% 10:12 ~ "KO_HSV",
          TRUE                ~ NA_character_
        ),
        treatment = case_when(
          well_row == "B" ~ "Co",
          well_row == "C" ~ "DG",
          well_row == "D" ~ "O",
          well_row == "E" ~ "DGO",
          well_row == "G" ~ "UNST",
          well_row == "A" ~ "DMSO_25uL",
          TRUE            ~ NA_character_
        )
      )

    # ---- Gates (same geometry as your Rmd)
    sqrcut <- matrix(
      c(8000,  20000,
        5000,      0,
       120000,  68000,
       120000,  90000),
      ncol = 2,
      byrow = TRUE
    )
    colnames(sqrcut) <- c("FSC.A", "FSC.H")
    pg <- polygonGate(filterId = "Singlets", gate = sqrcut)

    singleCell <- Subset(df, pg)

    fitc_threshold <- input$fitc_threshold
    rg_fitc <- rectangleGate("FITC.A" = c(0, fitc_threshold))
    singleCell_live <- Subset(singleCell, rg_fitc)

    puro_threshold <- input$puro_threshold
    gate_puro <- rectangleGate("APC.A" = c(puro_threshold, Inf))
    puro_cells <- Subset(singleCell_live, gate_puro)

    # ---- Cell counts
    samples <- sampleNames(df)
    raw_counts      <- purrr::map_int(seq_along(df),              ~ nrow(exprs(df[[.x]])))
    singlet_counts  <- purrr::map_int(seq_along(singleCell),      ~ nrow(exprs(singleCell[[.x]])))
    live_counts     <- purrr::map_int(seq_along(singleCell_live), ~ nrow(exprs(singleCell_live[[.x]])))

    cell_counts <- tibble(
      sample         = samples,
      n_raw          = raw_counts,
      n_singlets     = singlet_counts,
      n_live_singlet = live_counts
    ) %>% left_join(sample_map, by = "sample")

    # ---- Puromycin summary per sample
    sample_ids <- sampleNames(singleCell_live)
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

    # ---- Cell-level dataframe (live singlets)
    sample_ids_live <- sampleNames(singleCell_live)
    cell_level <- purrr::map_df(seq_along(sample_ids_live), function(i) {
      ff <- singleCell_live[[i]]
      as_tibble(exprs(ff)) %>% mutate(sample = sample_ids_live[i])
    }) %>% left_join(sample_map, by = "sample")

    cell_level_filtered <- cell_level %>%
      mutate(genotype = factor(genotype, levels = c("WT","WT_HSV","KO","KO_HSV"))) %>%
      filter(APC.A >= puro_threshold)

    geo_means <- cell_level_filtered %>%
      mutate(genotype = factor(genotype, levels = c("WT","WT_HSV","KO","KO_HSV"))) %>%
      filter(treatment %in% c("Co", "DG", "DGO", "O")) %>%
      group_by(genotype, treatment) %>%
      summarise(geo_mean = geomfi(APC.A), .groups = "drop")

    # ---- Scenith parameters
    scenith_dg <- geo_means %>%
      filter(treatment %in% c("Co", "DG", "DGO")) %>%
      select(genotype, treatment, geo_mean) %>%
      pivot_wider(names_from = treatment, values_from = geo_mean) %>%
      mutate(
        glucose_dependence = 100 * ((Co - DG) / (Co - DGO)),
        fao_aao_capacity   = 100 - glucose_dependence
      )

    scenith_o <- geo_means %>%
      filter(treatment %in% c("Co", "O", "DGO")) %>%
      select(genotype, treatment, geo_mean) %>%
      pivot_wider(names_from = treatment, values_from = geo_mean) %>%
      mutate(
        mito_dependence      = 100 * ((Co - O) / (Co - DGO)),
        glycolytic_capacity  = 100 - mito_dependence
      )

    scenith_summary <- scenith_dg %>%
      full_join(scenith_o, by = "genotype", suffix = c("_dg", "_o")) %>%
      select(
        genotype,
        Co_dg = Co_dg,
        DG,
        DGO_dg = DGO_dg,
        glucose_dependence,
        fao_aao_capacity,
        Co_o = Co_o,
        O,
        DGO_o = DGO_o,
        mito_dependence,
        glycolytic_capacity
      )

    # ---- QC plots
    # Representative indices (bounded)
    n_samp <- length(df)
    rep_singlet <- max(1, min(input$rep_singlet_idx, n_samp))
    rep_live    <- max(1, min(input$rep_live_idx,   length(singleCell)))

    p_qc_singlets <- ggcyto(df[rep_singlet], aes(x = FSC.A, y = FSC.H)) +
      geom_hex(bins = 60) +
      geom_gate(pg, colour = "red", size = 0.6) +
      geom_stats() +
      scale_fill_viridis_c(option = "magma") +
      theme_alba() +
      labs(
        title    = "Step 1: Singlet gate on FSC-A vs FSC-H",
        subtitle = paste0("Representative sample: ", sampleNames(df)[rep_singlet]),
        x        = "FSC-A", y = "FSC-H", fill = "Cell density"
      )

    # All-samples FSC plot (optional)
    p_fsc_all <- NULL
    if (isTRUE(input$make_all_fsc_plot)) {

      all_fsc <- purrr::map_df(seq_along(df), function(i){
        tibble(
          FSC.A  = exprs(df[[i]])[, "FSC.A"],
          FSC.H  = exprs(df[[i]])[, "FSC.H"],
          sample = sampleNames(df)[i]
        )
      }) %>% left_join(sample_map, by = "sample")

      inside_gate <- function(df_points, gate_mat) {
        sp::point.in.polygon(
          df_points$FSC.A,
          df_points$FSC.H,
          gate_mat[, 1],
          gate_mat[, 2]
        ) > 0
      }

      percent_singlets <- purrr::map_df(seq_along(df), function(i) {
        dat <- as_tibble(exprs(df[[i]])[, c("FSC.A", "FSC.H")])
        inside <- inside_gate(dat, sqrcut)
        tibble(
          sample       = sampleNames(df)[i],
          pct_singlets = 100 * mean(inside),
          n_raw        = nrow(dat)
        )
      })

      all_fsc2 <- all_fsc %>%
        left_join(percent_singlets, by = "sample") %>%
        mutate(sample_label = paste0(sample, " | ", genotype, " | ", treatment, " | n=", n_raw))

      gate_df <- as.data.frame(sqrcut)
      colnames(gate_df) <- c("FSC.A", "FSC.H")
      gate_df <- rbind(gate_df, gate_df[1, ])

      p_fsc_all <- ggplot(all_fsc2, aes(x = FSC.A, y = FSC.H)) +
        geom_hex(bins = 35) +
        geom_polygon(
          data        = gate_df,
          aes(x = FSC.A, y = FSC.H),
          inherit.aes = FALSE,
          color       = "red",
          linewidth   = 0.4,
          fill        = NA
        ) +
        scale_fill_viridis_c(option = "magma") +
        facet_wrap(~ sample_label, scales = "free", ncol = 8) +
        geom_text(
          data = dplyr::distinct(all_fsc2, sample_label, pct_singlets),
          aes(x = -Inf, y = Inf, label = paste0(round(pct_singlets, 1), "% singlets")),
          inherit.aes = FALSE,
          hjust = -0.1, vjust = 1.2,
          size = 2.5,
          color = "red",
          fontface = "bold"
        ) +
        theme_alba(8) +
        labs(
          title    = "FSC-A vs FSC-H for all samples (raw data)",
          subtitle = "Hex-binned FSC-A/FSC-H before gating; red polygon = singlet gate",
          x        = "FSC-A", y = "FSC-H", fill = "Cell density"
        )
    }

    fitc_max <- max(exprs(singleCell[[rep_live]])[, "FITC.A"], na.rm = TRUE)

    p_qc_live <- ggcyto(singleCell[rep_live], aes(x = FITC.A)) +
      annotate("rect", xmin = 1, xmax = fitc_threshold, ymin = -Inf, ymax = Inf,
               fill = "#1976D2", alpha = 0.15) +
      annotate("rect", xmin = fitc_threshold, xmax = fitc_max, ymin = -Inf, ymax = Inf,
               fill = "red", alpha = 0.10) +
      annotate(
        "segment",
        x = seq(fitc_threshold, fitc_max, length.out = 20),
        xend = seq(fitc_threshold, fitc_max, length.out = 20) + 500,
        y = -Inf, yend = Inf,
        colour = "red", alpha = 0.15, linewidth = 0.5
      ) +
      geom_density(alpha = 0.6, fill = "#1976D2") +
      geom_vline(xintercept = fitc_threshold, colour = "red", linetype = "dashed", linewidth = 0.7) +
      scale_x_log10() +
      theme_alba() +
      labs(
        title    = "Step 2: Live-cell gate on FITC-A",
        subtitle = paste0("Representative sample: ", sampleNames(singleCell)[rep_live]),
        x        = "FITC-A (log10)",
        y        = "Density"
      )

    # Puromycin QC: parse indices
    parse_idx <- function(txt, max_n) {
      x <- suppressWarnings(as.integer(str_split(txt, ",", simplify = TRUE)))
      x <- x[!is.na(x)]
      x <- x[x >= 1 & x <= max_n]
      if (length(x) == 0) x <- 1
      unique(x)
    }
    plot_list <- parse_idx(input$rep_puro_idx, length(singleCell_live))

    p_qc_puro <- ggcyto(singleCell_live[plot_list], aes(x = FSC.A, y = APC.A)) +
      geom_hex(bins = 60) +
      scale_y_log10(limits = c(1, NA)) +
      geom_hline(yintercept = puro_threshold, color = "red", linetype = "dashed") +
      scale_fill_viridis_c(option = "magma") +
      theme_alba() +
      facet_wrap(~ name) +
      labs(
        title    = "Step 3: Puromycin gate on APC-A",
        subtitle = paste0("Dashed line = APC-A threshold (", puro_threshold, ")"),
        x        = "FSC-A",
        y        = "APC-A (log10)",
        fill     = "Cell density"
      )

    p_cell_counts <- cell_counts %>%
      ggplot(aes(x = reorder(sample, n_live_singlet), y = n_live_singlet, fill = genotype)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = genotype_cols, na.value = "grey80") +
      theme_alba() +
      labs(
        title    = "QC: Number of cells remaining after singlet + live gating",
        subtitle = "Bars colored by genotype; samples ordered by retained live singlets",
        x        = "Sample",
        y        = "# Live singlet cells",
        fill     = "Genotype"
      )

    p_mean_live <- ggplot(puro_summary, aes(x = reorder(sample, mean_puro_live), y = mean_puro_live, fill = treatment)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = treatment_cols, na.value = "grey80") +
      theme_alba() +
      labs(
        title    = "Mean puromycin signal per sample (all live cells)",
        x        = "Sample",
        y        = "Mean APC-A",
        fill     = "Treatment"
      )

    p_mean_pos <- ggplot(puro_summary, aes(x = reorder(sample, mean_puro_pos), y = mean_puro_pos, fill = treatment)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = treatment_cols, na.value = "grey80") +
      theme_alba() +
      labs(
        title    = "Mean puromycin signal per sample (puromycin-positive cells)",
        x        = "Sample",
        y        = "Mean APC-A (puro+)",
        fill     = "Treatment"
      )

    p_mean_both <- ggpubr::ggarrange(p_mean_live, p_mean_pos, ncol = 1)

    p_co_dg_dgo_params <- cell_level_filtered %>%
      filter(treatment %in% c("Co", "DG", "DGO")) %>%
      ggplot(aes(x = APC.A, fill = treatment)) +
      geom_density(alpha = 0.6, aes(color = treatment)) +
      scale_x_log10(limits = c(1, 10e6)) +
      facet_wrap(~ genotype, ncol = 2) +
      scale_fill_manual(values = treatment_cols) +
      scale_colour_manual(values = treatment_cols, guide = "none") +
      geom_vline(
        data        = geo_means %>% filter(treatment %in% c("Co", "DG", "DGO")),
        aes(xintercept = geo_mean, colour = treatment),
        linetype    = "dashed",
        linewidth   = 0.6,
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      geom_segment(
        data        = scenith_dg,
        aes(x = DG, xend = Co, y = 0.5, yend = 0.5),
        arrow       = arrow(length = unit(0.15, "cm")),
        inherit.aes = FALSE
      ) +
      geom_label(
        data        = scenith_dg,
        aes(x = (DG + Co) / 0.8, y = 0.2,
            label = paste0("1. Glc dep = ", round(glucose_dependence, 1), "%")),
        size        = 4,
        inherit.aes = FALSE
      ) +
      geom_segment(
        data        = scenith_dg,
        aes(x = DGO, xend = DG, y = 1, yend = 1),
        arrow       = arrow(length = unit(0.15, "cm")),
        inherit.aes = FALSE,
        color       = "purple4"
      ) +
      geom_label(
        data        = scenith_dg,
        aes(x = (DGO + DG) / 8, y = 0.8,
            label = paste0("4. FAO/AAO cap = ", round(fao_aao_capacity, 1), "%")),
        size        = 4,
        inherit.aes = FALSE,
        color       = "purple4"
      ) +
      theme_alba(14) +
      labs(
        title    = "Scenith parameters – Glucose dependence and FAO/AAO capacity",
        subtitle = "Densities pooled across technical replicates; arrows use geometric means",
        x        = "APC-A (log10)",
        y        = "Density",
        fill     = "Treatment"
      )

    p_co_o_dgo_params <- cell_level_filtered %>%
      filter(treatment %in% c("Co", "O", "DGO")) %>%
      ggplot(aes(x = APC.A, fill = treatment)) +
      geom_density(alpha = 0.6, aes(color = treatment)) +
      scale_x_log10(limits = c(1, 10e6)) +
      facet_wrap(~ genotype, ncol = 2) +
      scale_fill_manual(values = treatment_cols) +
      scale_colour_manual(values = treatment_cols, guide = "none") +
      geom_vline(
        data        = geo_means %>% filter(treatment %in% c("Co", "O", "DGO")),
        aes(xintercept = geo_mean, colour = treatment),
        linetype    = "dashed",
        linewidth   = 0.6,
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      geom_segment(
        data        = scenith_o,
        aes(x = O, xend = Co, y = 0.5, yend = 0.5),
        arrow       = arrow(length = unit(0.15, "cm")),
        inherit.aes = FALSE,
        color       = "blue"
      ) +
      geom_label(
        data        = scenith_o,
        aes(x = (O + Co) / 0.8, y = 0.2,
            label = paste0("2. Mito dep = ", round(mito_dependence, 1), "%")),
        size        = 4,
        inherit.aes = FALSE,
        color       = "blue"
      ) +
      geom_segment(
        data        = scenith_o,
        aes(x = DGO, xend = O, y = 1, yend = 1),
        arrow       = arrow(length = unit(0.15, "cm")),
        inherit.aes = FALSE
      ) +
      geom_label(
        data        = scenith_o,
        aes(x = (DGO + O) / 5, y = 0.8,
            label = paste0("3. Glyc cap = ", round(glycolytic_capacity, 1), "%")),
        size        = 4,
        inherit.aes = FALSE
      ) +
      theme_alba(14) +
      labs(
        title    = "Scenith parameters – Mitochondrial dependence and glycolytic capacity",
        subtitle = "Densities pooled across technical replicates; arrows use geometric means",
        x        = "APC-A (log10)",
        y        = "Density",
        fill     = "Treatment"
      )

    bar_puro <- geo_means %>%
      filter(treatment %in% c("Co", "DG", "O", "DGO")) %>%
      mutate(treatment = factor(treatment, levels = c("Co", "DG", "O", "DGO")))

    p_puro_bar <- ggplot(bar_puro, aes(x = treatment, y = geo_mean, fill = treatment)) +
      geom_col(width = 0.7) +
      facet_wrap(~ genotype, ncol = 4) +
      scale_fill_manual(values = treatment_cols) +
      theme_alba(14) +
      labs(
        title    = "Translation per condition and genotype",
        subtitle = "Bars show geometric mean puromycin signal (APC-A) in live singlet cells",
        x        = "Condition",
        y        = "Geometric mean APC-A (puromycin MFI)",
        fill     = "Treatment"
      )

    list(
      df = df,
      sample_map = sample_map,
      sqrcut = sqrcut,
      pg = pg,
      singleCell = singleCell,
      singleCell_live = singleCell_live,
      puro_cells = puro_cells,
      cell_counts = cell_counts,
      puro_summary = puro_summary,
      geo_means = geo_means,
      scenith_summary = scenith_summary,
      plots = list(
        p_qc_singlets = p_qc_singlets,
        p_fsc_all = p_fsc_all,
        p_qc_live = p_qc_live,
        p_qc_puro = p_qc_puro,
        p_cell_counts = p_cell_counts,
        p_mean_both = p_mean_both,
        p_co_dg_dgo_params = p_co_dg_dgo_params,
        p_co_o_dgo_params = p_co_o_dgo_params,
        p_puro_bar = p_puro_bar
      )
    )
  }, ignoreInit = TRUE)

  output$status_box <- renderUI({
    if (is.null(input$fcs_files) || nrow(input$fcs_files) == 0) {
      return(tags$p(class = "busy", "No files loaded yet. Upload .fcs files in the sidebar."))
    }
    tags$p(class = "small-note",
           paste0("Uploaded files: ", nrow(input$fcs_files),
                  ". Click 'Run / Recompute analysis' to process.")
    )
  })

  # ---- Tables
  output$sample_map_tbl <- renderDT({
    req(results())
    datatable(
      results()$sample_map,
      options = list(pageLength = 10),
      caption = "Mapping of FCS files to plate wells, genotypes and treatments"
    )
  })

  output$cell_counts_tbl <- renderDT({
    req(results())
    datatable(
      results()$cell_counts %>%
        mutate(genotype = factor(genotype, levels = c("WT","WT_HSV","KO","KO_HSV"))) %>%
        arrange(genotype, treatment, sample),
      options = list(pageLength = 10),
      caption = "Cell counts at each gating step per sample"
    )
  })

  output$puro_summary_tbl <- renderDT({
    req(results())
    datatable(
      results()$puro_summary %>%
        mutate(genotype = factor(genotype, levels = c("WT","WT_HSV","KO","KO_HSV"))) %>%
        arrange(genotype, treatment, sample),
      options = list(pageLength = 10),
      caption = "Puromycin summary statistics per sample"
    )
  })

  output$scenith_tbl <- renderDT({
    req(results())
    datatable(
      results()$scenith_summary,
      options = list(pageLength = 10),
      caption = "Scenith-derived parameters per genotype (based on geometric mean APC-A)"
    )
  })

  # ---- Plots
  output$p_qc_singlets <- renderPlot({ req(results()); results()$plots$p_qc_singlets })
  output$p_qc_live     <- renderPlot({ req(results()); results()$plots$p_qc_live })
  output$p_qc_puro     <- renderPlot({ req(results()); results()$plots$p_qc_puro })

  output$all_fsc_notice <- renderUI({
    req(results())
    if (!isTRUE(input$make_all_fsc_plot)) {
      tags$p(class = "small-note",
             "Disabled. Enable it in the sidebar if you want this plot (can be slow and memory-heavy).")
    } else {
      tags$p(class = "small-note", "Rendering enabled.")
    }
  })

  output$p_fsc_all <- renderPlot({
    req(results())
    if (!isTRUE(input$make_all_fsc_plot)) return(invisible(NULL))
    results()$plots$p_fsc_all
  })

  output$p_cell_counts <- renderPlot({ req(results()); results()$plots$p_cell_counts })
  output$p_mean_both   <- renderPlot({ req(results()); results()$plots$p_mean_both })
  output$p_co_dg_dgo_params <- renderPlot({ req(results()); results()$plots$p_co_dg_dgo_params })
  output$p_co_o_dgo_params  <- renderPlot({ req(results()); results()$plots$p_co_o_dgo_params })
  output$p_puro_bar <- renderPlot({ req(results()); results()$plots$p_puro_bar })

  # ---- Downloads
  output$dl_qc_singlets <- downloadHandler(
    filename = function() "qc_singlets_plot.png",
    content = function(file) save_plot_png(results()$plots$p_qc_singlets, file, width = 7, height = 5)
  )
  output$dl_fsc_all <- downloadHandler(
    filename = function() "fsc_all_samples.png",
    content = function(file) {
      req(isTRUE(input$make_all_fsc_plot))
      save_plot_png(results()$plots$p_fsc_all, file, width = 18, height = 12)
    }
  )
  output$dl_qc_live <- downloadHandler(
    filename = function() "qc_live_plot.png",
    content = function(file) save_plot_png(results()$plots$p_qc_live, file, width = 7, height = 5)
  )
  output$dl_qc_puro <- downloadHandler(
    filename = function() "qc_puro_plot.png",
    content = function(file) save_plot_png(results()$plots$p_qc_puro, file, width = 10, height = 6)
  )
  output$dl_cell_counts <- downloadHandler(
    filename = function() "cell_counts_plot.png",
    content = function(file) save_plot_png(results()$plots$p_cell_counts, file, width = 10, height = 7)
  )
  output$dl_mean_puro <- downloadHandler(
    filename = function() "mean_puro_plots.png",
    content = function(file) save_plot_png(results()$plots$p_mean_both, file, width = 9, height = 8)
  )
  output$dl_co_dg_dgo <- downloadHandler(
    filename = function() "dist_co_dg_dgo_params.png",
    content = function(file) save_plot_png(results()$plots$p_co_dg_dgo_params, file, width = 12, height = 8)
  )
  output$dl_co_o_dgo <- downloadHandler(
    filename = function() "dist_co_o_dgo_params.png",
    content = function(file) save_plot_png(results()$plots$p_co_o_dgo_params, file, width = 12, height = 8)
  )
  output$dl_puro_bar <- downloadHandler(
    filename = function() "puro_bar_genotype.png",
    content = function(file) save_plot_png(results()$plots$p_puro_bar, file, width = 12, height = 6)
  )
}

shinyApp(ui, server)

