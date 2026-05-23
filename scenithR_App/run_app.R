# run_app.R

cran_pkgs <- c(
  "shiny", "dplyr", "tidyr", "ggplot2", "DT", "stringr", "purrr",
  "scales", "readr", "readxl", "ggridges", "ggpubr", "viridis",
  "ggbeeswarm", "sp", "magick", "imager"
)

bioc_pkgs <- c(
  "flowCore", "flowViz", "ggcyto", "openCyto"
)

need_cran <- setdiff(cran_pkgs, rownames(installed.packages()))

if (length(need_cran)) {
  install.packages(need_cran, repos = "https://cloud.r-project.org")
}

need_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))

if (length(need_bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }

  BiocManager::install(need_bioc, ask = FALSE, update = FALSE)
}

shiny::runApp(".")
