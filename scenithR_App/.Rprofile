# .Rprofile

cran_pkgs <- c(
  "shiny", "dplyr", "tidyr", "ggplot2", "DT", "stringr", "purrr",
  "scales", "readr", "readxl", "ggridges", "ggpubr", "viridis",
  "ggbeeswarm", "sp", "magick", "imager"
)

bioc_pkgs <- c(
  "flowCore", "flowViz", "ggcyto", "openCyto"
)

need_cran <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]

if (length(need_cran)) {
  message("Installing missing CRAN packages...")
  install.packages(need_cran, repos = "https://cloud.r-project.org")
}

need_bioc <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]

if (length(need_bioc)) {
  message("Installing missing Bioconductor packages...")

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }

  BiocManager::install(need_bioc, ask = FALSE, update = FALSE)
}
