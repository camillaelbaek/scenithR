# run_app.R
pkgs <- c(
  "shiny","dplyr","tidyr","ggplot2","DT","stringr","purrr","scales",
  "flowCore","flowViz","ggcyto","openCyto","ggridges",
  "ggpubr","viridis","ggbeeswarm","sp"
)

need <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(need)) {
  install.packages(need, repos = "https://cloud.r-project.org")
}

shiny::runApp(".", launch.browser = TRUE)

