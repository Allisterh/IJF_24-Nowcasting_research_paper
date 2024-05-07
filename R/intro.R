if (!"pacman" %in% installed.packages()[, "Package"]) install.packages("pacman", repos = "http://cran.r-project.org")
pacman::p_load(
    dplyr, tidyr, ggplot2, zoo, magrittr, readxl, stargazer, ucminf, goftest, lmtest,
    wesanderson, lubridate, Rcpp, RcppEigen, readr, RcppScoreDrivenDFM, ggh4x, rugarch,
    stringr, scoringRules, insight, HDInterval, isodistrreg, ggpubr, unikn, paletteer, binom
)

theme <- theme_bw() + theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    title = element_text(size = 11),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.box.margin = margin(-10, 0, -10, 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    axis.text = element_text(size = 11),
    axis.line.y = element_line(colour = "grey", size = 0.5, linetype = "solid"),
    axis.line.x = element_line(colour = "grey", size = 0.5, linetype = "solid"),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
)

theme_nolegend <- theme + theme(legend.position = "none")

select <- dplyr::select
