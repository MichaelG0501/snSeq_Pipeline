####################
# Analysis registry:
#   Status: active shared library; not a standalone workflow script
#   Script: analysis/lib/plot_helpers.R
#   Methodology: analysis/methodology/lib/shared_analysis_library_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Input: none
#   Output: shared plotting helper functions in the calling R session
####################

####################
# Shared plotting defaults for presentation-facing figures.
####################

pptx_theme_classic <- function(base_size = PLOT_DEFAULTS$pptx_base_size) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = PLOT_DEFAULTS$pptx_title_size, hjust = 0.5),
      axis.title = ggplot2::element_text(size = PLOT_DEFAULTS$pptx_axis_size, face = "bold"),
      axis.text = ggplot2::element_text(size = PLOT_DEFAULTS$pptx_axis_size),
      legend.title = ggplot2::element_text(size = PLOT_DEFAULTS$pptx_legend_size, face = "bold"),
      legend.text = ggplot2::element_text(size = PLOT_DEFAULTS$pptx_legend_size)
    )
}

pptx_theme_minimal <- function(base_size = PLOT_DEFAULTS$pptx_base_size) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = PLOT_DEFAULTS$pptx_title_size, hjust = 0.5),
      axis.title = ggplot2::element_text(size = PLOT_DEFAULTS$pptx_axis_size, face = "bold"),
      axis.text = ggplot2::element_text(size = PLOT_DEFAULTS$pptx_axis_size),
      legend.title = ggplot2::element_text(size = PLOT_DEFAULTS$pptx_legend_size, face = "bold"),
      legend.text = ggplot2::element_text(size = PLOT_DEFAULTS$pptx_legend_size),
      panel.grid.minor = ggplot2::element_blank()
    )
}

heatmap_text_gp <- function(size = PLOT_DEFAULTS$heatmap_row_font) {
  grid::gpar(fontsize = size)
}
