#' nichevol: Assessment of Species’ Ecological Niche Evolution Considering
#' Uncertainty in Reconstructions
#'
#' nichevol is a collection of tools that allow users to perform critical steps
#' in the process of assessing ecological niche evolution over phylogenies, with
#' uncertainty incorporated explicitly in reconstructions. The method proposed
#' here for ancestral reconstruction of ecological niches characterizes species'
#' niches using a bin-based approach that incorporates uncertainty in estimations.
#' Compared to other existing methods, the approaches presented here reduce risk
#' of overestimation of amounts and rates of ecological niche evolution. The
#' main analyses include: initial exploration of environmental data in occurrence
#' records and accessible areas, preparation of data for phylogenetic analyses,
#' executing comparative phylogenetic analyses of ecological niches, and plotting
#' for interpretations.
#'
#' @section Main functions in nichevol:
#' \code{\link{bin_ml_rec}}, \code{\link{bin_par_rec}}, \code{\link{bin_table}},
#' \code{\link{bin_tables}}, \code{\link{bin_tables0}}, \code{\link{hist_evalues}},
#' \code{\link{histograms_env}}, \code{\link{map_nichevol}},
#' \code{\link{niche_bars}}, \code{\link{nichevol_bars}},
#' \code{\link{niche_labels}}, \code{\link{nichevol_labels}},
#' \code{\link{niche_legend}}, \code{\link{nichevol_legend}},
#' \code{\link{set_uncertainty}}, \code{\link{smooth_rec}},
#' \code{\link{stats_eval}}, \code{\link{stats_evalues}}
#'
#' Other functions (important helpers)
#'
#' \code{\link{bin_env}}, \code{\link{pdf_histograms}},
#' \code{\link{rename_tips}}, \code{\link{score_tip}},
#' \code{\link{score_tree}}, \code{\link{sig_sq}}
#'
#' @docType package
#' @name nichevol
NULL
