#' nichevol: An R package for assessment of the evolution of ecological niches
#' considering uncertainty in niche reconstructions
#'
#' nichevol implements multiple tools to help in performing critical steps of
#' the process of assessment of evolution of species ecological niches with
#' uncertainty incorporated in reconstructions. The main functionalities of the
#' package include: intial exploration of environmental data in species records
#' and accessible areas, preparation of data for phylogenetic analyses,
#' phylogenetic analyses of ecological niches, and plotting for interpretations.
#'
#' @section Main functions in nichevol:
#' \code{\link{bin_ml_rec}}, \code{\link{bin_par_rec}},
#' \code{\link{bin_tables}}, \code{\link{bin_tables0}},
#' \code{\link{bin_tables_virtual}}, \code{\link{histograms_env}},
#' \code{\link{niche_bars}}, \code{\link{nichevol_bars}},
#' \code{\link{niche_labels}}, \code{\link{nichevol_labels}},
#' \code{\link{niche_legend}}, \code{\link{nichevol_legend}},
#' \code{\link{random_polygons}}, \code{\link{smooth_rec}},
#' \code{\link{stats_evalues}}, \code{\link{stats_evalues_virtual}}
#'
#' Other functions (important helpers)
#'
#' \code{\link{bin_env}}, \code{\link{bin_env_null}},
#' \code{\link{make_9blocks}}, \code{\link{pdf_histograms}},
#' \code{\link{rename_tips}}, \code{\link{rformat_type}},
#' \code{\link{score_tip}}, \code{\link{score_tree}}, \code{\link{sig_sq}}
#'
#' @docType package
#' @name nichevol
NULL
