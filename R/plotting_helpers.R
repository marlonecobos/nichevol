#' Labels to represent niches of tips and ancestors
#'
#' @description niche_labels helps in adding bar-type labels that represent
#' species ecological niches in one environmental variable.
#'
#' @param tree an object of class "phylo".
#' @param whole_rec_table matrix of environmental bins for all tips and nodes
#' derived from functions \code{\link{bin_par_rec}} or \code{\link{bin_ml_rec}}.
#' @param label_type (character) type of label; options are: "tip", "node", and
#' "tip_node". Default = "tip_node".
#' @param tip_offset (numeric) space between tips and the labels. Default = 0.015.
#' @param present (character) code indicating environmental bins in which the
#' species is present. Default = "1".
#' @param unknown (character) code indicating environmental bins in which the
#' species presence is unknown (uncertain). Default = "?".
#' @param present_col color for area of the bar representing environments where
#' the species is present. Default = "#e41a1c".
#' @param unknown_col color for area of the bar representing environments where
#' the species presence is unknown (uncertain). Default = "#969696".
#' @param absent_col color for area of the bar representing environments where
#' no change has been detected. Default = "#377eb8".
#' @param width value defining the width of niche bars; default = 1.
#' @param height value defining the height of niche bars; default = 1.
#'
#' @details
#' For the moment, only plots of type "phylogram" with "rightwards" or "leftwards"
#' directions, created with the function \code{\link[ape]{plot.phylo}} from the
#' package \code{ape} are supported.
#'
#' Ecological niches are represented in one environmental dimension with vertical
#' bars that indicate if the species is present, absent, or if its presence is
#' uncertain in the range of environmental conditions. Lower values of
#' environmental variables are represented in the lower part of the bar, and the
#' opposite part of the bar represents higher values.
#'
#' @importFrom graphics plot polygon
#'
#' @export
#'
#' @usage
#' niche_labels(tree, whole_rec_table, label_type = "tip_node",
#'              tip_offset = 0.015, present = "1", unknown = "?",
#'              present_col = "#e41a1c", unknown_col = "#969696",
#'              absent_col = "#377eb8", width = 1, height = 1)
#'
#' @examples
#' # a simple tree
#' data("tree5", package = "nichevol")
#'
#' # a matrix of niche charactes (1 = present, 0 = absent, ? = unknown)
#' dataTable <- cbind("241" = rep("1", length(tree5$tip.label)),
#'                    "242" = rep("1", length(tree5$tip.label)),
#'                    "243" = c("1", "1", "0", "0", "0"),
#'                    "244" = c("1", "1", "0", "0", "0"),
#'                    "245" = c("1", "?", "0", "0", "0"))
#' rownames(dataTable) <- tree5$tip.label
#'
#' # list with two objects (tree and character table)
#' treeWdata <- geiger::treedata(tree5, dataTable)
#'
#' # Maximum parsimony reconstruction
#' rec_tab <- smooth_rec(bin_par_rec(treeWdata))
#'
#' # plotting and adding labels
#' ape::plot.phylo(tree5, label.offset = 0.04)
#' niche_labels(tree5, rec_tab, height = 0.6)

niche_labels <- function(tree, whole_rec_table, label_type = "tip_node",
                         tip_offset = 0.015, present = "1", unknown = "?",
                         present_col = "#e41a1c", unknown_col = "#969696",
                         absent_col = "#377eb8", width = 1, height = 1) {
  if (missing(tree)) {stop("Argument 'tree' needs to be defined.")}
  if (missing(whole_rec_table)) {stop("Argument 'whole_rec_table' needs to be defined.")}
  if ("LogLik" %in% rownames(whole_rec_table)) {
    whole_rec_table <- whole_rec_table[1:(nrow(whole_rec_table) - 3), ]
  }

  # reorganizing character table
  tlab <- tree$tip.label
  nrt <- nrow(whole_rec_table)
  rns <- c(tlab, rownames(whole_rec_table)[(length(tlab) + 1):nrt])
  whole_rec_table <- rbind(whole_rec_table[tlab, ],
                           whole_rec_table[(length(tlab) + 1):nrt, ])
  rownames(whole_rec_table) <- rns

  # getting info from plot
  tp_info <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  if (tp_info$type != "phylogram") {
    stop("'niche_labels' can be used only for plots of type phylogram.")
  }
  if (!tp_info$direction %in% c("rightwards", "leftwards")) {
    stop("'niche_labels' can be used only for rightwards or leftwards phylograms.")
  }
  if (tp_info$direction == "leftwards") {tip_offset <- -tip_offset}

  xx <- tp_info$xx
  yy <- tp_info$yy
  edges <- tp_info$edge
  tpos <- 1:tp_info$Ntip
  npos <- (tp_info$Ntip + 1):(tp_info$Ntip + tp_info$Nnode)
  otips <- xx[tpos]
  rtip <- range(otips)
  if ((rtip[2] - rtip[1]) <= 0.00001) {
    xx[tpos] <- rep(max(otips), length(otips))
  }

  # organizing data
  tpol <- ncol(whole_rec_table)
  wt <- ((max(yy) / tp_info$Ntip) / 10) * (height * 6)
  wpol <- wt / tpol

  h_vertices <- seq(0, wt, wpol)

  # plotting bars
  if (label_type %in% c("tip", "node", "tip_node")) {
    if (label_type %in% c("tip", "tip_node")) {
      barss <- sapply(tpos, function(j) {
        ys <- yy[j] - (wt / 2)
        hver <- ys + h_vertices
        wdt <- 0.01 * width
        xs <- xx[j] + tip_offset; xs1 <- xs - (wdt / 2); xs2 <- xs + (wdt / 2)
        wver <- rep(c(xs1, xs2), each = 2)

        polys <- sapply(1:(length(h_vertices) - 1), function(x) {
          if (as.character(whole_rec_table[j, x]) == unknown) {
            pcolor <- unknown_col
          } else {
            pcolor <- ifelse(as.character(whole_rec_table[j, x]) == present,
                             present_col, absent_col)
          }

          yss <- c(hver[x], hver[x + 1], hver[x + 1], hver[x])

          polygon(x = wver, y = yss, col = pcolor, border = NA)
        })
      })
    }

    if (label_type %in% c("node", "tip_node")) {
      barss <- sapply(npos, function(j) {
        ys <- yy[j] - (wt / 2)
        hver <- ys + h_vertices
        wdt <- 0.01 * width
        xs <- xx[j]; xs1 <- xs - (wdt / 2); xs2 <- xs + (wdt / 2)
        wver <- rep(c(xs1, xs2), each = 2)

        polys <- sapply(1:(length(h_vertices) - 1), function(x) {
          if (as.character(whole_rec_table[j, x]) == unknown) {
            pcolor <- unknown_col
          } else {
            pcolor <- ifelse(as.character(whole_rec_table[j, x]) == present,
                             present_col, absent_col)
          }

          yss <- c(hver[x], hver[x + 1], hver[x + 1], hver[x])

          polygon(x = wver, y = yss, col = pcolor, border = NA)
        })
      })
    }

  } else {
    stop("Argument 'label_type' is not correct, see help(niche_labels).")
  }
}


#' Labels to represent changes of niche characteristics between ancestors and descendants
#'
#' @description nichevol_labels helps in adding bar-type labels that represent how
#' species' niches changed from ancestors to descendants.
#'
#' @param tree an object of class "phylo".
#' @param whole_rec_table matrix of reconstructed bins for nodes and species
#' derived from a process of maximum parsimony or maximum likelihood reconstruction.
#' See functions \code{\link{bin_par_rec}} or \code{\link{bin_ml_rec}}.
#' @param ancestor_line controls whether ancestor line is plotted.
#' Default = FALSE.
#' @param present (character) code indicating environmental bins in which the
#' species is present. Default = "1".
#' @param absent (character) code indicating environmental bins in which the
#' species is absent. Default = "0".
#' @param unknown (character) code indicating environmental bins in which the
#' species presence is unknown (uncertain). Default = "?".
#' @param present_col color for line representing environments where the species
#' is present. Default = "#252525".
#' @param unknown_col color for line representing environments where the species
#' presence is unknown (uncertain). Default = "#d9d9d9".
#' @param no_change_col color for area of the bar representing environments where
#' no change has been detected. Default = "#b2df8a".
#' @param retraction_col color for area of the bar representing environments where
#' niche retraction has been detected. Default = "#984ea3".
#' @param expansion_col color for area of the bar representing environments where
#' niche expansion has been detected. Default = "#4daf4a".
#' @param width value defining the width of bars representing changes in niches;
#' default = 1.
#' @param height value defining the height of bars representing changes in niches;
#' default = 1.
#'
#' @details
#' For the moment, only plots of type "phylogram" with "rightwards" or "leftwards"
#' directions, created with the function \code{\link[ape]{plot.phylo}} from the
#' package \code{ape} are supported.
#'
#' Evolution of ecological niches is represented in one environmental dimension,
#' with vertical bars indicating if the niche of the descendant has expanded,
#' retracted, or has not changed compared to its ancestor's niche. Lower values of
#' environmental variables are represented in the lower part of the bar, and the
#' opposite part of the bar represents higher values.
#'
#' Changes in niches (evolution) are defined as follows:
#' - if (ancestor == present & descendant == absent) {change <- "retraction"}
#' - if (ancestor == present & descendant == present) {change <- "no_change"}
#' - if (ancestor == present & descendant == unknown) {change <- "no_change"}
#' - if (ancestor == absent & descendant == present) {change <- "expansion"}
#' - if (ancestor == absent & descendant == absent) {change <- "no_change"}
#' - if (ancestor == absent & descendant == unknown) {change <- "no_change"}
#' - if (ancestor == unknown & descendant == absent) {change <- "no_change"}
#' - if (ancestor == unknown & descendant == present) {change <- "no_change"}
#' - if (ancestor == unknown & descendant == unknown) {change <- "no_change"}
#'
#' If \code{ancestor_line} is TRUE, the ancestor line will be plotted on the bar
#' representing niche evolution. The line will represent where, in the range of
#' environmental conditions, the ancestor was present, and where its presence is
#' uncertain (unknown).
#'
#' @importFrom graphics plot polygon lines
#'
#' @export
#'
#' @usage
#' nichevol_labels(tree, whole_rec_table, ancestor_line = FALSE,
#'   present = "1", absent = "0", unknown = "?",
#'   present_col = "#252525", unknown_col = "#d9d9d9",
#'   no_change_col = "#b2df8a", retraction_col = "#984ea3",
#'   expansion_col = "#4daf4a", width = 1, height = 1)
#'
#' @examples
#' # a simple tree
#' data("tree5", package = "nichevol")
#'
#' # a matrix of niche charactes (1 = present, 0 = absent, ? = unknown)
#' dataTable <- cbind("241" = rep("1", length(tree5$tip.label)),
#'                    "242" = rep("1", length(tree5$tip.label)),
#'                    "243" = c("1", "1", "0", "0", "0"),
#'                    "244" = c("1", "1", "0", "0", "0"),
#'                    "245" = c("1", "?", "0", "0", "0"))
#' rownames(dataTable) <- tree5$tip.label
#'
#' # list with two objects (tree and character table)
#' treeWdata <- geiger::treedata(tree5, dataTable)
#'
#' # Maximum parsimony reconstruction
#' rec_tab <- smooth_rec(bin_par_rec(treeWdata))
#'
#' # plotting and adding labels
#' ape::plot.phylo(tree5, label.offset = 0.04)
#' nichevol_labels(tree5, rec_tab, height = 0.6)

nichevol_labels <- function(tree, whole_rec_table, ancestor_line = FALSE,
                            present = "1", absent = "0", unknown = "?",
                            present_col = "#252525", unknown_col = "#d9d9d9",
                            no_change_col = "#b2df8a", retraction_col = "#984ea3",
                            expansion_col = "#4daf4a", width = 1, height = 1) {
  if (missing(tree)) {stop("Argument 'tree' needs to be defined.")}
  if (missing(whole_rec_table)) {stop("Argument 'whole_rec_table' needs to be defined.")}
  if ("LogLik" %in% rownames(whole_rec_table)) {
    whole_rec_table <- whole_rec_table[1:(nrow(whole_rec_table) - 3), ]
  }

  # reorganizing character table
  tlab <- tree$tip.label
  nrt <- nrow(whole_rec_table)
  rns <- c(tlab, rownames(whole_rec_table)[(length(tlab) + 1):nrt])
  whole_rec_table <- rbind(whole_rec_table[tlab, ],
                           whole_rec_table[(length(tlab) + 1):nrt, ])
  rownames(whole_rec_table) <- rns

  # getting info from plot
  tp_info <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  if (tp_info$type != "phylogram") {
    stop("'nichevol_labels' can be used only for plots of type phylogram.")
  }
  if (!tp_info$direction %in% c("rightwards", "leftwards")) {
    stop("'nichevol_labels' can be used only for rightwards or leftwards phylograms.")
  }

  xx <- tp_info$xx
  yy <- tp_info$yy
  edges <- tp_info$edge

  # preparing for niche evolution
  ## positions
  xb <- sapply(1:nrow(edges), function(x) {
    ds <- xx[edges[x, 1]] + ((xx[edges[x, 2]] - xx[edges[x, 1]]) / 2)
  })

  ## organizing data
  tpol <- ncol(whole_rec_table)
  wt <- ((max(yy) / tp_info$Ntip) / 10) * (height * 6)
  wpol <- wt / tpol

  h_vertices <- seq(0, wt, wpol)

  ## comparisons and plots
  comp_list <- lapply(1:nrow(edges), function(x) {
    comp <- sapply(1:ncol(whole_rec_table), function(z) {
      from <- whole_rec_table[edges[x, 1], z]
      to <- whole_rec_table[edges[x, 2], z]
      if (from == present & to == absent) {comp <- "loss"}
      if (from == present & to == present) {comp <- "nc"}
      if (from == present & to == unknown) {comp <- "nc"}
      if (from == absent & to == present) {comp <- "gain"}
      if (from == absent & to == absent) {comp <- "nc"}
      if (from == absent & to == unknown) {comp <- "nc"}
      if (from == unknown & to == absent) {comp <- "nc"}
      if (from == unknown & to == present) {comp <- "nc"}
      if (from == unknown & to == unknown) {comp <- "nc"}
      return(comp)
    })

    ## infor for lines and polygons per each site
    linesp <- whole_rec_table[edges[x, 1], ]

    ys <- yy[edges[x, 2]] - (wt / 2)
    hver <- ys + h_vertices

    wdt <- 0.01 * width
    xs <- xb[x]; xs1 <- xs - (wdt / 2); xs2 <- xs + (wdt / 2)
    wver <- rep(c(xs1, xs2), each = 2)
    x_line <- rep(mean(c(xs1, xs2)) - 0.002, 2)

    ## bar creation
    poly_lines <- sapply(1:(length(hver) - 1), function(y) {
      ### polygons
      if (comp[y] == "nc") {
        pcolor <- no_change_col
      } else {
        pcolor <- ifelse(comp[y] == "gain", expansion_col, retraction_col)
      }

      yss <- c(hver[y], hver[y + 1], hver[y + 1], hver[y])
      polygon(x = wver, y = yss, col = pcolor, border = NA)

      ### lines
      if (ancestor_line == TRUE) {
        if (linesp[y] == unknown) {
          pcolor <- unknown_col
        } else {
          pcolor <- ifelse(linesp[y] == present, present_col, "transparent")
        }

        yss <- c(hver[y], hver[y + 1])
        lines(x = x_line, y = yss, col = pcolor, lty = 1, lwd = 1.7)
      }
    })
  })
}


#' Legends for niche labels in phylogenetic trees
#'
#' @param position (character or numeric) position of legend. If character,
#' part of the plot (e.g., "topleft"), see \code{\link[graphics]{legend}}. If
#' numeric, vector of two values indicating x and y position (e.g., c(0.1, 6)).
#' @param legend (character) vector of length = three indicating the text to
#' identify environments with uncertain presence, presence, and absence of the
#' species. Default = c("Uncertain", "Present", "Not present").
#' @param pch point type as in \code{\link[graphics]{points}}. Default = 22.
#' @param pt.bg colors to represent what is in \code{legend}.
#' Default = c("#969696", "#e41a1c", "#377eb8").
#' @param col border of symbol (points). Default = "transparent".
#' @param pt.cex size of symbol (points). Default = 2.2.
#' @param bty legend border type. Default = "n".
#' @param ... Other arguments from function \code{\link[graphics]{legend}} other
#' than the ones described above.
#'
#' @importFrom graphics legend
#'
#' @export
#'
#' @usage
#' niche_legend(position, legend = c("Uncertain", "Present", "Not present"),
#'   pch = 22, pt.bg = c("#969696", "#e41a1c", "#377eb8"),
#'   col = "transparent", pt.cex = 2.2, bty = "n", ...)
#'
#' @examples
#' # a simple tree
#' data("tree5", package = "nichevol")
#'
#' # a matrix of niche charactes (1 = present, 0 = absent, ? = unknown)
#' dataTable <- cbind("241" = rep("1", length(tree5$tip.label)),
#'                    "242" = rep("1", length(tree5$tip.label)),
#'                    "243" = c("1", "1", "0", "0", "0"),
#'                    "244" = c("1", "1", "0", "0", "0"),
#'                    "245" = c("1", "?", "0", "0", "0"))
#' rownames(dataTable) <- tree5$tip.label
#'
#' # list with two objects (tree and character table)
#' treeWdata <- geiger::treedata(tree5, dataTable)
#'
#' # Maximum parsimony reconstruction
#' rec_tab <- smooth_rec(bin_par_rec(treeWdata))
#'
#' # plotting and adding labels and legend
#' ape::plot.phylo(tree5, label.offset = 0.04)
#' niche_labels(tree5, rec_tab, height = 0.6)
#' niche_legend(position = "topleft", cex = 0.7)

niche_legend <- function(position, legend = c("Uncertain", "Present", "Not present"),
                         pch = 22, pt.bg = c("#969696", "#e41a1c", "#377eb8"),
                         col = "transparent", pt.cex = 2.2, bty = "n", ...) {
  if (missing(position)) {stop("Argument position needs to be defined")}
  cp <- class(position)[1]
  if (!cp %in% c("character", "numeric")) {
    stop("Argument 'position' needs to be of class character or numeric.")
  }

  # legend
  if (cp == "character") {
    legend(position, legend = legend, bty = bty, pch = pch,
           pt.cex = pt.cex, col = col, pt.bg = pt.bg, ...)
  } else {
    legend(x = position[1], y = position[2], legend = legend, bty = bty,
           pch = pch, pt.cex = pt.cex, col = col, pt.bg = pt.bg, ...)
  }
}


#' Legends for niche evolution labels in phylogenetic trees
#'
#' @param position (character or numeric) position of legend. If character,
#' part of the plot (e.g., "topleft"), see \code{\link[graphics]{legend}}. If
#' numeric, vector of two values indicating x and y position (e.g., c(0.1, 6)).
#' @param ancestor_line whether or not ancestor line was plotted.
#' Default = FALSE.
#' @param ancestor_legend (character) vector of length = two indicating the text
#' to identify environments with uncertain presence and true presence of the
#' species. Default = c("Uncertain", "Present").
#' @param evol_legend (character) vector of length = three indicating the text
#' to identify environments where niches have not changed, have retracted or
#' expanded. Default = c("No change", "Retraction", "Expansion").
#' @param ancestor_col vector of two colors to represent what is indicated in
#' \code{ancestor_legend}. Default = c("#d9d9d9", "#252525").
#' @param evol_col vector of three colors to represent what is indicated in
#' \code{evol_legend}. Default = c("#b2df8a", "#984ea3", "#4daf4a").
#' @param pch point type as in \code{\link[graphics]{points}}. Default = 22.
#' @param pt.cex size of symbol (points). Default = 2.2.
#' @param lty line type see \code{\link[graphics]{par}}. Default = 1.
#' @param lwd line width see \code{\link[graphics]{par}}. Default = 1.
#' @param cex size of all elements in legend see \code{\link[graphics]{par}}.
#' Default = 1.
#' @param bty legend border type. Default = "n".
#' @param ... Other arguments from function \code{\link[graphics]{legend}} other
#' than the ones described above.
#'
#' @importFrom graphics legend
#'
#' @export
#'
#' @usage
#' nichevol_legend(position, ancestor_line = FALSE,
#'   ancestor_legend = c("Uncertain", "Present"),
#'   evol_legend = c("No change", "Retraction", "Expansion"),
#'   ancestor_col = c("#d9d9d9", "#252525"),
#'   evol_col = c("#b2df8a", "#984ea3", "#4daf4a"),
#'   pch = 22, pt.cex = 2.2, lty = 1, lwd = 1, cex = 1, bty = "n", ...)
#'
#' @examples
#' # a simple tree
#' data("tree5", package = "nichevol")
#'
#' # a matrix of niche charactes (1 = present, 0 = absent, ? = unknown)
#' dataTable <- cbind("241" = rep("1", length(tree5$tip.label)),
#'                    "242" = rep("1", length(tree5$tip.label)),
#'                    "243" = c("1", "1", "0", "0", "0"),
#'                    "244" = c("1", "1", "0", "0", "0"),
#'                    "245" = c("1", "?", "0", "0", "0"))
#' rownames(dataTable) <- tree5$tip.label
#'
#' # list with two objects (tree and character table)
#' treeWdata <- geiger::treedata(tree5, dataTable)
#'
#' # Maximum parsimony reconstruction
#' rec_tab <- smooth_rec(bin_par_rec(treeWdata))
#'
#' # plotting and adding labels and legend
#' ape::plot.phylo(tree5, label.offset = 0.04)
#' nichevol_labels(tree5, rec_tab, height = 0.6)
#' nichevol_legend(position = "bottomleft", cex = 0.7)

nichevol_legend <- function(position, ancestor_line = FALSE,
                            ancestor_legend = c("Uncertain", "Present"),
                            evol_legend = c("No change", "Retraction", "Expansion"),
                            ancestor_col = c("#d9d9d9", "#252525"),
                            evol_col = c("#b2df8a", "#984ea3", "#4daf4a"),
                            pch = 22, pt.cex = 2.2, lty = 1, lwd = 1,
                            cex = 1, bty = "n", ...) {
  if (missing(position)) {stop("Argument 'position' needs to be defined")}
  cp <- class(position)[1]
  if (!cp %in% c("character", "numeric")) {
    stop("Argument 'position' needs to be of class 'character' or 'numeric'.")
  }

  # legend
  if (ancestor_line == TRUE) {
    if (cp == "character") {
      legend(position, legend = c(ancestor_legend, evol_legend),
             cex = cex, lty = c(lty, NA, NA, NA, NA),
             lwd = lwd, col = c("transparent", NA, NA, NA, NA), bty = bty, ...)

      legend(position, legend = c("                ", "", "", "", ""), bty = "n",
             pch = pch, pt.bg = c(NA, NA, evol_col), pt.cex = pt.cex, lty = lty,
             col = "transparent", cex = cex)

      legend(position, legend = c("                  ", "", "", "", ""),
             bty = "n", lty = c(lty, lty, NA, NA, NA), lwd = lwd, cex = cex,
             col = c(ancestor_col, NA, NA, NA))
    } else {
      legend(x = position[1], y = position[2], cex = cex, bty = bty,
             legend = c(ancestor_legend, evol_legend),
             lty = c(lty, NA, NA, NA, NA), lwd = lwd,
             col = c("transparent", NA, NA, NA, NA), ...)

      legend(x = position[1], y = position[2], cex = cex,
             legend = c("                ", "", "", "", ""), bty = "n",
             pch = pch, pt.bg = c(NA, NA, evol_col), pt.cex = pt.cex, lty = lty,
             col = "transparent")

      legend(x = position[1], y = position[2], cex = cex,
             legend = c("                  ", "", "", "", ""),
             bty = "n", lty = c(lty, lty, NA, NA, NA), lwd = lwd,
             col = c(ancestor_col, NA, NA, NA))
    }
  } else {
    if (cp == "character") {
      legend(position, legend = evol_legend, pch = pch, pt.bg = evol_col,
             col = "transparent", pt.cex = pt.cex, cex = cex, bty = bty, ...)
    } else {
      legend(x = position[1], y = position[2], legend = evol_legend, pch = pch,
             pt.bg = evol_col, col = "transparent", pt.cex = pt.cex, cex = cex,
             bty = bty, ...)
    }
  }
}
