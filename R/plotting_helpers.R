#' Labels to represent niches of tips and ancestors
#'
#' @description niche_labels helps in adding bar-type labels that represent how
#' species niches are compared to others.
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
#' the species is present. Default = "red4".
#' @param unknown_col color for area of the bar representing environments where
#' the species presence is unknown (uncertain). Default = "lightblue".
#' @param absent_col color for area of the bar representing environments where
#' no change has been detected. Default = "royalblue1".
#'
#' @importFrom graphics plot polygon
#'
#' @export

niche_labels <- function(tree, whole_rec_table, label_type = "tip_node",
                         tip_offset = 0.015, present = "1", unknown = "?",
                         present_col = "red4", unknown_col = "lightblue",
                         absent_col = "royalblue1") {
  if (missing(tree)) {stop("Argument tree needs to be defined.")}
  if (missing(whole_rec_table)) {stop("Argument whole_rec_table needs to be defined.")}
  if ("LogLik" %in% rownames(whole_rec_table)) {
    whole_rec_table <- whole_rec_table[1:(nrow(whole_rec_table) - 3), ]
  }

  # reorganizing character table
  tlab <- tree$tip.label
  nrt <- nrow(whole_rec_table)
  rns <- c(tlab, rownames(whole_rec_table)[length(tlab):nrt])
  whole_rec_table <- rbind(whole_rec_table[tlab, ], whole_rec_table[length(tlab):nrt, ])
  rownames(whole_rec_table) <- rns

  # getting info from plot
  tp_info <- get("last_plot.phylo", envir = .PlotPhyloEnv)
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
  wt <- ((max(yy) / tp_info$Ntip) / 10) * 4
  wpol <- wt / tpol

  h_vertices <- seq(0, wt, wpol)

  # plotting bars
  if (label_type %in% c("tip", "node", "tip_node")) {
    if (label_type %in% c("tip", "tip_node")) {
      barss <- sapply(tpos, function(j) {
        ys <- yy[j] - (wt / 2)
        hver <- ys + h_vertices
        xs <- xx[j] + tip_offset; xs1 <- xs - 0.005; xs2 <- xs + 0.005
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
        xs <- xx[j]; xs1 <- xs - 0.005; xs2 <- xs + 0.005
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
    stop("Argument label_type is not correct, see function's help.")
  }
}


#' Labels to represent changes of between ancestors and descendants
#'
#' @description nichevol_labels helps in adding bar-type labels that represent how
#' species niches changed from ancestors to descendants.
#'
#' @param tree an object of class "phylo".
#' @param whole_rec_table matrix of reconstructed bins for nodes and species
#' derived from a process of maximum parsimony reconstruction.
#' @param ancestor_line controls whether ancestor line is plotted.
#' Default = FALSE.
#' @param present (character) code indicating environmental bins in which the
#' species is present. Default = "1".
#' @param absent (character) code indicating environmental bins in which the
#' species is absent. Default = "0".
#' @param unknown (character) code indicating environmental bins in which the
#' species presence is unknown (uncertain). Default = "?".
#' @param present_col color for line representing environments where the species
#' is present. Default = "grey10".
#' @param unknown_col color for line representing environments where the species
#' presence is unknown (uncertain). Default = "white".
#' @param no_change_col color for area of the bar representing environments where
#' no change has been detected. Default = "grey90".
#' @param retraction_col color for area of the bar representing environments where
#' niche retraction has been detected. Default = "dodgerblue3".
#' @param expansion_col color for area of the bar representing environments where
#' niche expansion has been detected. Default = "green1".
#'
#' @importFrom graphics plot polygon lines
#'
#' @export

nichevol_labels <- function(tree, whole_rec_table, ancestor_line = FALSE,
                            present = "1", absent = "0", unknown = "?",
                            present_col = "grey10", unknown_col = "orange",
                            no_change_col = "grey90", retraction_col = "dodgerblue3",
                            expansion_col = "green1") {
  if (missing(tree)) {stop("Argument tree needs to be defined.")}
  if (missing(whole_rec_table)) {stop("Argument whole_rec_table needs to be defined.")}
  if ("LogLik" %in% rownames(whole_rec_table)) {
    whole_rec_table <- whole_rec_table[1:(nrow(whole_rec_table) - 3), ]
  }

  # reorganizing character table
  tlab <- tree$tip.label
  nrt <- nrow(whole_rec_table)
  rns <- c(tlab, rownames(whole_rec_table)[length(tlab):nrt])
  whole_rec_table <- rbind(whole_rec_table[tlab, ], whole_rec_table[length(tlab):nrt, ])
  rownames(whole_rec_table) <- rns

  # getting info from plot
  tp_info <- get("last_plot.phylo", envir = .PlotPhyloEnv)
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
  wt <- ((max(yy) / tp_info$Ntip) / 10) * 4
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

    xs <- xb[x]; xs1 <- xs - 0.005; xs2 <- xs + 0.005
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
#' numeric, vector of two values indicating x and y postion (e.g., c(0.1, 6)).
#' @param legend (character) vector of length = three indicating the text to
#' identify environments with uncertain presence, presence, and absence of the
#' species. Default = c("Uncertain", "Present", "Not present").
#' @param pch point type as in \code{\link[graphics]{points}}. Default = 22.
#' @param pt.bg colors to represent what is in \code{legend}.
#' Default = c("lightblue", "red4", "royalblue1").
#' @param col border of symbol (points). Default = "transparent".
#' @param pt.cex size of symbol (points). Default = 2.2.
#' @param bty legend border type. Default = "n".
#' @param ... Other arguments from function \code{\link[graphics]{legend}} other
#' than the ones described above.
#'
#' @importFrom graphics legend
#'
#' @export

niche_legend <- function(position, legend = c("Uncertain", "Present", "Not present"),
                         pch = 22, pt.bg = c("lightblue", "red4", "royalblue1"),
                         col = "transparent", pt.cex = 2.2, bty = "n", ...) {
  if (missing(position)) {stop("Argument position needs to be defined")}
  cp <- class(position)[1]
  if (!cp %in% c("character", "numeric")) {
    stop("Argument position needs to be of class character or numeric.")
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
#' numeric, vector of two values indicating x and y postion (e.g., c(0.1, 6)).
#' @param legend (character) vector of length = five indicating the text to
#' identify environments with uncertain presence and presence of the species,
#' as well as, areas where niches have not changed, have retracted or expanded.
#' Default = c("Uncertain", "Present", "No change", "Retraction", "Expansion").
#' Order must be mantained.
#' @param pch point type as in \code{\link[graphics]{points}}. Default = 22.
#' @param pt.cex size of symbol (points). Default = 2.2.
#' @param col vector of five colors to represent what is in legend.
#' Default = c("orange", "grey10", "grey90", "dodgerblue3", "green1").
#' @param lty line type see \code{\link[graphics]{par}}. Default = 1.
#' @param lwd line width see \code{\link[graphics]{par}}. Default = 1.
#' @param cex size of all elements in legend see \code{\link[graphics]{par}}.
#' Deafult = 1.
#' @param bty legend border type. Default = "n".
#' @param ... Other arguments from function \code{\link[graphics]{legend}} other
#' than the ones described above.
#'
#' @importFrom graphics legend
#'
#' @export

nichevol_legend <- function(position, legend = c("Uncertain", "Present", "No change",
                                                 "Retraction", "Expansion"),
                            pch = 22, pt.cex = 2.2, lty = 1, lwd = 1,
                            col = c("orange", "grey10", "grey90", "dodgerblue3",
                                    "green1"), cex = 1, bty = "n", ...) {
  if (missing(position)) {stop("Argument position needs to be defined")}
  cp <- class(position)[1]
  if (!cp %in% c("character", "numeric")) {
    stop("Argument position needs to be of class character or numeric.")
  }

  # legend
  if (cp == "character") {
    legend(position, legend = legend, cex = cex, lty = c(lty, NA, NA, NA, NA),
           lwd = lwd, col = c("transparent", NA, NA, NA, NA), bty = bty, ...)

    legend(position, legend = c("                ", "", "", "", ""), bty = "n",
           pch = pch, pt.bg = c(NA, NA, col[3], col[4], col[5]),
           pt.cex = pt.cex, lty = lty, col = "transparent", cex = cex)

    legend(position, legend = c("                  ", "", "", "", ""), bty = "n",
           lty = c(lty, lty, NA, NA, NA), lwd = lwd, cex = cex,
           col = c(col[1], col[2], NA, NA, NA))
  } else {
    legend(x = position[1], y = position[2], cex = cex, legend = legend,
           lty = c(lty, NA, NA, NA, NA), lwd = lwd,
           col = c("transparent", NA, NA, NA, NA), bty = bty, ...)

    legend(x = position[1], y = position[2], cex = cex,
           legend = c("                ", "", "", "", ""), bty = "n",
           pch = pch, pt.bg = c(NA, NA, col[3], col[4], col[5]),
           pt.cex = pt.cex, lty = lty, col = "transparent")

    legend(x = position[1], y = position[2], cex = cex,
           legend = c("                  ", "", "", "", ""), bty = "n",
           lty = c(lty, lty, NA, NA, NA), lwd = lwd,
           col = c(col[1], col[2], NA, NA, NA))
  }
}
