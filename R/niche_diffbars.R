#' PNG bar figures for representing niche evolution
#'
#' @description nichevol_bars produces bar plots that represent how
#' species' niches (considering one environmental variable at a time) have evolved.
#' Bars are exported as png figures to an output directory for posterior use.
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
#' is present. Default = "#252525".
#' @param unknown_col color for line representing environments where the species
#' presence is unknown (uncertain). Default = "#d9d9d9".
#' @param no_change_col color for area of the bar representing environments where
#' no change has been detected. Default = "#b2df8a".
#' @param retraction_col color for area of the bar representing environments where
#' niche retraction has been detected. Default = "#984ea3".
#' @param expansion_col color for area of the bar representing environments where
#' niche expansion has been detected. Default = "#4daf4a".
#' @param width (numeric) width of the device in mm to be passed to the
#' \code{\link[grDevices]{png}} function. Default = 50.
#' @param height (numeric) height of the device in mm to be passed to the
#' \code{\link[grDevices]{png}} function. Default = 5.
#' @param res (numeric) nominal resolution in ppi to be passed to the
#' \code{\link[grDevices]{png}} function. Default = 300.
#' @param output_directory (character) name of the folder in which results will
#' be written. The directory will be created as part of the process.
#' @param overwrite (logical) whether or not to overwrite existing results in
#' \code{output_directory}. Default = FALSE.
#'
#' @details
#' Evolution of ecological niches is represented in one environmental dimension
#' with horizontal bars indicating if the niche of the descendant has expanded,
#' retracted, or has not changed compared to its ancestor. Lower values of
#' environmental variables are represented in the left part of the bar, higher
#' values at the right.
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
#' @return
#' A folder named as in \code{output_directory} containing all bar figures
#' produced, as well as a legend to describe what is plotted.
#'
#' @importFrom utils combn
#' @importFrom graphics par plot polygon legend lines
#' @importFrom grDevices dev.off png
#'
#' @export
#'
#' @usage
#' nichevol_bars(tree, whole_rec_table, ancestor_line = FALSE,
#'               present = "1", absent = "0", unknown = "?",
#'               present_col = "#252525", unknown_col = "#d9d9d9",
#'               no_change_col = "#b2df8a", retraction_col = "#984ea3",
#'               expansion_col = "#4daf4a", width = 50, height = 5,
#'               res = 300, output_directory, overwrite = FALSE)
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
#' # the running (before running, define a working directory)
#' nichevol_bars(tree5, rec_tab, output_directory = file.path(tempdir(), "evolbars"))

nichevol_bars <- function(tree, whole_rec_table, ancestor_line = FALSE,
                          present = "1", absent = "0", unknown = "?",
                          present_col = "#252525", unknown_col = "#d9d9d9",
                          no_change_col = "#b2df8a", retraction_col = "#984ea3",
                          expansion_col = "#4daf4a", width = 50, height = 5,
                          res = 300, output_directory, overwrite = FALSE) {

  # testing for potential errors
  if (missing(tree)) {stop("Argument 'tree' is needed to perform the analyses.")}
  if (missing(whole_rec_table)) {stop("Argument 'whole_rec_table' needs to be defined.")}
  if ("LogLik" %in% rownames(whole_rec_table)) {
    whole_rec_table <- whole_rec_table[1:(nrow(whole_rec_table) - 3), ]
  }
  if (missing(output_directory)) {
    stop("Argument 'output_directory' is missing.")
  } else {
    if (overwrite == FALSE & dir.exists(output_directory)) {
      stop("'output_directory' already exists, to replace it use 'overwrite' = TRUE.")
    }
    if (overwrite == TRUE & dir.exists(output_directory)) {
      unlink(x = output_directory, recursive = TRUE, force = TRUE)
    }
  }

  # par settings
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  # reorganizing character table
  tlab <- tree$tip.label
  nrt <- nrow(whole_rec_table)
  rns <- c(tlab, rownames(whole_rec_table)[(length(tlab) + 1):nrt])
  whole_rec_table <- rbind(whole_rec_table[tlab, ],
                           whole_rec_table[(length(tlab) + 1):nrt, ])
  rownames(whole_rec_table) <- rns

  edges <- tree$edge

  # preparing plotting parameters
  spp <- rownames(whole_rec_table)

  tpol <- ncol(whole_rec_table)
  wpol <- 1 / tpol

  h_vertices <- seq(0, 1, wpol)
  v_vertices <- rep(c(0, 0.05), each = 2)
  y_line <- rep(mean(c(0, 0.05)), 2)

  # combinations for comparisons
  all_comb <- cbind(spp[edges[, 1]], spp[edges[, 2]])
  all_nam <- paste("From_node", all_comb[, 1], "to", all_comb[, 2], sep = "_")

  # plotting data
  dir.create(output_directory) # folder

  # comparisons and plots
  comp_list <- lapply(1:nrow(all_comb), function(x){
    comp <- sapply(1:ncol(whole_rec_table), function(z) {
      from <- whole_rec_table[all_comb[x, 1], z]
      to <- whole_rec_table[all_comb[x, 2], z]
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

    # infor for lines
    linesp <- whole_rec_table[all_comb[x, 1], ]

    # bar creation
    bar_name <- paste0(output_directory, "/", all_nam[x], "_bar.png")

    png(filename = bar_name, width = width, height = height, units = "mm",
        bg = "transparent", res = res)
    par(mar = rep(0, 4))
    plot(x = c(0, 1), y = c(0, 0.05), col = "transparent", axes = FALSE)

    poly_lines <- sapply(1:(length(h_vertices) - 1), function(y) {
      # polygons
      if (comp[y] == "nc") {
        pcolor <- no_change_col
      } else {
        pcolor <- ifelse(comp[y] == "gain", expansion_col, retraction_col)
      }

      xs <- c(h_vertices[y], h_vertices[y + 1], h_vertices[y + 1], h_vertices[y])
      polygon(x = xs, y = v_vertices, col = pcolor, border = NA)

      # lines
      if (ancestor_line == TRUE) {
        if (linesp[y] == unknown) {
          pcolor <- unknown_col
        } else {
          pcolor <- ifelse(linesp[y] == present, present_col, "transparent")
        }

        xs <- c(h_vertices[y], h_vertices[y + 1])
        lines(x = xs, y = y_line, col = pcolor, lty = 1, lwd = 1.7)
      }
    })
    dev.off()

  })

  # legend
  png(filename = paste0(output_directory, "/0_Legend.png"), width = 50, height = 50,
      units = "mm", bg = "transparent", res = res)
  par(mar = rep(0, 4), cex = 1.2)
  plot(x = c(0, 0.5), y = c(0, 0.5), col = "transparent", axes = FALSE)

  if (ancestor_line == TRUE) {
    legend("center",
           legend = c("Uncertain", "Present", "No change", "Retraction", "Expansion"),
           lty = c(1, NA, NA, NA, NA), lwd = 2,
           col = c("transparent", NA, NA, NA, NA), bty = "n")

    legend("center", legend = c("                ", "", "", "", ""), bty = "n",
           pch = 22, pt.bg = c(NA, NA, no_change_col, retraction_col, expansion_col),
           pt.cex = 2.2, lty = 1, col = "transparent")

    legend("center", legend = c("                  ", "", "", "", ""),
           bty = "n", lty = c(1, 1, NA, NA, NA), lwd = 2,
           col = c(unknown_col, present_col, NA, NA, NA))
  } else {
    legend("center", legend = c("No change", "Retraction", "Expansion"),
           pch = 22, pt.bg = c(no_change_col, retraction_col, expansion_col),
           col = "transparent", pt.cex = 2.2, bty = "n")
  }

  invisible(dev.off())
}
