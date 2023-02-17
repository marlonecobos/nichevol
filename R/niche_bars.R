#' PNG bar figures to represent ecological niches of distinct taxa
#'
#' @description niche_bars produces bar plots that represent species
#' ecological niches in one environmental variable. Bars are exported as png
#' figures to an output directory for posterior use.
#'
#' @param tree an object of class "phylo".
#' @param whole_rec_table matrix of environmental bins for all tips and nodes
#' derived from functions \code{\link{bin_par_rec}} or \code{\link{bin_ml_rec}}.
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
#' Ecological niches are represented in one environmental dimension with vertical
#' bars that indicate if the species is present, absent, or if its presence is
#' uncertain in the range of environmental conditions. Lower values of
#' environmental variables are represented in the left part of the bar, and the
#' opposite part of the bar represents higher values.
#'
#' @return
#' A folder named as in \code{output_directory} containing all bar figures
#' produced, as well as a legend to describe what is plotted.
#'
#' @importFrom graphics par plot polygon legend
#' @importFrom grDevices dev.off png
#'
#' @export
#'
#' @usage
#' niche_bars(tree, whole_rec_table, present = "1", unknown = "?",
#'            present_col = "#e41a1c", unknown_col = "#969696",
#'            absent_col = "#377eb8", width = 50, height = 5, res = 300,
#'            output_directory, overwrite = FALSE)
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
#' niche_bars(tree5, rec_tab, output_directory = file.path(tempdir(), "nichebars"))


niche_bars <- function(tree, whole_rec_table, present = "1", unknown = "?",
                       present_col = "#e41a1c", unknown_col = "#969696",
                       absent_col = "#377eb8", width = 50, height = 5, res = 300,
                       output_directory, overwrite = FALSE) {

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

  # organizing data
  nnames <- rownames(whole_rec_table); nnames <- nnames[!nnames %in% tlab]
  spnames <- c(tlab, nnames)
  bnames <- c(tlab, paste0("Node", nnames))

  tpol <- ncol(whole_rec_table)
  wpol <- 1 / tpol

  h_vertices <- seq(0, 1, wpol)
  v_vertices <- rep(c(0, 0.05), each = 2)

  dir.create(output_directory)

  barss <- sapply(1:nrow(whole_rec_table), function(j) {
    bar_name <- paste0(output_directory, "/", bnames[j],
                       "_bar.png")

    png(filename = bar_name, width = width, height = height, units = "mm",
        bg = "transparent", res = res)
    par(mar = rep(0, 4))
    plot(x = c(0, 1), y = c(0, 0.05), col = "transparent", axes = FALSE)

    polys <- sapply(1:(length(h_vertices) - 1), function(x) {
      if (as.character(whole_rec_table[j, x]) == unknown) {
        pcolor <- unknown_col
      } else {
        pcolor <- ifelse(as.character(whole_rec_table[j, x]) == present,
                         present_col, absent_col)
      }

      xs <- c(h_vertices[x], h_vertices[x + 1], h_vertices[x + 1], h_vertices[x])

      polygon(x = xs, y = v_vertices, col = pcolor, border = NA)
    })
    dev.off()
  })

  # legend
  png(filename = paste0(output_directory, "/0_Legend.png"), width = 50,
      height = 30, units = "mm", bg = "transparent", res = res)
  par(mar = rep(0, 4), cex = 1.2)
  plot(x = c(0, 0.5), y = c(0, 0.5), col = "transparent", axes = FALSE)
  legend("center", legend = c("Uncertain", "Present", "Not present"),
         bty = "n", pch = 22, pt.cex = 2.2, col = "transparent",
         pt.bg = c(unknown_col, present_col, absent_col))
  invisible(dev.off())
}



