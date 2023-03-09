#' Maps of niche reconstructions and changes detected
#'
#' @description map_nichevol produces a SpatRaster layer representing
#' geographic areas corresponding to environmental bins of niche or events of
#' niche evolution detected in reconstructions.
#'
#' @param whole_rec_table matrix of environmental bins for all tips and nodes
#' derived from functions \code{\link{bin_par_rec}} or \code{\link{bin_ml_rec}}.
#' @param variable a SpatRaster layer corresponding to the variable for which
#' the reconstruction was performed (represented in \code{whole_rec_table}).
#' @param return (character) type of result to return. Options are: "niche",
#' "evolution", or "nichevol" (a combination of both). Default = "niche". If
#' "niche", values correspond to that defined in \code{from}. See Value.
#' @param from (character) if \code{return} = "niche" tip or node for which layer
#' will be prepared, otherwise, initial node from which niche comparison will be
#' performed. See example.
#' @param to (character) valid if \code{return} = "evolution" or "nichevol".
#' Tip or node to compare against \code{from} to detected changes.
#' Default = NULL. See example.
#' @param id_unknown (logical) whether to identify areas of unknown or uncertain
#' change. Default = TRUE. See details.
#' @param verbose (logical) whether messages should be printed. Default = TRUE.
#'
#' @details
#' Mapping is done following Cobos et al. (2021) <doi:10.1111/jav.02868>.
#' This allows to represent geographic areas with environments where
#'  niche expanded, retracted, or stayed stable (evolution). Niche is
#' represented as presence, absence, or unknown.
#'
#' Defining \code{id_unknown} = TRUE allows to map areas where niche or niche
#' change are uncertain. \code{id_unknown} = FALSE returns NA in areas with these
#' characteristics, hence they will not be visible when plotting the resulting
#' map.
#'
#' @return
#' A SpatRaster object classified according to values of niche in
#' \code{whole_rec_table}, and/or according to niche changes detected in
#' comparisons between an ancestor and a tip, or another more recent ancestor.
#'
#' Options of values resulting from classifications are as follow:
#'
#' If \code{return} = "niche":
#' | ID|        category|
#' |:--|----------------|
#' |  0|          Absent|
#' | 10|         Unknown|
#' |100|         Present|
#'
#' If \code{return} = "evolution":
#' | ID|        category|
#' |:--|----------------|
#' |  0|          Stable|
#' |  1|   Expansion low|
#' |  3|  Expansion high|
#' |  2| Retraction high|
#' |  4|  Retraction low|
#' | 10|         Unknown|
#'
#' If \code{return} = "nichevol":
#' | ID|        category|
#' |:--|----------------|
#' |  0|          Stable|
#' |  1|   Expansion low|
#' |  3|  Expansion high|
#' | 10|         Unknown|
#' |100|         Present|
#' |102| Retraction high|
#' |104|  Retraction low|
#'
#'
#' @usage
#' map_nichevol(whole_rec_table, variable, return = "niche", from, to = NULL,
#'              id_unknown = TRUE, verbose = TRUE)
#'
#' @export
#'
#' @importFrom terra nlyr NAflag minmax classify unique
#'
#' @examples
#' # a tree
#' data("tree", package = "nichevol")
#'
#' # raster variable
#' temp <- terra::rast(system.file("extdata", "temp.tif", package = "nichevol"))
#'
#' # results from reconstruction
#' data("par_rec_table", package = "nichevol")
#'
#' # rename tree tips
#' tree$tip.label <- rownames(par_rec_table)[1:6]
#'
#' # check in plot
#' plot.phylo(tree, label.offset = 0.02)
#' nodelabels()
#' nichevol_labels(tree, par_rec_table)
#'
#' # mapping nichevol
#' nevol_map <- map_nichevol(whole_rec_table = par_rec_table, variable = temp,
#'                           return = "nichevol", from = "9", to = "RD 6933")
#'
#' terra::plot(nevol_map)

map_nichevol <- function(whole_rec_table, variable, return = "niche",
                         from, to = NULL, id_unknown = TRUE, verbose = TRUE) {
  if (missing(whole_rec_table)) {stop("Argument 'whole_rec_table' must be defined.")}
  if (missing(variable)) {stop("Argument 'variable' must be defined.")}
  if (missing(from)) {stop("Argument 'from' must be defined.")}
  if (!return %in% c("niche", "evolution", "nichevol")) {
    stop("Argument 'return' is not valid.")
  }
  if (return != "niche" & is.null(to)) {stop("Argument 'to' must be defined.")}

  n_list <- find_niche_lims(whole_rec_table, tip_or_node = from,
                            return = "values")

  if (return == "niche") {
    e_list <- NULL
  } else {
    e_list <- find_evol_lims(whole_rec_table, from = from, to = to,
                             return = "values", verbose = verbose)
  }


  return(nichevol_layer(variable, niche_list = n_list, evol_list = e_list,
                        return = return, id_unknown = id_unknown))
}



# test if a value is continuous
is_continuous <- function(x) {
  if (length(x) == 1) {
    return(TRUE)
  } else {
    x <- as.numeric(x)
    return(identical(x, seq(min(x), max(x), 1)))
  }
}

# find limits of niche or places where niche evolution occurs
all_limits <- function(value_vector, split_list, type = c("niche", "change"),
                       changes_list = NULL) {
  if (missing(value_vector)) {stop("Argument 'value_vector' must be defined.")}
  if (missing(split_list)) {stop("Argument 'split_list' must be defined.")}
  if (type[1] == "change" & is.null(changes_list)) {
    stop("Argument 'changes_list' must be defined to detect changes.")
  }

  l <- lapply(1:length(split_list), function(x) {
    if (type[1] == "change") {
      value_vector[] <- changes_list[[x]]
    }

    j <- split_list[[x]]
    if (length(j) > 0) {
      if (is_continuous(j)) {
        c(j[j == min(j)], j[j == max(j)])
      } else {
        ls <- vector()
        ad <- 1
        while (ad < (length(value_vector) - 1)) {
          ocon <- sum(c(ad, ad + 1) %in% j) >= 1
          if (value_vector[ad] != value_vector[ad + 1] & ocon) {
            if (sum(j == ad) >= 1) {
              ls[ad] <- j[j == ad]
            } else {
              ls[ad] <- j[j == (ad + 1)]
            }
          }
          ad <- ad + 1
        }

        lsnona <- na.omit(ls)
        gt <- j[j %in% lsnona]

        if (max(gt) != max(j)) {
          gt <- c(gt, j[j == max(j)])

          if (length(gt) == 2) {
            gt <- c(gt, gt[2])
          }
        }
        if (1 %in% j) {
          gt <- c(j[1], gt)
        }
        return(gt)
      }
    } else {
      return(integer())
    }
  })

  if (type[1] == "change") {
    names(l) <- names(changes_list)
  } else {
    names(l) <- names(split_list)
  }

  return(l)
}

# find values of previous limits
critical_limits <- function(all_limits, return = c("values", "bins")) {
  if (missing(all_limits)) {stop("Argument 'all_limits' must be defined.")}

  if (return[1] == "values") {
    vals <- lapply(all_limits, function(x) {
      if (length(x) > 0) {
        v <- rep(1:2, ceiling(length(x) / 2))
        vals <- sapply(1:length(x), function (y) {
          strsplit(names(x)[y], " to ")[[1]][v[y]]
        })

        mat <- matrix(as.numeric(vals), ncol = 2, byrow = TRUE)
        colnames(mat) <- c("From", "To")

        return(mat)
      } else {
        mat <- matrix(ncol = 2); colnames(mat) <- c("From", "To")
        return(mat[!is.na(mat[, 1]), ])
      }
    })
  } else {
    vals <- lapply(all_limits, function(x) {
      mat <- matrix(x, ncol = 2, byrow = TRUE)
      colnames(mat) <- c("From", "To")
      return(mat)
    })
  }

  return(vals)
}

# find values of niche evolution using previous functions
find_evol_lims <- function(whole_rec_table, from, to, return = c("values", "bins"),
                           present = "1", absent = "0", unknown = "?",
                           verbose = TRUE) {
  if (missing(whole_rec_table)) {
    stop("Argument 'whole_rec_table' must be defined.")
  }
  if (missing(from)) {stop("Argument 'from' must be defined.")}
  if (missing(to)) {stop("Argument 'to' must be defined.")}

  # finding
  to_find <- whole_rec_table[c(from, to), ]
  tn <- whole_rec_table[1, ]

  ## changes
  expansion <- sapply(1:ncol(to_find), function(z) {
    from <- to_find[1, z]
    to <- to_find[2, z]
    return(ifelse(from == absent & to == present, TRUE, FALSE))
  })

  retraction <- sapply(1:ncol(to_find), function(z) {
    from <- to_find[1, z]
    to <- to_find[2, z]
    return(ifelse(from == present & to == absent, TRUE, FALSE))
  })

  stable <- !expansion & !retraction

  uns <- list(expansion = expansion, retraction = retraction, stable = stable)

  ## limits of changes
  lims <- lapply(uns, function(x) {
    tn[] <- x
    return(which(tn == TRUE))
  })

  alims <- all_limits(value_vector = tn, split_list = lims, type = "change",
                      changes_list = uns)

  ## critical values and bin numbers
  avals <- critical_limits(alims, return = "values")

  abins <- critical_limits(alims, return = "bins")

  # results
  if (sum(stable) >= length(tn) & verbose == TRUE) {
    message("No changes detected in comparisons.")
  }

  if (return[1] == "values") {
    return(avals)
  } else {
    return(abins)
  }
}

# find values of niche limits using previous functions
find_niche_lims <- function(whole_rec_table, tip_or_node, return = c("values", "bins"),
                            present = "1", absent = "0", unknown = "?") {
  if (missing(whole_rec_table)) {
    stop("Argument 'whole_rec_table' must be defined.")
  }
  if (missing(tip_or_node)) {stop("Argument 'tip_or_node' must be defined.")}
  if (length(tip_or_node) > 1) {stop("Argument 'tip_or_node' must be of length = 1.")}

  # finding
  tn <- whole_rec_table[tip_or_node, ]
  uns <- unique(tn)

  ns <- ifelse(uns == present, "present",
               ifelse(uns == absent, "absent", "unknown"))

  ## limits
  lims <- lapply(uns, function(x) {which(tn == x)})
  names(lims) <- ns

  alims <- all_limits(value_vector = tn, split_list = lims, type = "niche")

  ## critical values and bin numbers
  avals <- critical_limits(all_limits = alims, return = "values")

  abins <- critical_limits(all_limits = alims, return = "bins")

  # results
  if (return[1] == "values") {
    return(avals)
  } else {
    return(abins)
  }
}

# creates a raster layer to show areas of niche and niche evolution
nichevol_layer <- function(variable, niche_list,
                           return = c("niche", "evolution", "nichevol"),
                           evol_list = NULL, id_unknown = TRUE) {
  if (missing(variable)) {stop("Argument 'variable' must be defined.")}
  if (class(variable)[1] != "SpatRaster") {
    stop("Argument 'variable' must be an object of class 'SpatRaster'.")
  }
  if (missing(niche_list)) {stop("Argument 'niche_list' must be defined.")}
  if (!return[1] %in% c("nichevol", "niche", "evolution")) {
    stop("Arguments 'return' is not valid.")
  }
  if (return[1] %in% c("nichevol", "evolution") & is.null(evol_list)) {
    stop("Arguments 'evol_list' must be defined to identify changes in niche.")
  }
  if (terra::nlyr(variable) != 1) {
    stop("'variable' must contain only one layer.")
  }

  # preparing objects
  vmm <- terra::minmax(variable)
  nav <- terra::NAflag(variable)

  unk <- ifelse(id_unknown == TRUE, 10, nav)

  # reclassification and details
  if (return[1] %in% c("niche", "nichevol")) {
    ns <- names(niche_list)
    nvs <- ifelse(ns == "present", 100, ifelse(ns == "absent", 0, unk))
    rec_l <- lapply(1:length(niche_list), function(x) {
      cbind(niche_list[[x]], becomes = nvs[x])
    })
    recm <- do.call(rbind, rec_l)
    recm <- recm[order(recm[, 1]), ]

    for (i in 2:nrow(recm)) {
      recm[i, 1] <- recm[i - 1, 2]
    }

    # NAs for values outside ranges
    nnas <- range(c(recm[, 1:2]))

    recm <- rbind(c(From = vmm[1], To = nnas[1], becomes = nav), recm,
                  c(From = nnas[2], To = vmm[2], becomes = nav))

    ni <- terra::classify(variable, recm, others = nav)
  }

  if (return[1] %in% c("evolution", "nichevol")) {
    pres <- niche_list$present
    ns <- names(evol_list)
    nvs <- ifelse(ns == "expansion", 1, ifelse(ns == "retraction", 2, 0))
    rec_l <- lapply(1:length(evol_list), function(x) {
      etab <- evol_list[[x]]
      if (nrow(etab) > 0) {
        if(nrow(etab) > 1) {
          vs <- rep(nvs[x], nrow(etab))
          if (nvs[x] == 1) {
            vs <- seq(1, 10, 2)[1:length(vs)]
          }
          if (nvs[x] == 2) {
            vs <- seq(2, 10, 2)[1:length(vs)]
          }
          cbind(etab, becomes = vs)
        } else {
          nvst <- ifelse(etab[1, 2] <= pres[1, 1], 0, ifelse(nvs[x] == 0, 0, 2))
          vs <- nvs[x] + nvst
        }
        cbind(etab, becomes = vs)
      } else {
        cbind(etab, becomes = numeric())
      }
    })

    recm <- do.call(rbind, rec_l)
    recm <- recm[order(recm[, 1]), ]

    for (i in 2:nrow(recm)) {
      recm[i, 1] <- recm[i - 1, 2]
    }

    # NAs for values outside ranges
    nnas <- range(c(recm[, 1:2]))

    recm <- rbind(c(From = vmm[1], To = nnas[1], becomes = nav), recm,
                  c(From = nnas[2], To = vmm[2], becomes = nav))

    ev <- terra::classify(variable, recm, others = nav)

    if (return[1] == "nichevol") {
      ni <- ni + ev
    } else {
      ni <- ev
    }
  }

  # results with categories
  if (return[1] == "evolution") {
    cates <- data.frame(ID = c(0, 1, 3, 2, 4, 10),
                        category = c("Stable", "Expansion low",
                                     "Expansion high", "Retraction high",
                                     "Retraction low", "Unknown"))
  } else {
    cates <- data.frame(ID = c(100, 0, 1, 3, 102, 104, 10),
                        category = c("Present", "Absent", "Expansion low",
                                     "Expansion high", "Retraction high",
                                     "Retraction low", "Unknown"))
  }

  cates_ni <- cates[cates$ID %in% terra::unique(ni)[, 1], ]

  levels(ni) <- cates_ni

  return(ni)
}
