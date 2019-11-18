#' Generation of random polygons in a given area
#'
#' @description random_polygons creates polygons of random size and complexity within
#' a given SpatialPolygonsDataFrame, trying to fill the area with the resultant polygons
#' in at least nine quadrants. This is designed to simulate virtual species' Ms (aka
#' "training" or "background" regions).
#'
#' @param polygon SpatialPolygonsDataFrame object. CRS WGS84 is required.
#' @param style (character) algorithm to be used when creating polygons. Options
#' are: "TR" for vertices randomly located across the entire area; and "BR" for
#' vertices placed randomly across the entire area and in nine blocks derived
#' from dividing the area in equal number of random points. Default = "TR".
#' @param n_polygons (numeric) number of polygons to be created; default = 100.
#' @param n_vertices (numeric) maximum number of vertices for polygons.
#' @param minimum_distance (numeric) approximate minimum distance in km for
#' separation among vertices. Default = 10.
#' @param length_threshold (numeric) approximate distance in km for producing
#' concavity in polygons. Default = 5.
#' @param buffer_distance (numeric) approximate distance in km to buffer
#' resultant polygons. Default = 0.
#' @param overwrite (logical) whether or not to overwrite previous results.
#' @param output_directory (character) name of the folder in which results will
#' be written. Default = "Random_polygons".
#'
#' @details
#' Distances are approximate because 1 decimal degree is assumed to equal 111.32
#' km.
#'
#' Style for random polygons "BR" may help to get smaller and more uniformly
#' distributed across the area.
#'
#' @importFrom grDevices chull
#'
#' @export
#'
#' @return
#' A folder named as in \code{output_directory} containing all resultant
#' shapefiles of random polygons. Results will also be returned as a list of
#' spatial objects.

random_polygons <- function(polygon, style = "TR", n_polygons = 100, n_vertices = 25,
                            minimum_distance = 10, length_threshold = 5,
                            buffer_distance = 0, overwrite = FALSE,
                            output_directory = "Random_polygons") {
  if (missing(polygon)) {
    stop("Argument 'polygon' is necessary to perform the analysis.")
  }
  if (overwrite == FALSE & dir.exists(output_directory)) {
    stop("'output_directory' already exists. To replace it use overwrite = TRUE.")
  }

  cat("\nPreparing data...\n")
  m_dist <- minimum_distance / 111.32
  l_thres <- length_threshold / 111.32
  ext <- raster::extent(polygon)
  crsxy <- polygon@proj4string
  buff_dist <- buffer_distance / 111.32

  numx <- seq(ext[1], ext[2], m_dist)
  numy <- seq(ext[3], ext[4], m_dist)

  xyp <- cbind(sample(numx, 10000, replace = TRUE),
               sample(numy, 10000, replace = TRUE))

  xysp <- sp::SpatialPointsDataFrame(xyp, data = as.data.frame(xyp),
                                     proj4string = crsxy)
  xyin <- xysp[polygon, ]@data

  seeds <- sample(1:1000000, size = n_polygons)

  cat("\nCreating polygons:\n")
  pols <- lapply(1:n_polygons, function (x) {
    id <- sample(1:10, 1)
    xyg <- make_9blocks(xyin)
    set.seed(seeds[x])

    if (style == "TR") {
      xyi <- xyg[, 1:2]
    } else {
      if (id == 10) {xyi <- xyg[, 1:2]} else {xyi <- xyg[xyg[, 3] == id, 1:2]}
    }

    np <- sample(3:n_vertices, 1)
    xysam <- xyi[sample(nrow(xyi), np), ]
    xysam <- sp::SpatialPointsDataFrame(xysam, data = as.data.frame(xysam),
                                        proj4string = crsxy)
    type <- sample(1:2, 1)

    if (type == 1) {
      sppoints <- sf::st_multipoint(as.matrix(xysam@coords))
      sf_points <- sf::st_sf(sf::st_sfc(sppoints, crs = crsxy@projargs))
      concavehull <- concaveman::concaveman(sf_points, length_threshold = l_thres)
      hulls <- sf::as_Spatial(sf::st_zm(concavehull$polygons))
    } else {
      covexhull <- chull(xysam@coords)
      coord_pol <- xysam@coords[c(covexhull, covexhull[1]), ]
      poly_list <- list(sp::Polygons(list(sp::Polygon(coord_pol)),
                                     ID = 1))
      hulls <- sp::SpatialPolygons(poly_list, proj4string = crsxy)
    }

    if (buffer_distance > 0) {
      hulls <- suppressWarnings(rgeos::gBuffer(hulls, width = buff_dist))
    }

    pfin <- sp::SpatialPolygonsDataFrame(hulls, data = data.frame(RD = x),
                                         match.ID = FALSE)

    cat("\t", x, "of", n_polygons, "polygons\n")
    return(pfin)
  })

  r_polygons <- do.call(sp::rbind.SpatialPolygonsDataFrame, pols)
  rp_data <- r_polygons@data
  r_polygons <- rgeos::gIntersection(r_polygons, polygon, byid = TRUE,
                                     drop_lower_td = TRUE)
  r_polygons <- sp::SpatialPolygonsDataFrame(r_polygons, data = rp_data,
                                             match.ID = FALSE)

  # writing results
  cat("\nwritting results:\n")
  if (overwrite == TRUE & dir.exists(output_directory)) {
    unlink(output_directory, recursive = TRUE)
  }
  dir.create(output_directory)
  f_polygons <- lapply(1:nrow(rp_data), function(x) { ### change
    n_nam <- paste0("r_polygon", x)
    rgdal::writeOGR(obj = r_polygons[x, ], dsn = output_directory,
                    layer = n_nam, driver = "ESRI Shapefile")
    cat("\t", x, "of", nrow(rp_data), "polygons\n")
    return(r_polygons[x, ])
  })

  # returning results
  return(f_polygons)
}
