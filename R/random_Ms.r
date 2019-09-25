random_polygons <- function(polygon, n_polygons = 5, n_vertices = 25,
                            minimum_distance = 10, length_threshold = 5,
                            buffer_distance = 0, name = "Random_polygons",
                            overwrite = FALSE) {
  if (missing(polygon)) {
    stop("Argument polygon is necessary to perform the analysis.")
  }
  if (overwrite == FALSE & file.exists(paste0(name, ".shp"))) {
    stop("shapefile already exists, to replace it use overwrite = TRUE.")
  }

  m_dist <- minimum_distance / 111.32
  l_thres <- length_threshold / 111.32
  ext <- raster::extent(polygon)
  crsxy <- polygon@proj4string
  buff_dist <- buffer_distance / 111.32

  numx <- seq(ext[1], ext[2], m_dist)
  numy <- seq(ext[3], ext[4], m_dist)

  xyp <- cbind(sample(numx, 10000, replace = TRUE),
               sample(numy, 10000, replace = TRUE))

  xysp <- sp::SpatialPointsDataFrame(xyp, data = as.data.frame(xyp), proj4string = crsxy)
  xyin <- xysp[polygon, ]

  seeds <- sample(1:1000000, size = n_polygons)

  pols <- lapply(1:n_polygons, function (x) {
    set.seed(seeds[x])
    np <- sample(3:n_vertices, 1)
    xysam <- xyin[sample(nrow(xyin), np), ]
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

    hb <- suppressWarnings(rgeos::gBuffer(hulls, width = buff_dist))
    pfin <- sp::SpatialPolygonsDataFrame(hb, data = data.frame(RD = rep(x, length(hb))),
                                         match.ID = FALSE)

    pfin <- sp::spTransform(pfin, CRSobj = crsxy)
    cat(x, "of", n_polygons, "\n")
    return(pfin)
  })

  r_polygons <- do.call(sp::rbind.SpatialPolygonsDataFrame, pols)
  rp_data <- r_polygons@data
  r_polygons <- rgeos::gIntersection(r_polygons, polygon, byid = TRUE,
                                     drop_lower_td = TRUE)
  r_polygons <- sp::SpatialPolygonsDataFrame(r_polygons, data = rp_data,
                                             match.ID = FALSE)
  rgdal::writeOGR(obj = r_polygons, dsn = ".", layer = name,
                  driver = "ESRI Shapefile", overwrite_layer = overwrite)

  return(r_polygons)
}
