histograms_env <- function(M_folder, M_format = "shp", occ_folder, longitude_col, latitude_col,
                           vars_folder, vars_format = "GTiff", CL_lines = c(95, 99), col, round = FALSE,
                           round_names, multiplication_factor = 1, output_folder = "Histogram_check"){
   # checking for potential errors
  if (missing(col)) {
    col <- sort(gray.colors(length(CL_lines) + 1), decreasing = TRUE)[1:length(CL_lines)]
  }

  # formats
  if (M_format == "shp") {
    M_patt <- ".shp$"
    subs <- ".shp"
  }
  if (M_format == "ascii") {
    M_patt <- ".asc$"
    subs <- ".asc"
  }
  if (M_format == "GTiff") {
    M_patt <- ".tif$"
    subs <- ".tif"
  }
  if (M_format == "EHdr") {
    M_patt <- ".bil$"
    subs <- ".tif"
  }
  if (vars_format == "ascii") {v_patt <- ".asc$"}
  if (vars_format == "GTiff") {v_patt <- ".tif$"}
  if (vars_format == "EHdr") {v_patt <- ".bil$"}

  # needed data to start
  mlist <- gsub(subs, "",list.files(path = M_folder, pattern = paste0(".", M_format, "$")))
  occlist <- list.files(path = occ_folder, pattern = ".csv$", full.names = TRUE)
  spnames <- gsub(".csv", "", list.files(path = occ_folder, pattern = ".csv$"))
  variables <- raster::stack(list.files(path = vars_folder, pattern = v_patt, full.names = TRUE))

  if (round == TRUE) {
    rounds <- round(variables[[round_names]] * multiplication_factor)
    var_names <- names(variables)
    noround <- var_names[!var_names %in% round_names]
    variables <- raster::stack(variables[[noround]], rounds)
  }

  # directory for results
  dir.create(output_folder)

  for (i in 1:dim(variables)[3]) {
    # data
    df_layer <- list()
    occ_dfs <- list()
    yocc <- list()

    for (j in 1:length(occlist)) {
      ## M shapefiles
      M <- rgdal::readOGR(dsn = M_folder, layer = mlist[j], verbose = FALSE)

      ## occurrences
      occ <- read.csv(occlist[j])

      # processing
      ## get values of variables in M
      mvar <- raster::mask(raster::crop(variables[[i]], M), M)
      mval <- na.omit(raster::values(mvar))

      ## distance of each absolute value to median value
      medians <- median(mval)
      df_layer[[j]] <- abs(mval - medians)

      occval <- na.omit(raster::extract(mvar, occ[, c(longitude_col, latitude_col)]))
      occ_dfs[[j]] <- abs(occval - medians)

      yocc[[j]] <- max(table(df_layer[[j]]))

      cat("\t", j, "of", length(occlist), "species finished\n")
    }

    ## frecuency of each value in M
    pdf(file = paste0(output_folder, "/Histograms_", names(variables)[i], ".pdf"),
        width = 6, height = 9)
    par(mfrow = c(5, 2), cex = 0.6, mar = (c(4.5, 4, 3.5, 1) + 0.1))
    plots <- lapply(1:length(df_layer), function(x){
      hist(df_layer[[x]], main = gsub("_", " ", spnames[x]), xlab = "Values")
      points(occ_dfs[[x]], rep(yocc[[x]], length(occ_dfs[[x]])), col = "darkgreen", pch = 1)
      for (j in 1:length(CL_lines)) {
        limit <- floor(CL_lines[j] * length(df_layer[[x]]) / 100)
        abline(v = sort(df_layer[[x]])[limit], col = col[j])
      }
    })
    invisible(dev.off())

    cat(i, "OF", dim(variables)[3], "VARIABLES FINISHED\n")
  }
}
