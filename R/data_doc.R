#' Example of character table for six species
#'
#' A character table representing species ecological niches derived from
#' previous preparation processes. Each row represents a species and each column
#' a binary character in which one or more values of the environmental variable
#' are categorized as used "1", non used "0", or uncertain "?".
#'
#' @format A character matrix with 6 rows and 28 columns.
#'
#' @examples
#' data("character_table", package = "nichevol")
#'
#' head(character_table)
"character_table"


#' Example of table with results from parsimony reconstructions
#'
#' A character table representing species ecological niches derived from
#' previous preparation processes and reconstructed niches for ancestors.
#' Each row represents a species or a node and each column a binary
#' character in which one or more values of the environmental variable
#' are categorized as used "1", non used "0", or uncertain "?".
#'
#' @format A character matrix with 11 rows and 20 columns.
#'
#' @examples
#' data("par_rec_table", package = "nichevol")
#'
#' head(par_rec_table)
"par_rec_table"



#' Example of a phylogenetic tree for six species
#'
#' A phylogenetic tree with 6 species and their relationships.
#'
#' @format An object of class phylo for 6 species.
#'
#' @examples
#' data("tree", package = "nichevol")
#'
#' str(tree)
"tree"


#' Example of a phylogenetic tree for five species
#'
#' A phylogenetic tree with 5 species and their relationships.
#'
#' @format An object of class phylo for 5 species.
#'
#' @examples
#' data("tree5", package = "nichevol")
#'
#' str(tree5)
"tree5"


#' Example of a list containing a tree and a table of characters for six species
#'
#' A list of 2 elements (phy and data) resulting from using the function
#' \code{\link[geiger]{treedata}}.
#'
#' @format A list of 2 elements:
#' \describe{
#'   \item{phy}{object of class phylo for 6 species}
#'   \item{data}{matrix of 6 rows and 28 columns}
#' }
#'
#' @examples
#' data("tree_data", package = "nichevol")
#'
#' str(tree_data)
"tree_data"


#' Example of occurrence records for six species
#'
#' A list of 6 data.frames containing name and geographic coordinates for
#' 6 species.
#'
#' @format A list of 6 data.frames:
#' \describe{
#'   \item{species}{species name, a code in this example}
#'   \item{x}{longitude, longitude value}
#'   \item{y}{latitude, latitude value}
#' }
#'
#' @examples
#' data("occ_list", package = "nichevol")
#'
#' str(occ_list)
"occ_list"


#' Example of accessible areas for a species
#'
#' A SpatVector object representing the accessible area for a species.
#'
#' @format A SpatVector object.
#'
#' @name m
#' @aliases m1 m2 m3 m4 m5 m6
#'
#' @return No return value, used with function \code{\link[terra]{vect}} to
#' bring an example of an accessible area for a species.
#'
#' @examples
#' m1 <- terra::vect(system.file("extdata", "m1.gpkg", package = "nichevol"))
#'
#' terra::plot(m1)
NULL


#' Example of an environmental variable used in analysis
#'
#' A SpatRaster object representing the variable temperature.
#'
#' @format A SpatRaster object.
#'
#' @name temp
#'
#' @return No return value, used with function \code{\link[terra]{rast}} to
#' bring an example of an environmental variable used in analysis.
#'
#' @examples
#' temp <- terra::rast(system.file("extdata", "temp.tif", package = "nichevol"))
#'
#' terra::plot(temp)
NULL
