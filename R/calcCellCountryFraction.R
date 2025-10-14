#' @title calcCellCountryFraction
#' @description cell fraction belonging to a country based on LanduseInitialisation
#' @param cells lpjcell for 67420 cells or magpiecell for 59199 cells
#'
#' @return Clustered MAgPIE object on requested resolution
#' @author Florian Humpenoeder
#'
#' @examples
#' \dontrun{
#' calcOutput("calcCellCountryFraction", aggregate = FALSE)
#' }
calcCellCountryFraction <- function(cells = "lpjcell") {

  weight <- calcOutput("LandArea", cells = cells, aggregate = FALSE)

  if (cells == "lpjcell") {

    x <- new.magpie(cells_and_regions = getItems(weight, dim = 1),
                    years = NULL,
                    names = getItems(weight, dim = 1.3),
                    fill = 0)

    for (r in getItems(weight, dim = 1.3)) {
      x[r, , r] <- 1
    }

    getSets(x) <- c("x", "y", "iso", "year", "data")

  } else if (cells == "magpiecell") {

    x <- new.magpie(cells_and_regions = getItems(weight, dim = 1),
                    years = NULL,
                    names = getItems(weight, dim = 1.1),
                    fill = 0)

    for (r in getItems(x, dim = 1.1)) {
      x[r, , r] <- 1
    }
  } else {
    stop("Please select cells argument: lpjcell for 67420 or magpiecell for 59199")
  }

  return(list(x = x,
              weight = weight,
              unit = "share",
              description = "cell fraction belonging to a country",
              isocountries = FALSE))
}
