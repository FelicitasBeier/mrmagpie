#' @title calcAreaActuallyIrrigated
#' @description retrieves irrigated crop area from croparea intialization
#'
#' @param aggregationlevel default is iso
#' @param selectyears select years
#'
#' @return magpie object with results on cellular or iso country level
#' @author Felicitas Beier
#' @examples
#' \dontrun{
#' calcOutput("AreaActuallyIrrigated")
#' }
#'
#' @importFrom magclass dimSums

calcAreaActuallyIrrigated <- function(aggregationlevel = "iso", selectyears = "y1995") {
  # Read in data: crop- and water supply type specific crop area (in Mha):
  x <- calcOutput("Croparea", physical = TRUE, cellular = TRUE, irrigation = TRUE,
                  round = 6, aggregate = FALSE, years = selectyears, fallow = TRUE)

  # extract irrigated area:
  x <- dimSums(x[, , "irrigated"], dim = 3.1)
  # total irrigated area per cell (sum over different crop types)
  x <- dimSums(x, dim = 3)

  if (aggregationlevel == "iso") {
    # country to cell mapping
    countryToCell <- toolGetMapping(type = "cell", name = "CountryToCellMapping.csv", where = "mappingfolder")
    # aggregate data
    toolAggregate(x, dim = 1, rel = countryToCell)

  } else {
    out <- x
  }

  return(list(x = out,
              weight = NULL,
              unit = "Million ha",
              description = "Irrigated cropland area"))
}
