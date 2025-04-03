#' @title calcPeatland2
#' @description This function calculates degraded and intact peatland area at cell level.
#' The function takes degraded and intact peatland area from the Global Peatland Database 2022 (GPD2022) at the
#' national level and downscales the peatland area to grid cell level using gridded peatland area
#' from the Global Peatland Map 2.0 (GPM2)
#' The data has been provided by Alexandra Barthelmes.
#' @param countryLevel    Whether output shall be at country level.
#'                         Requires aggregate=FALSE in calcOutput.
#' @return magpie object in cellular resolution
#' @author Florian Humpenoeder
#'
#' @examples
#' \dontrun{
#' calcOutput("Peatland2", aggregate = FALSE)
#' }
#'
#' @importFrom madrat toolAggregate
#' @importFrom mstools toolGetMappingCoord2Country

calcPeatland2 <- function(countryLevel = FALSE) {

  # Country-level data on intact and degraded peatland from Global Peatland Database for 2022 (GPD2022)
  gpd2022 <- readSource("GPD2022", convert = TRUE)

  # Global Peatland Map 2.0; peatland location and extent; but no information on status peatlands (intact / degraded)
  gpm2 <- readSource("GPM2", convert = "onlycorrect")

  # Dissag. GPD2022 from country to cell with GPM2 as weight
  map <- toolGetMappingCoord2Country(pretty = TRUE)
  outCell   <- toolAggregate(x = toolIso2CellCountries(gpd2022, cells = "lpjcell"), rel = map,
                             weight = gpm2, dim = 1, from = "iso", to = "coords", zeroWeight = "allow")
  names(dimnames(outCell)) <- c("coords", "t", "d3")

  # calculate grid cell area of the earths sphere
  calArea <- function(ix, iy, res = 0.5, mha = 1) { # pixelarea in m2, mha as factor
    mha * (111.263 * 1000 * res) * (111.263 * 1000 * res) * cos(iy * pi / 180.)
  }
  grarea <- new.magpie(cells_and_regions = map$coords,
                       fill = calArea(map$lon, map$lat, mha = 10^-10))
  names(dimnames(grarea)) <- c("coords", "t", "d3")

  # reduce peatland area proportionally if larger than grid cell area
  outCellShr <- outCell / dimSums(outCell, dim = 3)
  outCell[dimSums(outCell, dim = 3) > grarea] <- (outCellShr * grarea)[dimSums(outCell, dim = 3) > grarea]


  dimnames(outCell) <- list("x.y.iso" = paste(map$coords, map$iso, sep = "."), "t" = NULL, "d3" = getNames(outCell))

  if (countryLevel) {
    outCell <- toolCountryFill(dimSums(outCell, dim = c("x", "y")), fill = 0)
  }

  description <- "Intact and degraded peatland area (Mha) by land-use type, based GPD 2022 and GPM2.0"

  return(list(
    x = outCell,
    weight = NULL,
    unit = "Mha",
    description = description,
    isocountries = FALSE
  ))
}
