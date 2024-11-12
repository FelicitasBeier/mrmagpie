#' @title calcIrrigation
#' @description This function extracts irrigation water (airrig: water applied additionally to rainfall)
#'              from LPJmL for MAgPIE
#'
#' @param lpjml         Defines LPJmL version for crop/grass and natveg specific inputs
#' @param climatetype   Switch between different climate scenarios
#' @param rainfedweight For clustering airrig is weighted with
#'                      cropland_irrigated + rainfedweight * cropland_rainfed (default: 0.01)
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier, Abhijeet Mishra
#'
#' @examples
#' \dontrun{
#' calcOutput("Irrigation", aggregate = FALSE)
#' }
#'
#' @importFrom magpiesets findset
#' @importFrom mrlandcore toolLPJmLVersion
#' @importFrom mstools toolGetMappingCoord2Country toolCoord2Isocell
#' @importFrom madrat toolGetMapping calcOutput toolAggregate
#' @importFrom magclass dimSums getItems getSets collapseNames
#' @importFrom withr local_options

calcIrrigation <- function(lpjml = "lpjml5.9.5-m1",
                           climatetype = "MRI-ESM2-0:ssp370",
                           rainfedweight = 0.01) {
  # Set size limit
  local_options(magclass_sizeLimit = 1e+12)

  # Read in airrig (irrigation water applied additionally to rainfall where irrigation takes place):
  airrigLPJ <- collapseNames(calcOutput("LPJmLharmonize", subtype = "cft_airrig",
                                        lpjmlversion = lpjml, climatetype = climatetype,
                                        aggregate = FALSE))

  # Load LPJmL to MAgPIE mapping to aggregate to MAgPIE crops
  mapping <- toolGetMapping("MAgPIE_LPJmL.csv", type = "sectoral", where = "mrlandcore")
  # Aggregate to MAgPIE crops
  airrigMAG <- toolAggregate(x = airrigLPJ, rel = mapping,
                             dim = 3.1, partrel = TRUE,
                             from = "LPJmL5", to = "MAgPIE")
  # Remove pasture (pasture is not irrigated in MAgPIE)
  airrigMAG <- airrigMAG[, , "pasture", invert = TRUE]
  # Remove negative airrig
  airrigMAG[airrigMAG < 0] <- 0

  # Check for NAs
  if (any(is.na(airrigMAG))) {
    stop("produced NA airrig")
  }

  # Clustering weight:
  totalCropland <- dimSums(calcOutput("Croparea",
                                      sectoral = "kcr", physical = TRUE,
                                      cellular = TRUE, cells = "lpjcell",
                                      years = "y1995", round = 6,
                                      irrigation = TRUE, aggregate = FALSE),
                           dim = 3.2)
  map <- toolGetMappingCoord2Country()
  getItems(totalCropland, dim = 1, raw = TRUE) <- paste(map$coords, map$iso, sep = ".")
  getSets(totalCropland) <- c("x", "y", "iso", "year", "irrigation")

  weightCropArea <- collapseNames(totalCropland[, , "irrigated"]) +
    rainfedweight * collapseNames(totalCropland[, , "rainfed"]) + 10e-10

  return(list(x = airrigMAG,
              weight = weightCropArea,
              unit = "m^3 per ha per yr",
              description = paste0("Irrigation water (water applied in addition to rainfall) for ",
                                   "different crop types following LPJmL irrigation system assumptions"),
              isocountries = FALSE))
}
