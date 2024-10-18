#' @title calcIrrigation
#' @description This function extracts irrigation water (airrig: water applied additionally to rainfall)
#'              from LPJmL for MAgPIE
#'
#' @param lpjml         Defines LPJmL version for crop/grass and natveg specific inputs
#' @param climatetype   Switch between different climate scenarios
#' @param cells         Number of cells to be returned:
#'                      "magpiecell" for 59199 cells or "lpjcell" for 67420 cells
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

calcIrrigation <- function(lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de", crop = "ggcmi_phase3_nchecks_9ca735cb"),
                           climatetype = "GSWP3-W5E5:historical", cells = "lpjcell", rainfedweight = 0.01) {
  # Extract arguments
  if (grepl("GSWP3-W5E5", climatetype)) {
    stage       <- "smoothed"
  } else {
    stage       <- "harmonized2020"
  }

  # Set size limit
  local_options(magclass_sizeLimit = 1e+12)

  # Read in airrig (irrigation water applied additionally to rainfall where irrigation takes place):
  airrigLPJ   <- collapseNames(calcOutput("LPJmL_new", version = lpjml["crop"],
                                          climatetype = climatetype, subtype = "irrig",
                                          stage = stage,
                                          aggregate = FALSE))

  # Load LPJmL to MAgPIE mapping to aggregate to MAgPIE crops
  mapping   <- toolGetMapping("MAgPIE_LPJmL.csv",
                              type = "sectoral", where = "mrlandcore")
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
  totalCropland <- dimSums(calcOutput("Croparea", sectoral = "kcr", physical = TRUE,
                                      cellular = TRUE, cells = "lpjcell",
                                      years = "y1995", round = 6,
                                      irrigation = TRUE, aggregate = FALSE),
                           dim = 3.2)
  map <- toolGetMappingCoord2Country()
  getItems(totalCropland, dim = 1, raw = TRUE) <- paste(map$coords, map$iso, sep = ".")
  getSets(totalCropland) <- c("x", "y", "iso", "year", "irrigation")

  weightCropArea <- collapseNames(totalCropland[, , "irrigated"]) +
    rainfedweight * collapseNames(totalCropland[, , "rainfed"]) + 10e-10

  # Reduce to 59199 cells
  if (cells == "magpiecell") {
    airrigMAG      <- toolCoord2Isocell(airrigMAG)
    weightCropArea <- toolCoord2Isocell(weightCropArea)
  }

  return(list(x = airrigMAG,
              weight = weightCropArea,
              unit = "m^3 per ha per yr",
              description = "Irrigation water (water applied in addition to rainfall) for
                             different crop types following LPJmL irrigation system assumptions",
              isocountries = FALSE))
}
