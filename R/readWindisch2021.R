#' @title readWindisch2021
#' @description Reads in data to calculate BphEffect, BphTCRE or BphMask depending on the chosen subtype.
#'              BphEffect: Biogeophysical temperature change of afforestation (degree C).
#'              (File is based on observation datasets of Bright et al. 2017 and Duveiller et al. 2018).
#'              BphMask: Mask of Datapoints of biogeophysical temperature change of afforestation (degree C)
#'              to be used as weight.
#'              (File is based on observation datasets of Bright et al. 2017 and Duveiller et al. 2018).
#'              BphTCRE: Transient Climate Response to accumulated doubling of CO2.
#'              (File is based on CMIP5 +1perc CO2 per year experiment.
#'              To be used in the translation to carbon equivalents of BphEffect)
#' @param subtype refordefor_BPHonly_05_new, annmean_pertCha_05_EW1, annstd_diff_pertCha_05_EW1
#'
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @author Felicitas Beier, Michael Windisch, Patrick v. Jeetze
#' @importFrom terra rast ext
#' @examples
#'
#' \dontrun{
#'   readSource("Windisch2021", convert="onlycorrect")
#' }
#'

readWindisch2021 <- function(subtype) {

  if (subtype == "refordefor_dT_ANN_nonzeromask_05") {
    # read in bph mask
    x <- read.magpie(paste0(subtype, ".nc"))
    getNames(x) <- "mask"
  } else {
    # read in raster objects
    x <- terra::rast(paste0(subtype, ".nc"))
    # replacing NA's by zero
    x[is.na(x[])] <- 0

    if (!grepl("refordefor", subtype)) {
      storage.mode(x[]) <- "double"
      x[] <- as.double(x[])
      storage.mode(x[])
    }

    # set the to longlat projection
    terra::crs(x) <- terra::crs(terra::rast(resolution = 0.5))
    terra::ext(x) <- terra::ext(-180, 180, -90, 90)

    # transform to magpie object
    x <- as.magpie(x)
  }

  # reduce to 67420 cells
  mapping <- toolGetMappingCoord2Country(pretty = TRUE)
  x       <- x[mapping$coords, , ]
  # clean naming
  getItems(x, dim = 1, raw = TRUE) <- paste(mapping$coords, mapping$iso, sep = ".")
  getSets(x)                       <- c("x", "y", "iso", "year", "data")

  # return magpie object
  return(x)
}
