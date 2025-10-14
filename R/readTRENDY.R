#' @title readTRENDY
#' @description read a TRENDY dataset
#' @author Michael Crawford
#'
#' @param subtype name of the TRENDY dataset to read. Current options include:
#'     "CABLEPOP", "CLASSIC", "CLM5.0", "DLEM", "ED", "ELM",
#'     "IBIS", "ISAM", "ISBACTRIP", "JSBACH", "JULES",
#'     "LPJ-GUESS", "LPJml", "LPJwsl", "lpxqs", "OCN",
#'     "ORCHIDEE", "SDGVM", "VISIT", "YIBS"
#'
#' @return a named list of terra::rast objects
#'
#' @examples
#' \dontrun{
#' readSource("TRENDY", subtype = "JSBACH")
#' }
readTRENDY <- function(subtype) {

  trendyYear <- 2022

  ncFiles <- list.files(file.path(".", subtype), pattern = "\\.nc$", full.names = TRUE)

  # Load cLitter, cSoil, cVeg rasters for a given subtype
  rasterList <- lapply(ncFiles, function(f) {
    raster <- terra::rast(f)
    nLayers <- terra::nlyr(raster)
    startYear <- trendyYear - nLayers + 1
    trendyYears <- seq(startYear, trendyYear)
    names(raster) <- paste0("y", trendyYears)
    return(raster)
  })

  names(rasterList) <- lapply(rasterList, terra::varnames)

  return(list(x = rasterList, class = "list"))
}
