#' @title readTRENDY
#' @description read a TRENDY dataset
#' @author Michael Crawford
#'
#' @param subtype name of the TRENDY dataset to read. Current options includes:
#'     "CABLEPOP", "CARDAMOM", "CLASSIC", "CLM5.0", "DLEM",
#'     "ELM", "IBIS", "ISBACTRIP", "JSBACH", "LPJ-GUESS",
#'     "LPJml", "LPJwsl", "lpxqs", "OCN", "ORCHIDEE",
#'     "SDGVM", "VISIT"
#'
#' @return a named list of terra::rast objects
#'
#' @examples
#' \dontrun{
#' readSource("TRENDY", subtype = "JSBACH")
#' }
readTRENDY <- function(subtype) {

  trendyYears  <- seq(1700, 2022)

  ncFiles <- list.files(file.path(".", subtype), pattern = "\\.nc$", full.names = TRUE)

  # Load cLitter, cSoil, cVeg rasters for a given subtype
  rasterList <- lapply(ncFiles, function(f) {
    raster <- terra::rast(f)
    names(raster) <- paste0("y", trendyYears)
    return(raster)
  })

  names(rasterList) <- lapply(rasterList, terra::varnames)

  return(list(x = rasterList, class = "list"))
}
