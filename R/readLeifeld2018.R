#' @title readLeifeld2018
#' @description read potential peatland area from Leifeld2018
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @author Florian Humpenoeder
#' @examples
#' \dontrun{
#' readSource("Leifeld2018", convert = "onlycorrect")
#' }
#' @importFrom magclass as.magpie
#' @importFrom mstools toolGetMappingCoord2Country
readLeifeld2018 <- function() {
  # projection is +proj=igh
  x <- terra::rast("Degradation_raster_homolosine_hires_rev4.tif")
  # re-project to regular grid
  r <- terra::rast(resolution = 0.5)
  rp2 <- suppressWarnings(terra::project(x, r))
  # get cell area
  cellArea <- terra::cellSize(rp2[[1]], unit = "ha", mask = FALSE) * 1e-6
  cellArea <- terra::mask(cellArea, rp2[[1]])
  # get spatial mapping
  map <- toolGetMappingCoord2Country(pretty = TRUE)
  # transform raster to magpie object
  x <- as.magpie(terra::extract(cellArea, map[c("lon", "lat")])[, -1], spatial = 1)
  # set dimension names
  dimnames(x) <- list("coords" = map$coords, "t" = NULL, "d3" = NULL)

  return(x)
}
