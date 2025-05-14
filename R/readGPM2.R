#' @title readGPM2
#' @description read peatland area from GPM2
#' @param subtype resolution ("1km" or "500m")
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @author Florian Humpenoeder
#' @examples
#' \dontrun{
#' readSource("GPM2", convert = "onlycorrect")
#' }
#' @importFrom magclass as.magpie
#' @importFrom mstools toolGetMappingCoord2Country
readGPM2 <- function(subtype = "1km") {
  previousOptions <- terra::terraOptions(print = FALSE)
  previousOptions$metadata <- NULL # prevent unknown options warning when resetting
  terra::terraOptions(tempdir = withr::local_tempdir(tmpdir = getConfig("tmpfolder")),
                      todisk = FALSE, memfrac = 0.5)
  withr::defer(do.call(terra::terraOptions, previousOptions))

  if (subtype == "1km") {
    # read-in file
    r <- terra::rast("peatGPA22WGS_2cl.tif")

    # choose both layers
    # 1 = peat dominated, 2 = peat in soil mosaic
    r1 <- terra::segregate(r, other = NA)
    r1 <- r1[[c(1, 2)]]

    # project r to lon/lat
    r2 <- terra::project(r1, "+proj=longlat +datum=WGS84", method = "near")

    # get cell area
    a <- terra::cellSize(r2, unit = "ha", mask = TRUE)
    a <- a * 1e-6

    # aggregate to 0.5 degree
    r3 <- terra::project(a, terra::rast(resolution = 0.5), method = "sum")

  } else if (subtype == "500m") {
    # read-in file
    r <- terra::vrt(list.files("500m", full.names = TRUE, pattern = ".tif$"), "500m/GPM2.0_500m.vrt", overwrite = TRUE)
    crs(r) <- "+proj=moll"

    # choose only 1st layer
    r1 <- terra::segregate(r, other = NA)
    r1 <- r1[[1]]

    # project r to lon/lat
    r2 <- terra::project(r1, "+proj=longlat +datum=WGS84", method = "near")

    # get cell area
    a <- terra::cellSize(r2, unit = "ha", mask = TRUE)
    a <- a * 1e-6

    # project or aggregate to 0.5 degree
    # use terra::aggregate because terra::project(a, terra::rast(resolution = 0.5), method = "sum") is not working
    # on the cluster (method = "sum" is the problem)
    r3 <- terra::aggregate(a, fact = 96, fun = "sum", na.rm = TRUE)

  }

  # get spatial mapping
  map <- toolGetMappingCoord2Country(pretty = TRUE)

  # transform raster to magpie object
  x <- as.magpie(terra::extract(r3, map[c("lon", "lat")])[, -1], spatial = 1)

  # set dimension names
  dimnames(x) <- list("coords" = map$coords, "t" = NULL, "d3" = NULL)

  return(x)
}
