#' @title readTravelTimeNelson2019
#' @description Read minimum travel time to cities and ports and
#' ports of various size, see metadata file in source folder
#' @param subtype currently only cities of 5, 20, or 50 thousand people ("cities5", "cities20", "cities50") or
#' ports of various sizes ("portsLarge|Medium|Small|VerySmall|Any")
#' @return gridded magpie object for 2015, minimum travel time to cities in minutes
#' @author David M Chen
#' @importFrom terra aggregate project rast classify focal
#' @importFrom raster brick extract

readTravelTimeNelson2019 <- function(subtype = "cities50") {
  terra::terraOptions(tempdir = withr::local_tempdir(tmpdir = getConfig("tmpfolder")),
                      todisk = TRUE, memfrac = 0.5)
  withr::defer(terra::terraOptions(tempdir = tempdir()))

  layers  <- c(cities5             = "travel_time_to_cities_12.tif",
               cities20            = "travel_time_to_cities_10.tif",
               cities50            = "travel_time_to_cities_11.tif",
               portsLarge          = "travel_time_to_ports_1.tif",
               portsMedium         = "travel_time_to_ports_2.tif",
               portsSmall          = "travel_time_to_ports_3.tif",
               portsVerySmall      = "travel_time_to_ports_4.tif",
               portsAny            = "travel_time_to_ports_5.tif")

  file <- toolSubtypeSelect(subtype, layers)
  file <- list.files(path = "./TravelTimeNelson_unzipped", pattern = file, full.names = TRUE)

  r <- rast(resolution = 0.5)

  x <- rast(file)
  x <- terra::classify(x, cbind(65534, 65536, NA), right = FALSE)
  x <- terra::aggregate(x, fact = 60, fun = "mean", na.rm = TRUE)
  x <- terra::project(x, r)

  # fill NAs with focal function (neighbour mean)
  x <- terra::focal(x, w = 9, mean, na.policy = "only", na.rm = TRUE)

  x <- raster::brick(x)

  # get spatial mapping
  map <- toolGetMappingCoord2Country(pretty = TRUE)
  # transform raster to magpie object
  out <- as.magpie(raster::extract(x, map[c("lon", "lat")]), spatial = 1)
  # set dimension names
  dimnames(out) <- list("x.y.iso" = paste(map$coords, map$iso, sep = "."),
                        "year" = NULL,
                        "data" = subtype)

  getItems(out, dim = 2) <- "y2015"

  # fill remote island NAs with high transport time
  out[is.na(out)] <- quantile(out, na.rm = TRUE, 0.90)

  return(out)
}
