#' @title readGFAD
#' @description Read GLobal Forest Age Dataset derived from MODIS and COPENICUS satellite data
#'
#' @return magpie object in cellular resolution
#' @author Abhijeet Mishra, Felicitas Beier
#' @examples
#' \dontrun{
#' readSource("GFAD", convert = "onlycorrect")
#' }
readGFAD <- function() {

  poulterData   <- "GFAD_V1-1.nc"

  mapping       <- mstools::toolGetMappingCoord2Country()

  ageClass      <- paste0("X", 1:15)
  forestPoulter <- c("NEEV", "NEDE", "BREV", "BRDC")

  temp <- NULL

  for (forestType in seq_along(forestPoulter)) {
    for (i in ageClass) {

      b <- suppressWarnings(raster::brick(poulterData, level = forestType))
      b <- raster::subset(b, i)

      mag <- as.magpie(b)
      mag <- mag[mapping$coords, , ]
      getItems(mag, dim = 2) <- "y2010"
      getItems(mag, dim = 3) <- i

      getSets(mag) <- c("x", "y", "t", "ac")
      mag <- add_dimension(x = mag, dim = 3.1,
                           add = "type",
                           nm = forestPoulter[forestType])

      temp <- mbind(temp, mag)
    }
  }

  out <- setYears(temp, NULL)

  return(out)
}
