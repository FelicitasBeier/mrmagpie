#' @title readGridPopGao
#' @description Read gridded population, by urban and rural, from Gao O'Neill and JOnes dataset,
#' see https://www.cgd.ucar.edu/iam/modeling/spatial-population-scenarios.html https://doi.org/10.7927/m30p-j498
#' @param subtype only "future" post-2000 available for this source
#' @author David Chen, Felicitas Beier
readGridPopGao <- function(subtype = "future") {

  if (subtype == "future") {
    ssps   <- c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5")
    years  <- seq(2010, 2100, 10)
    urbans <- c("rural", "urban")

    # mapping
    coordMapping <- toolGetMappingCoord2Country()

    .read <- function(year, urbans, ssps, mapping) {
      x <- NULL
      for (urb in urbans) {
        for (ssp in ssps) {

          if (urb == "rural") {
            ur <- "rur"
          } else if (urb == "urban") {
            ur <- "urb"
          }

          if (year == 2000) {
            t <- suppressWarnings(raster::raster(paste0("Base_year_data_2000_NetCDF/2000", urb, ".nc")))
          } else {
            t <- suppressWarnings(raster::raster(paste0(ssp, "_NetCDF/", urb, "/NetCDF/", tolower(ssp),
                                                        ur, year, ".nc")))
          }

          # aggregate and reproject
          r <- suppressWarnings(raster::raster(resolution =  0.5))
          t <- raster::aggregate(t, fact = 4, fun = sum, na.rm = TRUE)
          t <- suppressWarnings(raster::projectRaster(t, r, over = TRUE))

          # transform to magpie object
          tmp <- as.magpie(t)
          t   <- new.magpie(cells_and_regions = mapping$coords,
                            years = getYears(t),
                            names = getNames(t),
                            fill = 0)
          sameCells <- intersect(getItems(tmp, dim = 1), mapping$coords)
          t[sameCells, , ] <- tmp[sameCells, , ]

          # proper names
          getYears(t)          <- year
          getNames(t, dim = 1) <- ssp
          getNames(t, dim = 2) <- urb

          x <- mbind(x, t)
        }
      }
      return(x)
    }

    # read base year
    base <- .read(year = 2000, urbans = urbans, ssps = "SSP1", mapping = coordMapping)
    base <- add_columns(base, addnm = c("SSP2", "SSP3", "SSP4", "SSP5"), dim = 3.1)
    base[, , c("SSP2", "SSP3", "SSP4", "SSP5")] <- base[, , "SSP1"]

    # do the big loop with years as list (can't mbind both list and 3rd dim internally)
    out <- lapply(years, FUN = .read, urbans = urbans, ssps = ssps, mapping = coordMapping)
    out <- mbind(out)

    # bind base
    x <- mbind(base, out)

    # cell and set naming
    getCells(x) <- paste(getItems(x, dim = 1),
                         coordMapping$iso[coordMapping$coords == getItems(x, dim = 1)],
                         sep = ".")
    getSets(x)  <- c("x", "y", "iso", "year", "data", "newdim")

    return(x)
  } else if (subtype == "past") {
    stop("Data only available from 2000 onwards, use ISIMIP source instead")
  }
}
