#' @title readGridPopIsimip
#' @description Reads in past and future (SSP1-5) gridded population data, from
#' ISIMIP database, Past data is based on HYDE3.2, while future SSPs are based
#' on projections from Jones & O'Neill 2016
#'
#' @return A MAgPIE object, cellular 0.5deg resolution, of population (millions)
#'
#' @param subtype past (1965-2005) or future (2010-2100)
#'
#' @author David Chen, Marcos Alves, Felicitas Beier
#'
#' @importFrom raster extent brick
#' @importFrom magpiesets findset
#' @importFrom mstools toolGetMappingCoord2Country
readGridPopIsimip <- function(subtype) {

  coordMapping <- toolGetMappingCoord2Country()

  if (subtype == "past") {
    # extract information
    iniPast   <- as.numeric(gsub("y", "", findset("past")))[1]
    pastFile  <- "population_histsoc_0p5min_annual_1861-2005.nc"
    years     <- as.numeric(unlist(regmatches(pastFile,
                                              gregexpr("\\d{4}", pastFile))))
    yYears    <- paste0("y", seq(years[1], years[2], 1))

    # read in raster
    x         <- suppressWarnings(brick(pastFile))
    names(x)  <- yYears
    x         <- raster::subset(x, paste0("y", seq(iniPast, years[2], 1)))
    raster::extent(x) <- c(-180, 180, -90, 90)

    # transform to magpie object
    x <- as.magpie(x)
    x <- x[coordMapping$coords, , ]
  }

  if (subtype == "future") {

    files <- c("population_ssp1soc_0p5deg_annual_2006-2100.nc4",
               "population_ssp2soc_0p5deg_annual_2006-2100.nc4",
               "population_ssp3soc_0p5deg_annual_2006-2100.nc4",
               "population_ssp4soc_0p5deg_annual_2006-2100.nc4",
               "population_ssp5soc_0p5deg_annual_2006-2100.nc4")

    read <- function(file) {
      # extract information
      years     <- as.numeric(unlist(regmatches(file,
                                                gregexpr("\\d{4}", file))))
      yYears    <- paste0("y", seq(years[1], years[2], 1))

      # read in raster
      x         <- suppressWarnings(brick(file))
      names(x)  <- yYears

      # transform to magpie object
      x <- as.magpie(x)
      x <- x[coordMapping$coords, , ]

      return(x)
    }

    x         <- list()
    sspNaming <- NULL
    for (file in files) {
      sspNaming <- append(sspNaming,
                          paste0(toupper(unlist(regmatches(file,
                                                           gregexpr("ssp[0-9]", file))))))
      x[[file]] <- read(file)
    }
    x           <- mbind(x)
    getNames(x) <- sspNaming
  }

  getCells(x) <- paste(getItems(x, dim = 1),
                       coordMapping$iso[coordMapping$coords == getItems(x, dim = 1)],
                       sep = ".")
  getSets(x)  <- c("x", "y", "iso", "year", "data")

  return(x)
}
