#' @title calcGridPop
#'
#' @description Past and future (SSP1-5) population based on HYDE3.2 and Jones & O'Neill (2016)
#' Data is scaled to match WDI data from calcPopulation
#' NOTE that some scaling factors for the future (for small countries Gambia and Djibouti) are off,
#' data read in is 50% of WDI data, most likely due to large resolution
#'
#' @param subtype  time horizon to be returned.
#'                 Options: past (1965-2005),
#'                 future (2005-2010) or
#'                 all (divergence starts at year in harmonize_until)
#' @param source   default source (ISIMIP) or
#'                 Gao data (readGridPopGao) which is split into urban and rural.
#' @param cellular if true: half degree grid cell data returned
#' @param cells    number of halfdegree grid cells to be returned.
#'                 Options: "magpiecell" (59199), "lpjcell" (67420)
#' @param scale  if true: scales sum of gridded values to match country level totals
#' @param FiveYear TRUE for 5 year time steps, otherwise yearly from source
#' @param harmonize_until harmonization year until which SSPs diverge (default: 2015)
#' @param urban    TRUE to return only urban gridded population based on iso share
#'
#' @return Population in millions.
#'
#' @author David Chen, Felicitas Beier
#'
#' @importFrom magclass add_columns collapseNames where
#' @importFrom magpiesets findset
#' @importFrom madrat calcOutput toolGetMapping toolAggregate
#'
#' @examples
#' \dontrun{
#' calcOutput("GridPop", aggregate = FALSE)
#' }
#' @export
calcGridPop <- function(source = "ISIMIP", subtype = "all", # nolint
                        cellular = TRUE, cells = "lpjcell",
                        FiveYear = TRUE, scale = TRUE, # nolint
                        harmonize_until = 2015, urban = FALSE) { # nolint

  if (!cellular) {
    stop("Run calcPopulation instead")
  }

  # the linter assumes `source` refers to `base::source` which is undesirable
  # rename so we need to nolint only once
  src <- source # nolint

  # Gridded population data
  # past data
  if (subtype != "all") {
    if (src == "ISIMIP") {
      x <- readSource("GridPopIsimip", subtype = subtype, convert = FALSE)

      if (subtype == "past") {
        # Add scenario dimension and fill with same value
        x <- add_columns(x, dim = 3.1, addnm = c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5"))
        x[, , 2:6] <- x[, , 1]
        x <- (x[, , -1])
      }
    } else if (src == "Gao") {
      x <- readSource("GridPopGao", subtype = subtype, convert = FALSE)

      if (!urban && subtype == "future") {
        x <- collapseNames(dimSums(x, dim = 3.2))
      }
    }
  } else if (subtype == "all") {
    past <- calcOutput("GridPop", source = "ISIMIP", subtype = "past", cellular = cellular, cells = "lpjcell",
                       FiveYear = FiveYear, scale = FALSE, harmonize_until = 2015, urban = urban, aggregate = FALSE)
    future <- calcOutput("GridPop", source = src, subtype = "future", cellular = cellular, cells = "lpjcell",
                         FiveYear = FiveYear, scale = FALSE, harmonize_until = 2015, urban = urban, aggregate = FALSE)

    if (src == "Gao") {
      intYears <- seq(2005, 2095, 10)
      future <- time_interpolate(future, interpolated_year = intYears, integrate_interpolated_years = TRUE)

      if (urban) {
        # hold past rural urban share constant in each grid for now, based on year 2000

        ratio <- future[, 2000, ] / dimSums(future[, 2000, ], dim = 3.2)
        ratio[is.na(ratio)] <- 0
        ratio[is.infinite(ratio)] <- 0

        past <- setYears(ratio, NULL) * past
      }
    } else if (src == "ISIMIP") {
      if (urban) {
        message("This data source urban/rural is diaggregated using uniform country level value. ",
                "Use Gao source instead for cellular urban/rural population")

        # urban population at country-level
        urbanPop <- calcOutput("Urban", scenario = c("SSPs", "SDPs"), naming = "scenario", aggregate = FALSE)
        # disaggregate to cell level
        coordMapping <- toolGetMappingCoord2Country()
        urbanPop <- toolAggregate(urbanPop,
                                  rel = coordMapping,
                                  from = "iso", to = "coords", partrel = TRUE)
        getCells(urbanPop) <- paste(getItems(urbanPop, dim = 1),
                                    coordMapping$iso[coordMapping$coords == getItems(urbanPop, dim = 1)],
                                    sep = ".")
        getSets(urbanPop) <- c("x", "y", "iso", "year", "data")

        past <- past * urbanPop[, getYears(past), ]
        future <- future * urbanPop[, getYears(future), ]
      }
    }


    # harmonize future SSPs to divergence year by making them SSP2
    selectY <- 1:which(getYears(future, as.integer = TRUE) == harmonize_until)
    harmY   <- getYears(future, as.integer = TRUE)[selectY]
    future[, harmY, ] <- future[, harmY, "SSP2"]

    # take future years in case of overlap
    pYears <- setdiff(getYears(past), getYears(future))

    x <- mbind(past[, pYears, ], future)
    x <- toolHoldConstantBeyondEnd(x)
  }

  # Reduce number of grid cells to 59199
  if (cells == "magpiecell") {
    x <- toolCoord2Isocell(x, cells = cells)
  }

  # Add SDP, SDP_EI, SDP_RC and SDP_MC scenarios as copy of SSP1
  if ("SSP1" %in% getNames(x, dim = 1) && !("SDP" %in% getNames(x, dim = 1))) {
    combinedSDP <- x[, , "SSP1"]
    for  (i in c("SDP", "SDP_EI", "SDP_RC", "SDP_MC")) {
      getNames(combinedSDP) <- gsub("SSP1", i, getNames(x[, , "SSP1"]))
      x <- mbind(x, combinedSDP)
    }
  }

  if (scale) {
    ## Scale to match country-level data
    # Country-level population data (in million)
    pop <- calcOutput("Population", scenario = c("SSPs", "SDPs"), naming = "scenario", aggregate = FALSE)
    # aggregate to country-level and scale to match WDI country-level pop
    if (cells == "lpjcell") {
      agg <- collapseNames(dimSums(x, dim = c("x", "y")))
    } else if (cells == "magpiecell") {
      agg <- collapseNames(dimSums(x, dim = "cell"))
    }
    # fill missing years
    pop <- time_interpolate(pop, interpolated_year = getYears(agg))
    # scale from millions to units in agg
    pop <- pop * 1e6

    commonCountries <- intersect(getItems(agg, dim = 1), getItems(pop, dim = 1))
    # Note: 15 countries are dropped because they are not in the spatial
    #       data ("ABW" "AND" "ATA" "BES" "BLM" "BVT" "GIB" "LIE"
    #             "MAC" "MAF" "MCO" "SMR" "SXM" "VAT" "VGB")
    # scaling factor scF applied to every cell
    if (urban) {
      scF <- (dimSums(agg[commonCountries, , ], dim = 3.2)) / pop[commonCountries, , ]
    } else {
      scF <- agg[commonCountries, , ] / pop[commonCountries, , ]
    }

    x <- x[commonCountries, , ] / scF[commonCountries, , ]

    # Some countries have 0 population in grid, but the cells exist, so get filled.
    # even division of population across cells
    missing <- magclass::where(scF[commonCountries, , ] == 0)$true$regions
    for (i in missing) {
      x[i, , ] <- pop[i, , "SSP2"] / length(getCells(x[i, , ]))
    }
  }
  # unit conversion to million people
  x <- x / 1e6

  if (FiveYear) {
    years <- findset("time")
    x <- x[, intersect(years, getYears(x)), ]
  }

  # Checks
  if (any(is.na(x))) {
    stop("Function calcGridPop returned NAs.")
  }
  if (any(round(x, digits = 4) < 0)) {
    stop("Function calcGridPop returned negative values.")
  }

  return(list(x = x,
              weight = NULL,
              unit = "million",
              description = "Population in millions",
              isocountries = FALSE))
}
