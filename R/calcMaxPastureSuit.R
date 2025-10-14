#' @title calcMaxPastureSuit
#' @description Calculate maximum glassland suitable for pasture management based
#' on population and aridity criteria.
#'
#' @param lpjml Defines LPJmL version for crop/grass and natveg specific inputs
#' @param climatetype Switch between different climate scenarios
#' @param cells number of halfdegree grid cells to be returned.
#'              Options: "magpiecell" (59199), "lpjcell" (67420)
#'
#' @return List of magpie object with results on cluster level
#' @author Marcos Alves, Kristine Karstens, Alexandre Köberle
#' @examples
#' \dontrun{
#' calcOutput("MaxPastureSuit")
#' }
#' @importFrom raster area rasterFromXYZ


calcMaxPastureSuit <- function(climatetype = "MRI-ESM2-0:ssp126",
                               lpjml = "LPJmL4_for_MAgPIE_44ac93de",
                               cells = "lpjcell") {
  x <- toolSplitSubtype(climatetype, list(climatemodel = NULL, scenario = NULL))

  # Extract stage argument information
  if (grepl("GSWP3-W5E5", climatetype)) {
    stage <- "smoothed"
  } else {
    stage <- "harmonized2020"
  }

  # Drivers of managed pastures
  population <- calcOutput("GridPop",
    subtype = "all", cellular = TRUE, FiveYear = TRUE, harmonize_until = 2015,
    aggregate = FALSE
  )[, , toupper(substring(x$scenario, first = 0, last = 4))]

  cellPrep <- calcOutput("LPJmLClimateInput_new",
    climatetype = climatetype,
    variable = "precipitation:monthlySum",
    stage = stage,
    lpjmlVersion = lpjml,
    aggregate = FALSE
  )


  cellPet <- calcOutput(
    type = "LPJmL_new", climatetype = climatetype,
    subtype = "mpet",
    stage = stage,
    version = lpjml,
    aggregate = FALSE
  )

  yearsCellPet <- intersect(getYears(cellPet), findset("time"))
  yearsCellPrep <- intersect(findset("time"), getYears(cellPrep))
  years <- intersect(yearsCellPet, yearsCellPrep)
  cellPrep <- dimSums(cellPrep[, years, ], dim = 3)
  cellPet <- dimSums(cellPet[, years, ], dim = 3)

  # Cell area calculation
  land <- calcOutput("LanduseInitialisation", input_magpie = TRUE,
                     aggregate = FALSE, cellular = TRUE, cells = cells,
                     years = "y1995", round = 6)
  landArea <- setYears(dimSums(land, dim = 3), NULL)

  # Calculate population density in people per km2: (Million people / Mha) = (1e6 / 10000)
  popDensity <- population / landArea * 1e6 / 10000
  popDensity[is.infinite(popDensity)] <- 0
  popDensity[is.nan(popDensity)] <- 0

  years <- intersect(getYears(popDensity), getYears(cellPrep))

  # Following the methodology presented in HYDE3.2
  # (https://essd.copernicus.org/articles/9/927/2017/essd-9-927-2017.pdf),
  # to distinguish between these types of grasslands in a simple and transparent
  # way, also historically, a population density and an aridity index are applied,
  # based on expert judgement. Low animal densities can be related to either low
  # population density or to low productivity of the natural vegetation, which is
  # approximated via the aridity index. When the aridity index (defined as annual
  # precipitation divided by annual evapotranspiration) of a grid cell is greater than 0.5,
  # and population density is greater than 5 inhabitants km^2, then it is defined
  # as suitable for managed pastures, otherwise rangelands.

  aridity <- cellPrep[, years, ] / cellPet[, years, ]
  aridity[is.infinite(aridity) | is.nan(aridity)] <- 0

  aridity[aridity < 0.5] <- 0
  aridity[aridity >= 0.5] <- 1

  pastureSuit <- aridity
  popDensity <- popDensity[, getYears(pastureSuit), ]

  popDensity <- popDensity > 5
  pastureSuit <- pastureSuit * popDensity

  pastureSuitArea <- (pastureSuit * landArea * 100) / 1e6 # (from km2 (x100) to mha (/1e6))

  pastureSuitArea <- mstools::toolHoldConstantBeyondEnd(pastureSuitArea)
  getItems(pastureSuitArea, dim = 3)  <- NULL

  return(list(
    x = pastureSuitArea,
    weight = NULL,
    unit = "Mha",
    description = "Area suitable for pasture management in mha",
    isocountries = FALSE
  ))
}
