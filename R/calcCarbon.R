#' @title calcCarbon
#' @description Extract the carbon densities used in MAgPIE from a set of LPJmL runs
#'
#' @param lpjml       Defines LPJmL version for crop/grass and natveg specific inputs
#' @param climatetype Switch between different GCM climate scenarios
#'
#' @return magpie object in cellular resolution
#' @author Kristine Karstens, Patrick v. Jeetze, Michael Crawford
#'
#' @examples
#' \dontrun{
#' calcOutput("Carbon", aggregate = FALSE)
#' }
#'
#' @importFrom magpiesets findset

calcCarbon <- function(lpjml       = "lpjml5.9.5-m1",
                       climatetype = "MRI-ESM2-0:ssp370") {

  # HACKATHON – remove `years` once the grass run is extended back to 1930
  years <- paste0("y", seq(1950, 2100))

  # -----------------------------------------------------------------------------------------------
  # Generate carbon datasets

  # Collect `pnv` and `grass` carbon pools from LPJmL
  .getLPJmLCPools <- function(run, pools) {

    poolData <- lapply(pools, function(pool) {
      out <- calcOutput("LPJmLharmonize",
                        lpjmlversion = lpjml,
                        climatetype  = climatetype,
                        subtype      = paste(run, pool, sep = ":"),
                        aggregate    = FALSE)[, years, ]
      setNames(out, pool)
    })

    Reduce(mbind, poolData)
  }

  natveg <- .getLPJmLCPools(run = "pnv",   pools = c("vegc", "soilc", "litc"))
  grass  <- .getLPJmLCPools(run = "grass", pools = c("vegc", "soilc", "litc"))

  # Collect topsoil carbon
  topsoilc <- calcOutput("TopsoilCarbon",
                         lpjml       = lpjml,
                         climatetype = climatetype,
                         aggregate   = FALSE)[, years, ]

  # Calculate soil carbon loss share
  cshare <- calcOutput("SOCLossShare", aggregate = FALSE)

  # -----------------------------------------------------------------------------------------------
  # Format output object

  carbonStocks <- new.magpie(
    cells_and_regions = getCells(natveg),
    years = getYears(natveg),
    names = getNames(natveg),
    sets  = c("x.y.iso", "year", "data")
  )

  carbonStocks <- add_dimension(carbonStocks,
                                dim = 3.1,
                                add = "landtype",
                                nm = c("crop", "past", "forestry", "primforest", "secdforest", "urban", "other"))

  # -----------------------------------------------------------------------------------------------
  # Calculate the appropriate values for all land types and carbon types

  landIni <- calcOutput("LanduseInitialisation",
                        aggregate    = FALSE,
                        cellular     = TRUE,
                        nclasses     = "seven",
                        input_magpie = TRUE,
                        years        = "y1995",
                        round        = 6)

  # Factor 0.012 is based on the script subversion/svn/tools/carbon_cropland, executed at 30.07.2013
  carbonStocks[, , "crop.vegc"]  <- 0.012 * natveg[, , "vegc"]
  carbonStocks[, , "crop.litc"]  <- 0 # crops do not have litter
  carbonStocks[, , "crop.soilc"] <- cshare * topsoilc + (natveg[, , "soilc"] - topsoilc)

  carbonStocks[, , "past"]  <- grass

  # HACKATHON – The grass run is still FALSE
  # If soilc stocks are too large in ALLCROP grass run, use natveg stock for soilc in grasslands
  grasssoil   <- TRUE
  grassSoilc  <- dimSums(grass[, , "soilc"]  * landIni[, , "past"], dim = c(1, 2))
  natvegSoilc <- dimSums(natveg[, , "soilc"] * landIni[, , "past"], dim = c(1, 2))
  if (grassSoilc > natvegSoilc) {
    carbonStocks[, , "past.soilc"] <- natveg[, , "soilc"]
    grasssoil <- FALSE
  }

  carbonStocks[, , "forestry"]    <- natveg
  carbonStocks[, , "primforest"]  <- natveg
  carbonStocks[, , "secdforest"]  <- natveg
  carbonStocks[, , "urban"]       <- 0
  carbonStocks[, , "urban.soilc"] <- natveg[, , "soilc"]
  carbonStocks[, , "other"]       <- natveg

  mstools::toolExpectTrue(
    !any(is.na(carbonStocks)),
    "no NAs in carbonStocks",
    falseStatus = "warn"
  )

  # -----------------------------------------------------------------------------------------------
  # Calculate aggregation weight based on the potential forest area

  potForestArea <- calcOutput("PotentialForestArea",
                              refData     = "lpj",
                              lpjml       = lpjml,
                              climatetype = climatetype,
                              aggregate   = FALSE)[, years, ]

  weight <- new.magpie(
    cells_and_regions = getCells(carbonStocks),
    years = getYears(carbonStocks),
    names = getNames(carbonStocks),
    sets  = c("x.y.iso", "year", "data")
  )

  landArea <- dimSums(landIni, dim = 3)
  weight[, , ] <- landArea
  weight[, , c("primforest", "secdforest", "forestry")] <- potForestArea

  return(list(
    x           = carbonStocks,
    weight      = weight + 10^-10,
    unit        = "t per ha",
    description = "Carbon in tons per hectar for different land use types",
    note = ifelse(grasssoil,
                  "Pasture soil carbon stocks are based on allcrop run",
                  "Pasture soil carbon stocks are based on natveg run"),
    isocountries = FALSE
  ))
}
