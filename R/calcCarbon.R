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

calcCarbon <- function(lpjml       = "lpjml5.9.5-m1",
                       climatetype = "MRI-ESM2-0:ssp370") {

  # HACKATHON: remove `years` once the grass run is extended back to 1930
  years <- paste0("y", seq(1950, 2100))

  # -----------------------------------------------------------------------------------------------
  # Calculate carbon datasets

  # Calculate carbon pools from LPJmL
  .collectCarbonPools <- function(run, pools) {

    .processPool <- function(pool) {
      .y <- calcOutput("LPJmLHarmonize",
                       lpjmlversion = lpjml,
                       climatetype  = climatetype,
                       subtype      = paste0(run, ":", pool),
                       aggregate    = FALSE)[, years, ]

      setNames(.y, pool)
    }

    .x <- lapply(pools, .processPool)
    .x <- Reduce(mbind, .x)

    return(.x)
  }

  pnv   <- .collectCarbonPools(run = "pnv",   pools = c("vegc", "soilc", "litc"))

  # HACKATHON: For now we use vegc_avg for comparison purposes, but will switch to vegc before final PRs
  grass <- .collectCarbonPools(run = "grass", pools = c("vegc_avg", "soilc", "litc"))
  getItems(grass, dim = 3) <- c("vegc", "soilc", "litc")
  # nolint grass <- .collectCarbonPools(run = "grass", pools = c("vegc", "soilc", "litc"))

  # Calculate topsoil carbon
  topsoilc <- calcOutput("TopsoilCarbon",
                         lpjml       = lpjml,
                         climatetype = climatetype,
                         aggregate   = FALSE)[, years, ]

  # Calculate soil carbon loss share
  cshare <- calcOutput("SOCLossShare", aggregate = FALSE)

  # -----------------------------------------------------------------------------------------------
  # Format output object

  carbon <- new.magpie(
    cells_and_regions = getCells(pnv),
    years             = getYears(pnv),
    names             = getNames(pnv),
    sets              = c("x.y.iso", "year", "data")
  )

  carbon <- add_dimension(carbon,
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

  # Cropland
  # HACKATHON: We maintain the old version for testing, but will set crop.vegc to 0 for simplicity soon.
  # Factor 0.012 is based on the script subversion/svn/tools/carbon_cropland, executed at 30.07.2013
  carbon[, , "crop.vegc"]  <- 0.012 * pnv[, , "vegc"]

  # nolint carbon[, , "crop.vegc"]  <- 0 # crops have marginal vegetation carbon
  carbon[, , "crop.litc"]  <- 0 # crops do not have litter
  carbon[, , "crop.soilc"] <- cshare * topsoilc + (pnv[, , "soilc"] - topsoilc)

  # Pastures
  carbon[, , "past"]       <- grass
  carbon[, , "past.soilc"] <- pnv[, , "soilc"] # past.soilc is better captured by the pnv run

  # Forest type
  carbon[, , "primforest"] <- pnv
  carbon[, , "secdforest"] <- pnv
  carbon[, , "other"]      <- pnv
  carbon[, , "forestry"]   <- pnv

  # Urban
  carbon[, , "urban"]       <- 0
  carbon[, , "urban.soilc"] <- pnv[, , "soilc"]

  mstools::toolExpectTrue(
    !any(is.na(carbon)),
    "carbon densities contain no NAs",
    falseStatus = "warn"
  )

  mstools::toolExpectTrue(
    all(carbon > 0),
    "carbon densities are all positive",
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
    cells_and_regions = getCells(carbon),
    years = getYears(carbon),
    names = getNames(carbon),
    sets  = c("x.y.iso", "year", "data")
  )

  landArea <- dimSums(landIni, dim = 3)
  weight[, , ] <- landArea
  weight[, , c("primforest", "secdforest", "forestry")] <- potForestArea

  return(
    list(
      x            = carbon,
      weight       = weight + 10^-10,
      min          = 0,
      unit         = "t C per ha",
      description  = "Carbon in tons per hectar for different land use types and carbon pools",
      isocountries = FALSE
    )
  )
}
