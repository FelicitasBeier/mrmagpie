#' @title calcEFRRockstroem
#' @description This function calculates environmental flow requirements (EFR) for MAgPIE
#'              retrieved from LPJmL monthly discharge and water availability
#'              following the definition of the planetary boundary in Rockström et al. 2023
#'
#' @param lpjml       Defines LPJmL version for crop/grass and natveg specific inputs
#' @param climatetype Switch between different climate scenarios
#' @param stage       Degree of processing: raw, smoothed, harmonized, harmonized2020
#' @param seasonality grper (default): EFR in growing period per year;
#'                    total: EFR throughout the year;
#'                    monthly: monthly EFRs
#'
#' @import magclass
#' @import madrat
#' @importFrom stats quantile
#' @importFrom mstools toolHarmonize2Baseline
#' @importFrom mrlandcore toolLPJmLHarmonization
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier, Jens Heinke
#'
#' @examples
#' \dontrun{
#' calcOutput("EFRRockstroem", aggregate = FALSE)
#' }
#'

calcEFRRockstroem <- function(lpjml = "lpjml5.9.5-m1", climatetype = "MRI-ESM2-0:ssp370",
                              stage = "harmonized2020", seasonality = "grper") {
  # extract LPJmL config information
  cfg <- toolLPJmLHarmonization(lpjmlversion = lpjml,
                                climatetype = climatetype)

  #############################################################################
  # Definition of planetary boundary (PB) according to Rockström et al. 2023: #
  # less than 20% magnitude monthly surface flow alteration                   #
  # in all grid cells                                                         #
  # Translation to EFR:                                                       #
  # only 20% of monthly water flow can be withdrawn,                          #
  # i.e. 80% need to stay in the river in each grid cell                      #
  #############################################################################

  if (stage %in% c("raw", "smoothed")) {
    # Available water per month (smoothed)
    avlWaterMonth <- calcOutput("AvlWater", lpjml = lpjml, climatetype = climatetype,
                                seasonality = "monthly", stage = "smoothed",
                                aggregate = FALSE)

    # Monthly EFR: 80% of monthly available water
    efr <- 0.8 * avlWaterMonth

    ### aggregation to grper, total
    ### efr per cell per month
    if (seasonality == "monthly") {
      # Check for NAs
      if (any(is.na(efr))) {
        stop("calcEFRRockstroem produced NA monthly EFR")
      }
      out <- efr

      ### Total water available per cell per year
    } else if (seasonality == "total") {
      # Sum up over all month:
      efrTotal <- dimSums(efr, dim = 3)

      # Read in available water (for Smakthin calculation)
      avlWaterTotal <- calcOutput("AvlWater", lpjml = lpjml, climatetype = climatetype,
                                  seasonality = "total", stage = "smoothed",
                                  aggregate = FALSE)

      # Reduce EFR to 80% of available water where it exceeds this threshold
      efrTotal[which(efrTotal / avlWaterTotal > 0.8)] <-
        0.8 * avlWaterTotal[which(efrTotal / avlWaterTotal > 0.8)]

      # Check for NAs
      if (any(is.na(efrTotal))) {
        stop("calcEFRRockstroem produced NA efrTotal")
      }
      out <- efrTotal

      ### Water available in growing period per cell per year
    } else if (seasonality == "grper") {
      # magpie object with days per month with same dimension as EFR
      tmp <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      monthDays <- new.magpie(names = dimnames(efr)[[3]])
      monthDays[, , ] <- tmp
      monthDayMagpie <- as.magpie(efr)
      monthDayMagpie[, , ] <- 1
      monthDayMagpie <- monthDayMagpie * monthDays

      # Daily environmental flow requirements
      efrDay   <- efr / monthDayMagpie

      # Growing days per month
      growDays <- calcOutput("GrowingPeriod", lpjml = lpjml, climatetype = climatetype,
                             stage = "smoothed", yield_ratio = 0.1,
                             aggregate = FALSE)

      # Available water in growing period
      efrGrper <- efrDay * growDays
      # Available water in growing period per year
      efrGrper <- dimSums(efrGrper, dim = 3)
      # Read in available water (for Smakthin calculation)
      avlWaterGrper <- calcOutput("AvlWater", lpjml = lpjml, climatetype = climatetype,
                                  seasonality = "grper", stage = "smoothed",
                                  aggregate = FALSE)

      # Reduce EFR to 80% of available water where it exceeds this threshold
      efrGrper[which(efrGrper / avlWaterGrper > 0.8)] <-
        0.8 * avlWaterGrper[which(efrGrper / avlWaterGrper > 0.8)]

      # Check for NAs
      if (any(is.na(efrGrper))) {
        stop("calcEFRRockstroem produced NA efrGrper")
      }
      out <- efrGrper
    } else {
      stop("Specify seasonality: monthly, grper or total")
    }

  } else if (stage == "harmonized") {
    # Load baseline and climate EFR:
    baseline <- calcOutput("EFRRockstroem", lpjml = lpjml, climatetype = cfg$baselineHist,
                           seasonality = seasonality, stage = "smoothed",
                           aggregate = FALSE)

    if (climatetype == cfg$baselineHist) {

      out <- baseline

    } else {

      x   <- calcOutput("EFRRockstroem", lpjml = lpjml, climatetype = climatetype,
                        seasonality = seasonality, stage = "smoothed",
                        aggregate = FALSE)
      # Harmonize to baseline
      out <- toolHarmonize2Baseline(x = x, base = baseline, ref_year = cfg$baselineHist)
    }

  } else if (stage == "harmonized2020") {

    baseline2020 <- calcOutput("EFRRockstroem", lpjml = lpjml, climatetype = cfg$baselineGcm,
                               seasonality = seasonality, stage = "harmonized",
                               aggregate = FALSE)

    if (climatetype == cfg$baselineGcm) {

      out <- baseline2020

    } else {

      x        <- calcOutput("EFRRockstroem", lpjml = lpjml, climatetype = climatetype,
                             seasonality = seasonality, stage = "smoothed",
                             aggregate = FALSE)
      out <- toolHarmonize2Baseline(x, baseline2020, ref_year = cfg$refYearGcm)
    }

  } else {
    stop("Stage argument not supported!")
  }

  description <- paste0("EFRs according to water planetary boundary in ", seasonality)

  return(list(x = out,
              weight = NULL,
              unit = "mio. m^3",
              description = description,
              isocountries = FALSE))
}
