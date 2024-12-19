#' @title calcEFRSmakthin
#' @description This function calculates environmental flow requirements (EFR) for MAgPIE
#'              retrieved from LPJmL monthly discharge and water availability using the
#'              method of Smakthin et al. (2004)
#'
#' @param lpjml Defines LPJmL version for crop/grass and natveg specific inputs
#' @param climatetype Switch between different climate scenarios
#' @param stage Degree of processing: raw, smoothed, harmonized, harmonized2020
#' @param LFR_val Strictness of environmental flow requirements
#' @param HFR_LFR_less10 High flow requirements (share of total water for cells)
#'                       with LFR<10percent of total water
#' @param HFR_LFR_10_20 High flow requirements (share of total water for cells)
#'                      with 10percent < LFR < 20percent of total water
#' @param HFR_LFR_20_30 High flow requirements (share of total water for cells)
#'                      with 20percent < LFR < 30percent of total water
#' @param HFR_LFR_more30 High flow requirements (share of total water for cells)
#'                       with LFR>30percent of total water
#' @param seasonality grper (default): EFR in growing period per year; total:
#'                    EFR throughout the year; monthly: monthly EFRs
#'
#' @import magclass
#' @import madrat
#' @importFrom stats quantile
#' @importFrom mstools toolHarmonize2Baseline
#' @importFrom mrlandcore toolLPJmLHarmonize
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier, Abhijeet Mishra
#'
#' @examples
#' \dontrun{
#' calcOutput("EFRSmakthin", aggregate = FALSE)
#' }
#'
calcEFRSmakthin <- function(lpjml = "lpjml5.9.5-m1", climatetype = "MRI-ESM2-0:ssp370",
                            stage = "harmonized2020",
                            LFR_val = 0.1, HFR_LFR_less10 = 0.2, HFR_LFR_10_20 = 0.15, # nolint
                            HFR_LFR_20_30 = 0.07, HFR_LFR_more30 = 0.00,               # nolint
                            seasonality = "grper") {
  # extract LPJmL config information
  cfg <- toolLPJmLHarmonize(lpjmlversion = lpjml,
                            climatetype = climatetype)

  if (stage %in% c("raw", "smoothed")) {
    ############################################################
    # Step 1 Determine monthly discharge low flow requirements #
    #        (lfrMonthlyDischarge)                           #
    ############################################################

    ### Monthly Discharge
    monthlyDischargeMagpie <- calcOutput("LPJmLTransform", subtype = "pnv:discharge",
                                         lpjmlversion = lpjml, climatetype = climatetype,
                                         stage = "raw:cut", aggregate = FALSE)
    # Extract years for quantile calculation
    years <- getYears(monthlyDischargeMagpie, as.integer = TRUE)
    years <- seq(years[1] + 7, years[length(years)], by = 1)
    # Transform to array (faster calculation)
    monthlyDischargeMagpie <-  as.array(collapseNames(monthlyDischargeMagpie))
    # Empty array with magpie object names
    lfrQuant <- array(NA, dim = c(dim(monthlyDischargeMagpie)[1], length(years)),
                      dimnames = list(dimnames(monthlyDischargeMagpie)[[1]], paste("y", years, sep = "")))

    ### Calculate lfrQuant
    ## Note: LFRs correspond to the Q90-value (fair condition)
    ## (i.e. to the discharge that is exceeded in nine out of ten months)
    ## Q75-value (good condition), or Q50 (natural condition) (Bonsch et al. 2015).
    ## This is calculated via the 10%- (25%-, 50%-) quantile of monthly discharge.

    # Quantile calculation: Yearly LFR quantile value
    for (year in years) {
      # get the LFR_val quantile in range of 8 years for each year for all cells
      neededYears <- seq(year - 7, year, by = 1)
      lfrQuant[, paste("y", year, sep = "")] <-
        apply(monthlyDischargeMagpie[, paste("y", neededYears, sep = ""), ],
              MARGIN = c(1), quantile, probs = LFR_val)
    }
    # Time-smooth lfrQuant
    lfrQuant <- as.magpie(lfrQuant, spatial = 1)
    lfrQuant <- toolFillYears(lfrQuant, getYears(monthlyDischargeMagpie, as.integer = TRUE))

    if (stage == "smoothed") lfrQuant <- toolSmooth(lfrQuant)

    # Raw monthly discharge no longer needed at this point
    rm(monthlyDischargeMagpie)

    ### Read in smoothed monthly discharge
    monthlyDischargeMagpie <- calcOutput("LPJmLTransform", subtype = "pnv:discharge",
                                         lpjmlversion = lpjml, climatetype = climatetype,
                                         stage = "smoothed:cut", aggregate = FALSE)

    # Transform to array (faster calculation)
    lfrQuant <- as.array(collapseNames(lfrQuant))
    monthlyDischargeMagpie <- as.array(collapseNames(monthlyDischargeMagpie))

    ### Calculate LFR discharge values for each month
    # If lfrQuant < magpie_discharge: take lfrQuant
    # Else: take magpie_discharge
    lfrMonthlyDischarge <- monthlyDischargeMagpie
    for (month in 1:12) {
      tmp1 <- as.vector(lfrQuant)
      tmp2 <- as.vector(monthlyDischargeMagpie[, , month])
      lfrMonthlyDischarge[, , month] <- pmin(tmp1, tmp2)
    }
    # Remove no longer needed objects
    rm(lfrQuant)


    ################################################
    # Step 2 Determine low flow requirements (LFR) #
    #        from available water per month        #
    ################################################
    ### Available water per month (smoothed)
    avlWaterMonth <- calcOutput("AvlWater", lpjml = lpjml, climatetype = climatetype,
                                seasonality = "monthly", stage = "smoothed",
                                aggregate = FALSE)

    # Transform to array for faster calculation
    avlWaterMonth <- as.array(collapseNames(avlWaterMonth))

    # Empty array
    lfr       <- avlWaterMonth
    lfr[, , ] <- NA

    ### Calculate LFRs
    lfr <- avlWaterMonth * (lfrMonthlyDischarge / monthlyDischargeMagpie)
    # There are NA's where monthlyDischargeMagpie was 0, replace with 0:
    lfr[is.nan(lfr)] <- 0

    ###################################################################
    # Step 3 Determine monthly high flow requirements (HFR)           #
    #        based on the ratio between LFR_month and avlWaterMonth   #
    ###################################################################
    ## Note: "For rivers with low Q90 values, high-flow events are important
    ## for river channel maintenance, wetland flooding, and riparian vegetation.
    ## HFRs of 20% of available water are therefore assigned to rivers with a
    ## low fraction of Q90 in total discharge. Rivers with a more stable flow
    ## regime receive a lower HFR." (Bonsch et al. 2015)
    hfr       <- lfr
    hfr[, , ] <- NA

    hfr[lfr < 0.1 * avlWaterMonth]  <- HFR_LFR_less10 * avlWaterMonth[lfr < 0.1 * avlWaterMonth]
    hfr[lfr >= 0.1 * avlWaterMonth] <- HFR_LFR_10_20  * avlWaterMonth[lfr >= 0.1 * avlWaterMonth]
    hfr[lfr >= 0.2 * avlWaterMonth] <- HFR_LFR_20_30  * avlWaterMonth[lfr >= 0.2 * avlWaterMonth]
    hfr[lfr >= 0.3 * avlWaterMonth] <- HFR_LFR_more30 * avlWaterMonth[lfr >= 0.3 * avlWaterMonth]
    hfr[avlWaterMonth <= 0]       <- 0

    efr <- lfr + hfr
    efr <- as.magpie(efr, spatial = 1)

    ### aggregation to grper, total
    ### efr per cell per month
    if (seasonality == "monthly") {
      # Check for NAs
      if (any(is.na(efr))) {
        stop("calcEFRSmakthin produced NA efr")
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

      # Reduce EFR to 50% of available water where it exceeds this threshold (according to Smakhtin 2004)
      efrTotal[which(efrTotal / avlWaterTotal > 0.5)] <-
        0.5 * avlWaterTotal[which(efrTotal / avlWaterTotal > 0.5)]

      # Check for NAs
      if (any(is.na(efrTotal))) {
        stop("calcEFRSmakthin produced NA efrTotal")
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
      getItems(growDays, dim = 3) <- c(1:12)
      getSets(growDays) <- c("x", "y", "iso", "year", "month")
      yrs <- intersect(getItems(growDays, dim = 2), getItems(efr, dim = 2))

      # Available water in growing period
      efrGrper <- efrDay[, yrs, ] * growDays[, yrs, ]
      # Available water in growing period per year
      efrGrper <- dimSums(efrGrper, dim = 3)
      # Read in available water (for Smakthin calculation)
      avlWaterGrper <- calcOutput("AvlWater", lpjml = lpjml, climatetype = climatetype,
                                  seasonality = "grper", stage = "smoothed",
                                  aggregate = FALSE)

      # Reduce EFR to 50% of available water where it exceeds this threshold (according to Smakhtin 2004)
      efrGrper[which(efrGrper / avlWaterGrper > 0.5)] <-
        0.5 * avlWaterGrper[which(efrGrper / avlWaterGrper > 0.5)]

      # Check for NAs
      if (any(is.na(efrGrper))) {
        stop("calcEFRSmakthin produced NA efrGrper")
      }
      out <- efrGrper
    } else {
      stop("Specify seasonality: monthly, grper or total")
    }

  } else if (stage == "harmonized") {
    # load historical baseline
    baseline <- calcOutput("EFRSmakthin", stage = "smoothed",
                           lpjml = lpjml, climatetype = cfg$baselineHist,
                           seasonality = seasonality,
                           LFR_val = LFR_val,
                           HFR_LFR_less10 = HFR_LFR_less10, HFR_LFR_10_20 = HFR_LFR_10_20,
                           HFR_LFR_20_30 = HFR_LFR_20_30, HFR_LFR_more30 = HFR_LFR_more30,
                           aggregate = FALSE)

    if (climatetype == cfg$baselineHist) {
      # no further harmonization required when climatetype is historical baseline
      out <- baseline

    } else {
      # load smoothed future scenario
      x   <- calcOutput("EFRSmakthin", lpjml = lpjml, climatetype = climatetype,
                        seasonality = seasonality, stage = "smoothed",
                        LFR_val = LFR_val,
                        HFR_LFR_less10 = HFR_LFR_less10, HFR_LFR_10_20 = HFR_LFR_10_20,
                        HFR_LFR_20_30 = HFR_LFR_20_30, HFR_LFR_more30 = HFR_LFR_more30,
                        aggregate = FALSE)
      # harmonize future scenario to baseline
      out <- toolHarmonize2Baseline(x = x, base = baseline, ref_year = cfg$refYearHist)
    }

  } else if (stage == "harmonized2020") {
    # load harmonized baseline GCM scenario
    baseline2020 <- calcOutput("EFRSmakthin", stage = "harmonized",
                               lpjml = lpjml, climatetype = cfg$baselineGcm,
                               seasonality = seasonality,
                               LFR_val = LFR_val,
                               HFR_LFR_less10 = HFR_LFR_less10, HFR_LFR_10_20 = HFR_LFR_10_20,
                               HFR_LFR_20_30 = HFR_LFR_20_30, HFR_LFR_more30 = HFR_LFR_more30,
                               aggregate = FALSE)

    if (climatetype == cfg$baselineGcm) {
      # no further harmonization required if climatetype is baseline GCM
      out <- baseline2020

    } else {
      # load smoothed future scenario
      x        <- calcOutput("EFRSmakthin", stage = "smoothed",
                             lpjml = lpjml, climatetype = climatetype,
                             seasonality = seasonality,
                             LFR_val = LFR_val,
                             HFR_LFR_less10 = HFR_LFR_less10, HFR_LFR_10_20 = HFR_LFR_10_20,
                             HFR_LFR_20_30 = HFR_LFR_20_30, HFR_LFR_more30 = HFR_LFR_more30,
                             aggregate = FALSE)
      # harmonize future scenario to baseline
      out <- toolHarmonize2Baseline(x, baseline2020, ref_year = cfg$refYearGcm)
    }

  } else {
    stop("Stage argument not supported!")
  }

  description <- paste0("EFR according to Smakthin method in ", seasonality)

  return(list(x = out,
              weight = NULL,
              unit = "mio. m^3",
              description = description,
              isocountries = FALSE))
}
