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
#' @param cells       lpjcell for 67420 cells or magpiecell for 59199 cells
#'
#' @importFrom mstools toolHarmonize2Baseline
#' @importFrom mrlandcore toolLPJmLVersion
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier, Abhijeet Mishra
#'
#' @examples
#' \dontrun{
#' calcOutput("EFRSmakthin", aggregate = FALSE)
#' }
#'
calcEFRSmakthin <- function(lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de",
                                      crop = "ggcmi_phase3_nchecks_9ca735cb"),
                            climatetype = "GSWP3-W5E5:historical", stage = "harmonized2020",
                            LFR_val = 0.1, HFR_LFR_less10 = 0.2, HFR_LFR_10_20 = 0.15, # nolint
                            HFR_LFR_20_30 = 0.07, HFR_LFR_more30 = 0.00,               # nolint
                            seasonality = "grper", cells = "lpjcell") {
  # Create settings for LPJmL from version and climatetype argument
  cfgNatveg <- toolLPJmLVersion(version = lpjml["natveg"], climatetype = climatetype)
  cfgCrop   <- toolLPJmLVersion(version = lpjml["crop"],   climatetype = climatetype)

  lpjmlReadin  <- c(natveg = unname(cfgNatveg$readin_version),
                    crop   = unname(cfgCrop$readin_version))

  lpjmlBaseline <- c(natveg = unname(cfgNatveg$baseline_version),
                     crop   = unname(cfgCrop$baseline_version))

  if (stage %in% c("raw", "smoothed")) {
    ############################################################
    # Step 1 Determine monthly discharge low flow requirements #
    #        (lfrMonthlyDischarge)                           #
    ############################################################

    ### Monthly Discharge
    monthlyDischargeMagpie <- calcOutput("LPJmL_new", version = lpjmlReadin["natveg"],
                                         climatetype = climatetype, subtype = "mdischarge",
                                         aggregate = FALSE, stage = "raw")
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
    # The window length is constant (8 years * 12 months), so the type-7
    # interpolation index is the same for every cell and every year. Precompute
    # it once and use an inline single-probability type-7 estimator (with a
    # partial sort) instead of stats::quantile, which is numerically identical
    # but avoids per-call validation/multi-prob overhead (millions of calls).
    windowLength      <- 8L * 12L
    quantileIndex     <- 1 + (windowLength - 1) * LFR_val
    lowerRank         <- floor(quantileIndex)
    interpolationFrac <- quantileIndex - lowerRank
    lfrQuantile <- function(dischargeWindow) {
      sorted <- sort.int(dischargeWindow, partial = c(lowerRank, lowerRank + 1L))
      (1 - interpolationFrac) * sorted[lowerRank] + interpolationFrac * sorted[lowerRank + 1L]
    }
    for (year in years) {
      # get the LFR_val quantile in range of 8 years for each year for all cells
      neededYears <- seq(year - 7, year, by = 1)
      lfrQuant[, paste("y", year, sep = "")] <-
        apply(monthlyDischargeMagpie[, paste("y", neededYears, sep = ""), ],
              MARGIN = c(1), lfrQuantile)
    }
    # Time-smooth lfrQuant
    lfrQuant <- as.magpie(lfrQuant, spatial = 1)
    lfrQuant <- toolFillYears(lfrQuant, getYears(monthlyDischargeMagpie, as.integer = TRUE))

    if (stage == "smoothed") lfrQuant <- mstools::toolSmooth(lfrQuant)

    # Raw monthly discharge no longer needed at this point
    rm(monthlyDischargeMagpie)

    ### Read in smoothed monthly discharge
    monthlyDischargeMagpie <- calcOutput("LPJmL_new", version = lpjmlReadin["natveg"],
                                         climatetype = climatetype, subtype = "mdischarge",
                                         aggregate = FALSE, stage = "smoothed")

    # Transform to array (faster calculation)
    lfrQuant <- as.array(collapseNames(lfrQuant))
    monthlyDischargeMagpie <- as.array(collapseNames(monthlyDischargeMagpie))

    ### Calculate LFR discharge values for each month
    # If lfrQuant < magpie_discharge: take lfrQuant
    # Else: take magpie_discharge
    # lfrQuant is cells x years; monthlyDischargeMagpie is cells x years x 12.
    # Column-major recycling of the cells*years vector across each month block
    # makes this equivalent to the per-month pmin loop, with fewer allocations.
    lfrMonthlyDischarge    <- monthlyDischargeMagpie
    lfrMonthlyDischarge[]  <- pmin(as.vector(lfrQuant), monthlyDischargeMagpie)
    # Remove no longer needed objects
    rm(lfrQuant)


    ################################################
    # Step 2 Determine low flow requirements (LFR) #
    #        from available water per month        #
    ################################################
    ### Available water per month (smoothed)
    avlWaterMonth <- calcOutput("AvlWater", lpjml = lpjmlReadin, climatetype = climatetype,
                                seasonality = "monthly", stage = "smoothed",
                                aggregate = FALSE, cells = "lpjcell")

    # Transform to array for faster calculation
    avlWaterMonth <- as.array(collapseNames(avlWaterMonth))

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
    # hfr allocated with the right shape; every cell is assigned below.
    hfr <- lfr

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
      avlWaterTotal <- calcOutput("AvlWater", lpjml = lpjmlReadin, climatetype = climatetype,
                                  seasonality = "total", stage = "smoothed",
                                  aggregate = FALSE, cells = "lpjcell")

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
      growDays <- calcOutput("GrowingPeriod", lpjml = lpjmlReadin, climatetype = climatetype,
                             stage = "smoothed", yield_ratio = 0.1,
                             aggregate = FALSE, cells = "lpjcell")

      # Available water in growing period
      efrGrper <- efrDay * growDays
      # Available water in growing period per year
      efrGrper <- dimSums(efrGrper, dim = 3)
      # Read in available water (for Smakthin calculation)
      avlWaterGrper <- calcOutput("AvlWater", lpjml = lpjmlReadin, climatetype = climatetype,
                                  seasonality = "grper", stage = "smoothed",
                                  aggregate = FALSE, cells = "lpjcell")

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
    # Load baseline and climate EFR:
    baseline <- calcOutput("EFRSmakthin", lpjml = lpjmlBaseline, climatetype = cfgNatveg$baseline_hist,
                           seasonality = seasonality, stage = "smoothed",
                           LFR_val = LFR_val,
                           HFR_LFR_less10 = HFR_LFR_less10, HFR_LFR_10_20 = HFR_LFR_10_20,
                           HFR_LFR_20_30 = HFR_LFR_20_30, HFR_LFR_more30 = HFR_LFR_more30,
                           aggregate = FALSE, cells = "lpjcell")

    if (climatetype == cfgNatveg$baseline_hist) {

      out <- baseline

    } else {

      x   <- calcOutput("EFRSmakthin", lpjml = lpjmlReadin, climatetype = climatetype,
                        seasonality = seasonality, stage = "smoothed",
                        LFR_val = LFR_val,
                        HFR_LFR_less10 = HFR_LFR_less10, HFR_LFR_10_20 = HFR_LFR_10_20,
                        HFR_LFR_20_30 = HFR_LFR_20_30, HFR_LFR_more30 = HFR_LFR_more30,
                        aggregate = FALSE, cells = "lpjcell")
      # Harmonize to baseline
      out <- toolHarmonize2Baseline(x = x, base = baseline, ref_year = cfgNatveg$ref_year_hist)
    }

  } else if (stage == "harmonized2020") {

    baseline2020 <- calcOutput("EFRSmakthin", lpjml = lpjmlBaseline, climatetype = cfgNatveg$baseline_gcm,
                               seasonality = seasonality, stage = "harmonized",
                               LFR_val = LFR_val,
                               HFR_LFR_less10 = HFR_LFR_less10, HFR_LFR_10_20 = HFR_LFR_10_20,
                               HFR_LFR_20_30 = HFR_LFR_20_30, HFR_LFR_more30 = HFR_LFR_more30,
                               aggregate = FALSE, cells = "lpjcell")

    if (climatetype == cfgNatveg$baseline_gcm) {

      out <- baseline2020

    } else {

      x        <- calcOutput("EFRSmakthin", lpjml = lpjmlReadin, climatetype = climatetype,
                             seasonality = seasonality, stage = "smoothed",
                             LFR_val = LFR_val,
                             HFR_LFR_less10 = HFR_LFR_less10, HFR_LFR_10_20 = HFR_LFR_10_20,
                             HFR_LFR_20_30 = HFR_LFR_20_30, HFR_LFR_more30 = HFR_LFR_more30,
                             aggregate = FALSE, cells = "lpjcell")
      out      <- toolHarmonize2Baseline(x, baseline2020, ref_year = cfgNatveg$ref_year_gcm)
    }

  } else {
    stop("Stage argument not supported!")
  }

  description <- paste0("EFR according to Smakthin method in ", seasonality)

  if (cells == "magpiecell") {
    out <- toolCoord2Isocell(out, cells = cells)
  }

  return(list(x = out,
              weight = NULL,
              unit = "mio. m^3",
              description = description,
              isocountries = FALSE))
}
