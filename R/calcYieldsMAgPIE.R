#' @title calcYieldsMAgPIE
#'
#' @description
#' Orchestrates MAgPIE yields by combining LPJmL-based yields, optional ISIMIP
#' replacement, India yield scaling, and aggregation weights.
#'
#' @param datasource    Choose LPJmL version for main crop inputs and optionally ISIMIP version
#'                      For ISIMIP choose crop model/gcm/rcp/co2 combination formatted like this:
#'                      "yields:EPIC-IIASA:ukesm1-0-ll:ssp585:default:3b"
#' @param climatetype   Climate scenario or historical baseline
#' @param selectyears   Years to be returned
#' @param multicropping Multicropping activated (TRUE) or not (FALSE) and
#'                      Multiple Cropping Suitability mask selected
#' @param calibration   NULL for uncalibrated yields or list for FAO calibration
#' @param weighting     area-based weighting options for yield aggregation
#' @param marginal_land marginal land setting used for weighting
#' @param indiaYields   If TRUE returns scaled yields for rainfed crops in India
#' @param scaleFactor   Integer value by which indiaYields will be scaled
#'
#' @return magpie object in cellular resolution
#'
#' @author Kristine Karstens, Felicitas Beier
#'
#' @examples
#' \dontrun{
#' calcOutput("YieldsMAgPIE", aggregate = FALSE)
#' }
#'
#' @importFrom magclass getYears getNames as.magpie
#' @importFrom mrlandcore toolLPJmLHarmonize
#' @importFrom mstools toolHarmonize2Baseline

calcYieldsMAgPIE <- function(datasource = list(lpjml = "ggcmi_phase3_nchecks_9ca735cb", isimip = NULL),
                             climatetype = "GSWP3-W5E5:historical",
                             selectyears = seq(1965, 2100, by = 5),
                             multicropping = FALSE,
                             calibration = NULL,
                             weighting = "totalCrop",
                             marginal_land = "q33_marginal:rainfed_and_irrigated", # nolint [object_name_liner]
                             indiaYields = FALSE,
                             scaleFactor = 0.3) { # nolint

  # LPJmL-based yields (MAgPIE crops)
  yields <- calcOutput("YieldsMAgPIEcrops",
                       lpjml = datasource[["lpjml"]],
                       climatetype = climatetype,
                       selectyears = selectyears,
                       multicropping = multicropping,
                       calibration = calibration,
                       aggregate = FALSE)

  # Optional: use ISIMIP yields instead of LPJmL yields
  if (!is.null(datasource[["isimip"]])) {
    isimipYields <- calcOutput("ISIMIP3bYields", subtype = datasource[["isimip"]],
                               cells = "lpjcell", aggregate = FALSE)
    commonVars   <- intersect(getNames(yields), getNames(isimipYields))
    commonYears  <- intersect(getYears(yields), getYears(isimipYields))

    # harmonize to LPJml
    cfg       <- toolLPJmLHarmonize(lpjmlversion = datasource[["lpjml"]],
                                    climatetype  = climatetype)
    repHarmon <- toolHarmonize2Baseline(x = isimipYields[, commonYears, commonVars],
                                        base = yields[, commonYears, commonVars],
                                        ref_year = cfg$refYearGcm)
    gc()
    # convert to array for memory
    yields    <- as.array(yields)
    repHarmon <- as.array(repHarmon)
    yields[, commonYears, commonVars] <- repHarmon[, commonYears, commonVars]
    yields    <- as.magpie(yields, spatial = 1)
    repHarmon <- as.magpie(repHarmon, spatial = 1)

    #### Question (Jens, Kristine, Jan): Would we want to be able to apply
    #### the multiple cropping logic on ISIMIP yields as well?
    #### If so: I would suggest to use the same yield increase factor
    #### (as calculated from LPJmL grass GPP)
    #### and move the multicropping yield increase calculation from calcYieldsLPJmL
    #### to a separate tool-/calc-Function
  }

  # Special case for India case study
  if (indiaYields) {
    yields["IND", , "rainfed"] <- yields["IND", , "rainfed"] * scaleFactor
  }

  # Weights
  weights <- calcOutput("YieldsAggregationWeights",
                        weighting = weighting,
                        marginal_land = marginal_land,
                        lpjml = datasource[["lpjml"]],
                        climatetype = climatetype,
                        multicropping = multicropping,
                        aggregate = FALSE)

  return(list(x            = yields,
              weight       = weights,
              unit         = "tDM per ha",
              description  = "Yields for different MAgPIE crop types",
              isocountries = FALSE))
}
