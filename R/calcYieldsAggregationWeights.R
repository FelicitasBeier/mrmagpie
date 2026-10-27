#' @title calcYieldsAggregationWeights
#'
#' @description Calculates the crop area weightings to use for yield aggregation.
#'
#' @param weighting     area-based weighting options for yield aggregation
#'                      ("totalCrop" (default),
#'                      "totalLUspecific", "cropSpecific", "crop+irrigSpecific",
#'                      "avlCropland", "avlCropland+potentiallyIrrigatedAreas",
#'                      "avlCropland+avlPasture")
#' @param lpjml         lpjml version, only required if potentially irrigated areas
#'                      are used as aggregation weight
#' @param climatetype   different climate scenarios, only required if potentially
#'                      irrigated areas are used as aggregation weight
#' @param multicropping Multicropping activated (TRUE) or not (FALSE) and
#'                      Multiple Cropping Suitability mask selected
#'                      (e.g. TRUE:actual:total; TRUE:none; FALSE)
#' @param marginal_land  Defines which share of marginal land should be included
#'                       if avlCropland is chosen as weighting option.
#'
#' @return magpie object in cellular resolution
#'
#' @author Kristine Karstens, Felicitas Beier
#'
#' @examples
#' \dontrun{
#' calcOutput("YieldsAggregationWeights", aggregate = FALSE)
#' }
#'
#' @importFrom magpiesets findset
#' @importFrom magclass dimSums new.magpie add_dimension setYears mbind
#' @importFrom madrat toolGetMapping
#' @importFrom mstools toolGetMappingCoord2Country

calcYieldsAggregationWeights <- function(weighting = "totalCrop",
                                         lpjml = NULL, climatetype = NULL,
                                         multicropping = NULL,
                                         marginal_land = "q33_marginal:rainfed_and_irrigated") { # nolint

  # extract dimension information
  yieldNames <- toolGetMapping("MAgPIE_LPJmL.csv", type = "sectoral",
                               where = "mappingfolder")$MAgPIE
  isos       <- toolGetMappingCoord2Country()
  yieldCells <- paste(isos$coords, isos$iso, sep = ".")

  # magpie crop types
  kcr <- magpiesets::findset("kcr")

  # object with correct dimensions
  cropAreaWeight <- add_dimension(new.magpie(cells_and_regions = yieldCells,
                                             years = NULL,
                                             names = yieldNames,
                                             fill = NA),
                                  dim = 3.2,
                                  add = "irrigation",
                                  nm = c("rainfed", "irrigated"))

  #########################################################
  ############ Weight for spatial aggregation #############
  #########################################################
  #### Current cropland as basis for aggregation ####
  if (weighting == "totalCrop") {
    # croparea in initialization year as weight for all irrigated/rainfed crops and pasture areas
    totalCroparea <- dimSums(calcOutput("Croparea", sectoral = "kcr", physical = TRUE, irrigation = FALSE,
                                        cellular = TRUE, aggregate = FALSE,
                                        years = "y1995", round = 6),
                             dim = 3)
    cropAreaWeight[, , ] <- totalCroparea

  } else if (weighting %in% c("totalLUspecific", "cropSpecific", "crop+irrigSpecific")) {
    # crop area in initialization year
    crop <- calcOutput("Croparea", sectoral = "kcr", physical = TRUE, irrigation = TRUE,
                       cellular = TRUE, aggregate = FALSE, years = "y1995", round = 6)
    # pasture area in initialization year
    past <- calcOutput("LanduseInitialisation", aggregate = FALSE, cellular = TRUE, nclasses = "seven",
                       input_magpie = TRUE, years = "y1995", round = 6)[, , "past"]

    if (weighting == "crop+irrigSpecific") {
      # every irrigated/rainfed crop is weighted with its specific crop area in the initialization year
      cropAreaWeight[, , kcr] <- crop
      cropAreaWeight[, , "pasture"]      <- mbind(setNames(past, "irrigated"),
                                                  setNames(past, "rainfed"))

    } else if (weighting == "cropSpecific") {
      # every crop is weighted with its specific crop area in the initialization year
      cropAreaWeight[, , kcr] <- dimSums(crop, dim = "irrigation")
      cropAreaWeight[, , "pasture"] <- past

    } else {
      # total croparea as weight for crops
      cropAreaWeight[, , kcr] <- dimSums(crop, dim = 3)
      # pasture area as weight for pasture
      cropAreaWeight[, , "pasture"] <- past

    }

    #### Available cropland as basis for aggregation ####
  } else if (weighting == "avlCropland") {
    # available cropland as weight for all irrigated/rainfed crops and pasture areas
    avlCrop <- setNames(calcOutput("AvlCropland", marginal_land = marginal_land,
                                   country_level = FALSE, aggregate = FALSE),
                        NULL)

    cropAreaWeight[, , ] <- avlCrop

  } else if (weighting == "avlCropland+potIrrig+avlPasture") {

    # Available cropland as weight for all rainfed crops
    avlCrop <- setNames(calcOutput("AvlCropland", marginal_land = marginal_land,
                                   country_level = FALSE, aggregate = FALSE),
                        NULL)

    # Potentially irrigated areas in the initialization year are used as weight
    # for irrigated crops and pastures.
    # It is a static (non-scenario dependent) weight that ensures that
    # biophysical potential yields from LPJmL at 0.5 degree spatial resolution
    # are aggregated to MAgPIE cluster-level resolution under consideration
    # of not only available cropland (as in the case of rainfed yields),
    # but also criteria that determine whether irrigation is possible in the
    # specific location. This includes water availability and yield gain through
    # irrigation.
    # It is not equivalent to currently irrigated areas in 1995. It considers
    # where irrigation would be possible under the climatic conditions in 1995
    # and the given scenario assumptions.
    # The weight has to be the same for the different yield scenarios (gsadapt vs.
    # constgsadapt) that enter MAgPIE because their difference should merely
    # capture the different assumptions in terms of growing period, not differences
    # in irrigated areas that affect yields. This is a separate management decision
    # that is handled via potentially irrigated areas within MAgPIE and should
    # not enter the weighting.
    pia <- calcOutput("PotIrrigAreas", cropAggregation = TRUE,
                      lpjml = gsub("+scen_constgsadapt_crops", "", lpjml, fixed = TRUE),
                      climatetype = climatetype,
                      multicropping = multicropping,
                      # To Do: (as part of mrwater integration)
                      #        create new function (calcPotIrrigAreasMAgPIE)
                      #        that returns PIA as required by MAgPIE
                      #        and set standard settings there. Then call here and select
                      #        BAU and iniyear.
                      selectyears = seq(1995, 2100, by = 5), iniyear = 1995,
                      efrMethod = "VMF:fair", irrigationsystem = "initialization",
                      accessibilityrule = "CV:2", rankmethod = "USD_m3:GLO:TRUE",
                      gainthreshold = 10, allocationrule = "optimization",
                      yieldcalib = FALSE, comAg = TRUE, ### double-check: shouldn't yield calib be activated?
                      fossilGW = TRUE, transDist = 100,
                      landScen = "potCropland:NULL", cropmix = "hist_total",
                      aggregate = FALSE)[, "y1995", "off"][, , "ssp2"]

    # Read in land use classes in 1995 for pasture area weight
    lu1995  <- setYears(calcOutput("LanduseInitialisation", aggregate = FALSE,
                                   cellular = TRUE, nclasses = "seven",
                                   input_magpie = TRUE, years = "y1995", round = 6),
                        NULL)

    # Use available cropland for all rainfed crops
    cropAreaWeight[, , kcr][, , "rainfed"] <- avlCrop
    # Use available cropland or other land use classes in 1995 if these exceed
    # available cropland for pasture yields
    cropAreaWeight[, , "pasture"] <- pmax(avlCrop,
                                          dimSums(lu1995[, , c("crop",
                                                               "primforest",
                                                               "secdforest",
                                                               "forestry",
                                                               "past")],
                                                  dim = 3))
    # Use potentially irrigated areas for irrigated crop and pasutre yields
    cropAreaWeight[, , "irrigated"] <- pia

    ## To Do: (as part of cropare downscaling improvement)
    # add additional if here that excludes protected areas
    # and use new croparea suitability rules
    # calcOutput("ProtectedAreaBaseline", nclasses = "seven",
    # magpie_input = TRUE,
    # aggregate = FALSE, round = NULL

  } else if (weighting == "avlCropland+avlPasture") {
    # available cropland as weight for all irrigated/rainfed crops
    avlCrop <- setNames(calcOutput("AvlCropland", marginal_land = marginal_land,
                                   country_level = FALSE, aggregate = FALSE),
                        "avlCrop")
    # land use classes in 1995
    lu1995  <- setYears(calcOutput("LanduseInitialisation", aggregate = FALSE,
                                   cellular = TRUE, nclasses = "seven",
                                   input_magpie = TRUE, years = "y1995", round = 6),
                        NULL)

    cropAreaWeight[, , kcr] <- avlCrop

    # weight for pasture is available cropland or other land use classes in 1995 if these exceed available cropland
    cropAreaWeight[, , "pasture"] <- pmax(avlCrop,
                                          dimSums(lu1995[, , c("crop",
                                                               "primforest",
                                                               "secdforest",
                                                               "forestry",
                                                               "past")],
                                                  dim = 3))

  } else {
    warning("No area-based weight for yield aggregation selected or ",
            "weighting setting is not available. ",
            "Using uniform weight.")
    cropAreaWeight <- add_dimension(new.magpie(cells_and_regions = yieldCells,
                                               years = NULL,
                                               names = yieldNames,
                                               fill = 1),
                                    dim = 3.2,
                                    add = "irrigation",
                                    nm = c("rainfed", "irrigated"))
  }

  if (any(is.na(cropAreaWeight))) stop("NAs in weights.")

  # add small number in case of 0s in cropAreaWeight
  x <- cropAreaWeight + 10e-10

  return(list(x            = x,
              weight       = NULL,
              unit         = "Mha",
              description  = "Area-based weight for yields",
              isocountries = FALSE))
}
