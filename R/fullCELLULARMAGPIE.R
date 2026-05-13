#' @title fullCELLULARMAGPIE
#' @description Function that produces the complete cellular data set required
#'              for running the MAgPIE model.
#'
#' @param rev data revision which should be used as input (numeric_version).
#' @param ctype aggregation clustering type, which is a combination of a single letter,
#'              indicating the cluster methodology, and a number, indicating the number
#'              of resulting clusters. Available methodologies are
#'              - hierarchical clustering (h),
#'              - normalized k-means clustering (n) and
#'              - combined hierarchical/normalized k-means clustering (c).
#'              In the latter hierarchical clustering is used to determine the cluster
#'              distribution among regions whereas normalized k-means is used for the
#'              clustering within a region.
#' @param dev development suffix to distinguish development versions for the same data revision.
#'            This can be useful to distinguish parallel lines of development.
#' @param climatetype Global Circulation Model to be used
#' @param lpjml Defines LPJmL version for crop/grass and natveg specific inputs
#' @param isimip Defines isimip crop model input which replace maiz, tece, rice_pro and soybean
#' @param emu_id Pasture Soil carbon emulator ID
#' @param clusterweight Should specific regions be resolved with more or less detail?
#'                      Values > 1 mean higher share, < 1 lower share
#'                      e.g. cfg$clusterweight <- c(LAM=2) means that
#'                      a higher level of detail for region LAM if set to NULL
#'                      all weights will be assumed to be 1. Examples:
#'                      c(LAM=1.5,SSA=1.5,OAS=1.5) or c(LAM=2,SSA=2,OAS=2)
#' \code{\link[madrat]{setConfig}} (e.g. for setting the mainfolder if not already set properly).
#'
#' @author Kristine Karstens, Jan Philipp Dietrich
#' @seealso
#' \code{\link[madrat]{readSource}},\code{\link[madrat]{getCalculations}},\code{\link[madrat]{calcOutput}},
#' \code{\link[madrat]{setConfig}}
#' @examples
#' \dontrun{
#' retrieveData("CELLULARMAGPIE", rev = numeric_version("12"),
#'              mainfolder = "pathtowhereallfilesarestored")
#' }
#' @importFrom madrat setConfig getConfig
#' @importFrom magpiesets findset
#' @importFrom digest digest
#' @importFrom luplot plotregionscluster
#' @importFrom ggplot2 ggsave
#' @importFrom withr local_options

fullCELLULARMAGPIE <- function(rev = numeric_version("0.1"), dev = "",
                               ctype = "c200",
                               climatetype = "MRI-ESM2-0:ssp370",
                               lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de",
                                         crop = "ggcmi_phase3_nchecks_9ca735cb",
                                         grass = "lpjml5p2_pasture"),
                               isimip = NULL,
                               clusterweight = NULL,
                               emu_id = NULL) { # nolint

  "!# @pucArguments ctype clusterweight"

  withr::local_options(magclass_sizeLimit = 1e+12)

  ### Version settings ###
  if (rev < numeric_version("4.94")) {
    stop("mrmagpie(>= 1.35.2) does not support revision below 4.94 anymore. ",
         "Please use an older snapshot/version of the library, if you need older revisions.")
  }
  cells       <- "lpjcell"

  message(paste0("Start preprocessing for \n climatescenario: ", climatetype,
                 "\n LPJmL-Versions: ", paste(names(lpjml), lpjml, sep = "->", collapse = ", "),
                 "\n clusterweight: ", paste(names(clusterweight), clusterweight, sep = ":", collapse = ", "),
                 "\n isimip yield subtype: ", paste(names(isimip), isimip, sep = ":", collapse = ", ")))

  # Create version tag (will be returned at the very end of this function)
  versionTag <- paste(ctype,
                      gsub(":", "-", climatetype),
                      paste0("lpjml-", digest::digest(lpjml, algo = getConfig("hash"))),
                      sep = "_")
  versionTag <- ifelse(is.null(isimip), versionTag,
                       paste0(versionTag, "_isimip-",
                              digest::digest(isimip, algo = getConfig("hash"))))
  versionTag <- ifelse(is.null(clusterweight), versionTag,
                       paste0(versionTag, "_clusterweight-",
                              digest::digest(clusterweight, algo = getConfig("hash"))))
  versionTag <- ifelse(is.null(emu_id), versionTag,
                       paste0(versionTag, "_gsoilc-", emu_id))


  magYearsPastLong <- c("y1995", "y2000", "y2005", "y2010", "y2015")
  magYears         <- findset("time")
  shortYears       <- findset("t_all")
  lpjYears         <- seq(1995, 2100, by = 5)
  roundArea        <- 5
  stats            <- c("summary", "sum")

  # Clustering based on 67420 cells
  map      <- calcOutput("Cluster", ctype = ctype, weight = clusterweight, lpjml = lpjml,
                         clusterdata = "yield_airrig", aggregate = FALSE)
  weightID <- ifelse(is.null(clusterweight), "", paste0("_", names(clusterweight), clusterweight, collapse = ""))
  clustermapname <- sub("\\.[^.]*$", ".rds",
                        paste0("clustermap_rev", rev, dev, "_", ctype, "_67420",
                               weightID, "_", getConfig("regionmapping")))
  addMapping(clustermapname, map)

  # plot map with regions and clusters
  clustermap <- readRDS(clustermapname) # nolint
  p <- plotregionscluster(clustermap, cells = "lpjcell") # nolint
  suppressWarnings(ggsave(sub(".rds", ".pdf", sub("clustermap", "spamplot", clustermapname)),
                          p, height = 6, width = 10, scale = 1))

  # distinguish between region and superregion if mapping provides this distinction
  mapReg      <- toolGetMapping(getConfig("regionmapping"), type = "regional", where = "mappingfolder")
  superregion <- ifelse("superregion" %in% colnames(mapReg), "superregion", "region")

  # 09 drivers
  calcOutput("GridPop", source = "Gao", subtype = "all", harmonize_until = 2015, urban = FALSE,
             cellular = TRUE, cells = cells, aggregate = "cluster",
             years = magYears, round = 6, outputStatistics = stats, file = "f09_pop_grid.cs3")

  calcOutput("GridPop", source = "Gao", subtype = "all", harmonize_until = 2015, urban = TRUE,
             cellular = TRUE, cells = cells, aggregate = "cluster",
             years = magYears, round = 6, outputStatistics = stats, file = "f09_urbanpop_grid.cs3")

  # 10 land
  ## seven land classes
  calcOutput("LanduseInitialisation", nclasses = "seven",
             aggregate = FALSE, cellular = TRUE, cells = cells,
             input_magpie = TRUE, selectyears = magYearsPastLong,
             round = NULL, outputStatistics = stats, file = "avl_land_t_0.5.mz")
  calcOutput("LanduseInitialisation", nclasses = "seven",
             aggregate = "cluster", cellular = TRUE, cells = cells,
             input_magpie = TRUE, selectyears = magYearsPastLong,
             round = roundArea, outputStatistics = stats, file = paste0("avl_land_t_", ctype, ".mz"))
  calcOutput("LanduseInitialisation", nclasses = "seven",
             aggregate = FALSE, cellular = FALSE, cells = cells,
             input_magpie = TRUE, selectyears = magYearsPastLong,
             round = roundArea, outputStatistics = stats, file = "avl_land_t_iso.cs3")

  ## nine land classes
  calcOutput("LanduseInitialisation", nclasses = "nine",
             aggregate = FALSE, cellular = TRUE, cells = cells,
             input_magpie = TRUE, selectyears = magYearsPastLong,
             round = NULL, outputStatistics = stats, file = "avl_land_full_t_0.5.mz")
  calcOutput("LanduseInitialisation", nclasses = "nine",
             aggregate = "cluster", cellular = TRUE, cells = cells,
             input_magpie = TRUE, selectyears = magYearsPastLong,
             round = roundArea, outputStatistics = stats, file = paste0("avl_land_full_t_", ctype, ".mz"))
  calcOutput("LanduseInitialisation", nclasses = "nine",
             aggregate = FALSE, cellular = FALSE, cells = cells,
             input_magpie = TRUE, selectyears = magYearsPastLong,
             round = roundArea, outputStatistics = stats, file = "avl_land_full_t_iso.cs3")

  calcOutput("AvlLandSi", cells = cells, aggregate = FALSE,
             round = NULL, outputStatistics = stats, file = "avl_land_si_0.5.mz")
  calcOutput("AvlLandSi", cells = cells, aggregate = "cluster",
             round = roundArea, outputStatistics = stats, file = paste0("avl_land_si_", ctype, ".mz"))

  # 13 TC
  calcOutput("PastrTauHist", round = 2, past_mngmt = "mdef",
             outputStatistics = stats, file = "f13_pastr_tau_hist.csv",
             cells = cells, aggregate = superregion)

  # 14 yields
  if (dev == "+calibYield") {

    calcOutput("YieldsCalibrated", aggregate = "cluster", cells = cells,
               source = c(lpjml = lpjml[["crop"]], isimip = isimip),
               climatetype = climatetype, round = 2, years = lpjYears,
               outputStatistics = stats, file = paste0("lpj_yields_", ctype, ".mz"))

  } else if (grepl("india", dev)) {

    calcOutput("Yields", source = c(lpjml = lpjml[["crop"]], isimip = isimip),
               cells = cells, aggregate = FALSE,
               climatetype = climatetype, round = NULL, years = lpjYears,
               outputStatistics = stats, file = "lpj_yields_0.5.mz",
               weighting = "crop+irrigSpecific", indiaYields = TRUE, scaleFactor = 0.5)


    calcOutput("Yields", aggregate = "cluster", cells = cells,
               source = c(lpjml = lpjml[["crop"]], isimip = isimip),
               climatetype = climatetype, round = 2, years = lpjYears,
               outputStatistics = stats, file = paste0("lpj_yields_", ctype, ".mz"),
               weighting = "crop+irrigSpecific", indiaYields = TRUE, scaleFactor = 0.5)

  } else {

    calcOutput("Yields", source = c(lpjml = lpjml[["crop"]], isimip = isimip),
               aggregate = FALSE, cells = cells,
               climatetype = climatetype, round = NULL, years = lpjYears,
               outputStatistics = stats, file = "lpj_yields_0.5.mz",
               weighting = ifelse(grepl("YieldWeights_", dev), gsub("YieldWeights_", "", dev), "totalCrop"))

    calcOutput("Yields", aggregate = "cluster", cells = cells,
               source = c(lpjml = lpjml[["crop"]], isimip = isimip),
               climatetype = climatetype, round = 2, years = lpjYears,
               outputStatistics = stats, file = paste0("lpj_yields_", ctype, ".mz"),
               weighting = ifelse(grepl("YieldWeights_", dev), gsub("YieldWeights_", "", dev), "totalCrop"))

  }

  calcOutput("DegradationYieldReduction", aggregate = "cluster", round = 6, cells = cells,
             outputStatistics = stats, file = paste0("f14_degradation_yld_reduc_", ctype, ".mz"))


  # 22 land conservation
  calcOutput("ProtectedAreaBaseline", nclasses = "seven",
             cells = cells, magpie_input = TRUE,
             aggregate = FALSE, round = NULL, outputStatistics = stats, file = "wdpa_baseline_0.5.mz")
  calcOutput("ProtectedAreaBaseline", nclasses = "seven",
             cells = cells, magpie_input = TRUE,
             aggregate = "cluster", round = roundArea,
             outputStatistics = stats, file = paste0("wdpa_baseline_", ctype, ".mz"))

  calcOutput("ConservationPriorities", nclasses = "seven", cells = cells,
             aggregate = FALSE, round = NULL, outputStatistics = stats, file = "consv_prio_areas_0.5.mz")
  calcOutput("ConservationPriorities", nclasses = "seven", cells = cells,
             aggregate = "cluster", round = roundArea,
             outputStatistics = stats, file = paste0("consv_prio_areas_", ctype, ".mz"))

  calcOutput("ProtectArea", bhifl = TRUE,
             cells = cells, aggregate = "cluster", round = roundArea,
             outputStatistics = stats, file = paste0("protect_area_", ctype, ".mz"))


  # 30 crop
  calcOutput("Croparea", sectoral = "kcr", physical = TRUE,
             cellular = TRUE, irrigation = FALSE, round = roundArea,
             aggregate = "cluster", outputStatistics = stats,
             file = paste0("f30_croparea_initialisation_", ctype, ".mz"))
  calcOutput("Croparea", sectoral = "kcr", physical = TRUE,
             cellular = TRUE, irrigation = TRUE, round = roundArea,
             aggregate = "cluster", outputStatistics = stats,
             file = paste0("f30_croparea_w_initialisation_", ctype, ".mz"))

  calcOutput("MulticroppingIntensity", selectyears = magYears, round = 4,
             scenario = "crop", sectoral = "kcr",
             outputStatistics = stats, file = "fm_multicropping.cs3",
             aggregate = "cluster")

  ## For cellular comparison
  calcOutput("MAPSPAM", subtype = "physical",  aggregate = FALSE, round = NULL,
             outputStatistics = stats, file = paste0("MAPSPAM_croparea_0.5.mz"))
  calcOutput("Croparea", sectoral = "kcr", physical = TRUE, cellular = TRUE,
             irrigation = TRUE, round = NULL,
             aggregate = FALSE, outputStatistics = stats, file = "LUH3_croparea_0.5.mz")

  calcOutput("AvlCropland", marginal_land = "magpie", cell_upper_bound = 0.9,
             aggregate = FALSE, cells = cells,
             round = NULL, outputStatistics = stats, file = "avl_cropland_0.5.mz")
  calcOutput("AvlCropland", marginal_land = "magpie", cell_upper_bound = 0.9,
             aggregate = "cluster", cells = cells,
             round = roundArea, outputStatistics = stats, file = paste0("avl_cropland_", ctype, ".mz"))
  calcOutput("AvlCropland", marginal_land = "magpie", cell_upper_bound = 0.9,
             aggregate = FALSE, cells = cells, country_level = TRUE,
             round = roundArea, outputStatistics = stats, file = "avl_cropland_iso.cs3")

  calcOutput("CroplandTreecover",
             aggregate = FALSE, cells = cells,
             round = NULL, outputStatistics = stats, file = "CroplandTreecover_0.5.mz")
  calcOutput("CroplandTreecover",
             aggregate = "cluster", cells = cells,
             round = roundArea, outputStatistics = stats, file = paste0("CroplandTreecover_", ctype, ".mz"))
  calcOutput("CroplandTreecover",
             aggregate = FALSE, cells = cells, countryLevel = TRUE,
             round = roundArea, outputStatistics = stats, file = "CroplandTreecover_iso.cs2")

  calcOutput("SNVTargetCropland",
             aggregate = FALSE, cells = cells,
             round = NULL, outputStatistics = stats, file = "SNVTargetCropland_0.5.mz")
  calcOutput("SNVTargetCropland",
             aggregate = "cluster", cells = cells,
             round = roundArea, outputStatistics = stats, file = paste0("SNVTargetCropland_", ctype, ".mz"))

  # 31 past
  if (grepl("+grasslandRealization", dev)) {
    calcOutput("GrasslandBiomass",  round = 3, outputStatistics = stats, file = "f31_grass_bio_hist.cs3",
               cells = cells, aggregate = "region")
    calcOutput("LUH3", aggregate = "cluster", landuseTypes = "LUH3",
               cellular = TRUE,
               outputStatistics = stats, file = paste0("f31_LUH3_", ctype, ".mz"))
    # hard coded climate scenario for harmonization of data
    calcOutput("GrasslandsYields", lpjml = lpjml[["grass"]], climatetype = "MRI-ESM2-0:ssp126",
               subtype = "/co2/Nreturn0p5", # nolint
               lsu_levels = c(seq(0, 2.2, 0.2), 2.5), past_mngmt = "mdef",
               outputStatistics = stats, file = paste0("f31_grassl_yld_", ctype, ".mz"), years = magYears,
               cells = cells, aggregate = "cluster")
    calcOutput("GrasslandsYields", lpjml = lpjml[["grass"]], climatetype = "MRI-ESM2-0:ssp126",
               subtype = "/co2/Nreturn0p5", # nolint
               lsu_levels = c(seq(0, 2.2, 0.2), 2.5), past_mngmt = "mdef",
               outputStatistics = stats, file = paste0("f31_grassl_yld.mz"), years = magYears,
               cells = cells, aggregate = FALSE)
  }

  calcOutput("ClimateClass", aggregate = "cluster", datasource = "koeppen", cells = cells,
             years = "y2001", outputStatistics = stats,
             file = paste0("koeppen_geiger_", ctype, ".mz"))  # years available: 1951, 1976, 2001
  calcOutput("ClimateClass", aggregate = "cluster", datasource = "ipcc", cells = cells,
             outputStatistics = stats, file = paste0("ipcc_climate_zones_", ctype, ".mz"))
  calcOutput("CellCountryFraction", aggregate = "cluster", cells = cells,
             outputStatistics = stats, file = paste0("cell_country_fraction_", ctype, ".mz"))

  # 32 forestry
  calcOutput("AfforestationMask", subtype = "noboreal", cells = cells, aggregate = "cluster", round = roundArea,
             outputStatistics = stats, file = paste0("aff_noboreal_", ctype, ".mz"))
  calcOutput("AfforestationMask", subtype = "onlytropical", cells = cells, aggregate = "cluster", round = roundArea,
             outputStatistics = stats, file = paste0("aff_onlytropical_", ctype, ".mz"))
  calcOutput("AfforestationMask", subtype = "unrestricted", cells = cells, aggregate = "cluster", round = roundArea,
             outputStatistics = stats, file = paste0("aff_unrestricted_", ctype, ".mz"))

  calcOutput("NpiNdcAdAolcPol", aggregate = "cluster", cells = cells,
             round = roundArea, outputStatistics = stats, file = paste0("npi_ndc_ad_aolc_pol_", ctype, ".mz"))
  calcOutput("NpiNdcAffPol",    aggregate = "cluster", cells = cells,
             round = roundArea, outputStatistics = stats, file = paste0("npi_ndc_aff_pol_", ctype, ".mz"))

  calcOutput("BphEffect", aggregate = "cluster", cells = cells,
             outputStatistics = stats, file = paste0("f32_bph_effect_noTCRE_", ctype, ".mz"))
  calcOutput("BphTCRE",   aggregate = "cluster", cells = cells,
             outputStatistics = stats, file = paste0("f32_localTCRE_", ctype, ".mz"))
  calcOutput("BphMask",   aggregate = "cluster", cells = cells,
             outputStatistics = stats, file = paste0("f32_bph_mask_", ctype, ".mz"))

  # 34 urban land
  if (dev == "+GaoUrbanLand") {

    calcOutput("UrbanLandFuture", subtype = "Gao",
               aggregate = FALSE, cells = cells,
               round = NULL, years = shortYears,
               outputStatistics = stats, file = "f34_urbanland_0.5.mz")
    calcOutput("UrbanLandFuture", subtype = "Gao",
               aggregate = "cluster", cells = cells,
               round = roundArea, years = shortYears,
               outputStatistics = stats, file = paste0("f34_urbanland_", ctype, ".mz"))
  } else {

    calcOutput("UrbanLandFuture", subtype = "LUH3",
               aggregate = FALSE, cells = cells,
               round = NULL, years = shortYears,
               outputStatistics = stats, file = "f34_urbanland_0.5.mz")
    calcOutput("UrbanLandFuture", subtype = "LUH3",
               aggregate = "cluster", cells = cells,
               round = roundArea, years = shortYears,
               outputStatistics = stats, file = paste0("f34_urbanland_", ctype, ".mz"))
  }

  # 35 natveg
  calcOutput("AgeClassDistribution", round = 6,
             aggregate = "cluster", cells = cells,
             outputStatistics = stats, file = paste0("forestageclasses_", ctype, ".mz"))

  calcOutput("PotentialForestArea",
             refData = "lpj", cells = cells, lpjml = lpjml, climatetype = climatetype, years = lpjYears,
             aggregate = FALSE, round = NULL, outputStatistics = stats, file = "pot_forest_area_0.5.mz")

  calcOutput("PotentialForestArea",
             refData = "lpj", cells = cells, lpjml = lpjml, climatetype = climatetype, years = lpjYears,
             aggregate = "cluster", round = roundArea,
             outputStatistics = stats, file = paste0("pot_forest_area_", ctype, ".mz"))

  # 37 labour prod
  calcOutput("LabourProdImpactEmu", aggregate = "cluster", cells = cells, subtype = "impact",
             round = 6, outputStatistics = stats, file = paste0("f37_labourprodimpact_", ctype, ".mz"))
  calcOutput("LabourProdImpactEmu", aggregate = "cluster", cells = cells, subtype = "relief",
             round = 6, outputStatistics = stats, file = paste0("f37_labourprodrelief_", ctype, ".mz"))

  # 40
  calcOutput("TransportTime", aggregate = "cluster", cells = cells,
             round = 6, outputStatistics = stats, file = paste0("transport_distance_", ctype, ".mz"))
  calcOutput("TransportTime", aggregate = FALSE, cells = cells,
             round = 6, outputStatistics = stats, file = "transport_distance.mz")
  calcOutput("TransportCosts", aggregate = "GLO", round = 4, outputStatistics = stats, file = "f40_transport_costs.csv")

  # 41 area equipped for irrigation
  calcOutput("AreaEquippedForIrrigation",
             aggregate = "cluster", cellular = TRUE,
             selectyears = magYearsPastLong, round = roundArea,
             outputStatistics = stats, file = paste0("avl_irrig_", ctype, ".mz"))

  # 42 water demand
  calcOutput("Irrigation", lpjml = lpjml, years = lpjYears, climatetype = climatetype,
             cells = cells, aggregate = "cluster", round = 6,
             outputStatistics = stats, file = paste0("lpj_airrig_", ctype, ".mz"))
  if (grepl("+griddedL2Mcomp", dev)) {
    calcOutput("Irrigation", lpjml = lpjml, years = lpjYears, climatetype = climatetype,
               aggregate = FALSE, round = NULL,
               outputStatistics = stats, file = "lpj_airrig_0.5.mz")
  }

  # dummy Growing Period
  calcOutput("GrowingPeriod", lpjml = lpjml, years = lpjYears,
             climatetype = climatetype, yield_ratio = 0.1,
             aggregate = FALSE, cells = cells,
             round = 2, outputStatistics = stats, file = "lpj_grper_0.5.mz")

  # 43 water availability
  calcOutput("AvlWater", lpjml = lpjml, years = lpjYears,
             climatetype = climatetype, seasonality = "grper",
             aggregate = "cluster",  cells = cells,
             round = 6, outputStatistics = stats, file = paste0("lpj_watavail_grper_", ctype, ".mz"))
  calcOutput("AvlWater", lpjml = lpjml, years = lpjYears,
             climatetype = climatetype, seasonality = "total",
             aggregate = "cluster", cells = cells,
             round = 6, outputStatistics = stats, file = paste0("lpj_watavail_total_", ctype, ".mz"))

  calcOutput("EFRSmakthin", lpjml = lpjml, years = lpjYears, climatetype = climatetype,
             aggregate = "cluster", cells = cells,
             round = 6, seasonality = "grper",
             outputStatistics = stats, file = paste0("lpj_envflow_grper_", ctype, ".mz"))
  calcOutput("EFRSmakthin", lpjml = lpjml, years = lpjYears, climatetype = climatetype,
             aggregate = "cluster", cells = cells,
             round = 6, seasonality = "total",
             outputStatistics = stats, file = paste0("lpj_envflow_total_", ctype, ".mz"))

  if (dev == "EFRtest") {
    calcOutput("EnvmtlFlow", lpjml = lpjml, years = lpjYears, climatetype = climatetype,
               aggregate = "cluster",
               round = 6, seasonality = "grper",
               outputStatistics = stats, file = paste0("envflow_grper_", ctype, ".cs3"))
    calcOutput("EnvmtlFlow", lpjml = lpjml, years = lpjYears, climatetype = climatetype,
               aggregate = "cluster",
               round = 6, seasonality = "total",
               outputStatistics = stats, file = paste0("envflow_total_", ctype, ".cs3"))
  }

  calcOutput("WaterUseNonAg", datasource = "WATERGAP_ISIMIP", usetype = "all:all",
             selectyears = lpjYears, seasonality = "grper", lpjml = lpjml, climatetype = climatetype,
             aggregate = "cluster", cells = cells,
             outputStatistics = stats, file = paste0("watdem_nonagr_grper_", ctype, ".mz"))

  calcOutput("WaterUseNonAg", datasource = "WATERGAP_ISIMIP", usetype = "all:all",
             selectyears = lpjYears, seasonality = "total", lpjml = lpjml, climatetype = climatetype,
             aggregate = "cluster", cells = cells,
             outputStatistics = stats, file = paste0("watdem_nonagr_total_", ctype, ".mz"))

  # 44 biodiversity
  calcOutput("BiomeType", aggregate = "cluster", cells = cells, round = roundArea,
             outputStatistics = stats, file = paste0("biorealm_biome_", ctype, ".mz"))
  calcOutput("Luh2SideLayers", aggregate = "cluster", cells = cells,
             round = roundArea, outputStatistics = stats, file = paste0("luh2_side_layers_", ctype, ".mz"))
  calcOutput("Luh2SideLayers", aggregate = FALSE, cells = cells,
             round = NULL, outputStatistics = stats, file = "luh2_side_layers_0.5.mz")
  calcOutput("RRLayer", aggregate = "cluster", cells = cells,
             round = roundArea, outputStatistics = stats, file = paste0("rr_layer_", ctype, ".mz"))

  # 50 nitrogen
  calcOutput("AtmosphericDepositionRates", cellular = TRUE, aggregate = FALSE, round = NULL, cells = cells,
             outputStatistics = stats, file = "f50_AtmosphericDepositionRates_0.5.mz")
  calcOutput("NitrogenFixationRateNatural", aggregate = FALSE, round = NULL, cells = cells,
             outputStatistics = stats, file = "f50_NitrogenFixationRateNatural_0.5.mz")

  calcOutput("AtmosphericDepositionRates", cellular = TRUE, aggregate = "cluster", round = 6, cells = cells,
             outputStatistics = stats, file = paste0("f50_AtmosphericDepositionRates_", ctype, ".mz"))
  calcOutput("NitrogenFixationRateNatural", cells = cells, aggregate = "cluster", round = 6,
             outputStatistics = stats, file = paste0("f50_NitrogenFixationRateNatural_", ctype, ".mz"))

  calcOutput("SchulteUebbing", aggregate = FALSE, round = NULL,
             outputStatistics = stats, file = "criticalNitrogenSurplus_0.5.mz")

  # 52 carbon
  calcOutput("Carbon", aggregate = FALSE, lpjml = lpjml, climatetype = climatetype, cells = cells,
             round = NULL, years = lpjYears, outputStatistics = stats, file = "lpj_carbon_stocks_0.5.mz")
  calcOutput("TopsoilCarbon", aggregate = FALSE, lpjml = lpjml, climatetype = climatetype, cells = cells,
             round = NULL, years = lpjYears, outputStatistics = stats, file = "lpj_carbon_topsoil_0.5.mz")

  calcOutput("Carbon", aggregate = "cluster", lpjml = lpjml, climatetype = climatetype,  cells = cells,
             round = 6, years = lpjYears, outputStatistics = stats, file = paste0("lpj_carbon_stocks_", ctype, ".mz"))
  calcOutput("TopsoilCarbon", aggregate = "cluster", lpjml = lpjml, climatetype = climatetype, cells = cells,
             round = 6, years = lpjYears, outputStatistics = stats, file = paste0("lpj_carbon_topsoil_", ctype, ".mz"))

  # 58 peatland
  calcOutput("Peatland", subtype = "degraded", cells = cells, aggregate = FALSE,
             round = NULL, outputStatistics = stats, file = "f58_peatland_degrad_0.5.mz")
  calcOutput("Peatland", subtype = "intact",   cells = cells, aggregate = FALSE,
             round = NULL, outputStatistics = stats, file = "f58_peatland_intact_0.5.mz")
  calcOutput("Peatland", subtype = "degraded", cells = cells, aggregate = "cluster", round = roundArea,
             outputStatistics = stats, file = paste0("f58_peatland_degrad_", ctype, ".mz"))
  calcOutput("Peatland", subtype = "intact",   cells = cells, aggregate = "cluster", round = roundArea,
             outputStatistics = stats, file = paste0("f58_peatland_intact_", ctype, ".mz"))

  calcOutput("Peatland2", aggregate = FALSE, cells = cells, round = NULL,
             outputStatistics = stats, file = "f58_peatland_area_0.5.mz")
  calcOutput("Peatland2", aggregate = "cluster", cells = cells, round = roundArea,
             outputStatistics = stats, file = paste0("f58_peatland_area_", ctype, ".mz"))
  calcOutput("Peatland2", aggregate = FALSE, cells = cells, round = roundArea,
             countryLevel = TRUE, outputStatistics = stats, file = "f58_peatland_area_iso.cs3")

  # 59 som
  calcOutput("SOMinitialsiationPools", aggregate = "cluster", round = 6, cells = cells,
             outputStatistics = stats, file = paste0("f59_som_initialisation_pools_", ctype, ".mz"))
  calcOutput("SOCLossShare", aggregate = "cluster", rate = "loss", round = 6, cells = cells,
             outputStatistics = stats, file = paste0("cshare_released_", ctype, ".mz"))

  if (grepl("newSOC", dev)) {

    if (grepl("region", dev)) {
      aggregateLevel <- "region"
      cellular <- FALSE
    } else if (grepl("cluster", dev)) {
      aggregateLevel <- "cluster"
      cellular <- TRUE
    } else {
      stop("Undefined aggregation level.")
    }
    histClimatetype <- toolLPJmLVersion(version     = lpjml[["natveg"]],
                                        climatetype = climatetype)$baseline_hist
    calcOutput("LanduseInitialisation", nclasses = "seven",
               aggregate = aggregateLevel, cellular = cellular, cells = cells,
               input_magpie = TRUE, selectyears = "y1990",
               round = roundArea, outputStatistics = stats, file = "f59_land_y1990.cs3")
    calcOutput("CarbonInputMultiplier", aggregate = aggregateLevel,
               round = 6, outputStatistics = stats, file = "f59_cinput_multiplier.cs3")
    calcOutput("CarbonInputMultiplier", inputType = "kcr", aggregate = aggregateLevel,
               round = 6, outputStatistics = stats, file = "f59_cinput_multiplier_residue.cs3")
    calcOutput("CarbonInputMultiplier", inputType = "kli", aggregate = aggregateLevel,
               round = 6, outputStatistics = stats, file = "f59_cinput_multiplier_manure.cs3")
    calcOutput("LitterSoilinput", aggregate = aggregateLevel, years = lpjYears,
               lpjmlNatveg = lpjml[["natveg"]], climatetype = climatetype,
               fixFpc = TRUE, round = 6, outputStatistics = stats, file = "f59_litter_input.cs3")
    calcOutput("DecayFuture", aggregate = aggregateLevel, years = lpjYears,
               lpjmlNatveg = lpjml[["natveg"]], climatetype = climatetype,
               round = 6, outputStatistics = stats, file = "f59_topsoilc_decay.cs3")
    calcOutput("SoilCarbon", aggregate = aggregateLevel, years = "y1990", output = "actualstate",
               lpjmlNatveg = lpjml[["natveg"]], climatetype = histClimatetype,
               round = 6, outputStatistics = stats, file = "f59_topsoilc_actualstate.cs3")
    calcOutput("SoilCarbon", aggregate = aggregateLevel, years = "y1990", output = "naturalstate",
               lpjmlNatveg = lpjml[["natveg"]], climatetype = histClimatetype,
               round = 6, outputStatistics = stats, file = "f59_topsoilc_naturalstate.cs3")
  }

  ##### AGGREGATION ######

  # create info file
  writeInfo <- function(file, lpjmlData, resHigh, resOut, rev, cluster) {
    functioncall <- paste(deparse(sys.call(-3)), collapse = "")

    map <- toolGetMapping(type = "regional", where = "mappingfolder", name = getConfig("regionmapping"))
    regionscode <- regionscode(map)

    info <- c("lpj2magpie settings:",
              paste("* LPJmL data:", lpjmlData),
              paste("* Revision:", rev),
              "", "aggregation settings:",
              paste("* Input resolution:", resHigh),
              paste("* Output resolution:", resOut),
              paste("* Regionscode:", regionscode),
              "* Number of clusters per region:",
              paste(format(names(cluster), width = 5, justify = "right"), collapse = ""),
              paste(format(cluster, width = 5, justify = "right"), collapse = ""),
              paste("* Call:", functioncall))

    base::cat(info, file = file, sep = "\n")
  }
  nrClusterPerRegion <- substr(attributes(p$data)$legend_text, 6,
                               nchar(attributes(p$data)$legend_text) - 1)

  writeInfo(file = "info.txt",
            lpjmlData = climatetype,
            resHigh = "0.5",
            resOut = ctype,
            rev = rev,
            cluster = nrClusterPerRegion)

  mstools::toolWriteMadratLog()

  return(list(tag = versionTag,
              pucTag = sub("^[^_]*_", "", versionTag)))
}
