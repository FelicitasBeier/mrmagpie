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
#' \code{\link{setConfig}} (e.g. for setting the mainfolder if not already set properly).
#'
#' @author Kristine Karstens, Jan Philipp Dietrich
#' @seealso
#' \code{\link{readSource}},\code{\link{getCalculations}},\code{\link{calcOutput}},\code{\link{setConfig}}
#' @examples
#' \dontrun{
#' retrieveData("CELLULARMAGPIE", rev = numeric_version("12"),
#'              mainfolder = "pathtowhereallfilesarestored")
#' }
#' @importFrom madrat setConfig getConfig
#' @importFrom magpiesets findset
#' @importFrom digest digest
#' @importFrom stringr str_split
#' @importFrom luplot plotregionscluster
#' @importFrom ggplot2 ggsave
#' @importFrom withr local_options

fullCELLULARMAGPIE <- function(rev = numeric_version("0.1"), dev = "",
                               ctype = "c200",
                               climatetype = "MRI-ESM2-0:ssp370",
                               lpjml = "lpjml5.9.5-m1",
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
  climatescen <- str_split(climatetype, ":")[[1]][2]

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
             cellular = TRUE,  aggregate = "cluster",
             years = magYears, round = 6, file = "f09_pop_grid.cs3")

  calcOutput("GridPop", source = "Gao", subtype = "all", harmonize_until = 2015, urban = TRUE,
             cellular = TRUE,  aggregate = "cluster",
             years = magYears, round = 6, file = "f09_urbanpop_grid.cs3")

  # 10 land
  ## seven land classes
  calcOutput("LanduseInitialisation", nclasses = "seven",
             aggregate = FALSE, cellular = TRUE,
             input_magpie = TRUE, selectyears = magYearsPastLong,
             round = roundArea, file = "avl_land_t_0.5.mz")
  calcOutput("LanduseInitialisation", nclasses = "seven",
             aggregate = "cluster", cellular = TRUE,
             input_magpie = TRUE, selectyears = magYearsPastLong,
             round = roundArea, file = paste0("avl_land_t_", ctype, ".mz"))
  calcOutput("LanduseInitialisation", nclasses = "seven",
             aggregate = FALSE, cellular = FALSE,
             input_magpie = TRUE, selectyears = magYearsPastLong,
             round = roundArea, file = paste0("avl_land_t_iso.cs3"))

  ## nine land classes
  calcOutput("LanduseInitialisation", nclasses = "nine",
             aggregate = FALSE, cellular = TRUE,
             input_magpie = TRUE, selectyears = magYearsPastLong,
             round = roundArea, file = "avl_land_full_t_0.5.mz")
  calcOutput("LanduseInitialisation", nclasses = "nine",
             aggregate = "cluster", cellular = TRUE,
             input_magpie = TRUE, selectyears = magYearsPastLong,
             round = roundArea, file = paste0("avl_land_full_t_", ctype, ".mz"))
  calcOutput("LanduseInitialisation", nclasses = "nine",
             aggregate = FALSE, cellular = FALSE,
             input_magpie = TRUE, selectyears = magYearsPastLong,
             round = roundArea, file = paste0("avl_land_full_t_iso.cs3"))

  calcOutput("AvlLandSi",  aggregate = FALSE,
             round = roundArea, file = "avl_land_si_0.5.mz")
  calcOutput("AvlLandSi",  aggregate = "cluster",
             round = roundArea, file = paste0("avl_land_si_", ctype, ".mz"))

  # 13 TC
  calcOutput("PastrTauHist", round = 2, past_mngmt = "mdef",
             file = "f13_pastr_tau_hist.csv",
             aggregate = superregion)

  # 14 yields
  if (dev == "+calibYield") {

    calcOutput("YieldsCalibrated", aggregate = "cluster",
               datasource = c(lpjml = lpjml, isimip = isimip),
               climatetype = climatetype, round = 2, years = lpjYears,
               file = paste0("lpj_yields_", ctype, ".mz"))

  } else if (grepl("india", dev)) {

    calcOutput("Yields", datasource = c(lpjml = lpjml, isimip = isimip), aggregate = FALSE,
               climatetype = climatetype, round = 2, years = lpjYears, file = paste0("lpj_yields_0.5.mz"),
               weighting = "crop+irrigSpecific", indiaYields = TRUE, scaleFactor = 0.5)


    calcOutput("Yields", aggregate = "cluster",
               datasource = c(lpjml = lpjml, isimip = isimip),
               climatetype = climatetype, round = 2, years = lpjYears, file = paste0("lpj_yields_", ctype, ".mz"),
               weighting = "crop+irrigSpecific", indiaYields = TRUE, scaleFactor = 0.5)

  } else {

    calcOutput("Yields", datasource = c(lpjml = lpjml, isimip = isimip),
               aggregate = FALSE,
               climatetype = climatetype, round = 2, years = lpjYears, file = paste0("lpj_yields_0.5.mz"),
               weighting = ifelse(grepl("YieldWeights_", dev), gsub("YieldWeights_", "", dev), "totalCrop"))

    calcOutput("Yields", aggregate = "cluster",
               datasource = c(lpjml = lpjml, isimip = isimip),
               climatetype = climatetype, round = 2, years = lpjYears, file = paste0("lpj_yields_", ctype, ".mz"),
               weighting = ifelse(grepl("YieldWeights_", dev), gsub("YieldWeights_", "", dev), "totalCrop"))

  }

  calcOutput("DegradationYieldReduction", aggregate = "cluster", round = 6,
             file = paste0("f14_degradation_yld_reduc_", ctype, ".mz"))


  # 22 land conservation
  calcOutput("ProtectedAreaBaseline", nclasses = "seven",
             magpie_input = TRUE,
             aggregate = FALSE, round = roundArea, file = "wdpa_baseline_0.5.mz")
  calcOutput("ProtectedAreaBaseline", nclasses = "seven",
             magpie_input = TRUE,
             aggregate = "cluster", round = roundArea,
             file = paste0("wdpa_baseline_", ctype, ".mz"))

  calcOutput("ConservationPriorities", nclasses = "seven",
             aggregate = FALSE, round = roundArea, file = "consv_prio_areas_0.5.mz")
  calcOutput("ConservationPriorities", nclasses = "seven",
             aggregate = "cluster", round = roundArea, file = paste0("consv_prio_areas_", ctype, ".mz"))

  calcOutput("ProtectArea", bhifl = TRUE,
             aggregate = "cluster", round = roundArea,
             file = paste0("protect_area_", ctype, ".mz"))


  # 30 crop
  calcOutput("Croparea", sectoral = "kcr", physical = TRUE,
             cellular = TRUE,  irrigation = FALSE, round = roundArea,
             aggregate = "cluster", file = paste0("f30_croparea_initialisation_", ctype, ".mz"))
  calcOutput("Croparea", sectoral = "kcr", physical = TRUE,
             cellular = TRUE,  irrigation = TRUE, round = roundArea,
             aggregate = "cluster", file = paste0("f30_croparea_w_initialisation_", ctype, ".mz"))
  ## For cellular comparison
  calcOutput("MAPSPAM", subtype = "physical",  aggregate = FALSE,
             file = paste0("MAPSPAM_croparea_0.5.mz"))
  calcOutput("Croparea", sectoral = "kcr", physical = TRUE, cellular = TRUE,
             irrigation = TRUE, round = roundArea,
             aggregate = FALSE, file = paste0("LUH2_croparea_0.5.mz"))

  calcOutput("AvlCropland", marginal_land = "magpie", cell_upper_bound = 0.9,
             aggregate = FALSE,
             round = roundArea, file = "avl_cropland_0.5.mz")
  calcOutput("AvlCropland", marginal_land = "magpie", cell_upper_bound = 0.9,
             aggregate = "cluster",
             round = roundArea, file = paste0("avl_cropland_", ctype, ".mz"))
  calcOutput("AvlCropland", marginal_land = "magpie", cell_upper_bound = 0.9,
             aggregate = FALSE,  country_level = TRUE,
             round = roundArea, file = paste0("avl_cropland_iso.cs3"))

  calcOutput("CroplandTreecover",
             aggregate = FALSE,
             round = roundArea, file = "CroplandTreecover_0.5.mz")
  calcOutput("CroplandTreecover",
             aggregate = "cluster",
             round = roundArea, file = paste0("CroplandTreecover_", ctype, ".mz"))
  calcOutput("CroplandTreecover",
             aggregate = FALSE,  countryLevel = TRUE,
             round = roundArea, file = "CroplandTreecover_iso.cs2")

  calcOutput("SNVTargetCropland",
             aggregate = FALSE,
             round = roundArea, file = "SNVTargetCropland_0.5.mz")
  calcOutput("SNVTargetCropland",
             aggregate = "cluster",
             round = roundArea, file = paste0("SNVTargetCropland_", ctype, ".mz"))

  # 31 past
  if (grepl("+grasslandRealization", dev)) {
    calcOutput("GrasslandBiomass",  round = 3, file = "f31_grass_bio_hist.cs3",
               aggregate = "region")
    calcOutput("LUH2v2", aggregate = "cluster", landuse_types = "LUH2v2",
               cellular = TRUE,
               file = paste0("f31_LUH2v2_", ctype, ".mz"))
    # hard coded climate scenario for harmonization of data
    calcOutput("GrasslandsYields", lpjml = lpjml[["grass"]], climatetype = "MRI-ESM2-0:ssp126",
               subtype = "/co2/Nreturn0p5", # nolint
               lsu_levels = c(seq(0, 2.2, 0.2), 2.5), past_mngmt = "mdef",
               file = paste0("f31_grassl_yld_", ctype, ".mz"), years = magYears,
               aggregate = "cluster")
    calcOutput("GrasslandsYields", lpjml = lpjml[["grass"]], climatetype = "MRI-ESM2-0:ssp126",
               subtype = "/co2/Nreturn0p5", # nolint
               lsu_levels = c(seq(0, 2.2, 0.2), 2.5), past_mngmt = "mdef",
               file = paste0("f31_grassl_yld.mz"), years = magYears,
               aggregate = FALSE)
    calcOutput("MaxPastureSuit", climatetype = climatetype, lpjml =  lpjml[["natveg"]],
               file = paste0("f31_max_managed_pasture_", ctype, ".mz"), years = magYears, aggregate = "cluster")
    calcOutput("MaxPastureSuit", climatetype = climatetype, lpjml =  lpjml[["natveg"]],
               file = "f31_max_managed_pasture.mz", years = magYears, aggregate = FALSE)
  }

  if (grepl("+PastrMngtLevels", dev)) {
    calcOutput("PastrMngtLevels", climatetype = paste0("MRI-ESM2-0", ":", climatescen),
               options = c("brazil_1", "brazil_2", "brazil_4"), cost_level = c(1, 2, 3),
               file = "PastrMngtLevels.mz", aggregate = FALSE)
  }

  calcOutput("ClimateClass", aggregate = "cluster", datasource = "koeppen",  years = "y2001",
             file = paste0("koeppen_geiger_", ctype, ".mz"))          # years available: 1951, 1976, 2001
  calcOutput("ClimateClass", aggregate = "cluster", datasource = "ipcc",
             file = paste0("ipcc_climate_zones_", ctype, ".mz"))
  calcOutput("CellCountryFraction", aggregate = "cluster",
             file = paste0("cell_country_fraction_", ctype, ".mz"))

  # 32 forestry
  calcOutput("AfforestationMask", subtype = "noboreal",  aggregate = "cluster", round = roundArea,
             file = paste0("aff_noboreal_", ctype, ".mz"))
  calcOutput("AfforestationMask", subtype = "onlytropical",  aggregate = "cluster", round = roundArea,
             file = paste0("aff_onlytropical_", ctype, ".mz"))
  calcOutput("AfforestationMask", subtype = "unrestricted",  aggregate = "cluster", round = roundArea,
             file = paste0("aff_unrestricted_", ctype, ".mz"))

  calcOutput("NpiNdcAdAolcPol", aggregate = "cluster",
             round = roundArea, file = paste0("npi_ndc_ad_aolc_pol_", ctype, ".mz"))
  calcOutput("NpiNdcAffPol",    aggregate = "cluster",
             round = roundArea, file = paste0("npi_ndc_aff_pol_", ctype, ".mz"))

  calcOutput("BphEffect", aggregate = "cluster",
             file = paste0("f32_bph_effect_noTCRE_", ctype, ".mz"))
  calcOutput("BphTCRE",   aggregate = "cluster",
             file = paste0("f32_localTCRE_", ctype, ".mz"))
  calcOutput("BphMask",   aggregate = "cluster",
             file = paste0("f32_bph_mask_", ctype, ".mz"))

  # 34 urban land
  if (dev == "+GaoUrbanLand") {

    calcOutput("UrbanLandFuture", subtype = "Gao",
               aggregate = FALSE,
               round = roundArea, years = shortYears,
               file = "f34_urbanland_0.5.mz")
    calcOutput("UrbanLandFuture", subtype = "Gao",
               aggregate = "cluster",
               round = roundArea, years = shortYears,
               file = paste0("f34_urbanland_", ctype, ".mz"))
  } else {

    calcOutput("UrbanLandFuture", subtype = "LUH2v2",
               aggregate = FALSE,
               round = roundArea, years = shortYears,
               file = "f34_urbanland_0.5.mz")
    calcOutput("UrbanLandFuture", subtype = "LUH2v2",
               aggregate = "cluster",
               round = roundArea, years = shortYears,
               file = paste0("f34_urbanland_", ctype, ".mz"))
  }

  # 35 natveg
  calcOutput("AgeClassDistribution", round = 6,
             aggregate = "cluster",
             file = paste0("forestageclasses_", ctype, ".mz"))

  calcOutput("PotentialForestArea",
             refData = "lpj", lpjml = lpjml, climatetype = climatetype, years = lpjYears,
             aggregate = FALSE, round = roundArea, file = "pot_forest_area_0.5.mz")

  calcOutput("PotentialForestArea",
             refData = "lpj", lpjml = lpjml, climatetype = climatetype, years = lpjYears,
             aggregate = "cluster", round = roundArea, file = paste0("pot_forest_area_", ctype, ".mz"))

  # 37 labour prod
  calcOutput("LabourProdImpactEmu", aggregate = "cluster",  subtype = "impact",
             round = 6, file = paste0("f37_labourprodimpact_", ctype, ".mz"))
  calcOutput("LabourProdImpactEmu", aggregate = "cluster",  subtype = "relief",
             round = 6, file = paste0("f37_labourprodrelief_", ctype, ".mz"))

  # 40
  calcOutput("TransportTime", aggregate = "cluster",
             round = 6, file = paste0("transport_distance_", ctype, ".mz"))
  calcOutput("TransportTime", aggregate = FALSE,
             round = 6, file = "transport_distance.mz")
  calcOutput("TransportCosts", aggregate = "GLO", round = 4, file = "f40_transport_costs.csv")
  calcOutput("TransportCosts", transport = "nonlocal", aggregate = "GLO", round = 4,
             file = "f40_transport_costs_nonlocal.csv")

  # 41 area equipped for irrigation
  calcOutput("AreaEquippedForIrrigation",
             aggregate = "cluster", cellular = TRUE,
             selectyears = magYearsPastLong, round = roundArea,
             file = paste0("avl_irrig_", ctype, ".mz"))

  # 42 water demand
  calcOutput("Irrigation", lpjml = lpjml, years = lpjYears, climatetype = climatetype,
             aggregate = "cluster", round = 6,
             file = paste0("lpj_airrig_", ctype, ".mz"))

  # dummy Growing Period
  calcOutput("GrowingPeriod", lpjml = lpjml, years = lpjYears,
             climatetype = climatetype, yield_ratio = 0.1,
             aggregate = FALSE,
             round = 2, file = "lpj_grper_0.5.mz")

  # 43 water availability
  calcOutput("AvlWater", lpjml = lpjml, years = lpjYears,
             climatetype = climatetype, seasonality = "grper",
             aggregate = "cluster",
             round = 6, file = paste0("lpj_watavail_grper_", ctype, ".mz"))
  calcOutput("AvlWater", lpjml = lpjml, years = lpjYears,
             climatetype = climatetype, seasonality = "total",
             aggregate = "cluster",
             round = 6, file = paste0("lpj_watavail_total_", ctype, ".mz"))

  calcOutput("EFRSmakthin", lpjml = lpjml, years = lpjYears, climatetype = climatetype,
             aggregate = "cluster",
             round = 6, seasonality = "grper",
             file = paste0("lpj_envflow_grper_", ctype, ".mz"))
  calcOutput("EFRSmakthin", lpjml = lpjml, years = lpjYears, climatetype = climatetype,
             aggregate = "cluster",
             round = 6, seasonality = "total",
             file = paste0("lpj_envflow_total_", ctype, ".mz"))

  if (dev == "EFRtest") {
    calcOutput("EnvmtlFlow", lpjml = lpjml, years = lpjYears, climatetype = climatetype,
               aggregate = "cluster",
               round = 6, seasonality = "grper",
               file = paste0("envflow_grper_", ctype, ".cs3"))
    calcOutput("EnvmtlFlow", lpjml = lpjml, years = lpjYears, climatetype = climatetype,
               aggregate = "cluster",
               round = 6, seasonality = "total",
               file = paste0("envflow_total_", ctype, ".cs3"))
  }

  calcOutput("WaterUseNonAg", datasource = "WATERGAP_ISIMIP", usetype = "all:all",
             selectyears = lpjYears, seasonality = "grper", lpjml = lpjml, climatetype = climatetype,
             aggregate = "cluster",
             file = paste0("watdem_nonagr_grper_", ctype, ".mz"))

  calcOutput("WaterUseNonAg", datasource = "WATERGAP_ISIMIP", usetype = "all:all",
             selectyears = lpjYears, seasonality = "total", lpjml = lpjml, climatetype = climatetype,
             aggregate = "cluster",
             file = paste0("watdem_nonagr_total_", ctype, ".mz"))

  # 44 biodiversity
  calcOutput("BiomeType", aggregate = "cluster",  round = roundArea,
             file = paste0("biorealm_biome_", ctype, ".mz"))
  calcOutput("Luh2SideLayers", aggregate = "cluster",
             round = roundArea, file = paste0("luh2_side_layers_", ctype, ".mz"))
  calcOutput("Luh2SideLayers", aggregate = FALSE,
             round = roundArea, file = "luh2_side_layers_0.5.mz")
  calcOutput("RRLayer", aggregate = "cluster",
             round = roundArea, file = paste0("rr_layer_", ctype, ".mz"))

  # 50 nitrogen
  calcOutput("AtmosphericDepositionRates", cellular = TRUE, aggregate = FALSE, round = 6,
             file = "f50_AtmosphericDepositionRates_0.5.mz")
  calcOutput("NitrogenFixationRateNatural", aggregate = FALSE, round = 6,
             file = "f50_NitrogenFixationRateNatural_0.5.mz")

  calcOutput("AtmosphericDepositionRates", cellular = TRUE, aggregate = "cluster", round = 6,
             file = paste0("f50_AtmosphericDepositionRates_", ctype, ".mz"))
  calcOutput("NitrogenFixationRateNatural",  aggregate = "cluster", round = 6,
             file = paste0("f50_NitrogenFixationRateNatural_", ctype, ".mz"))
  calcOutput("SchulteUebbing", aggregate = FALSE, file = "criticalNitrogenSurplus_0.5.mz")

  # 52 carbon
  calcOutput("Carbon", aggregate = FALSE, lpjml = lpjml, climatetype = climatetype,
             round = 6, years = "y1995", file = "lpj_carbon_stocks_0.5.mz")
  calcOutput("TopsoilCarbon", aggregate = FALSE, lpjml = lpjml, climatetype = climatetype,
             round = 6, years = "y1995", file = "lpj_carbon_topsoil_0.5.mz")

  calcOutput("Carbon", aggregate = "cluster", lpjml = lpjml, climatetype = climatetype,
             round = 6, years = lpjYears, file = paste0("lpj_carbon_stocks_", ctype, ".mz"))
  calcOutput("TopsoilCarbon", aggregate = "cluster", lpjml = lpjml, climatetype = climatetype,
             round = 6, years = lpjYears, file = paste0("lpj_carbon_topsoil_", ctype, ".mz"))

  # 58 peatland
  calcOutput("Peatland", subtype = "degraded",  aggregate = FALSE,
             round = roundArea, file = "f58_peatland_degrad_0.5.mz")
  calcOutput("Peatland", subtype = "intact",    aggregate = FALSE,
             round = roundArea, file = "f58_peatland_intact_0.5.mz")
  calcOutput("Peatland", subtype = "degraded",  aggregate = "cluster", round = roundArea,
             file = paste0("f58_peatland_degrad_", ctype, ".mz"))
  calcOutput("Peatland", subtype = "intact",    aggregate = "cluster", round = roundArea,
             file = paste0("f58_peatland_intact_", ctype, ".mz"))

  calcOutput("Peatland2", aggregate = FALSE,  round = roundArea,
             file = "f58_peatland_area_0.5.mz")
  calcOutput("Peatland2", aggregate = "cluster",  round = roundArea,
             file = paste0("f58_peatland_area_", ctype, ".mz"))
  calcOutput("Peatland2", aggregate = FALSE,  round = roundArea,
             countryLevel = TRUE, file = "f58_peatland_area_iso.cs3")

  # 59 som
  calcOutput("SOMinitialsiationPools", aggregate = "cluster", round = 6,
             file = paste0("f59_som_initialisation_pools_", ctype, ".mz"))
  calcOutput("SOCLossShare", aggregate = "cluster", rate = "loss", round = 6,
             file = paste0("cshare_released_", ctype, ".mz"))

  if (grepl("newSOC", dev)) {

    histClimatetype <- toolLPJmLVersion(version     = lpjml[["natveg"]],
                                        climatetype = climatetype)$baseline_hist
    calcOutput("LanduseInitialisation", nclasses = "seven",
               aggregate = "region", cellular = FALSE,
               input_magpie = TRUE, selectyears = "y1990",
               round = roundArea, file = "f59_land_y1990.cs3")
    calcOutput("CarbonInputMultiplier", aggregate = "region",
               round = 6, file = "f59_cinput_multiplier.cs3")
    calcOutput("CarbonInputMultiplier", inputType = "kcr", aggregate = "region",
               round = 6, file = "f59_cinput_multiplier_residue.cs3")
    calcOutput("CarbonInputMultiplier", inputType = "kli", aggregate = "region",
               round = 6, file = "f59_cinput_multiplier_manure.cs3")
    calcOutput("LitterSoilinput", aggregate = "region", years = lpjYears,
               lpjmlNatveg = lpjml[["natveg"]], climatetype = climatetype,
               fixFpc = TRUE, round = 6, file = "f59_litter_input.cs3")
    calcOutput("DecayFuture", aggregate = "region", years = lpjYears,
               lpjmlNatveg = lpjml[["natveg"]], climatetype = climatetype,
               round = 6, file = "f59_topsoilc_decay.cs3")
    calcOutput("SoilCarbon", aggregate = "region", years = "y1990", output = "actualstate",
               lpjmlNatveg = lpjml[["natveg"]], climatetype = histClimatetype,
               round = 6, file = "f59_topsoilc_actualstate.cs3")
    calcOutput("SoilCarbon", aggregate = "region", years = "y1990", output = "naturalstate",
               lpjmlNatveg = lpjml[["natveg"]], climatetype = histClimatetype,
               round = 6, file = "f59_topsoilc_naturalstate.cs3")
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

  return(list(tag = versionTag,
              pucTag = sub("^[^_]*_", "", versionTag)))
}
