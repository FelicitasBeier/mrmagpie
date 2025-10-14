#' @title calcGrasslandBiomass
#' @description Calculates pasture biomass demand for the historical period split
#' between rangelands andmanaged pastures.
#' @param cells "magpiecell" for 59199 cells or "lpjcell" for 67420 cells
#' @return Regional biomass demand
#' @author Marcos Alves
#' @seealso \code{\link[madrat]{calcOutput}}, \code{\link[mrcommons]{calcFAOmassbalance}},
#' \code{\link[madrat]{readSource}}
#' @examples
#' \dontrun{
#' calcOutput("GrasslandBiomass")
#' }
#' @importFrom mstools toolCoord2Isocell
#' @importFrom magpiesets findset

calcGrasslandBiomass <- function(cells = "lpjcell") {

  # biomass production
  biomass <- calcOutput("FAOmassbalance", aggregate = FALSE)[, , "production.dm"][, , "pasture"]
  biomass <- collapseNames(biomass)

  land <- calcOutput("LanduseInitialisation", cellular = TRUE, cells = cells, selectyears = seq(1965, 2015, 5),
                     nclasses = "nine", aggregate = FALSE)
  grasslLand <- land[, , c("past", "range")]
  grasslLand <- setNames(grasslLand, c("pastr", "range"))

  # select years
  years <- intersect(getYears(biomass), getYears(land))
  biomass <- biomass[, years, ]
  land <- land[, years, ]
  grasslLand <- grasslLand[, years, ]

  # India, Bangladesh, and Pakistan have livestock numbers unusually high, in comparison to their
  # estimated grazing lands leading to very high grassland biomass demand. A share of this increased
  # demand can be fulfilled by alternative feed sources that are not accounted for in MAgPIE, such as
  # roadside grazing. To reduce demand/area ratio, I'm adding to the areas of managed pastures the category
  # "other land" that is assumed to be used for grazing in those three countries.

  grasslLand["IND", , "pastr"] <- grasslLand["IND", , "pastr"] +
    setNames(dimSums(land["IND", , c("primother", "secdother")], dim = 3), "pastr")
  grasslLand["BGD", , "pastr"] <- grasslLand["BGD", , "pastr"] +
    setNames(dimSums(land["BGD", , c("primother", "secdother")], dim = 3), "pastr")
  grasslLand["PAK", , "pastr"] <- grasslLand["PAK", , "pastr"] +
    setNames(dimSums(land["PAK", , c("primother", "secdother")], dim = 3), "pastr")

  grassYld <- calcOutput("GrasslandsYields",
                         lpjml = "lpjml5p2_pasture",
                         climatetype = paste0("MRI-ESM2-0", ":", "ssp245"),
                         subtype = "/co2/Nreturn0p5", # nolint
                         lsu_levels = c(seq(0, 2.2, 0.2), 2.5), past_mngmt = "mdef", # nolint
                         cells = cells,
                         aggregate = FALSE)[, , "rainfed"]
  grassYld <- collapseNames(grassYld)
  grassYld[grassYld == 0.00088]  <- 0

  potBioMass <- grasslLand * grassYld[, getYears(grasslLand), ]

  potBioMassShare <- setNames(potBioMass[, , "pastr"] / dimSums(potBioMass, dim = 3), "pastr")
  potBioMassShare <- add_columns(potBioMassShare, addnm = "range", dim = 3.1)
  potBioMassShare[, , "range"] <- potBioMass[, , "range"] / dimSums(potBioMass, dim = 3)
  potBioMassShare[is.nan(potBioMassShare)]       <- 0
  potBioMassShare[is.infinite(potBioMassShare)]  <- 1

  livestock <- setNames(readSource("GLW3", subtype = "Aw"),
                        "liv_numb")
  if (cells == "magpiecell") {
    livestock <- toolCoord2Isocell(livestock)
  }
  livestock[livestock < 1] <- 0

  # I am working with the assumption that in most places, the proportional distribution of
  # livestock heads was kept the same between 1965 and 2010. One notable exception is Brasil,
  # which experienced a shift in its agricultural frontier towards the west with increases in
  # rangelands used for cattle production. A "fader" is used to correct the proportional share
  # of livestock being in rangelands and managed pastures, and avoid distortions on the amount
  # of grass biomass production assigned to each system.

  livstSplit <- livestock * potBioMassShare
  livstSplit <- collapseNames(livstSplit)
  if (cells == "lpjcell") {
    livstSplitCtry <- dimSums(livstSplit, dim = c("x", "y"))
  } else if (cells == "magpiecell") {
    livstSplitCtry <- dimSums(livstSplit, dim = 1.2)
  } else {
    stop("Please select cells argument in calcGrasslandBiomass")
  }
  livstShareCtry <- livstSplitCtry[, , "pastr"] / dimSums(livstSplitCtry, dim = 3)
  livstShareCtry[is.nan(livstShareCtry) | is.infinite(livstShareCtry)] <- 0
  livstShareCtry <- add_columns(livstShareCtry, addnm = "range", dim = 3.1)
  livstShareCtry[, , "range"] <- 1 - livstShareCtry[, , "pastr"]

  # I am splitting biomass consumption assuming the share
  # between animals reared on rangelands and pastures correlates linearly
  # with the production of grass in pastures and rangelands in a country. That can be
  # derived by the fact that the feedbaskets assume the same feed ingredients shares
  # within a country.
  commonCtrs <- intersect(getItems(biomass, dim = 1),
                          unique(getItems(livstShareCtry, dim = 1)))
  biomassSplit <- biomass[commonCtrs, , ] * livstShareCtry[commonCtrs, , ]
  biomassSplit <- toolCountryFill(biomassSplit, fill = 0)

  # OBS: we could add a diet correction factor for animals being reared in different system
  # (as to say: even thought the LSUs number might be the same, animals on
  # rangelands might eat less grass than animals on managed pastures)

  return(list(x = biomassSplit,
              weight = NULL,
              isocountries = FALSE,
              unit = "tDM",
              description = "Pasture biomass demand"))
}
