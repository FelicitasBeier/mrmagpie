#' @title calcPotentialForestArea
#'
#' @description Calculates the area than can be potentially covered by forests,
#'              based on environmental conditions.
#'
#' @param refData  Determines the reference data that the estimated potential
#'                 forest area is derived from (currently only "lpj")
#' @param cells            magpiecell (59199 cells) or lpjcell (67420 cells)
#' @param countryLevel    Whether output shall be at country level.
#'                        Requires aggregate=FALSE in calcOutput.
#' @param lpjml Defines LPJmL version for crop/grass and natveg specific inputs.
#'              Only relevant, if refData = "lpj".
#' @param climatetype Switch between different GCM climate scenarios.
#'                    Only relevant, if refData = "lpj".
#' @param grassyCorrection If TRUE, the potential forest area is reduced by the grassland-ecoregion
#'                    cover fraction (RESOLVE 2017 biomes 7-10, calcGrassyEcoregions), correcting the
#'                    LPJmL overestimation of forest in open grassy ecosystems.
#'
#' @return magpie object in cellular resolution
#' @author Patrick v. Jeetze, Florian Humpenoeder
#'
#' @examples
#' \dontrun{
#' calcOutput("PotentialForestArea", aggregate = FALSE)
#' }
#'
#' @importFrom madrat readSource calcOutput toolCountryFill
#' @importFrom magclass dimSums getYears setYears
#' @importFrom mstools toolCoord2Isocell

calcPotentialForestArea <- function(refData = "lpj", countryLevel = FALSE, cells = "lpjcell",
                                    lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de"),
                                    climatetype = "MRI-ESM2-0:ssp370", grassyCorrection = TRUE) {
  if (refData == "lpj") {
    vegc <- calcOutput("LPJmL_new",
      version = lpjml["natveg"],
      climatetype = climatetype,
      subtype = "vegc", stage = "harmonized2020",
      aggregate = FALSE
    )

    potForest <- toolConditionalReplace(vegc, c("<20", ">=20"), c(0, 1))

    landIni <- calcOutput("LanduseInitialisation",
      aggregate = FALSE, cellular = TRUE, nclasses = "seven",
      input_magpie = TRUE, cells = "lpjcell", years = "y1995", round = 6
    )
    landArea <- dimSums(landIni, dim = 3)

    potForestArea <- potForest * landArea

    if (grassyCorrection) {
      # reduce the potential forest area by the grassland-ecoregion cover fraction, where LPJmL
      # overestimates forest in open grassy ecosystems
      grassy <- calcOutput("GrassyEcoregions", aggregate = FALSE)
      potForestArea <- potForestArea * (1 - setYears(grassy, NULL))
    }

  }

  getNames(potForestArea) <- NULL

  if (countryLevel) {
    out <- toolCountryFill(dimSums(potForestArea, dim = c("x", "y")), fill = 0)
  } else {
    if (cells == "magpiecell") {
      out <- toolCoord2Isocell(potForestArea)
    } else if (cells == "lpjcell") {
      out <- potForestArea
    } else {
      stop("Please specify cells argument")
    }
  }

  return(list(
    x = out,
    weight = NULL,
    unit = "Mha",
    description = "Potential forest area",
    isocountries = countryLevel
  ))
}
