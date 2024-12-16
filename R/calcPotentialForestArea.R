#' @title calcPotentialForestArea
#'
#' @description Calculates the area than can be potentially covered by forests,
#'              based on environmental conditions.
#'
#' @param refData Determines the reference data that the estimated potential
#'                forest area is derived from (currently only "lpj")
#' @param lpjml Defines LPJmL version for crop/grass and natveg specific inputs.
#'              Only relevant, if refData = "lpj".
#' @param climatetype Switch between different GCM climate scenarios.
#'                    Only relevant, if refData = "lpj".
#' @param countryLevel Whether output shall be at country level.
#'                     Requires aggregate=FALSE in calcOutput.
#'
#' @return magpie object in cellular resolution
#' @author Patrick v. Jeetze, Michael Crawford
#'
#' @examples
#' \dontrun{
#' calcOutput("PotentialForestArea", aggregate = FALSE)
#' }

calcPotentialForestArea <- function(refData      = "lpj",
                                    lpjml        = "lpjml5.9.5-m1",
                                    climatetype  = "MRI-ESM2-0:ssp370",
                                    countryLevel = FALSE) {
  if (refData == "lpj") {
    vegc <- calcOutput(
      "LPJmLharmonize",
      lpjmlversion = lpjml,
      climatetype  = climatetype,
      subtype      = "pnv:vegc",
      aggregate    = FALSE
    )

    potForest <- toolConditionalReplace(vegc, c("<20", ">=20"), c(0, 1))

    landArea <- calcOutput("LandArea", aggregate = FALSE)

    potForestArea <- potForest * landArea

  } else {
    stop("calcPotentialForestArea can only use LPJmL data at this time")
  }

  getNames(potForestArea) <- NULL

  if (countryLevel) {
    out <- toolCountryFill(dimSums(potForestArea, dim = c("x", "y")), fill = 0)
  } else {
    out <- potForestArea
  }

  return(
    list(
      x            = out,
      weight       = NULL,
      unit         = "Mha",
      description  = "Potential forest area",
      isocountries = countryLevel
    )
  )
}
