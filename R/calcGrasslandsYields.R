#' @title calcGrasslandsYields
#' @description Calculates rangelands maximum output and managed pastures yields
#' @param past_mngmt pasture areas management option
#' @param lsu_levels Livestock unit levels in the source folder
#' @param subtype Switch between different climate scenarios
#' @param lpjml Defines LPJmL version for crop/grass and natveg specific inputs
#' @param climatetype Global Circulation Model to be used
#' @param cells       "magpiecell" for 59199 cells or "lpjcell" for 67420 cells
#' @return magpie object in cellular resolution
#' @author Marcos Alves
#'
#' @examples
#' \dontrun{
#' calcOutput("GrasslandsYields")
#' }
calcGrasslandsYields <- function(lpjml = "lpjml5p2_pasture", climatetype = "MRI-ESM2-0:ssp370",
                                 cells = "lpjcell",
                                 subtype = "/co2/Nreturn0p5", # nolint
                                 lsu_levels = c(seq(0, 2, 0.2), 2.5), past_mngmt = "mdef") { # nolint

  gCm2yTotDMy <- (10000 * 2.21 / 1e6)
  x <- calcOutput("RangelandsMaxNew", lsuLevels = lsu_levels,
                  lpjml = lpjml, climatetype = climatetype,
                  scenario = paste0(subtype, "/limN"),
                  report = "harvest", aggregate = FALSE)

  if (past_mngmt == "mdef") {
    n <- "/unlimN"
  } else if (past_mngmt == "me2") {
    n <- "/limN"
  } else {
    stop("past_mngmt not available yet")
  }

  y <- calcOutput("Pastr_new", past_mngmt = past_mngmt,
                  lpjml = lpjml, climatetype = climatetype,
                  scenario = paste0(subtype, n), aggregate = FALSE)
  invalid <- (y - x) < 0
  invalid <- collapseNames(invalid)
  # substituting pastr yields that are smaller than rangelands by the value of rangeland yields
  y[invalid]  <- x[invalid]
  pasture <- mbind(x, y)
  pasture <- mstools::toolHoldConstantBeyondEnd(pasture)
  x <- pasture * gCm2yTotDMy

  # reduce to old grid cell format
  if (cells == "magpiecell") {
    x <- toolCoord2Isocell(x)
  }
  # Land area as weights
  landArea <- calcOutput("LandArea", cells = cells,
                         aggregate = FALSE)

  return(
    list(
      x = x,
      weight = landArea,
      unit = "t/DM/y",
      description = paste0("Maximum grasslands yields obtained with rangelands ",
                           "and managed pastures yields for",
                           past_mngmt),
      isocountries = FALSE
    )
  )
}
