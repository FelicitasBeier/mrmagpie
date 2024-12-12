#' @title calcTopsoilCarbon
#' @description This function extracts topsoil carbon densities from LPJ to MAgPIE
#'
#' @param lpjml       Defines LPJmL version for crop/grass and natveg specific inputs
#' @param climatetype Switch between different GCM climate scenarios
#'
#' @return magpie object in cellular resolution
#' @author Kristine Karstens
#'
#' @examples
#' \dontrun{
#' calcOutput("TopsoilCarbon", aggregate = FALSE)
#' }
#'
#' @importFrom magpiesets findset
#' @importFrom mstools toolCoord2Isocell

calcTopsoilCarbon <- function(lpjml       = "lpjml5.9.5-m1",
                              climatetype = "MRI-ESM2-0:ssp370") {

  soilcLayerNatveg <- calcOutput("LPJmLharmonize",
                                 lpjmlversion = lpjml,
                                 climatetype  = climatetype,
                                 subtype      = "grass:soilc_layer",
                                 aggregate    = FALSE)

  topsoilc           <- soilcLayerNatveg[, , 1] + 1 / 3 * soilcLayerNatveg[, , 2]
  getNames(topsoilc) <- "topsoilc"

  mstools::toolExpectTrue(
    !any(is.na(topsoilc)),
    "no NAs in topsoilc",
    falseStatus = "warn"
  )

  weight <- calcOutput("LandArea", aggregate = FALSE)

  return(list(x            = topsoilc,
              weight       = weight,
              unit         = "t per ha",
              description  = "Topsoil carbon in tons per hectar for natural vegetation.",
              isocountries = FALSE))
}
