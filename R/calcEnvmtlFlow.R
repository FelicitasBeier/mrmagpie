#' @title calcEnvmtlFlow
#' @description This function calculates environmental flow requirements (EFR) for MAgPIE
#'              retrieved from LPJmL monthly discharge and water availability
#'
#' @param lpjml Defines LPJmL version for crop/grass and natveg specific inputs
#' @param climatetype Switch between different climate scenarios
#' @param seasonality grper (default): EFR in growing period per year; total:
#'                    EFR throughout the year; monthly: monthly EFRs
#'
#' @import magclass
#' @import madrat
#' @importFrom stats quantile
#' @importFrom mstools toolHarmonize2Baseline
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier
#'
#' @examples
#' \dontrun{
#' calcOutput("EnvmtlFlow", aggregate = FALSE)
#' }
#'
calcEnvmtlFlow <- function(lpjml = "lpjml5.9.5-m1",  climatetype = "GSWP3-W5E5:historical", seasonality = "grper") {

  ### Environemtal flow requirements (EFR) based on Smakthin (2004) method ###
  # EFRs that maintain a "fair ecosystem status", i.e. with Q90 low flow requirements (LFR)
  fair <- calcOutput("EFRSmakthin", lpjml = lpjml, climatetype = climatetype,
                     seasonality = seasonality, LFR_val = 0.1, aggregate = FALSE)
  getItems(fair, dim = 3) <- "fair"

  # EFRs that maintain a "good ecosystem status", i.e. with Q75 low flow requirements (LFR)
  # Note: High flow requirements (HFRs) are the same across different conservation statuses.
  #       Only LFRs vary.
  good <- calcOutput("EFRSmakthin", lpjml = lpjml, climatetype = climatetype,
                     seasonality = seasonality, LFR_val = 0.25, aggregate = FALSE)
  getItems(good, dim = 3) <- "good"

  out <- mbind(fair, good)

  ### EFRs according to planetary boundary (PB) as of Rockström et al. (2023) ###
  # 80% of monthly flow is reserved for environment
  pb <- calcOutput("EFRRockstroem", lpjml = lpjml, climatetype = climatetype,
                   seasonality = seasonality, aggregate = FALSE)
  getItems(pb, dim = 3) <- "pb"
  out <- mbind(out, pb)

  description <- paste0("EFR in ", seasonality)

  return(list(x = out,
              weight = NULL,
              unit = "mio. m^3",
              description = description,
              isocountries = FALSE))
}
