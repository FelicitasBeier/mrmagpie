#' @title calcPeatland
#' @description This function calculates degraded and intact peatland area at cell level.
#' The function takes degraded and intact peatland area from the Global Peatland Database
#' (GPD) at the national level and downscales the peatland area to grid cell level using
#' gridded potential peatland area.
#' The GPD has been provided by Alexandra Barthelmes. The potential peatland area has been
#' provided by Leifeld_2018 (DOI 10.1038/s41467-018-03406-6).
#' @param subtype degraded (default) or intact
#'
#' @return magpie object in cellular resolution
#' @author Florian Humpenoeder
#'
#' @examples
#' \dontrun{
#' calcOutput("Peatland", aggregate = FALSE)
#' }
#'
#' @importFrom magclass clean_magpie
#' @importFrom madrat toolGetMapping toolAggregate
#' @importFrom mstools toolIso2CellCountries toolCoord2Isocell

calcPeatland <- function(subtype = "degraded") {

  gpd         <- readSource("GPD", convert = TRUE)
  potPeatArea <- readSource("Leifeld2018", convert = "onlycorrect") + 10^-10

  # Total and drained peatland area
  peatAreaTotal   <- collapseNames(gpd[, , "PeatAreaTotal"])
  names(dimnames(peatAreaTotal)) <- c("iso", "t", "d3")
  peatAreaDrained <- collapseNames(gpd[, , "PeatAreaDrained"])
  names(dimnames(peatAreaDrained)) <- c("iso", "t", "d3")

  # Dissag. from country to cell
  map <- toolGetMappingCoord2Country(pretty = TRUE)
  peatAreaTotal   <- toolAggregate(x = toolIso2CellCountries(peatAreaTotal), rel = map,
                                   weight = potPeatArea, dim = 1, from = "iso", to = "coords")
  names(dimnames(peatAreaTotal)) <- c("coords", "t", "d3")
  peatAreaDrained <- toolAggregate(x = toolIso2CellCountries(peatAreaDrained), rel = map,
                                   weight = potPeatArea, dim = 1, from = "iso", to = "coords")
  names(dimnames(peatAreaDrained)) <- c("coords", "t", "d3")

  # potPeatArea is the upper limit of peatland area in a cell; this will reduce the peatland area of GPD!
  peatAreaTotal[peatAreaTotal > potPeatArea] <- potPeatArea[peatAreaTotal > potPeatArea]
  peatAreaDrained[peatAreaDrained > potPeatArea] <- potPeatArea[peatAreaDrained > potPeatArea]

  dimnames(peatAreaTotal) <- list("x.y.iso" = paste(map$coords, map$iso, sep = "."), "t" = NULL, "d3" = NULL)
  dimnames(peatAreaDrained) <- list("x.y.iso" = paste(map$coords, map$iso, sep = "."), "t" = NULL, "d3" = NULL)
  peatAreaIntact <- peatAreaTotal - peatAreaDrained

  if (subtype == "degraded") {
    description <- "Degraded peatland area (Mha) in 0.5 degree resolution based on Humpenoeder
                    et al 2020 (DOI 10.1088/1748-9326/abae2a)"
    x <- peatAreaDrained
  } else if (subtype == "intact") {
    description <- "Intact peatland area (Mha) in 0.5 degree resolution based on Humpenoeder
                    et al 2020 (DOI 10.1088/1748-9326/abae2a)"
    x <- peatAreaIntact
  }

  return(list(x = x,
              weight = NULL,
              unit = "Mha",
              description = description,
              isocountries = FALSE))

}
