#' @title calcAgeClassDistribution
#' @description This function calculates forest area in 15 age classes based on the
#' Global Forest Age Dataset (GFAD) from Poulter et al. 2019
#' @param cells lpjcell for 67420 cells or magpiecell for 59199 cells
#'
#' @return magpie object in cluster resolution
#' @author Abhijeet Mishra, Felicitas Beier, Florian Humpenoeder
#'
#' @examples
#' \dontrun{
#' calcOutput("AgeClassDistribution", aggregate = FALSE)
#' }
#'
#' @importFrom magclass where
#' @importFrom mstools toolCoord2Isocell

calcAgeClassDistribution <- function(cells = "lpjcell") {
  # Cell fraction from poulter data set
  poulterDataset <- readSource("GFAD", convert = "onlycorrect")

  # Calculate cellarea
  mapping <- toolGetMappingCoord2Country()
  cb <- toolGetMapping("LPJ_CellBelongingsToCountries.csv",
                       type = "cell", where = "mrcommons")
  cellArea <- (111e3 * 0.5) * (111e3 * 0.5) * cos(cb$lat / 180 * pi)
  cellArea <- as.magpie(cellArea, spatial = 1) /  1e10 # convert to Mha
  getItems(cellArea, dim = 1, raw = TRUE) <- paste(mapping$coords, mapping$iso, sep = ".")
  getSets(cellArea) <- c("x", "y", "iso", "year", "data")

  ######################
  out <- poulterDataset * cellArea

  out <- dimSums(out, dim = 3.1)

  getNames(out) <- gsub(pattern = "X",
                        replacement = "class",
                        x = getNames(out))

  if (cells == "magpiecell") {
    out      <- toolCoord2Isocell(out, cells = cells)
    cellArea <- toolCoord2Isocell(cellArea, cells = cells)
  }

  return(list(x = out,
              weight = NULL,
              unit = "Mha",
              description = "Forest area in 15 age classes based the Global Forest Age
                Dataset (GFAD) from Poulter et al. 2019",
              isocountries = FALSE))
}
