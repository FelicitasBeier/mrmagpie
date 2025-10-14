#' @title calcAvlLandSi
#' @description Extracts si0 and nsi0 areas based on Ramankutty dataset
#'
#' @param cells magpiecell (59199 cells) or lpjcell (67420 cells)
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier
#'
#' @examples
#' \dontrun{
#' calcOutput("AvlLandSi", aggregate = FALSE)
#' }
#'
#' @importFrom madrat readSource calcOutput
#' @importFrom magclass dimSums getCells getYears getNames mbind collapseDim
#' @importFrom mstools toolCoord2Isocell toolGetMappingCoord2Country
#' @importFrom magpiesets addLocation
#'

calcAvlLandSi <- function(cells = "lpjcell") {
  # input data (Ramankutty)
  x        <- readSource("Ramankutty", convert = "onlycorrect")

  # land area according to LUH in initialization year (1995) [note: landarea is constant, so the year does not matter]
  landarea <- calcOutput("LUH3",
                         landuseTypes = "magpie",
                         aggregate = FALSE,
                         cellular = TRUE,
                         irrigation = FALSE,
                         years = "y1995")
  croparea <- landarea[, , "crop"]
  landarea <- dimSums(landarea, dim = 3)

  # add missing cells to Ramankutty data (fill with 0)
  coords67420 <- paste(getItems(landarea, dim = "x", full = TRUE),
                       getItems(landarea, dim = "y", full = TRUE),
                       sep = ".")
  coordsRamankutty <-  getItems(x, dim = 1)
  tmp <- new.magpie(cells_and_regions = setdiff(coords67420, coordsRamankutty),
                    years = getYears(x),
                    names = getNames(x),
                    fill = 0)
  x <- mbind(x, tmp)
  x <- x[coords67420, , ]
  getItems(x, dim = 1, raw = TRUE) <- getItems(landarea, dim = 1)

  #### Calculations
  # cell is either suitable or not suitable for cropland
  isSI0               <- x
  isSI0[isSI0 > 0.1]  <- 1
  isSI0[isSI0 <= 0.1] <- 0
  getItems(isSI0, dim = 3) <- NULL
  getItems(isSI0, dim = 2) <- NULL

  # suitable
  si0  <- landarea * isSI0
  # correction of suitable land to LUH croparea (land is declared as suitable where LUH reports cropland)
  si0  <- pmax(croparea, si0)
  # not suitable for cropland
  nsi0 <- landarea - si0

  out  <- mbind(setNames(si0, "si0"), setNames(nsi0, "nsi0"))

  if (cells == "magpiecell") {
    out <- toolCoord2Isocell(out)
  } else if (cells == "lpjcell") {
    out <- out
  } else {
    stop("Please specify cells argument")
  }

  return(list(x = out,
              weight = NULL,
              unit = "Mha",
              description = paste0("si and nsi0 areas based on Ramankutty",
                                   "suitability information and LUH area",
                                   "information from initialization year"),
              isocountries = FALSE))
}
