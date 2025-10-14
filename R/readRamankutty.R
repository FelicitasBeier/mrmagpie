#' @title readRamankutty
#'
#' @description Read in data of Ramankutty dataset
#' (Source: Ramankutty N, Foley JA, Norman J and McSweeney K (2002) The global
#' distribution of cultivable lands: current patterns and sensitivity to possible
#' climate change. Global Ecology and Biogeography, 11, 377-392.).
#' Link to data: https://www.nelson.wisc.edu/sage/data-and-models/global-land-use/grid.php
#'
#' @return magpie object
#'
#' @examples
#' \dontrun{
#' readSource("Ramankutty", convert = "onlycorrect")
#' }
#'
#' @author Felicitas Beier
#'
readRamankutty <- function() {
  x <- suppressWarnings(raster::brick("raw_data_land_suit_land_suit_0.50x0.50.nc"))
  return(as.magpie(x))
}
