#' @title readLabourProdImpactEmu
#' @description read in labour productivity impacts from climate change emulated by the LAMACLIMA project
#' @description based on method of Orlov et al. 2019. Economics of Disasters and Climate Change, 3(3), 191-211.
#' @return magpie object of gridded productivity loss in percent (0-100)
#' @author Michael Windisch, Florian Humpenöder, Felicitas Beier
#' @seealso \code{\link[madrat]{readSource}}
#' @importFrom magclass new.magpie read.magpie mbind
#' @importFrom mstools toolCoord2Isocell

readLabourProdImpactEmu <- function() {

  mapping <- toolGetMappingCoord2Country(pretty = TRUE)
  x       <- NULL
  for (exp in (c("CTL_rcp585", "CTL_rcp119", "FRST", "CROP", "HARV", "IRR"))) {
    for (int in c("300", "400")) {
      for (sta in c("ensmean", "ensstd")) {
        for (fct in c("ISO", "HOTHAPS")) {
          tmp <- read.magpie(paste0(exp, "/", exp, "_", sta, "_laborprod_wby_grdays_", fct, "_", int, "_0.5.nc"))
          # reduce to 67420 cells and fill missings with 0
          inter <- new.magpie(cells_and_regions = mapping$coords,
                              years = getItems(tmp, dim = 2),
                              names = getItems(tmp, dim = 3),
                              fill = 0)
          tmp <- tmp[intersect(getItems(inter, dim = 1), getItems(tmp, dim = 1)), , ]
          inter[intersect(getItems(inter, dim = 1), getItems(tmp, dim = 1)), , ] <- tmp
          inter           <- toolCoord2Isocoord(inter)
          # clean up object dimensions
          getYears(inter) <- seq(from = 1995, to = 2095, by = 1)
          inter           <- collapseNames(inter)
          getNames(inter) <- paste0(exp, ".", fct, ".", int, "W", ".", sta)
          getSets(inter)  <- c("x", "y", "iso", "year", "data")
          # bind to previous inputs
          x <- mbind(x, inter)
        }
      }
    }
  }

  return(x)

}
