#' @title readLabourProdImpactOrlov
#' @description read in labour productivity impacts from climate change from Orlov
#'              (see Orlov et al. 2019. Economic Losses of Heat-Induced Reductions 
#'              in Outdoor Worker Productivity: a Case Study of Europe.
#'              Economics of Disasters and Climate Change, 3(3), 191-211.)
#' @return magpie object of gridded productivity as share of 1 (full productivity)
#' @param subtype subtype of choice between indoor outdoor work, GCM, work intesnsity (300W medium, 400W high, rcp)
#' @author David Chen
#' @seealso \code{\link[madrat]{readSource}}
#' @importFrom magclass as.magpie
#' @importFrom ncdf4 nc_open ncvar_get

readLabourProdImpactOrlov <- function(subtype = "IPSL-CM5A-LR_rcp85_wbgtod_hothaps_400W.nc") {

  files <- vector()

  for (od in (c("od","id"))) {
    for (int in c("300W", "400W")) {
      for (gcm in c("IPSL-CM5A-LR", "GFDL-ESM2M", "HadGEM2-ES")) {
        for (rcp in c("rcp26", "rcp60", "rcp85")) {
          files <- append(files, paste0(gcm, "_", rcp, "_",
                                        "wbgt", od, "_", "hothaps",
                                        "_", int, ".nc" ))
        }
      }
    }
  }

  names(files) <- files
  file         <- toolSubtypeSelect(subtype, files)

  # read in selected file
  x <- read.magpie(file)
  # sort, rename and replace missing with 0
  mapping <- toolGetMappingCoord2Country()
  out <- new.magpie(cells_and_regions = mapping$coords,
                    years = getItems(x, dim = 2),
                    names = getItems(x, dim = 3),
                    fill = 0)
  out[intersect(getItems(x, dim = 1), mapping$coords), , ] <- x[intersect(getItems(x, dim = 1), mapping$coords), , ]
  getItems(out, dim = 1, raw = TRUE) <- paste(mapping$coords, mapping$iso, sep = ".")
  getSets(out) <- c("x", "y", "iso", "year", "data")

  return(x)

}
