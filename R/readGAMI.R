#' @title readGAMI
#' @description Read the Global Age Mapping Integration (GAMI) v2.1 forest-age dataset
#' (Besnard et al. 2024, GFZ Data Services, doi:10.5880/GFZ.1.4.2023.006), 0.5-degree
#' `class_fraction` product. GAMI provides the *within-forest* age distribution: for every
#' grid cell the fraction of the forested area in each of 12 age classes (sums to 1 where
#' forest exists), as a 20-member ensemble for 2010 and 2020. This read returns the
#' ensemble mean on GAMI's native 0.5-degree grid (720 x 360 cells) x 2 years x 12 age
#' classes. Cleaning (fill / NA to 0) is done in \code{correctGAMI}; the projection onto the
#' 67420 lpj-cell grid, the forest-area weighting and the mapping onto the 15 GFAD-style
#' MAgPIE age classes are done in \code{calcAgeClassDistribution}.
#'
#' GAMI's `forest_age` is a five-dimensional variable (time, longitude, latitude, age_class,
#' 20-member ensemble). terra/GDAL flattens the three non-spatial dimensions into unlabelled
#' layers - every age class is exposed as `age_class=0` and the member order is scrambled -
#' so the ensemble mean cannot be grouped reliably from the layer names. ncdf4 addresses the
#' named dimensions explicitly and is therefore used here instead of the usual terra reader.
#'
#' @return magpie object on GAMI's native 0.5-degree grid: cells x 2 years (2010, 2020) x 12 age classes
#' @author Florian Humpenoeder
#' @seealso \code{\link{calcAgeClassDistribution}}, \code{readGFAD}
#' @examples
#' \dontrun{
#' readSource("GAMI", convert = "onlycorrect")
#' }
readGAMI <- function() {

  file <- "GAMIv2-1_2010-2020_class_fraction_0deg50.nc"

  # using ncdf4, because terra::rast() flattens the non-spatial dims into 480 layers,
  # but labels every age class age_class=0 and scrambles member order
  nc <- ncdf4::nc_open(file)
  withr::defer(ncdf4::nc_close(nc))
  glon     <- ncdf4::ncvar_get(nc, "longitude")
  glat     <- ncdf4::ncvar_get(nc, "latitude")
  ageClass <- ncdf4::ncvar_get(nc, "age_class")           # "0-20", ..., "200-299", ">299"

  # forest_age dims (stored order): time(2), longitude(720), latitude(360), age_class(12), members(20)
  arr <- ncdf4::ncvar_get(nc, "forest_age")
  arr <- apply(arr, c(1, 2, 3, 4), mean, na.rm = TRUE)    # ensemble mean -> [time, lon, lat, age]

  # native 0.5-degree grid, cells named "<lon>.<lat>" (p-notation), longitude varying fastest
  cells   <- paste(gsub("\\.", "p", rep(glon, times = length(glat))),
                   gsub("\\.", "p", rep(glat, each  = length(glon))), sep = ".")
  years   <- c("y2010", "y2020")
  classNm <- paste0("gami", seq_along(ageClass))          # gami1..gami12
  out <- new.magpie(cells_and_regions = cells, years = years, names = classNm, fill = 0)
  for (tt in c(1, 2)) {
    for (k in seq_along(ageClass)) {
      out[, years[tt], classNm[k]] <- as.vector(arr[tt, , , k])
    }
  }
  getSets(out) <- c("x", "y", "year", "ageClass")
  # keep the human-readable class edges (0-20, ..., >299) for the calc's remap
  attr(out, "ageClassLabels") <- ageClass

  return(out)
}
