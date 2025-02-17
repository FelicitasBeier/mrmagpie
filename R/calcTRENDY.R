#' @title calcTRENDY
#' @description Resample and use inverse distance weighted interpolation to align
#'              natural vegetation carbon densities from TRENDY with MAgPIE's input requirements
#' @author Michael Crawford
#'
#' @param subtype name of the TRENDY dataset to read. Current options includes:
#'     "CABLEPOP", "CARDAMOM", "CLASSIC", "CLM5.0", "DLEM",
#'     "ELM", "IBIS", "ISBACTRIP", "JSBACH", "LPJ-GUESS",
#'     "LPJml", "LPJwsl", "lpxqs", "OCN", "ORCHIDEE",
#'     "SDGVM", "VISIT"
#'
#' @return a magclass object with model, carbon pool in third dimension
#'
#' @examples
#' \dontrun{
#' readSource("TRENDY", subtype = "JSBACH")
#' }
#'
#' @importFrom dplyr %>%

calcTRENDY <- function(subtype) {

  # ---------------------------------------------------------------------------------------
  # Generate rasterList, eliminate unnecessary years
  rasterList <- madrat::readSource("TRENDY", subtype = subtype)

  magpieYears <- paste0("y", seq(1995, 2020, 5))

  rasterList <- lapply(rasterList, function(raster) {
    raster[[magpieYears]]
  })

  # ---------------------------------------------------------------------------------------
  # Generate template for fine scale IDW interpolation
  template <- calcOutput("LandArea", aggregate = FALSE)

  template[, , ] <- 0
  template <- magclass::as.SpatRaster(template)

  # ---------------------------------------------------------------------------------------
  # Load cLitter, cSoil, cVeg rasters for a given subtype
  interpolateRaster <- function(raster, rasterName) {

    interpolateLayer <- function(layer) {

      # Convert the source layer to a SpatialPointsDataFrame with coordinates
      r <- as.data.frame(layer, xy = TRUE, na.rm = TRUE)
      names(r) <- c("x", "y", "value")
      sp::coordinates(r) <- ~ x + y

      # Create the prediction grid from the template
      t <- as.data.frame(template, xy = TRUE, na.rm = TRUE) %>% dplyr::select("x", "y")
      sp::coordinates(t) <- ~ x + y

      # Perform IDW interpolation
      idwResult <- gstat::idw(value ~ 1, r, newdata = t, idp = 2)

      # Convert the result back to a raster
      idwRaster <- terra::rast(
        nrow = nrow(template), ncol = ncol(template),
        xmin = terra::xmin(template), xmax = terra::xmax(template),
        ymin = terra::ymin(template), ymax = terra::ymax(template)
      )
      terra::crs(idwRaster) <- terra::crs(template)
      terra::values(idwRaster) <- NA

      # Remove non-MAgPIE cells with a mask using the template
      mask <- !is.na(terra::values(template))
      terra::values(idwRaster)[mask] <- idwResult$var1.pred

      return(idwRaster)
    }

    raster <- terra::app(raster, function(x) x * 10) # kg C m-2 to t C ha-1
    raster <- terra::resample(raster, template, method = "bilinear")

    layerNames <- names(raster)
    resultsLayers <- lapply(X = layerNames, FUN = function(year) {
      layer <- terra::subset(raster, year)
      interpolateLayer(layer)
    })

    idwRaster <- terra::rast(resultsLayers)
    names(idwRaster) <- layerNames

    rasterMag <- magclass::as.magpie(idwRaster)
    rasterMag <- magclass::collapseDim(rasterMag, dim = 3)
    rasterMag <- magclass::add_dimension(rasterMag, dim = 3.1, add = "c_pools", nm = rasterName)

    return(rasterMag)
  }

  idwList <- Map(interpolateRaster, rasterList, names(rasterList))

  idwMag <- magclass::mbind(idwList)

  idwMag <- magclass::add_dimension(idwMag, dim = 3.1, add = "model", nm = subtype)

  return(list(
    x            = idwMag,
    weight       = NULL,
    unit         = "t C ha-1",
    min          = 0,
    description  = paste0("TRENDY dataset for the ESM << ", subtype, " >>, interpolated to the 0.5 deg."),
    isocountries = FALSE
  ))
}
