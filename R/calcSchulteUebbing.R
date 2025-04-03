#' @title calcSchulteUebbing
#' @description use inverse distance weighted interpolation on critical nitrogen surplus data set of
#' Schulte-Uebbing et al. 2022
#' @author Michael Crawford
#'
#' @return a magpie object containing critical nitrogen surplus
#'
#' @examples
#'
#' \dontrun{
#'   calcOutput("SchulteUebbing")
#' }
#' @importFrom rlang .data

calcSchulteUebbing <- function() {

  # Read and filter out rows where 'value' is missing
  criticalNitrogenSurplus <- readSource("SchulteUebbing") %>%
    dplyr::as_tibble() %>%
    dplyr::select(.data$x, .data$y, .data$value) %>%
    dplyr::filter(!is.na(.data$value))

  sp::coordinates(criticalNitrogenSurplus) <- ~x + y

  # Use template as the fine prediction grid
  template <- madrat::calcOutput("LandArea", aggregate = FALSE)
  template[, , ] <- 0
  template <- magclass::as.SpatRaster(template)

  # Convert template to SpatialPointsDataFrame for IDW interpolation, using a descriptive name
  templatePoints <- as.data.frame(template, xy = TRUE, na.rm = TRUE) %>%
    dplyr::select(.data$x, .data$y)
  sp::coordinates(templatePoints) <- ~x + y

  # Perform IDW interpolation
  idwResult <- gstat::idw(value ~ 1, criticalNitrogenSurplus, newdata = templatePoints, idp = 2)

  # Convert the result back to a raster
  idwRaster <- terra::rast(nrows = nrow(template), ncols = ncol(template),
                           xmin = terra::xmin(template), xmax = terra::xmax(template),
                           ymin = terra::ymin(template), ymax = terra::ymax(template))
  terra::crs(idwRaster) <- terra::crs(template)
  terra::values(idwRaster) <- NA

  # Remove NAs (oceans, e.g.)
  mask <- !is.na(terra::values(template))
  terra::values(idwRaster)[mask] <- idwResult$var1.pred

  out <- as.magpie(idwRaster) |>
    magclass::setYears(nm = "y2010") |>
    magclass::collapseDim("data1")

  getNames(out) <- "Critical nitrogen surplus"

  list(x            = out,
       weight       = NULL,
       unit         = "kg N ha-1 yr-1",
       min          = -395.96,
       description  = "Critical nitrogen surplus from Schulte-Uebbing et al. 2022",
       isocountries = FALSE)
}
