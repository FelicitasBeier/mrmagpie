#' @title calcAfforestationMask
#' @description Afforestation mask for where afforestation possible
#' @param subtype afforestation mask sub type
#' @param cells "magpiecell" or "lpjcell"
#' @return magpie object in cellular resolution
#' @author David Chen, Florian Humpenoeder
#'
#' @examples
#' \dontrun{
#' calcOutput("AfforestationMask", subtype = "noboreal", aggregate = FALSE)
#' }
#'
#' @importFrom magclass as.magpie
#' @importFrom mstools toolGetMappingCoord2Country toolCoord2Isocell

calcAfforestationMask <- function(subtype, cells = "lpjcell") {

  if (subtype == "unrestricted") {
    r <- terra::rast(resolution = 0.5, vals = 1)
  } else if (subtype == "noboreal") { # Exclude boreal regions > 50deg N
    r <- terra::rast(resolution = 0.5, vals = 1)
    lat <- terra::init(r, "y") |> terra::mask(r)
    r <- terra::mask(r, lat > 50, maskvalues = TRUE, updatevalue = 0)
  } else if (subtype == "onlytropical") { # only tropical areas between 20deg S and 20deg N
    r <- terra::rast(resolution = 0.5, vals = 1)
    lat <- terra::init(r, "y") |> terra::mask(r)
    r <- terra::mask(r, lat > 20, maskvalues = TRUE, updatevalue = 0)
    r <- terra::mask(r, lat < -20, maskvalues = TRUE, updatevalue = 0)
  }


  # get spatial mapping
  map <- toolGetMappingCoord2Country(pretty = TRUE)
  # transform raster to magpie object
  x <- as.magpie(terra::extract(r, map[c("lon", "lat")])[, -1], spatial = 1)
  dimnames(x) <- list("x.y.iso" = paste(map$coords, map$iso, sep = "."), "t" = NULL, "d3" = NULL)

  if (cells == "magpiecell") {
    x <- toolCoord2Isocell(x, cells = cells)
  }

  weight <- calcOutput("LandArea", aggregate = FALSE, cells = cells)

  return(list(x = x,
              weight = weight,
              unit = "binary",
              description = "",
              isocountries = FALSE))
}
