#' @title readMehta2024
#' @description reads in Global Area Equipped for Irrigation for years 1900-2015 from Mehta et al. (2022)
#'
#' @param subtype data subtype to be downloaded
#'                Subtypes available:
#'                Combination of version (v3, v4) and data source:
#'                'GMIA': gridded base map for downscaling from Stefan et al. (2013).
#'                        Global Map of Irrigation Areas version 5.
#'                'Meier2018': gridded base map for downscaling from Meier, et al. (2018).
#'                             Global Irrigated Areas.
#'                Separated by "_" (e.g., subtype = "v4_GMIA")
#'
#' @author  Felicitas Beier
#' @seealso [correctMehta2024()]
#' @examples
#' \dontrun{
#' a <- readSource("Mehta2024")
#' }
#' @importFrom terra aggregate project rast global
#' @importFrom magclass as.magpie
#' @importFrom mstools toolGetMappingCoord2Country

readMehta2024 <- function(subtype = "v4_GMIA") {

  if (grepl("GMIA", subtype)) {
    dataname <- "G_AEI_"
    itemname <- "AEI_Mehta2024_Siebert2013"
  } else if (grepl("Meier2018", subtype)) {
    dataname <- "MEIER_G_AEI_"
    itemname <- "AEI_Mehta2024_Meier2018"
  } else {
    stop("The selected subtype is not available for downloadMehta2024.
         Please select 'GMIA' or 'Meier2018' in combination with version
         argument (v3 or v4).")
  }

  years  <- c(seq(1900, 1970, by = 10),
              seq(1980, 2015, by = 5))

  files <- c(paste0(dataname, years, ".ASC"))

  .transformObject <- function(x) {

    resolution <- terra::rast(resolution = 0.5)

    # global sum of AEI (in ha)
    checkSum <- terra::global(x, sum, na.rm = TRUE)
    # aggregate to 0.5 degree
    x <- suppressWarnings(terra::aggregate(x, fact = 6, fun = "sum", na.rm = TRUE))
    # Check whether sum before and after aggregation is the same.
    if (any(round(checkSum - terra::global(x, sum, na.rm = TRUE), digits = 0) != 0)) {
      warning(paste0("The global sum of AEI before and after aggregation differ: ",
                     "Deviation is: ",
                     round(checkSum - terra::global(x, sum, na.rm = TRUE), digits = 4)))
      stop("There is an issue with the aggregation. Please check mrmagpie::readMehta2024")
    }
    x <- suppressWarnings(terra::project(x, resolution))
    x <- as.magpie(x)

    return(x)
  }

  # read in data and transform to magpie object
  out <- NULL
  for (file in files) {

    aei <- terra::rast(file)
    print(paste0("Read in ", file))
    aei <- .transformObject(x = aei)

    getItems(aei, dim = 2) <- gsub(dataname, "y", getItems(aei, dim = 3))
    getItems(aei, dim = 3) <- itemname

    out <- mbind(out, aei)
  }

  # reduce number of cells
  map67420 <- readRDS(system.file("extdata", "mapLPJcells2Coords.rds",
                                  package = "magpiesets"))
  out      <- out[map67420$coords, , ]

  # rename cells
  map           <- toolGetMappingCoord2Country()
  out           <- out[map$coords, , ]
  getCells(out) <- paste(map$coords, map$iso, sep = ".")
  getSets(out)  <- c("x", "y", "iso", "year", "data")

  # transform units to Mha
  out <- out * 1e-06

  return(out)
}
