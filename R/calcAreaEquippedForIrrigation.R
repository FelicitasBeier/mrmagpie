#' @title calcAreaEquippedForIrrigation
#' @description Calculates the area equipped for irrigation based on LUH3 or
#'              Mehta data sets.
#'              For LUH3, it assumes, that all cropland irrigated in the last
#'              20 years at least once is equipped for irrigation.
#'              Mehta et al. (2022) directly report Global Area Equipped for
#'              Irrigation for the years 1900-2015
#'
#' @param cellular    if TRUE: 0.5 degree resolution returned
#' @param selectyears default on "past"
#'
#' @return List of magpie objects with results on country/cellular level,
#'         weight on country level, unit and description.
#'
#' @author Benjamin Leon Bodirsky, Kristine Karstens, Felicitas Beier
#'
#' @importFrom mstools toolHoldConstant
#'
#' @seealso
#' [calcLanduseInitialisation()]
#' @examples
#' \dontrun{
#' calcOutput("AreaEquippedForIrrigation", source = "LUH3", cellular = TRUE, aggregate = FALSE)
#' }
calcAreaEquippedForIrrigation <- function(cellular = FALSE,
                                          selectyears = "past_til2020") {

  selectyears <- sort(magpiesets::findset(selectyears, noset = "original"))

  ##########################################
  #### Read in LUH3 irrigated area data ####
  ##########################################
  yearsNeeded <- as.integer(substring(selectyears, 2))
  lastYear <- utils::tail(yearsNeeded, 1)
  yearsNeeded <- (yearsNeeded[1] - 20):lastYear

  # calcLUH3 materialises the LUH3 raster in long format, so its peak memory
  # scales with the number of years requested per call. Load in chunks and
  # reduce each one to irrigated area per cell before the next one is read.
  chunkSize <- 20
  chunks <- split(yearsNeeded, ceiling(seq_along(yearsNeeded) / chunkSize))

  x <- mbind(lapply(chunks, function(yrs) {
    luh3 <- calcOutput("LUH3",
                       landuseTypes = "magpie",
                       irrigation = TRUE,
                       cellular = TRUE,
                       yrs = yrs,
                       aggregate = FALSE)
    collapseNames(luh3[, , "irrigated"])
  }))

  years <- as.numeric(substring(selectyears, 2))

  # Cropland that it is irrigated at least once in the last 20 years
  # is defined as "equipped for irrigation".
  luh <- mbind(lapply(years, function(year) {
    span <- paste0("y", (year - 20):year)
    setYears(Reduce(pmax, lapply(span, function(y) x[, y, ])), paste0("y", year))
  }))
  getItems(luh, dim = 3) <- "LUH3"

  # Naming of first dimension:
  # Temporarily (until 67k preprocessing merged)
  mapping <- toolGetMappingCoord2Country()
  getItems(luh, dim = 1, raw = TRUE) <- paste(mapping$coords, mapping$iso, sep = ".")
  # Temporarily (until 67k preprocessing merged)

  # rename sets
  getSets(luh) <- c("x", "y", "iso", "year", "data")

  ########################################
  ### Read in Mehta et al. (2024) data ###
  ########################################
  .readMehta <- function(subtype, itemName) {
    m <- readSource("Mehta2024", subtype = subtype, convert = "onlycorrect")
    m <- time_interpolate(m, interpolated_year = selectyears)
    # remove negative values introduced by time interpolation
    m <- pmax(m, 0)
    m <- m[, intersect(getItems(m, dim = 2), selectyears), ]
    getItems(m, dim = 3) <- itemName
    return(m)
  }

  mehta <- mbind(.readMehta("v4_GMIA", "Mehta2024_Siebert2013"),
                 .readMehta("v4_Meier2018", "Mehta2024_Meier2018"))

  #########################
  ### Combine data sets ###
  #########################
  years <- intersect(getItems(luh, dim = 2), getItems(mehta, dim = 2))
  out   <- mbind(luh[, years, ], mehta[, years, ])

  ##############
  ### Output ###
  ##############
  # aggregate to iso level
  if (!cellular) {
    out <- dimSums(out, dim = c("x", "y"))
    out <- toolCountryFill(out, fill = 0)
  }

  return(list(x            = out,
              weight       = NULL,
              unit         = "million ha",
              description  = "Area equipped for irrigation",
              isocountries = !cellular))
}
