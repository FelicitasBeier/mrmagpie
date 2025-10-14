#' @title calcRangelandsMaxNew
#' @description Calculates rangelands maximum output
#' @param report Either 'harvest' or 'lsu/ha' controlling what values are output by the function.
#' @param lsuLevels Livestock unit levels in the source folder
#' @param lpjml Defines LPJmL version for crop/grass and natveg specific inputs
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param scenario specify ssp scenario
#' @param cells    "magpiecell" for 59199 cells or "lpjcell" for 67420 cells
#' @return magpie object in cellular resolution
#' @author Marcos Alves
calcRangelandsMaxNew <- function(lsuLevels = c(seq(0, 2.2, 0.2), 2.5),
                                 lpjml = "lpjml5p2_pasture",
                                 climatetype = "MRI-ESM2-0:ssp370",
                                 scenario = "/co2/Nreturn0p5/limN", report = "harvest", #nolint
                                 cells = "lpjcell") {
  # tidyr variables initiation to avoid problems with buildlibrary()
  water <- NULL
  year  <- NULL
  cell  <- NULL
  lsuHa <- NULL
  value <- NULL

  lsuLevels <- gsub("\\.", "p", lsuLevels)
  years <- seq(1965, 2100, 5)
  y <- list()
  for (lsu in lsuLevels) {
    .subtype <- paste0(lpjml, ":", climatetype, paste0(scenario, "/", lsu))
    message(.subtype)
    hist <- readSource("LPJmL_new", subtype = paste(.subtype, "grass_pft_hist", sep = ":"),
                       convert = FALSE)[, , "mgrass"]
    scen <- readSource("LPJmL_new", subtype = paste(.subtype, "grass_pft_scen", sep = ":"),
                       convert = FALSE)[, , "mgrass"]
    x <- mbind(hist, scen)
    x <- x[, years, ]
    getNames(x) <- gsub("mgrass", lsu, getNames(x))
    y[[lsu]] <- x
  }
  y <- mbind(y)
  w <- as.array(y)
  w <- as.data.frame(w)
  w$cell <- rownames(w)
  y <-
    tidyr::pivot_longer(
      w,
      cols = -dplyr::last_col(),
      names_to = c("year", ".value", "water"),
      names_sep = "\\."
    )
  y <- tidyr::pivot_longer(y, cols = dplyr::matches("^[0-9]*p*[0-9]*$"), names_to = "lsuHa")
  y <- dplyr::mutate(y, year = substr(year, 2, 5), lsuHa = as.numeric(gsub("p", ".", lsuHa)))
  y <- dplyr::group_by(y, cell, year, water)
  y <- dplyr::filter(y, value == max(value))
  y <- dplyr::filter(y, lsuHa == min(lsuHa))

  if (report == "harvest") {
    maxHarvest <- as.magpie(y[, -4], spatial = 1, tidy = TRUE, replacement = ".")
    maxHarvest <- add_dimension(maxHarvest, dim = 3.1, add = "CO2", nm = "range")

    x <- maxHarvest
    unit <- "gC/m2/y"
    description <- "Maximum pasture yields obtained with rangelands"
  }

  if (report == "lsu/ha") {
    optmLsu <- as.magpie(y[, -5], spatial = 1, tidy = TRUE, replacement = ".")
    optmLsu[optmLsu > 2.5] <- 2.5
    optmLsu <- add_dimension(optmLsu, dim = 3.1, add = "CO2", nm = "range")

    x <- optmLsu
    unit <- "lsu/ha"
    description <- "Optimal LSU density that corresponds to the maximum grass yields"
  }

  getSets(x, fulldim = FALSE)[1] <- c("x.y.iso")

  # reduce to old grid cell format
  if (cells == "magpiecell") {
    x <- toolCoord2Isocell(x)
  }

  return(list(x = x,
              weight = NULL,
              unit = unit,
              description = description,
              isocountries = FALSE))
}
