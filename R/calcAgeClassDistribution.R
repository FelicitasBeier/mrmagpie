#' @title calcAgeClassDistribution
#' @description Forest area in 15 age classes. With \code{dataset = "GFAD"} (default) the
#' distribution comes directly from the Global Forest Age Dataset (GFAD v1.1, Poulter et al.
#' 2019): per-cell class fractions times cell area, summed over the four forest types.
#' With \code{dataset = "GAMI"} the *age distribution* comes from GAMI v2.1 (Besnard et al.
#' 2024) while the per-cell forest *area* is taken from the GFAD product. GAMI supplies only
#' the within-forest age shares (12 classes, summing to 1); keeping the GFAD forest-area map
#' as the spatial weight conserves each cell's total forest area and isolates the pure
#' age-structure change. GAMI's 12 twenty-year classes are mapped onto the 15 GFAD ten-year
#' classes (each split 50/50 into two adjacent classes; all classes >=140 yr, i.e. GAMI
#' `140-160 ... >299`, collapse into `class15`, the mature bin routed to `acx` in module 28).
#' Where GAMI reports no forest but GFAD does, the GFAD age shape is retained (area-conserving).
#' @param cells lpjcell for 67420 cells
#' @param dataset "GFAD" (default) or "GAMI"
#' @param gamiYear GAMI reference year when \code{dataset = "GAMI"}: "y2010" (default, matches
#'   the GFAD vintage) or "y2020"
#'
#' @return magpie object in cluster resolution
#' @author Abhijeet Mishra, Felicitas Beier, Florian Humpenoeder
#'
#' @examples
#' \dontrun{
#' calcOutput("AgeClassDistribution", aggregate = FALSE)
#' calcOutput("AgeClassDistribution", dataset = "GAMI", aggregate = FALSE)
#' }

calcAgeClassDistribution <- function(cells = "lpjcell", dataset = "GFAD", gamiYear = "y2010") {

  if (cells == "magpiecell") {
    stop("magpiecell is no longer supported")
  }

  if (dataset == "GFAD") {

    # Cell fraction from Poulter data set
    poulterDataset <- readSource("GFAD", convert = "onlycorrect")

    # Calculate cellarea
    mapping <- mstools::toolGetMappingCoord2Country()
    cb <- toolGetMapping("LPJ_CellBelongingsToCountries.csv",
                         type = "cell", where = "mrcommons")
    # 111e3 m ~ length of one degree of latitude; (111e3 * 0.5)^2 = area of a 0.5-degree cell
    # at the equator (m2); cos(lat) corrects for the convergence of meridians towards the poles
    cellArea <- (111e3 * 0.5) * (111e3 * 0.5) * cos(cb$lat / 180 * pi)
    cellArea <- as.magpie(cellArea, spatial = 1) /  1e10 # convert to Mha
    getItems(cellArea, dim = 1, raw = TRUE) <- paste(mapping$coords, mapping$iso, sep = ".")
    getSets(cellArea) <- c("x", "y", "iso", "year", "data")

    out <- poulterDataset * cellArea
    out <- dimSums(out, dim = 3.1)
    getNames(out) <- gsub(pattern = "X", replacement = "class", x = getNames(out))

    descr <- "Forest area in 15 age classes based on the Global Forest Age Dataset (GFAD) from Poulter et al. 2019"

  } else if (dataset == "GAMI") {

    # within-forest age shares from GAMI v2.1 (native 0.5-degree grid, 12 classes, sum to 1 where forest exists)
    gami <- readSource("GAMI", convert = "onlycorrect")

    # project GAMI onto the 67420 lpjcell grid and attach ISO codes (done here, not in the read function)
    coordMap <- mstools::toolGetMappingCoord2Country()
    gami <- gami[coordMap$coords, , ]
    getItems(gami, dim = 1, raw = TRUE) <- paste(coordMap$coords, coordMap$iso, sep = ".")
    getSets(gami) <- c("x", "y", "iso", "year", "ageClass")

    if (!gamiYear %in% getItems(gami, 2)) {
      stop("calcAgeClassDistribution: gamiYear must be one of ", paste(getItems(gami, 2), collapse = ", "))
    }
    shares <- collapseNames(gami[, gamiYear, ])

    # per-cell forest AREA from the GFAD product (same forest-area map -> isolates the age change)
    gfad  <- calcOutput("AgeClassDistribution", dataset = "GFAD", aggregate = FALSE)
    cellIt <- getItems(gfad, 1)
    gfadM  <- matrix(as.numeric(gfad[cellIt, , ]), nrow = length(cellIt), ncol = 15,
                     dimnames = list(cellIt, paste0("class", 1:15)))
    forestArea <- rowSums(gfadM)
    sharesM <- matrix(as.numeric(shares[cellIt, , ]), nrow = length(cellIt), ncol = 12,
                      dimnames = list(cellIt, NULL))

    # remap 12 GAMI 20-yr classes -> 15 GFAD 10-yr classes
    remap <- matrix(0, nrow = 12, ncol = 15)
    for (g in 1:7) {
      remap[g, 2 * g - 1] <- 0.5
      remap[g, 2 * g]     <- 0.5
    }
    remap[8:12, 15] <- 1                                   # all GAMI classes >=140 yr -> class15 (acx)
    shape15 <- sharesM %*% remap                           # rowSums = 1 (forest) or 0

    # fallback: keep the GFAD age shape where GAMI has no forest but GFAD does (area-conserving)
    gfadShape <- gfadM / pmax(forestArea, 1e-12)
    useGami   <- rowSums(sharesM) > 1e-6
    shape15[!useGami, ] <- gfadShape[!useGami, ]

    areaByClass <- shape15 * forestArea                    # per-cell total forest area == GFAD (conserved)

    out <- gfad                                            # reuse the GFAD magpie skeleton (cells x 1 x 15)
    out[cellIt, , ] <- areaByClass

    descr <- paste0("Forest area in 15 age classes: age distribution from GAMI v2.1 (Besnard et al. 2024, ",
                    gamiYear, ") on the GFAD-derived per-cell forest-area weight")

  } else {
    stop("calcAgeClassDistribution: 'dataset' must be 'GFAD' or 'GAMI'")
  }

  return(list(x = out,
              weight = NULL,
              unit = "Mha",
              description = descr,
              isocountries = FALSE))
}
