#' @title calcClusterBase
#' @description Reads a series of MAgPIE files and combines them to a matrix
#'              which is then used for calculating a clustering.
#' @param clusterdata similarity data to be used to determine clusters:
#'                    yield_airrig (current default) or yield_increment
#' @param lpjml defines LPJmL version for crop/grass and natveg specific inputs
#' @return A matrix containing the data
#' @author Jan Philipp Dietrich, Felicitas Beier
#' @seealso \code{\link{calcCluster}}
#' @importFrom magclass wrap read.magpie ndata
#' @importFrom madrat toolGetMapping

calcClusterBase <- function(clusterdata = "yield_airrig", lpjml = "lpjml5.9.5-m1") {
  d <- list()
  # read in data which should be used to determine cluster
  if (clusterdata == "yield_airrig") {

    d$yld    <- calcOutput("Yields", datasource = c(lpjml = lpjml), selectyears = 1995, aggregate = FALSE)
    d$irrig  <- calcOutput("Irrigation", lpjml = lpjml, years = 1995, aggregate = FALSE)
    d$td     <- calcOutput("TransportTime", aggregate = FALSE)[, , rep(1, floor(ndata(d$yld) / 2))]
    gridpop  <- collapseNames(calcOutput("GridPop", source = "Gao", urban = FALSE,
                                         aggregate = FALSE, years = 1995)[, , "SSP2"])
    d$gp  <- gridpop[, , rep(1, floor(ndata(d$yld) / 2))]

  } else if (clusterdata == "yield_increment") {

    yield    <- calcOutput("Yields", datasource = c(lpjml = lpjml), selectyears = 1995, aggregate = FALSE)
    d$yld    <- collapseNames(yield[, , "rainfed"])
    d$irrig  <- (collapseNames(yield[, , "irrigated"][, , "pasture", invert = TRUE])
                 - collapseNames(yield[, , "rainfed"][, , "pasture", invert = TRUE]))
    d$td     <- calcOutput("TransportTime", aggregate = FALSE)[, , rep(1, floor(ndata(d$yld) / 2))]
    gridpop  <- collapseNames(calcOutput("GridPop", source = "Gao", urban = FALSE,
                                         aggregate = FALSE, years = 1995)[, , "SSP2"])
    d$gp     <- gridpop[, , rep(1, floor(ndata(d$yld) / 2))]

  } else {
    stop("Specify clusterdata argument: yield_airrig (old clustering data)
          or yield_increment (new clustering data)")
  }

  # Unify structures
  d$gp <- setNames(d$gp, paste0("pop", seq_len(ndata(d$gp))))
  d$yld <- setNames(d$yld, sub(".", "_", getNames(d$yld), fixed = TRUE))
  d <- lapply(d, function(m) {
    # Ensure all datasets have 1995 as the year (they should have been selected for that)
    m <- setYears(m, "y1995")
    getSets(m, fulldim = FALSE)[3] <- "data"
    return(m)
  })

  # Combine into one object
  cdata <- do.call(mbind, d)
  cdata <- cdata[, , dimSums(cdata, dim = c(1, 2)) != 0]
  cdata <- wrap(cdata, list(1, c(2, 3)))

  cdata <- scale(cdata) # Scale converts the magpie object to a matrix
  colnames(cdata) <- paste0("i", seq_len(ncol(cdata)))

  # Transform to magpie object
  cdata <- as.magpie(cdata, spatial = 1)

  return(list(x            = cdata,
              weight       = NULL,
              unit         = "1",
              description  = "Similarity matrix as basis for clustering",
              isocountries = FALSE,
              putInPUC     = TRUE))
}
