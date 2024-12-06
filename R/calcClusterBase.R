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

calcClusterBase <- function(clusterdata = "yield_airrig",
                            lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de",
                                      crop = "ggcmi_phase3_nchecks_9ca735cb")) {
  d <- list()
  # read in data which should be used to determine cluster
  if (clusterdata == "yield_airrig") {

    d$yld    <- calcOutput("Yields", datasource = c(lpjml = lpjml[["crop"]]),
                           selectyears = 1995, aggregate = FALSE)
    d$irrig  <- calcOutput("Irrigation", lpjml = lpjml, cells = "lpjcell",
                           years = 1995, aggregate = FALSE)
    d$td     <- calcOutput("TransportTime", cells = "lpjcell",
                           aggregate = FALSE)[, , rep(1, floor(ndata(d$yld) / 2))]
    gridpop <- collapseNames(calcOutput("GridPop", source = "Gao", urban = FALSE,
                                        aggregate = FALSE, years = 1995)[, , "SSP2"])
    d$gp  <- gridpop[, , rep(1, floor(ndata(d$yld) / 2))]

  } else if (clusterdata == "yield_increment") {

    yield    <- calcOutput("Yields", datasource = c(lpjml = lpjml[["crop"]]),
                           selectyears = 1995, aggregate = FALSE)
    d$yld    <- collapseNames(yield[, , "rainfed"])
    d$irrig  <- (collapseNames(yield[, , "irrigated"][, , "pasture", invert = TRUE])
                 - collapseNames(yield[, , "rainfed"][, , "pasture", invert = TRUE]))
    d$td     <- calcOutput("TransportTime", cells = "lpjcell",
                           aggregate = FALSE)[, , rep(1, ndata(d$yld))]
    gridpop <- collapseNames(calcOutput("GridPop", source = "Gao", urban = FALSE,
                                        aggregate = FALSE, years = 1995)[, , "SSP2"])
    d$gp  <- gridpop[, , rep(1, floor(ndata(d$yld)))]

  } else {
    stop("Specify clusterdata argument: yield_airrig (old clustering data)
          or yield_increment (new clustering data)")
  }

  cdata <- do.call(cbind, lapply(d, wrap, list(1, c(2, 3))))
  cdata <- cdata[, colSums(cdata) != 0]
  cdata <- scale(cdata)
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
