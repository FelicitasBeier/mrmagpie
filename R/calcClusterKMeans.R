#' @title calcClusterKMeans
#'
#' @description Performs MAgPIE kmeans clustering and
#'              calculates corresponding spam relation matrix
#'
#' @param regionscode regionscode of the regional mapping to be used.
#'                    Must agree with the regionscode of the mapping
#'                    mentioned in the madrat config!
#'                    Can be retrieved via \code{regionscode()}.
#' @param ncluster    The desired total number of clusters.
#' @param weight      named vector with weighting factors for each region
#'                    for the cluster distribution, e.g. weight=c(AFR=3,EUR=0.5).
#'                    weight > 1 will grant more cluster to a region and
#'                    weight < 1 less cluster than by default.
#' @param cpr         cells-per-region information as returned by toolClusterPerRegionManual.
#'                    Weight and ncluster are ignored in case that cpr is provided!
#' @param seed        a single value, interpreted as an integer, or NULL, to define seed for random calculations
#' @param clusterdata similarity data to be used to determine clusters: yield_airrig (current default)
#'                    or yield_increment
#' @param lpjml       defines LPJmL version for crop/grass and natveg specific inputs
#' @return A mapping between regions and clusters
#' @author Jan Philipp Dietrich
#' @seealso \code{\link{toolClusterPerRegionManual}}, \code{\link{calcClusterHierarchical}}
calcClusterKMeans <- function(regionscode, ncluster, weight = NULL, cpr = NULL, seed = 42,
                              lpjml = "lpjml5.9.5-m1", clusterdata = "yield_airrig") {
  # read in clustering base data
  cdata <- toolApplyRegionNames(cdata = calcOutput("ClusterBase", aggregate = FALSE,
                                                   lpjml = lpjml, clusterdata = clusterdata),
                                regionscode = regionscode)

  if (is.null(cpr)) {
    cpr <- toolClusterPerRegion(sub("^[^\\.]*\\.", "", getCells(cdata)), ncluster, weight)
  }

  out <- new.magpie(cells_and_regions = getItems(cdata, dim = 1),
                    fill = 1)
  getSets(out, fulldim = FALSE)[1] <- "country.region.cell.cluster"

  # Transform to array
  cdata <- as.array(cdata)[, , ]

  # Clustering
  ccount <- 0
  set.seed(seed)
  for (r in dimnames(cpr)[[1]]) {
    cells <- grep(paste0(".", r, "."), dimnames(cdata)[[1]], fixed = TRUE)
    # We use the MacQueen algorithm as the HartiganWong algorithm does not converge reliably, as some
    # points in cdata are very close.
    fit <- stats::kmeans(cdata[cells, ], cpr[r, "clusters"], iter.max = 50000, algorithm = "MacQueen")
    getCells(out)[cells] <- paste(getCells(out)[cells], fit$cluster + ccount, sep = ".")
    ccount <- ccount + cpr[r, "clusters"]
  }

  set.seed(NULL)

  return(list(x = out,
              weight = NULL,
              unit = "1",
              description = "KMeans cluster mapping between cells and clusters",
              isocountries = FALSE))
}
