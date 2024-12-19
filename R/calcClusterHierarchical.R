#' @title calcClusterHierarchical
#'
#' @description Performs MAgPIE hierarchical clustering and calculates corresponding spam
#'              relation matrix
#'
#' As the creation of a clustering tree is very time consuming the function
#' checks first in the input folder if the corresponding data already exists
#' and if not it stores the tree information in the input folder for later use
#' in the next execution of this function.
#'
#' @param regionscode regionscode of the regional mapping to be used.
#'                    Must agree with the regionscode of the mapping mentioned
#'                    in the madrat config!
#'                    Can be retrieved via \code{regionscode()}.
#' @param ncluster    The desired total number of clusters.
#' @param clusterdata similarity data to be used to determine clusters:
#'                    yield_airrig (current default) or yield_increment
#' @param mode        Clustering type. At the moment you can choose between
#'                    complete linkage clustering (h),
#'                    single linkage clustering (s) and
#'                    Ward clustering (w).
#' @param lpjml       defines LPJmL version for crop/grass and natveg specific inputs
#' @param weight      named vector with weighting factors for each region for the
#'                    cluster distribution, e.g. weight=c(AFR=3,EUR=0.5).
#'                    weight > 1 will grant more cluster to a region and
#'                    weight < 1 less cluster than by default.
#' @return A mapping between regions and clusters
#' @author Jan Philipp Dietrich
#' @importFrom magclass getCells ncells getRegions new.magpie getSets getSets<-
#' @importFrom stats hclust cutree
#' @seealso \code{\link{calcCluster}}, \code{\link{calcClusterKMeans}}
#' @export

calcClusterHierarchical <- function(regionscode, ncluster, lpjml = "lpjml5.9.5-m1",
                                    clusterdata = "yield_airrig", mode = "h", weight = NULL) {

  fullfit <- attributes(calcOutput("ClusterTreeHierarchical", regionscode = regionscode, mode = mode, weight = weight,
                                   lpjml = lpjml, clusterdata = clusterdata, aggregate = FALSE))$hclust

  clusters <- cutree(fullfit, k = ncluster)
  # sort clusters by regions
  cl <- NULL
  regions <- unique(sub("^.*\\.(.*)\\..*$", "\\1", fullfit$labels))
  for (r in regions) {
    cl <- c(cl, unique(clusters[grep(paste0(".", r, "."), names(clusters), fixed = TRUE)]))
  }
  if (length(cl) != ncluster) {
    stop("Some cluster seem to exist across region borders or there are less clusters than demanded!")
  }
  tmp <- order(cl)[clusters]
  names(tmp) <- names(clusters)
  clusters <- tmp
  out <- new.magpie(paste(names(clusters), clusters, sep = "."), fill = 1)
  getSets(out, fulldim = FALSE)[1] <- "country.region.cell.cluster"

  return(list(x = out,
              weight = NULL,
              unit = "1",
              description = "Hierarchical cluster mapping between cells and clusters",
              isocountries = FALSE))
}
