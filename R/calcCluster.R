#' @title calcCluster
#' @description This function calculates the aggregation mapping for a given cluster methodology
#' @param ctype aggregation clustering type, which is a combination of
#'              a single letter, indicating the cluster methodology, and
#'              a number, indicating the number of resulting clusters.
#'              Available methodologies are hierarchical clustering (h),
#'              normalized k-means clustering (n),
#'              combined hierarchical/normalized k-means clustering (c) and
#'              manual setting for clusters per region (m).
#'              In the combined clustering hierarchical clustering is used to
#'              determine the cluster distribution among regions whereasit is
#'              manually set for the m type. Both use normalized
#'              k-means for the clustering within a region.
#' @param regionscode regionscode of the regional mapping to be used.
#'                    Must agree with the regionscode of the mapping
#'                    mentioned in the madrat config!
#'                    Can be retrieved via \code{regionscode()}.
#' @param seed Seed for Random Number Generation. If set to NULL it is chosen automatically, if set to an integer it
#' will always return the same pseudo-random numbers (useful to get identical clusters under identical inputs for
#' n and c clustering)
#' @param weight Should specific regions be resolved with more or less detail? Values > 1 mean higher share, < 1 lower
#' share e.g. cfg$cluster_weight <- c(LAM=2) means that a higher level of detail for region LAM if set to NULL all
#' weights will be assumed to be 1 (examples: c(LAM=1.5,SSA=1.5,OAS=1.5), c(LAM=2,SSA=2,OAS=2))
#' @param lpjml defines LPJmL version for crop/grass and natveg specific inputs
#' @param clusterdata similarity data to be used to determine clusters: yield_airrig (current default)
#' or yield_increment
#' @return map from cells to clusters as data.frame
#' @author Jan Philipp Dietrich
#'
#' @examples
#' \dontrun{
#' calcOutput("Cluster", ctype = "c200", aggregate = FALSE)
#' }
#' @importFrom madrat calcOutput

calcCluster <- function(ctype, regionscode = madrat::regionscode(), seed = 42, weight = NULL,
                        lpjml = "lpjml5.9.5-m1", clusterdata = "yield_airrig") {

  mode <- substr(ctype, 0, 1)
  ncluster <- as.integer(substring(ctype, 2))

  if (mode == "n") {
    mapping <- calcOutput("ClusterKMeans", regionscode = regionscode, ncluster = ncluster, weight = weight,
                          seed = seed, lpjml = lpjml, clusterdata = clusterdata, aggregate = FALSE)
  } else if (mode == "h" || mode == "w" || mode == "s") {
    mapping <- calcOutput("ClusterHierarchical", regionscode = regionscode, ncluster = ncluster,
                          mode = mode, weight = weight, lpjml = lpjml, clusterdata = clusterdata, aggregate = FALSE)
  } else if (mode == "c") {
    calcCPR <- function(x) {
      clusters <- table(sub("\\..*$", "", unique(sub("\\..*\\.", ".", x))))
      cells <- table(sub("\\..*$", "", x))
      return(cbind(cells, clusters))
    }
    tmpmap  <- calcOutput("ClusterHierarchical", regionscode = regionscode, ncluster = ncluster,
                          mode = "h", weight = weight, lpjml = lpjml, clusterdata = clusterdata, aggregate = FALSE)
    mapping <- calcOutput("ClusterKMeans", regionscode = regionscode, ncluster = ncluster,
                          weight = weight, cpr = calcCPR(sub("^[^\\.]*\\.", "", getCells(tmpmap))),
                          seed = seed, lpjml = lpjml, clusterdata = clusterdata, aggregate = FALSE)
  } else if (mode == "m") {

    cdata <- toolApplyRegionNames(calcOutput("ClusterBase", aggregate = FALSE, lpjml = lpjml,
                                             clusterdata = clusterdata), regionscode)
    cpr   <- toolClusterPerRegionManual(sub("^[^\\.]*\\.", "", getCells(cdata)), ncluster = ncluster,
                                        ncluster2reg = weight)

    mapping <- calcOutput("ClusterKMeans", regionscode = regionscode, ncluster = ncluster, weight = weight,
                          cpr = cpr, seed = seed, lpjml = lpjml, clusterdata = clusterdata, aggregate = FALSE)
  } else {
    stop("Unkown clustering mode ", mode, "!")
  }

  pattern <- "^(.*)\\.(.*)\\.(.*)\\.(.*)$"
  mapping <- data.frame(cell    = sub(pattern, "\\1.\\3", getCells(mapping)),
                        cluster = sub(pattern, "\\2.\\4", getCells(mapping)),
                        region  = sub(pattern, "\\2", getCells(mapping)),
                        country = sub(pattern, "\\1", getCells(mapping)),
                        global  = "GLO")

  # Add coordinates to cell (x.y.iso)
  coord2country <- toolGetMappingCoord2Country()
  mapping$cell  <- paste(coord2country$coords, coord2country$iso, sep = ".")

  return(list(x = mapping,
              weight = NULL,
              unit = "1",
              class = "data.frame",
              description = "Mapping between cells and cluster",
              isocountries = FALSE,
              putInPUC = FALSE))
}
