#' @title calcClusterTreeHierarchical
#'
#' @description calculates hierarchical clustering tree
#'
#' @param regionscode regionscode of the regional mapping to be used.
#'                    Must agree with the regionscode of the mapping
#'                    mentioned in the madrat config! Can be retrieved via \code{regionscode()}.
#' @param mode        Clustering type. At the moment you can choose between
#'                    complete linkage clustering (h),
#'                    single linkage clustering (s) and
#'                    Ward clustering (w).
#' @param weight      named vector with weighting factors for each region for the cluster distribution,
#'                    e.g. weight = c(AFR = 3, EUR = 0.5).
#'                    weight > 1 will grant more cluster to a region and
#'                    weight < 1 less cluster than by default.
#' @param lpjml       defines LPJmL version for crop/grass and natveg specific inputs
#' @param clusterdata similarity data to be used to determine clusters:
#'                    yield_airrig (current default) or yield_increment
#' @return A spam relation matrix
#' @author Jan Philipp Dietrich
#' @importFrom magclass getCells ncells getRegions getItems
#' @importFrom stats hclust cutree
#' @export
calcClusterTreeHierarchical <- function(regionscode, mode = "h", weight = NULL,
                                        lpjml = "lpjml5.9.5-m1", clusterdata = "yield_airrig") {

  # read in cluster base data
  cdata <- toolApplyRegionNames(cdata = calcOutput("ClusterBase", lpjml = lpjml,
                                                   clusterdata = clusterdata, aggregate = FALSE),
                                regionscode = regionscode)
  # transform to array
  cdata <- as.array(cdata)[, , ]

  # extract regional weights
  calcw <- function(weight, regions) {
    w <- rep(1, length(regions))
    names(w) <- regions
    if (!is.null(weight)) w[names(weight)] <- weight
    return(w)
  }

  fullfit <- list()
  fullfit$labels <- getItems(cdata, dim = 1)
  fullfit$order  <- rep(NA, ncells(cdata))
  names(fullfit$order) <- fullfit$labels
  nrow   <- 0
  nrows  <- NULL
  weight <- calcw(weight, getItems(cdata, "region"))

  for (r in getItems(cdata, "region")) {

    cells <- grep(paste0(".", r, "."), dimnames(cdata)[[1]], fixed = TRUE)
    dist  <- dist(cdata[cells, ], method = "euclidean") * weight[r]

    if (mode == "h") {
      fit <- hclust(dist, method = "complete")
    } else if (mode == "w") {
      fit <- hclust(dist^2, method = "ward")
    } else if (mode == "s") {
      fit <- hclust(dist, method = "single")
    } else {
      stop("Unknown mode \"", mode, "\"")
    }
    cellnumber <- as.numeric(sub("^.*\\.", "", fit$labels))
    fit$merge[fit$merge < 0] <- cellnumber[fit$merge[fit$merge < 0] * (-1)] * (-1)
    fit$merge[fit$merge > 0] <- fit$merge[fit$merge > 0] + nrow
    fit$order <- cellnumber[fit$order]

    fullfit$merge  <- rbind(fullfit$merge, fit$merge)
    fullfit$height <- c(fullfit$height, fit$height)

    fullfit$order[fit$labels] <- fit$order

    nrow  <- nrow(fullfit$merge)
    nrows <- c(nrows, nrow)
  }

  # add links between regions with huge heights
  fullfit$height <- c(fullfit$height, rep(2 * max(fullfit$height), length(nrows) - 1))

  for (n in nrows[-length(nrows)]) {
    fullfit$merge <- rbind(fullfit$merge, c(n, dim(fullfit$merge)[1]))
  }

  # bring data in the right order, adapt row numbers accordingly
  o <- order(fullfit$height)
  fullfit$merge <- fullfit$merge[o, ]
  b <- order(fullfit$merge)
  b <- b[min(which(fullfit$merge[b] > 0)):length(b)]
  fullfit$merge[b] <- order(o)[-length(o)]
  fullfit$height <- fullfit$height[o]
  attr(fullfit, "class") <- "hclust"

  x <- as.magpie(1)
  attr(x, "hclust") <- fullfit

  return(list(x = x,
              weight = NULL,
              unit = "1",
              description = "Hierarchical cluster tree provided as attribute hclust",
              isocountries = FALSE))
}
