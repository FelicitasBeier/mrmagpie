# calcClusterHierarchical

Performs MAgPIE hierarchical clustering and calculates corresponding
spam relation matrix

As the creation of a clustering tree is very time consuming the function
checks first in the input folder if the corresponding data already
exists and if not it stores the tree information in the input folder for
later use in the next execution of this function.

## Usage

``` r
calcClusterHierarchical(
  regionscode,
  ncluster,
  lpjml = c(natveg = "LPJmL4", crop = "LPJmL5"),
  clusterdata = "yield_airrig",
  mode = "h",
  weight = NULL
)
```

## Arguments

- regionscode:

  regionscode of the regional mapping to be used. Must agree with the
  regionscode of the mapping mentioned in the madrat config! Can be
  retrieved via `regionscode()`.

- ncluster:

  The desired total number of clusters.

- lpjml:

  defines LPJmL version for crop/grass and natveg specific inputs

- clusterdata:

  similarity data to be used to determine clusters: yield_airrig
  (current default) or yield_increment

- mode:

  Clustering type. At the moment you can choose between complete linkage
  clustering (h), single linkage clustering (s) and Ward clustering (w).

- weight:

  named vector with weighting factors for each region for the cluster
  distribution, e.g. weight=c(AFR=3,EUR=0.5). weight \> 1 will grant
  more cluster to a region and weight \< 1 less cluster than by default.

## Value

A mapping between regions and clusters

## See also

[`calcCluster`](calcCluster.md),
[`calcClusterKMeans`](calcClusterKMeans.md)

## Author

Jan Philipp Dietrich
