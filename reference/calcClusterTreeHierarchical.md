# calcClusterTreeHierarchical

calculates hierarchical clustering tree

## Usage

``` r
calcClusterTreeHierarchical(
  regionscode,
  mode = "h",
  weight = NULL,
  lpjml = c(natveg = "LPJmL4", crop = "LPJmL5"),
  clusterdata = "yield_airrig"
)
```

## Arguments

- regionscode:

  regionscode of the regional mapping to be used. Must agree with the
  regionscode of the mapping mentioned in the madrat config! Can be
  retrieved via `regionscode()`.

- mode:

  Clustering type. At the moment you can choose between complete linkage
  clustering (h), single linkage clustering (s) and Ward clustering (w).

- weight:

  named vector with weighting factors for each region for the cluster
  distribution, e.g. weight = c(AFR = 3, EUR = 0.5). weight \> 1 will
  grant more cluster to a region and weight \< 1 less cluster than by
  default.

- lpjml:

  defines LPJmL version for crop/grass and natveg specific inputs

- clusterdata:

  similarity data to be used to determine clusters: yield_airrig
  (current default) or yield_increment

## Value

A magpie object

## Author

Jan Philipp Dietrich
