# calcCluster

This function calculates the aggregation mapping for a given cluster
methodology

## Usage

``` r
calcCluster(
  ctype,
  regionscode = madrat::regionscode(),
  seed = 42,
  weight = NULL,
  lpjml = c(natveg = "LPJmL4", crop = "LPJmL5"),
  clusterdata = "yield_airrig"
)
```

## Arguments

- ctype:

  aggregation clustering type, which is a combination of a single
  letter, indicating the cluster methodology, and a number, indicating
  the number of resulting clusters. Available methodologies are
  hierarchical clustering (h), normalized k-means clustering (n),
  combined hierarchical/normalized k-means clustering (c) and manual
  setting for clusters per region (m). In the combined clustering
  hierarchical clustering is used to determine the cluster distribution
  among regions whereasit is manually set for the m type. Both use
  normalized k-means for the clustering within a region.

- regionscode:

  regionscode of the regional mapping to be used. Must agree with the
  regionscode of the mapping mentioned in the madrat config! Can be
  retrieved via `regionscode()`.

- seed:

  Seed for Random Number Generation. If set to NULL it is chosen
  automatically, if set to an integer it will always return the same
  pseudo-random numbers (useful to get identical clusters under
  identical inputs for n and c clustering)

- weight:

  Should specific regions be resolved with more or less detail? Values
  \> 1 mean higher share, \< 1 lower share e.g. cfg\$cluster_weight \<-
  c(LAM=2) means that a higher level of detail for region LAM if set to
  NULL all weights will be assumed to be 1 (examples:
  c(LAM=1.5,SSA=1.5,OAS=1.5), c(LAM=2,SSA=2,OAS=2))

- lpjml:

  defines LPJmL version for crop/grass and natveg specific inputs

- clusterdata:

  similarity data to be used to determine clusters: yield_airrig
  (current default) or yield_increment

## Value

map from cells to clusters as data.frame

## Author

Jan Philipp Dietrich

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("Cluster", ctype = "c200", aggregate = FALSE)
} # }
```
