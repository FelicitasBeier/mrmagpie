# calcClusterBase

Reads a series of MAgPIE files and combines them to a matrix which is
then used for calculating a clustering.

## Usage

``` r
calcClusterBase(
  clusterdata = "yield_airrig",
  lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de", crop =
    "ggcmi_phase3_nchecks_9ca735cb")
)
```

## Arguments

- clusterdata:

  similarity data to be used to determine clusters: yield_airrig
  (current default) or yield_increment

- lpjml:

  defines LPJmL version for crop/grass and natveg specific inputs

## Value

A matrix containing the data

## See also

[`calcCluster`](calcCluster.md)

## Author

Jan Philipp Dietrich, Felicitas Beier
