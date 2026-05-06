# toolClusterPerRegion

This function calculates an appropriate number of clusters per region as
it is needed for ClusterKMeans

## Usage

``` r
toolClusterPerRegion(cells, ncluster, weight = NULL)
```

## Arguments

- cells:

  spatial names as returned by `getCells`

- ncluster:

  The desired total number of clusters.

- weight:

  named vector with weighting factors for each region for the cluster
  distribution, e.g. weight=c(AFR=3,EUR=0.5). weight \> 1 will grant
  more cluster to a region and weight \< 1 less cluster than by default.

## Value

A matrix with regions in rows and number of cells and clusters in
columns

## See also

[`calcClusterKMeans`](calcClusterKMeans.md),
[`calcClusterBase`](calcClusterBase.md)

## Author

Jan Philipp Dietrich
