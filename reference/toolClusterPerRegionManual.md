# toolClusterPerRegionManual

This function translates weights into number of clusters per region as
it is needed for ClusterKMeans. Weights have to sum up to total number
of clusters.

## Usage

``` r
toolClusterPerRegionManual(cells, ncluster, ncluster2reg)
```

## Arguments

- cells:

  spatial names as returned by `getCells`

- ncluster:

  The desired total number of clusters.

- ncluster2reg:

  named vector with numbers per region

## Value

A matrix with regions in rows and number of cells and clusters in
columns

## See also

[`calcClusterKMeans`](calcClusterKMeans.md),
[`calcClusterBase`](calcClusterBase.md)

## Author

Kristine Karstens
