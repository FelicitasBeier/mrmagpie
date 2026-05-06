# Apply region names

This tool function replaces country names with region names in the
spatial dimension of the object. To avoid mixing up of cache files with
different regional aggregation the regioncode needs to supplied and
checked as well. Only if the supplied regions code agrees with the
region mapping currently chosen the function will return the data.

## Usage

``` r
toolApplyRegionNames(cdata, regionscode)
```

## Arguments

- cdata:

  a cluster data file as produced by cluster_base

- regionscode:

  regionscode of the regional mapping to be used. Must agree with the
  regionscode of the mapping mentioned in the madrat config! Can be
  retrieved via `regionscode()`.

## Value

the cluster data file with region names in spatial dimension rather than
country names

## See also

[`calcClusterKMeans`](calcClusterKMeans.md),
[`calcClusterBase`](calcClusterBase.md)

## Author

Jan Philipp Dietrich, Felicitas Beier
