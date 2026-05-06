# readLabourProdImpactOrlov

read in labour productivity impacts from climate change from Orlov (see
Orlov et al. 2019. Economic Losses of Heat-Induced Reductions in Outdoor
Worker Productivity: a Case Study of Europe. Economics of Disasters and
Climate Change, 3(3), 191-211.)

## Usage

``` r
readLabourProdImpactOrlov(
  subtype = "IPSL-CM5A-LR_rcp85_wbgtod_hothaps_400W.nc"
)
```

## Arguments

- subtype:

  subtype of choice between indoor outdoor work, GCM, work intesnsity
  (300W medium, 400W high, rcp)

## Value

magpie object of gridded productivity as share of 1 (full productivity)

## See also

[`readSource`](https://rdrr.io/pkg/madrat/man/readSource.html)

## Author

David Chen
