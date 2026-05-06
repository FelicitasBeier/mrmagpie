# calcLabourProdImpact

Labour productivity impacts

## Usage

``` r
calcLabourProdImpact(
  timestep = "5year",
  subtype = "Orlov",
  cellular = TRUE,
  cells = "lpjcell"
)
```

## Arguments

- timestep:

  5year or yearly

- subtype:

  data source comes from

- cellular:

  cellular is true

- cells:

  "magpiecell" or "lpjcell"

## Value

List of magpie objects with results on 0.5deg grid level, weights based
on production value, unit (ratio) and description.

## Author

David Chen
