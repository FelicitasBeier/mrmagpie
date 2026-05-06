# calcLabourProdImpactEmu

Spatial and temporal aggr. of labour productivity impacts from climate
change emulated by LAMACLIMA

based on method of Orlov et al. 2019. Economics of Disasters and Climate
Change, 3(3), 191-211.

## Usage

``` r
calcLabourProdImpactEmu(
  timestep = "5year",
  cellular = TRUE,
  subtype = "impact",
  cells = "lpjcell"
)
```

## Arguments

- timestep:

  5-year or yearly

- cellular:

  cellular is true

- subtype:

  impact for rcp based laborprod decrease, relief for LCLM based relief
  of impact

- cells:

  "magpiecell" or "lpjcell"

## Value

List of magpie object of gridded (0.5) labour productivity as percentage
of full labour prod 1

## Author

Michael Windisch, Florian Humpenöder
