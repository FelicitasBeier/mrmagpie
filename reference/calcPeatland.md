# calcPeatland

This function calculates degraded and intact peatland area at cell
level. The function takes degraded and intact peatland area from the
Global Peatland Database (GPD) at the national level and downscales the
peatland area to grid cell level using gridded potential peatland area.
The GPD has been provided by Alexandra Barthelmes. The potential
peatland area has been provided by Leifeld_2018 (DOI
10.1038/s41467-018-03406-6).

## Usage

``` r
calcPeatland(subtype = "degraded", cells = "lpjcell")
```

## Arguments

- subtype:

  degraded (default) or intact

- cells:

  "magpiecell" or "lpjcell"

## Value

magpie object in cellular resolution

## Author

Florian Humpenoeder

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("Peatland", aggregate = FALSE)
} # }
```
