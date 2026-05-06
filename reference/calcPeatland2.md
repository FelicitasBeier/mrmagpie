# calcPeatland2

This function calculates degraded and intact peatland area at cell
level. The function takes degraded and intact peatland area from the
Global Peatland Database 2022 (GPD2022) at the national level and
downscales the peatland area to grid cell level using gridded peatland
area from the Global Peatland Map 2.0 (GPM2) The data has been provided
by Alexandra Barthelmes.

## Usage

``` r
calcPeatland2(cells = "magpiecell", countryLevel = FALSE)
```

## Arguments

- cells:

  number of cells to be returned: magpiecell (59199), lpjcell (67420)

- countryLevel:

  Whether output shall be at country level. Requires aggregate=FALSE in
  calcOutput.

## Value

magpie object in cellular resolution

## Author

Florian Humpenoeder

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("Peatland2", aggregate = FALSE)
} # }
```
