# calcBphTCRE

Transient Climate Response to accumulated doubling of CO2. File based on
CMIP5 +1perc CO2 per year experiment. To be used in the translation to
carbon equivalents of BphEffect

## Usage

``` r
calcBphTCRE(cells = "lpjcell")
```

## Arguments

- cells:

  lpjcell for 67420 cells or magpiecell for 59199 cells

## Value

magpie object in cellular resolution

## Author

Michael Windisch, Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("BphTCRE", aggregate = FALSE)
} # }
```
