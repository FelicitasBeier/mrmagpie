# calcBphMask

Mask of Datapoints of biogeophysical temperature change of afforestation
(degree C) to be used as weight. File is based on observation datasets
of Bright et al. 2017 and Duveiller et al. 2018

## Usage

``` r
calcBphMask(cells = "lpjcell")
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
calcOutput("BphMask", aggregate = FALSE)
} # }
```
