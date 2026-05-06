# calcBphEffect

Biogeophysical temperature change of afforestation (degree C). File is
based on observation datasets of Bright et al. 2017 and Duveiller et al.
2018

## Usage

``` r
calcBphEffect(cells = "lpjcell")
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
calcOutput("BphEffect", aggregate = FALSE)
} # }
```
