# calcAgeClassDistribution

This function calculates forest area in 15 age classes based on the
Global Forest Age Dataset (GFAD) from Poulter et al. 2019

## Usage

``` r
calcAgeClassDistribution(cells = "lpjcell")
```

## Arguments

- cells:

  lpjcell for 67420 cells or magpiecell for 59199 cells

## Value

magpie object in cluster resolution

## Author

Abhijeet Mishra, Felicitas Beier, Florian Humpenoeder

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("AgeClassDistribution", aggregate = FALSE)
} # }
```
