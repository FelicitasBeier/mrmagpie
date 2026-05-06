# calcDegradationYieldReduction

Function creates dummy file for including yield reduction coefficients
to represent land degradation

## Usage

``` r
calcDegradationYieldReduction(cells = "lpjcell")
```

## Arguments

- cells:

  number of halfdegree grid cells to be returned. Options: "magpiecell"
  (59199), "lpjcell" (67420)

## Value

magpie object in cellular resolution

## Author

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("DegradationYieldReduction", aggregate = FALSE)
} # }
```
