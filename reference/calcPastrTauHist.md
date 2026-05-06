# calcPastrTauHist

Calculates managed pastures Tau based on FAO yield trends for 1995.

## Usage

``` r
calcPastrTauHist(past_mngmt = "mdef", cells = "lpjcell")
```

## Arguments

- past_mngmt:

  Pasture management reference yield

- cells:

  "magpiecell" for 59199 cells or "lpjcell" for 67420 cells

## Value

List of magpie objects with results on country level, weight on country
level, unit and description.

## Author

Marcos Alves

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("PastrTauHist", past_mngmt)
} # }
```
