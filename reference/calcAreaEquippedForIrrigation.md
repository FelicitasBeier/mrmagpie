# calcAreaEquippedForIrrigation

Calculates the area equipped for irrigation based on LU2v2 or Mehta data
sets. For LUH3, it assumes, that all cropland irrigated in the last 20
years at least once is equipped for irrigation. Mehta et al. (2022)
directly report Global Area Equipped for Irrigation for the years
1900-2015

## Usage

``` r
calcAreaEquippedForIrrigation(cellular = FALSE, selectyears = "past_til2020")
```

## Arguments

- cellular:

  if TRUE: 0.5 degree resolution returned

- selectyears:

  default on "past"

## Value

List of magpie objects with results on country/cellular level, weight on
country level, unit and description.

## See also

\[calcLanduseInitialisation()\]

## Author

Benjamin Leon Bodirsky, Kristine Karstens, Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("AreaEquippedForIrrigation", source = "LUH3", cellular = TRUE, aggregate = FALSE)
} # }
```
