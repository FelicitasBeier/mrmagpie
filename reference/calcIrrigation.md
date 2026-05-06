# calcIrrigation

This function extracts irrigation water (airrig: water applied
additionally to rainfall) from LPJmL for MAgPIE

## Usage

``` r
calcIrrigation(
  lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de", crop =
    "ggcmi_phase3_nchecks_9ca735cb"),
  climatetype = "GSWP3-W5E5:historical",
  cells = "lpjcell",
  rainfedweight = 0.01
)
```

## Arguments

- lpjml:

  Defines LPJmL version for crop/grass and natveg specific inputs

- climatetype:

  Switch between different climate scenarios

- cells:

  Number of cells to be returned: "magpiecell" for 59199 cells or
  "lpjcell" for 67420 cells

- rainfedweight:

  For clustering airrig is weighted with cropland_irrigated +
  rainfedweight \* cropland_rainfed (default: 0.01)

## Value

magpie object in cellular resolution

## Author

Felicitas Beier, Abhijeet Mishra

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("Irrigation", aggregate = FALSE)
} # }
```
