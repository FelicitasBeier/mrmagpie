# calcCarbonTests

This function extracts carbon densities from LPJ to MAgPIE

## Usage

``` r
calcCarbonTests(
  lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de", crop =
    "ggcmi_phase3_nchecks_9ca735cb"),
  climatetype = "GSWP3-W5E5:historical",
  stage = "raw"
)
```

## Arguments

- lpjml:

  Defines LPJmL version for crop/grass and natveg specific inputs

- climatetype:

  Switch between different GCM climate scenarios

- stage:

  Switch for raw data or harmonization

## Value

magpie object in cellular resolution

## Author

Kristine Karstens, Florian Humpenoeder

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("CarbonTests", aggregate = FALSE)
} # }
```
