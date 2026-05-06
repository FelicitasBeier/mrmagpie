# calcEnvmtlFlow

This function calculates environmental flow requirements (EFR) for
MAgPIE retrieved from LPJmL monthly discharge and water availability

## Usage

``` r
calcEnvmtlFlow(
  lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de", crop =
    "ggcmi_phase3_nchecks_9ca735cb"),
  climatetype = "GSWP3-W5E5:historical",
  stage = "harmonized2020",
  seasonality = "grper"
)
```

## Arguments

- lpjml:

  Defines LPJmL version for crop/grass and natveg specific inputs

- climatetype:

  Switch between different climate scenarios

- stage:

  Degree of processing: raw, smoothed, harmonized, harmonized2020

- seasonality:

  grper (default): EFR in growing period per year; total: EFR throughout
  the year; monthly: monthly EFRs

## Value

magpie object in cellular resolution

## Author

Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("EnvmtlFlow", aggregate = FALSE)
} # }
```
