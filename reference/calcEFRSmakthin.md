# calcEFRSmakthin

This function calculates environmental flow requirements (EFR) for
MAgPIE retrieved from LPJmL monthly discharge and water availability
using the method of Smakthin et al. (2004)

## Usage

``` r
calcEFRSmakthin(
  lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de", crop =
    "ggcmi_phase3_nchecks_9ca735cb"),
  climatetype = "GSWP3-W5E5:historical",
  stage = "harmonized2020",
  LFR_val = 0.1,
  HFR_LFR_less10 = 0.2,
  HFR_LFR_10_20 = 0.15,
  HFR_LFR_20_30 = 0.07,
  HFR_LFR_more30 = 0,
  seasonality = "grper",
  cells = "lpjcell"
)
```

## Arguments

- lpjml:

  Defines LPJmL version for crop/grass and natveg specific inputs

- climatetype:

  Switch between different climate scenarios

- stage:

  Degree of processing: raw, smoothed, harmonized, harmonized2020

- LFR_val:

  Strictness of environmental flow requirements

- HFR_LFR_less10:

  High flow requirements (share of total water for cells) with
  LFR\<10percent of total water

- HFR_LFR_10_20:

  High flow requirements (share of total water for cells) with 10percent
  \< LFR \< 20percent of total water

- HFR_LFR_20_30:

  High flow requirements (share of total water for cells) with 20percent
  \< LFR \< 30percent of total water

- HFR_LFR_more30:

  High flow requirements (share of total water for cells) with
  LFR\>30percent of total water

- seasonality:

  grper (default): EFR in growing period per year; total: EFR throughout
  the year; monthly: monthly EFRs

- cells:

  lpjcell for 67420 cells or magpiecell for 59199 cells

## Value

magpie object in cellular resolution

## Author

Felicitas Beier, Abhijeet Mishra

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("EFRSmakthin", aggregate = FALSE)
} # }
```
