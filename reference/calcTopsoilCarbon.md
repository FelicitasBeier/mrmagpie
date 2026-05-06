# calcTopsoilCarbon

This function extracts topsoil carbon densities from LPJ to MAgPIE

## Usage

``` r
calcTopsoilCarbon(
  cells = "lpjcell",
  lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de", crop =
    "ggcmi_phase3_nchecks_9ca735cb"),
  climatetype = "GSWP3-W5E5:historical"
)
```

## Arguments

- cells:

  "magpiecell" for 59199 cells or "lpjcell" for 67420 cells

- lpjml:

  Defines LPJmL version for crop/grass and natveg specific inputs

- climatetype:

  Switch between different GCM climate scenarios

## Value

magpie object in cellular resolution

## Author

Kristine Karstens

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("TopsoilCarbon", aggregate = FALSE)
} # }
```
