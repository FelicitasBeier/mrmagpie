# calcCarbon

This function extracts carbon densities from LPJ to MAgPIE

## Usage

``` r
calcCarbon(
  lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de", crop =
    "ggcmi_phase3_nchecks_9ca735cb"),
  climatetype = "GSWP3-W5E5:historical",
  cells = "lpjcell"
)
```

## Arguments

- lpjml:

  Defines LPJmL version for crop/grass and natveg specific inputs

- climatetype:

  Switch between different GCM climate scenarios

- cells:

  "magpiecell" for 59199 cells or "lpjcell" for 67420 cells

## Value

magpie object in cellular resolution

## Author

Kristine Karstens, Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("Carbon", aggregate = FALSE)
} # }
```
