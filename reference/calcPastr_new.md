# calcPastr_new

Calculates managed pasture yields

## Usage

``` r
calcPastr_new(
  past_mngmt = "me2",
  lpjml = "lpjml5p2_pasture",
  climatetype = "MRI-ESM2-0:ssp370",
  scenario = "/co2/Nreturn0p5/limN",
  cells = "lpjcell"
)
```

## Arguments

- past_mngmt:

  pasture areas management option

- lpjml:

  Defines LPJmL version for crop/grass and natveg specific inputs

- climatetype:

  Switch between different climate scenarios (default: "CRU_4")

- scenario:

  specify ssp scenario

- cells:

  "magpiecell" for 59199 cells or "lpjcell" for 67420 cells

## Value

magpie object in cellular resolution

## Author

Marcos Alves

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("Pastr_new", past_mngmt = "me2", lpjml = "LPJml_pastr", climatetype)
} # }
```
