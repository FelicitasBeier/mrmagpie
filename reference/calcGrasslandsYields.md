# calcGrasslandsYields

Calculates rangelands maximum output and managed pastures yields

## Usage

``` r
calcGrasslandsYields(
  lpjml = "lpjml5p2_pasture",
  climatetype = "MRI-ESM2-0:ssp370",
  cells = "lpjcell",
  subtype = "/co2/Nreturn0p5",
  lsu_levels = c(seq(0, 2, 0.2), 2.5),
  past_mngmt = "mdef"
)
```

## Arguments

- lpjml:

  Defines LPJmL version for crop/grass and natveg specific inputs

- climatetype:

  Global Circulation Model to be used

- cells:

  "magpiecell" for 59199 cells or "lpjcell" for 67420 cells

- subtype:

  Switch between different climate scenarios

- lsu_levels:

  Livestock unit levels in the source folder

- past_mngmt:

  pasture areas management option

## Value

magpie object in cellular resolution

## Author

Marcos Alves

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("GrasslandsYields")
} # }
```
