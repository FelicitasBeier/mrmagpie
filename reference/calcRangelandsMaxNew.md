# calcRangelandsMaxNew

Calculates rangelands maximum output

## Usage

``` r
calcRangelandsMaxNew(
  lsuLevels = c(seq(0, 2.2, 0.2), 2.5),
  lpjml = "lpjml5p2_pasture",
  climatetype = "MRI-ESM2-0:ssp370",
  scenario = "/co2/Nreturn0p5/limN",
  report = "harvest",
  cells = "lpjcell"
)
```

## Arguments

- lsuLevels:

  Livestock unit levels in the source folder

- lpjml:

  Defines LPJmL version for crop/grass and natveg specific inputs

- climatetype:

  Switch between different climate scenarios (default: "CRU_4")

- scenario:

  specify ssp scenario

- report:

  Either 'harvest' or 'lsu/ha' controlling what values are output by the
  function.

- cells:

  "magpiecell" for 59199 cells or "lpjcell" for 67420 cells

## Value

magpie object in cellular resolution

## Author

Marcos Alves
