# calcPotentialForestArea

Calculates the area than can be potentially covered by forests, based on
environmental conditions.

## Usage

``` r
calcPotentialForestArea(
  refData = "lpj",
  countryLevel = FALSE,
  cells = "lpjcell",
  lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de"),
  climatetype = "MRI-ESM2-0:ssp370",
  grassyCorrection = TRUE
)
```

## Arguments

- refData:

  Determines the reference data that the estimated potential forest area
  is derived from (currently only "lpj")

- countryLevel:

  Whether output shall be at country level. Requires aggregate=FALSE in
  calcOutput.

- cells:

  magpiecell (59199 cells) or lpjcell (67420 cells)

- lpjml:

  Defines LPJmL version for crop/grass and natveg specific inputs. Only
  relevant, if refData = "lpj".

- climatetype:

  Switch between different GCM climate scenarios. Only relevant, if
  refData = "lpj".

- grassyCorrection:

  If TRUE, the potential forest area is reduced by the
  grassland-ecoregion cover fraction (RESOLVE 2017 biomes 7-10,
  calcGrassyEcoregions), correcting the LPJmL overestimation of forest
  in open grassy ecosystems.

## Value

magpie object in cellular resolution

## Author

Patrick v. Jeetze, Florian Humpenoeder

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("PotentialForestArea", aggregate = FALSE)
} # }
```
