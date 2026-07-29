# calcAgeClassDistribution

Forest area in 15 age classes. With `dataset = "GFAD"` (default) the
distribution comes directly from the Global Forest Age Dataset (GFAD
v1.1, Poulter et al. 2019): per-cell class fractions times cell area,
summed over the four forest types. With `dataset = "GAMI"` the \*age
distribution\* comes from GAMI v2.1 (Besnard et al. 2024) while the
per-cell forest \*area\* is taken from the GFAD product. GAMI supplies
only the within-forest age shares (12 classes, summing to 1); keeping
the GFAD forest-area map as the spatial weight conserves each cell's
total forest area and isolates the pure age-structure change. GAMI's 12
twenty-year classes are mapped onto the 15 GFAD ten-year classes (each
split 50/50 into two adjacent classes; all classes \>=140 yr, i.e. GAMI
\`140-160 ... \>299\`, collapse into \`class15\`, the mature bin routed
to \`acx\` in module 28). Where GAMI reports no forest but GFAD does,
the GFAD age shape is retained (area-conserving).

## Usage

``` r
calcAgeClassDistribution(
  cells = "lpjcell",
  dataset = "GFAD",
  gamiYear = "y2010"
)
```

## Arguments

- cells:

  lpjcell for 67420 cells

- dataset:

  "GFAD" (default) or "GAMI"

- gamiYear:

  GAMI reference year when `dataset = "GAMI"`: "y2010" (default, matches
  the GFAD vintage) or "y2020"

## Value

magpie object in cluster resolution

## Author

Abhijeet Mishra, Felicitas Beier, Florian Humpenoeder

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("AgeClassDistribution", aggregate = FALSE)
calcOutput("AgeClassDistribution", dataset = "GAMI", aggregate = FALSE)
} # }
```
