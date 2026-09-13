# calcGrassyEcoregions

Per-cell areal cover fraction of the RESOLVE Ecoregions 2017 grassy
biomes (Dinerstein et al. 2017): tropical and subtropical, temperate,
flooded, and montane grasslands, savannas and shrublands (biomes 7-10).
Used to correct the LPJmL-derived potential forest area, which
overestimates forest cover in these open ecosystems, where tree planting
and forest expansion are ecologically inappropriate (Veldman et al.
2015).

## Usage

``` r
calcGrassyEcoregions()
```

## Value

magpie object in cellular resolution (67420 lpjcells)

## Author

Florian Humpenoeder

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("GrassyEcoregions", aggregate = FALSE)
} # }
```
