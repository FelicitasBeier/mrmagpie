# readDinerstein2017

Reads the RESOLVE Terrestrial Ecoregions 2017 (Dinerstein et al. 2017)
and returns the areal cover fraction of a biome group per 0.5 degree
cell.

## Usage

``` r
readDinerstein2017(subtype = "grassland")
```

## Arguments

- subtype:

  Biome group whose areal cover fraction is returned. Currently only
  `"grassland"` = RESOLVE biomes 7 (tropical and subtropical grasslands,
  savannas and shrublands), 8 (temperate grasslands, savannas and
  shrublands), 9 (flooded grasslands and savannas) and 10 (montane
  grasslands and shrublands).

## Value

magpie object (67420 cells) with the cover fraction of the selected
biome group. Cells outside any ecoregion of that group are NA and set to
zero in [`correctDinerstein2017`](correctDinerstein2017.md).

## See also

[`correctDinerstein2017`](correctDinerstein2017.md)

## Author

Florian Humpenoeder

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("Dinerstein2017", subtype = "grassland", convert = "onlycorrect")
} # }
```
