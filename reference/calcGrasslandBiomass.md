# calcGrasslandBiomass

Calculates pasture biomass demand for the historical period split
between rangelands andmanaged pastures.

## Usage

``` r
calcGrasslandBiomass(cells = "lpjcell")
```

## Arguments

- cells:

  "magpiecell" for 59199 cells or "lpjcell" for 67420 cells

## Value

Regional biomass demand

## See also

[`calcOutput`](https://rdrr.io/pkg/madrat/man/calcOutput.html),
[`calcFAOmassbalance`](https://rdrr.io/pkg/mrcommons/man/calcFAOmassbalance.html),
[`readSource`](https://rdrr.io/pkg/madrat/man/readSource.html)

## Author

Marcos Alves

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("GrasslandBiomass")
} # }
```
