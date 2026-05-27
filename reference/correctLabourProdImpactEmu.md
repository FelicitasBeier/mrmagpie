# correctLabourProdImpactEmu

correct labour productivity impacts from climate change emulated by the
LAMACLIMA project

based on method of Orlov et al. 2019. Economics of Disasters and Climate
Change, 3(3), 191-211.

## Usage

``` r
correctLabourProdImpactEmu(x)
```

## Arguments

- x:

  magpie object provided by the read function

## Value

List of magpie objects with results on cellular level, weight, unit and
description.

## See also

[`readLabourProdImpactEmu`](readLabourProdImpactEmu.md)

## Author

Michael Windisch

## Examples

``` r
if (FALSE) { # \dontrun{
  readSource("LabourProdImpactEmu", convert="onlycorrect")
} # }
```
