# correctDinerstein2017

Sets missing values (cells outside any ecoregion of the selected biome
group) to zero.

## Usage

``` r
correctDinerstein2017(x)
```

## Arguments

- x:

  magpie object provided by
  [`readDinerstein2017`](readDinerstein2017.md)

## Value

magpie object on cellular level

## See also

[`readDinerstein2017`](readDinerstein2017.md)

## Author

Florian Humpenoeder

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("Dinerstein2017", subtype = "grassland", convert = "onlycorrect")
} # }
```
