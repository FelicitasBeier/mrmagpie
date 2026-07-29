# correctGAMI

Correct the GAMI v2.1 forest-age dataset: replace missing and negative
class-fraction values (ocean / non-forest fill) by 0. The object stays
on GAMI's native 0.5-degree grid; the projection onto the 67420 lpj-cell
grid happens in `calcAgeClassDistribution`.

## Usage

``` r
correctGAMI(x)
```

## Arguments

- x:

  magpie object provided by `readGAMI`

## Value

magpie object on GAMI's native 0.5-degree grid, cleaned

## See also

[`readGAMI`](readGAMI.md)

## Author

Florian Humpenoeder

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("GAMI", convert = "onlycorrect")
} # }
```
