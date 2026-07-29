# downloadGAMI

Downloads the Global Age Mapping Integration (GAMI) v2.1 forest-age
dataset (Besnard et al. 2024, GFZ Data Services,
doi:10.5880/GFZ.1.4.2023.006), the 0.5-degree \`class_fraction\` product
read by [`readGAMI`](readGAMI.md). GAMI provides the within-forest age
distribution (12 age classes, 20-member ensemble) for 2010 and 2020.

## Usage

``` r
downloadGAMI()
```

## Value

metadata list describing the downloaded GAMI source

## See also

\[downloadSource()\], [`readGAMI`](readGAMI.md)

## Author

Florian Humpenoeder

## Examples

``` r
if (FALSE) { # \dontrun{
downloadSource("GAMI")
} # }
```
