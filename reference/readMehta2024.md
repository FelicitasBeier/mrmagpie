# readMehta2024

reads in Global Area Equipped for Irrigation for years 1900-2015 from
Mehta et al. (2022)

## Usage

``` r
readMehta2024(subtype = "v4_GMIA")
```

## Arguments

- subtype:

  data subtype to be downloaded Subtypes available: Combination of
  version (v3, v4) and data source: 'GMIA': gridded base map for
  downscaling from Stefan et al. (2013). Global Map of Irrigation Areas
  version 5. 'Meier2018': gridded base map for downscaling from Meier,
  et al. (2018). Global Irrigated Areas. Separated by "\_" (e.g.,
  subtype = "v4_GMIA")

## See also

\[correctMehta2024()\]

## Author

Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
a <- readSource("Mehta2024")
} # }
```
