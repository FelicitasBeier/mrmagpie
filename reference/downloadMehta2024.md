# downloadMehta2024

download Global Area Equipped for Irrigation Dataset 1900-2015 from
Mehta et al. (2024). Gridded dataset is created based on (sub-)national
statistics from FAOSTAT, AQUASTAT, EUROSTAT and country's census data
downscaled using two alternative gridded irrigation maps (GMIA from
Siebert et al. 2013 and Meier et al. 2018)

## Usage

``` r
downloadMehta2024(subtype = "v4_GMIA")
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

\[downloadSource()\] \[readMehta2024()\]

## Author

Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
a <- downloadSource()
} # }
```
