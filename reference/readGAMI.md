# readGAMI

Read the Global Age Mapping Integration (GAMI) v2.1 forest-age dataset
(Besnard et al. 2024, GFZ Data Services, doi:10.5880/GFZ.1.4.2023.006),
0.5-degree \`class_fraction\` product. GAMI provides the
\*within-forest\* age distribution: for every grid cell the fraction of
the forested area in each of 12 age classes (sums to 1 where forest
exists), as a 20-member ensemble for 2010 and 2020. This read returns
the ensemble mean on GAMI's native 0.5-degree grid (720 x 360 cells) x 2
years x 12 age classes. Cleaning (fill / NA to 0) is done in
`correctGAMI`; the projection onto the 67420 lpj-cell grid, the
forest-area weighting and the mapping onto the 15 GFAD-style MAgPIE age
classes are done in `calcAgeClassDistribution`.

GAMI's \`forest_age\` is a five-dimensional variable (time, longitude,
latitude, age_class, 20-member ensemble). terra/GDAL flattens the three
non-spatial dimensions into unlabelled layers - every age class is
exposed as \`age_class=0\` and the member order is scrambled - so the
ensemble mean cannot be grouped reliably from the layer names. ncdf4
addresses the named dimensions explicitly and is therefore used here
instead of the usual terra reader.

## Usage

``` r
readGAMI()
```

## Value

magpie object on GAMI's native 0.5-degree grid: cells x 2 years (2010,
2020) x 12 age classes

## See also

[`calcAgeClassDistribution`](calcAgeClassDistribution.md), `readGFAD`

## Author

Florian Humpenoeder

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("GAMI", convert = "onlycorrect")
} # }
```
