# fullCELLULARMAGPIE

Function that produces the complete cellular data set required for
running the MAgPIE model.

## Usage

``` r
fullCELLULARMAGPIE(
  rev = numeric_version("0.1"),
  dev = "",
  ctype = "c200",
  climatetype = "MRI-ESM2-0:ssp370",
  lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de", crop =
    "ggcmi_phase3_nchecks_9ca735cb", grass = "lpjml5p2_pasture"),
  isimip = NULL,
  clusterweight = NULL,
  emu_id = NULL
)
```

## Arguments

- rev:

  data revision which should be used as input (numeric_version).

- dev:

  development suffix to distinguish development versions for the same
  data revision. This can be useful to distinguish parallel lines of
  development.

- ctype:

  aggregation clustering type, which is a combination of a single
  letter, indicating the cluster methodology, and a number, indicating
  the number of resulting clusters. Available methodologies are -
  hierarchical clustering (h), - normalized k-means clustering (n) and -
  combined hierarchical/normalized k-means clustering (c). In the latter
  hierarchical clustering is used to determine the cluster distribution
  among regions whereas normalized k-means is used for the clustering
  within a region.

- climatetype:

  Global Circulation Model to be used

- lpjml:

  Defines LPJmL version for crop/grass and natveg specific inputs

- isimip:

  Defines isimip crop model input which replace maiz, tece, rice_pro and
  soybean

- clusterweight:

  Should specific regions be resolved with more or less detail? Values
  \> 1 mean higher share, \< 1 lower share e.g. cfg\$clusterweight \<-
  c(LAM=2) means that a higher level of detail for region LAM if set to
  NULL all weights will be assumed to be 1. Examples:
  c(LAM=1.5,SSA=1.5,OAS=1.5) or c(LAM=2,SSA=2,OAS=2)
  [`setConfig`](https://rdrr.io/pkg/madrat/man/setConfig.html) (e.g. for
  setting the mainfolder if not already set properly).

- emu_id:

  Pasture Soil carbon emulator ID

## See also

[`readSource`](https://rdrr.io/pkg/madrat/man/readSource.html),[`getCalculations`](https://rdrr.io/pkg/madrat/man/getCalculations.html),[`calcOutput`](https://rdrr.io/pkg/madrat/man/calcOutput.html),
[`setConfig`](https://rdrr.io/pkg/madrat/man/setConfig.html)

## Author

Kristine Karstens, Jan Philipp Dietrich

## Examples

``` r
if (FALSE) { # \dontrun{
retrieveData("CELLULARMAGPIE", rev = numeric_version("12"),
             mainfolder = "pathtowhereallfilesarestored")
} # }
```
