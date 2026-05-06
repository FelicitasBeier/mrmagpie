# calcGridPop

Past and future (SSP1-5) population based on HYDE3.2 and Jones & O'Neill
(2016) Data is scaled to match WDI data from calcPopulation NOTE that
some scaling factors for the future (for small countries Gambia and
Djibouti) are off, data read in is 50

## Usage

``` r
calcGridPop(
  source = "ISIMIP",
  subtype = "all",
  cellular = TRUE,
  cells = "lpjcell",
  FiveYear = TRUE,
  scale = TRUE,
  harmonize_until = 2015,
  urban = FALSE
)
```

## Arguments

- source:

  default source (ISIMIP) or Gao data (readGridPopGao) which is split
  into urban and rural.

- subtype:

  time horizon to be returned. Options: past (1965-2005), future
  (2005-2010) or all (divergence starts at year in harmonize_until)

- cellular:

  if true: half degree grid cell data returned

- cells:

  number of halfdegree grid cells to be returned. Options: "magpiecell"
  (59199), "lpjcell" (67420)

- FiveYear:

  TRUE for 5 year time steps, otherwise yearly from source

- scale:

  if true: scales sum of gridded values to match country level totals

- harmonize_until:

  harmonization year until which SSPs diverge (default: 2015)

- urban:

  TRUE to return only urban gridded population based on iso share

## Value

Population in millions.

## Author

David Chen, Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("GridPop", aggregate = FALSE)
} # }
```
