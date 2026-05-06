# calcTransportTime

Function extracts travel time to major cities in minutes

## Usage

``` r
calcTransportTime(subtype = "cities50", cells = "lpjcell")
```

## Arguments

- subtype:

  currently only cities of 5, 20, or 50 thousand people ("cities5",
  "cities20", "cities50") or ports of various sizes
  ("portsLarge\|Medium\|Small\|VerySmall\|Any")

- cells:

  number of cells to be returned: magpiecell (59199), lpjcell (67420)

## Value

magpie object in cellular resolution

## Author

David Chen

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("TransportTime", aggregate = FALSE)
} # }
```
