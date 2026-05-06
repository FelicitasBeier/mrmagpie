# readTravelTimeNelson2019

Read minimum travel time to cities and ports and ports of various size,
see metadata file in source folder

## Usage

``` r
readTravelTimeNelson2019(subtype = "cities50")
```

## Arguments

- subtype:

  currently only cities of 5, 20, or 50 thousand people ("cities5",
  "cities20", "cities50") or ports of various sizes
  ("portsLarge\|Medium\|Small\|VerySmall\|Any")

## Value

gridded magpie object for 2015, minimum travel time to cities in minutes

## Author

David M Chen
