# calcTransportCosts

calculates country-level transport costs from GTAP total transport
costs, cellular production, and cellular travel time

## Usage

``` r
calcTransportCosts(transport = "all", gtapVersion = "9")
```

## Arguments

- transport:

  "all" or "nonlocal". "all" means all production incurs transport
  costs, while "nonlocal" sees only production greater than local rural
  consumption with transport costs

- gtapVersion:

  "9" or "81"

## Value

List of magpie objects with results on country level, weight on country
level, unit and description.

## See also

\[calcTransportTime()\], \[calcGTAPTotalTransportCosts()\]

## Author

David M Chen

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("TransportCosts_new")
} # }
```
