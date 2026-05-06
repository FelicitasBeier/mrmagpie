# calcFoodDemandGridded

Calculates grid-level food demand, note also includes food and feed

## Usage

``` r
calcFoodDemandGridded(
  attribute = "dm",
  prod = "k",
  feed = TRUE,
  cells = "lpjcell"
)
```

## Arguments

- attribute:

  dm or calories ("ge") or other massbalance attribute

- prod:

  for memory reasons

- feed:

  whether to include feed demand in the gridded demand

- cells:

  magpiecell or lpjcell (default)

## Value

Gridded magpie object of food demand disaggregated by rural urban pop

## Author

David M Chen

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("FoodDemandGridded")
} # }
```
