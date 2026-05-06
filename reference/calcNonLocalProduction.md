# calcNonLocalProduction

Calculates grid-level amount of food that would need to be transported,
assuming that food produced in the grid cell is first consumed by local
population i.e. amount of food greater than local rural demand, split
into that which feeds the local urban population, and that which exceeds
total local demand and is available to export

## Usage

``` r
calcNonLocalProduction(cells = "lpjcell")
```

## Arguments

- cells:

  magpiecell or lpjcell (default)

## Author

David M Chen

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("NonLocalTransport")
} # }
```
