# calcAfforestationMask

Afforestation mask for where afforestation possible

## Usage

``` r
calcAfforestationMask(subtype, cells = "lpjcell")
```

## Arguments

- subtype:

  afforestation mask sub type

- cells:

  "magpiecell" or "lpjcell"

## Value

magpie object in cellular resolution

## Author

David Chen, Florian Humpenoeder

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("AfforestationMask", subtype = "noboreal", aggregate = FALSE)
} # }
```
