# readBendingTheCurve

Read bending the curve data

## Usage

``` r
readBendingTheCurve(subtype)
```

## Arguments

- subtype:

  Data used in the Bending the Curve initiative. Type "rr_layer" for the
  range-size rarity layer and "luh2_side_layers" for the LUH2 Side
  Layers.

## Value

List of magpie objects with results on cellular level, weight, unit and
description.

## Author

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
  readSource("BendingTheCurve", subtype="rr_layer", convert="onlycorrect")
} # }
```
