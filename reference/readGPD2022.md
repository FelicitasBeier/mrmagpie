# readGPD2022

read x Data from the Global Peatland Database provided by Alexandra
Barthelmes. The original xls file has been clean-up manually (country
names). Turkey had two identical entries in the original xls file.
Sources: "Inventory Reports and National Communications UNFCC 2014",
"soil and peatland science", "European Mires Book" , "own estimates
(incl. GIS data)",

## Usage

``` r
readGPD2022()
```

## Value

List of magpie objects with results on cellular level, weight, unit and
description.

## Author

Florian Humpenoeder

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("x", convert = "onlycorrect")
} # }
```
