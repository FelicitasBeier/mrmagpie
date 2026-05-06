# readTRENDY

read a TRENDY dataset

## Usage

``` r
readTRENDY(subtype)
```

## Arguments

- subtype:

  name of the TRENDY dataset to read. Current options include:
  "CABLEPOP", "CLASSIC", "CLM5.0", "DLEM", "ED", "ELM", "IBIS", "ISAM",
  "ISBACTRIP", "JSBACH", "JULES", "LPJ-GUESS", "LPJml", "LPJwsl",
  "lpxqs", "OCN", "ORCHIDEE", "SDGVM", "VISIT", "YIBS"

## Value

a named list of terra::rast objects

## Author

Michael Crawford

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("TRENDY", subtype = "JSBACH")
} # }
```
