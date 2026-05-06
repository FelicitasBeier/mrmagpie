# calcTRENDY

Resample and use inverse distance weighted interpolation to align
natural vegetation carbon densities from TRENDY with MAgPIE's input
requirements

## Usage

``` r
calcTRENDY(subtype)
```

## Arguments

- subtype:

  name of the TRENDY dataset to read. Current options includes:
  "CABLEPOP", "CARDAMOM", "CLASSIC", "CLM5.0", "DLEM", "ELM", "IBIS",
  "ISBACTRIP", "JSBACH", "LPJ-GUESS", "LPJml", "LPJwsl", "lpxqs", "OCN",
  "ORCHIDEE", "SDGVM", "VISIT"

## Value

a magclass object with model, carbon pool in third dimension

## Author

Michael Crawford

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("TRENDY", subtype = "JSBACH")
} # }
```
