# The logo of multiCSN

The multiCSN logo, using ASCII or Unicode characters Use
[cli::ansi_strip](https://cli.r-lib.org/reference/ansi_strip.html) to
get rid of the colors.

## Usage

``` r
multicsn_logo(unicode = cli::is_utf8_output())
```

## Arguments

- unicode:

  Unicode symbols on UTF-8 platforms. Default is
  [cli::is_utf8_output](https://cli.r-lib.org/reference/is_utf8_output.html).

## Value

A character vector with class `multicsn_logo`.

## References

<https://github.com/tidyverse/tidyverse/blob/main/R/logo.R>

## Examples

``` r
multicsn_logo()
#>           ⬢          .        ⬡             ⬢     .
#>                            ____  _  ___________ _   __
#>           ____ ___  __  __/ / /_(_)/ ____/ ___// | / /
#>          / __ `__ ./ / / / / __/ // /    .__ ./  |/ /
#>         / / / / / / /_/ / / /_/ // /___ ___/ / /|  /
#>        /_/ /_/ /_/.__,_/_/.__/_/ .____//____/_/ |_/
#>       ⬡               ⬢      .        ⬡          ⬢
```
