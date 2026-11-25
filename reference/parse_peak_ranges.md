# Parse peak identifiers

Parse peak identifiers

## Usage

``` r
parse_peak_ranges(peaks, sep = NULL)
```

## Arguments

- peaks:

  Character vector of peak identifiers.

- sep:

  Two non-empty literal strings separating chromosome from start and
  start from end, respectively. The default, `NULL`, recognizes
  `chromosome:start-end`, `chromosome-start-end`, and
  `chromosome_start_end` automatically.

## Value

Named genomic ranges.
