# Internal helper to normalize regulators/targets specification into a per-celltype named list.

Internal helper to normalize regulators/targets specification into a
per-celltype named list.

## Usage

``` r
parse_regulators_targets(x, celltypes, arg_name = "targets")
```

## Arguments

- x:

  Regulators/targets specification: a character vector, or a named list
  of character vectors with names matching `celltypes`.

- celltypes:

  Cell types the specification is expanded to.

- arg_name:

  Argument name used in error messages.
