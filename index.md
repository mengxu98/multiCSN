# ***multiCSN***

## **Introduction**

[*`multiCSN`*](https://mengxu98.github.io/multiCSN/) is an R package for
***infer***ring ***C***ell-***S***pecific gene regulatory ***N***etwork
from single-cell omics data.

![multiCSN workflow
diagram](https://raw.githubusercontent.com/mengxu98/figures/main/multiCSN/multiCSN.svg#gh-light-mode-only)

![multiCSN workflow
diagram](https://raw.githubusercontent.com/mengxu98/figures/main/multiCSN/multiCSN-dark.svg#gh-dark-mode-only)

## **Installation**

Install the development version from
[*`GitHub`*](https://github.com/mengxu98/multiCSN) use
[*`pak`*](https://github.com/r-lib/pak):

``` r

if (!require("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pak("mengxu98/multiCSN")
```

## **Usage**

### **Examples**

#### For a `matrix` object.

``` r

library(multiCSN)
data("example_matrix", package = "inferCSN")

network <- inferCSN(
  example_matrix
)
```

#### For a `matrix` object, initiate it to `Network` object.

``` r

library(multiCSN)
data("example_matrix", package = "inferCSN")

object <- initiate_object(
  example_matrix
)
object
object <- inferCSN(
  object
)

network_table <- export_csn(
  object
)
head(network_table)
```

More functions and usages about
[*`multiCSN`*](https://mengxu98.github.io/multiCSN/)? Please reference
[*`here`*](https://mengxu98.github.io/multiCSN/reference/index.html).
