# Divides grnDF into epochs, filters interactions between genes not in same or consecutive epochs

Divides grnDF into epochs, filters interactions between genes not in
same or consecutive epochs

## Usage

``` r
epochGRN(grnDF, epochs, epoch_network = NULL)
```

## Arguments

- grnDF:

  result of GRN reconstruction

- epochs:

  result of running assign_epochs

- epoch_network:

  dataframe outlining higher level epoch connectivity (i.e. epoch
  transition network). If NULL, will assume epochs is ordered linear
  trajectory

## Value

list of GRNs across epochs and transitions
