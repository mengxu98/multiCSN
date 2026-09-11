# Adds interaction type to dynamic differential network

Adds interaction type to dynamic differential network

## Usage

``` r
add_type(diffnet, type, grnDF_on, grnDF_offlist)
```

## Arguments

- diffnet:

  diffnet

- type:

  "on" or "off" depending on the type of differential network. If "on"
  will assign type based on grnDF_on. Otherwise interaction assigned
  from grnDF_offlist.

- grnDF_on:

  the static network in which diffnet edges are active

- grnDF_offlist:

  list of static networks in which diffnet edges are inactive

## Value

community assignments of nodes in the dynamic network
