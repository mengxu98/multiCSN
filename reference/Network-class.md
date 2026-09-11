# The Network class

The Network object stores the inferred network itself, information about
the fitting process as well as graph representations of the network.

## Slots

- `data`:

  A matrix.

- `modules`:

  A list TF modules.

- `regulators`:

  A named list containing the transcription factors included in the
  network.

- `targets`:

  A named list containing the targets included in the network.

- `metrics`:

  A dataframe with goodness of fit measures.

- `coefficients`:

  A dataframe with the fitted coefficients.

- `graphs`:

  Graphical representations of the inferred network.

- `network`:

  A dataframe containing the network edges and their properties.

- `params`:

  A named list with GRN inference parameters.

- `perturbation`:

  A `Perturbation` object containing perturbation analysis results.
