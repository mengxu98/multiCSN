# Changelog

## multiCSN 0.0.1

- Trajectory partitions, gene assignments, subnetworks, and heatmap
  annotations now use state terminology consistently, with default IDs
  such as state_1. Existing saved assignments are normalized when read
  by the new interfaces.

- Initial release.

- Removed the unused dynamic_genes and dynamic_genes_new interfaces,
  their legacy GAM implementation, and the gam dependency. Use
  scop::RunDynamicFeatures directly; state inference already uses this
  path.

- XGBoost fitting uses xgb.train with DMatrix input and maps nthread =
  -1 to the native all-thread setting, avoiding removed high-level API
  arguments.
