# multiCSN 0.1.0

* The independent layered workflow has an explicit input and result contract.
  Candidate preparation and fitting require one joined `data` layer per assay;
  per-group normalization also requires one joined `counts` layer.
* Candidate preparation rebuilds the multiCSN catalogue. In-memory fitting
  returns the complete scientific tables without storing a Seurat object by default.
* Fitting returns valid empty endpoint tables when the peak detection threshold
  removes every candidate region.
* Constant TF expression produces empty TF layers with explicit accounting;
  the independent region-gene layer is retained.

* Checkpoint hashes use SHA-256 through `digest`, retaining support for R 4.1.
* Explicit peak separators are respected when normalizing chromatin assays,
  and normalized values remain aligned when peak names are canonicalized.
* Boolean endpoint validation returns `FALSE` for malformed result files and
  invalid summary counts.
* Require inferCSN 1.2.6 for the shared deletion-evidence ranking API.

* `prepare_candidates()` builds the unfiltered candidate catalogue and motif
  matches in one call; `find_motifs()` uses the bundled `motifs` and
  `motif2tf` by default.
* `fit_layered_network()` gains `preset = c("state", "timepoint")` for the two
  common fitting designs and `output_dir` to write the network tables after a
  complete fit, keeping resumable checkpoints beside the output until then.
* `validate_network_endpoints()` checks a saved result directory;
  `read_network_endpoints()` detects the file format and can return a single
  table; `classify_edge_transitions()` labels the state transition.
* New evaluation helpers: `evaluate_network()`, `ranking_metrics()`,
  `tie_aware_auroc()`, `tie_aware_average_precision()`,
  `expected_precision_at_k()`, `signed_ordinal_weight()`, `as_edge_table()`
  and `weight_matrix_to_edges()`.
* Peak coordinate parsing no longer calls the deprecated Signac
  `StringToGRanges()` and `GRangesToString()`.
* `create_seurat_object()` works with current Signac annotation checks.
* Removed `define_states_new()`, `assign_genes_to_states_new1()` and
  `assign_genes_to_states_simple()`, which duplicated `define_states()` and
  `assign_genes_to_states()`, and `JI_across_topregs()` and
  `biglist_compute_betweenness_degree()`, which failed on every call.
* Removed arguments that had no effect: `lag_fraction` and `lag_steps` from
  the `Network` and `Seurat` methods of `inferCSN()`, `pseudotime` from the
  `Network` method and `seed` from the `Seurat` method; the
  per-method arguments of the `fit_models()` generic; `label_points`,
  `max_overlaps` and `grid_n` from `plotPerturbation()`; `weight_cutoff` and
  `compare_random` from `plot_gene_rank()` for Seurat objects; `group.by`
  from `plot_rewiring()` and `plot_state_rewiring()`; `celltypes` from
  `metrics()` for `Network`; `method` from `compute_JI_topregs()`;
  `directed` from `compute_pagerank()`; `assay`
  from `get_pseudotime()` for Seurat.
* `fit_models()` passes `extend` to the GREAT peak-to-gene mapping.
* Fewer hard dependencies: `ComplexHeatmap`, `circlize`, `thisplot`,
  `tidygraph`, `ggraph`, `ggnetwork`, `patchwork` and `rlang` move to
  Suggests and are needed only by the plotting functions. Optional packages
  are installed on first use through `thisutils::check_r(strict = TRUE)`
  (thisutils >= 0.5.4), which stops only if installation fails; `purrr` is no longer used.
* `network_ordinal_weight()` uses `inferCSN::deletion_rank_weight()`.

# multiCSN 0.0.1

* Trajectory partitions, gene assignments, subnetworks, and heatmap annotations
  now use state terminology consistently, with default IDs such as state_1.
  Existing saved assignments are normalized when read by the new interfaces.

* Initial release.

* Removed the unused dynamic_genes and dynamic_genes_new interfaces, their
  legacy GAM implementation, and the gam dependency. Use
  scop::RunDynamicFeatures directly; state inference already uses this path.
* XGBoost fitting uses xgb.train with DMatrix input and maps nthread = -1 to
  the native all-thread setting, avoiding removed high-level API arguments.
