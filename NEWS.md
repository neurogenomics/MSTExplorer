# MSTExplorer 1.0.5

## Bug fixes

* Tests
  - Use `force_new=TRUE` where some tests occasionally fail with cached
  files.
  - `test-load_example_results`: Update test files.
  - `test-prioritise_targets`: Remove unused arguments and change input size.
  - `test-prioritise_targets_network`: Process `top_targets` to include effect
  - `test-plot_differential_outcomes`: Use non-specific plot name in
  `patchwork::wrap_plots`.
  - `test-plot_differential_outcomes`: Wrap p3 in `expect_error` to prevent
  test failure even if error was handled.
  variable.
* Vignettes
  - `MSTExplorer`: Update effect variable to `fold_change`.
* `add_logfc`
  - Return `results` with new column rather than directly modifying the original
  input.
  - Update references (`add_logfc(results)` -> `results <- add_logfc(results)`).
* `ttd_check`, `plot_differential_outcomes`
  - Add check for disease column names before executing, if required.
  `HPOExplorer::add_disease` on input.
* `plot_ttd`
  - Remove `fill` aesthetic for `geom_text` (doesn't exist anymore).
* `extract_help`
  - [DEVELOPEMENT ONLY] Look for help docs only in legitimate pkg installation
  paths.
* `subset_results`
  - Add new `effect_var` argument.
  - Adjust default `effect_threshold` to 0.1.
* `add_symptom_results`
  - Only merge `results` and `phenotypes_to_genes` if required (prevents
  column duplicates with altered names).
* Add missing import: simona

# MSTExplorer 1.0.4

## New features

* `run_phenomix`: can now take multiple `annotLevel`s at once (will iterate).

# MSTExplorer 1.0.3

## New features

* add function `run_congenital_enrichment`
* add function `plot_celltype_severity`

# MSTExplorer 1.0.2

## New features

* add function `prioritise_targets_grid`

# MSTExplorer 1.0.1

## Bug fixes

* `plot_bar_dendro_facets`
  - Add `hpo` arg to pass HPO directly.
  - Fix `data_summary` so that denominator is the number of on-target cell types only.

# MSTExplorer 1.0.0

## New features

* Rename `MSTExplorer`.
* Revamp to work with `KGExplorer`.

# MultiEWCE 0.1.10

## New features

* `gene_results`
  - Return merged versions of `$results` as well as `$gene_data` 
    (when available).

# MultiEWCE 0.1.9

## New features

* New function: `get_bg`:
  - Creates and caches background made with `gprofiler`.
* New function: `standardise_genes`
  - May move this to `orthogene` package as it's quite generally useful.
* Allow users to set min number of genes:
  - `min_genes` arg in `gen_results` / `ewce_para`.

## Bug fixes

* `gen_results` / `ewce_para`
  - `bg` was incorrectly set to use only genes in `gene_data`.
  - Now uses `get_bg` to create background using *gprofiler*.
* `get_valid_gene_lists`
  - Throw error when 0 valid gene lists found.
* `get_unfinished_list_names`
  - Fix example
* Use *rworkflows@dev* to avoid API limit.

# MultiEWCE 0.1.8

## New features

* `MultiEWCE` finally gets a hex sticker!
* `gen_results` / `gen_overlap`
  - Check for existing results and import if already there.
  - Name all results "gen_results.rds" or "gen_overlap.rds" to avoid 
    rerunning duplicate analyses on HPC.
* Update *rworkflows.yml*
    
## Bug fixes

* Fix unit tests and examples to use "hpo_id" instead of "hpo_name".
* `load_hpo_graph`: 
  - export
  - Regenerate and update "hpo_graph.rds" file.
* Drastically reduce time to run examples.
* `ontology_plot`
  - Fix function and add test.
* *DESCRIPTION*
  - `Depends: R (>= 2.10)` --> `Depends: R (>= 4.1)`, 
    to ensure `|>` function available.
  - Rewrite `Description` field to reflect `MultiEWCE`'s current purpose.

# MultiEWCE 0.1.7

## New features

* Update to coordinate with `HPOExplorer` updates.
* New funcs:
  - `plot_ontology_levels`
  - `ttd_check`
  - `ttd_plot`
  - `ttd_import`

# MultiEWCE 0.1.6

## New features

* Require `HPOExplorer (>= 0.99.10)`
* Switch terms:
  - "HPO_ID" --> "hpo_id"
  - "Phenotype" --> "hpo_name"
  - "FREQUENCY" --> "gene_freq"
  - "Onset" --> "onset"
  - "Modifier" --> "modifier"
  - "Aspect" --> "aspect"
  - "Gene" --> "gene_symbol"
  - "DatabaseID" --> "disease_id" 
  - "LinkID" --> "disease_id"
* Update all "data" objects.
* `get_data`
  - Add `tag` arg.

## Bug fixes

* `load_example_results`
  - Change `tag` to "latest".
  - Update colnames dynamically.
* `load_example_ctd`
  - Change `tag` to "latest".
* `load_hpo_graph`
  - Change `tag` to "latest".
* `map_tissues`
  - Fix docs.
* `agg_results`
  - Add "hpo_id" column before aggregation.
* `correlation_heatmap`
  - Ensure row annot order is correct.
* `add_ctd`
  - Fix example.

# MultiEWCE 0.1.5

## New features

* `prioritise_targets`
  - New arg `evidence_score_threshold` to utilise GenCC evidence scores provided via `HPOExplorer`.
  - Merge symptom-level driver genes after unlisting the "intersection" column
    instead of using indirect approach that on nested "intersection" column.
* `prioritise_targets_network`
  - Add feature to resize plot after double-click.
* New function: `add_tissues`
  - Add tissues that each cell type is found in, using a celltype-tissue mapping file.
  - Regenerate "DescartesHuman_celltype_mapping.csv" file as the one previously generated wrong
  (had every celltype in every tissue)
  
## Bug fixes

* `prioritise_targets_network`
  - Fill screen better in saved plots with new default args: 
    `width = "100%", height = "90vh"`

# MultiEWCE 0.1.4

## New features

* `prioritise_targets`
  - Move filtering steps and arg docs inside respective `HPOExplorer::add_*` functions.
* `gen_results`
  - New arg `parallel_boot` lets users choose which level to parallelise over.
* New exported func: `gen_overlap`
  - Compute simple overlap enrichment results. 
  - Overcomes requirement of >=4 genes.
  - Faster than EWCE (but less robust)
* New internal func: `save_results`
* `load_example_results`
  - New arg `force_new`

## Bug fixes

* `get_valid_gene_lists`
  - Can now handle any `list_name_column`
  
# MultiEWCE 0.1.3

## New features

* `ggnetwork_plot_full`:
  - Handles multiple celltypes by aggregation.
* New internal function: `agg_results`.
* New functions for generating/visualizing prioritised targets:
  - `prioritise_targets`
  - `prioritise_targets_heatmap`
  - `prioritise_targets_network`
  - `plot_report`
  - `correlation_heatmap`
  - `frequency_bar`
  - `plot_frequency_histogram` 
  
* Added `example_targets` data:
  - Speeds up run time of examples/tests.
  

## Bug fixes

* Remove unused function: `count_results`
* Make `get_unfinished_list_names` much more efficient.
* Add leeway to `gen_results` test, as the number of sig results 
  varies from run to run.


# MultiEWCE 0.1.2

## New features

* New functions:
  - `load_example_results`
  - `ontology_plot`
  - `terminal_celltypes`
  - various support functions
* Add `ontologyPlot` as new *Import*. 
* Was able to replicate Momoko's results!!!

## Bug fixes

* `load_example_ctd`
  - Pass `file` to `piggyback::pb_download`
* Remove redundant `get_gene_list` function 
  (now handled by `HPOExplorer::get_gene_lists`).
* Make `get_valid_gene_lists` much more efficient and consider 
  the intersect between ctd/gene_data genes.
* `merge_results`:
  - Handle `NULL` results in list.

# MultiEWCE 0.1.1

## New features

* Added a `NEWS.md` file to track changes to the package.
* Added `rworkflows`
* `ewce_para`
  - Harmonize arguments with `EWCE`
  - Set default `list_names` arg
  - Output named list of saved files instead of `TRUE`.
  - Add messages that appear in parallel. Make EWCE message silent. 
* `ewce_plot`:
  - Simply make this a shallow wrapper for `EWCE::ewce_plot` 
  as the latter function has since been fixed.
* Speed up with `lapply` throughout.
* `gen_results`:
  - Allow speedup by skipping writing to RDS and 
     returning merged results directly instead.
* Use `@inheritParams` throughout.
* Use `\link` throughout.
* `load_example_ctd`: 
  - Cache file  
* `merge_results`:
  - Now handles lists of file path AND list of results lists.
* Remove `RDA_assign_load`
  - Not used and supplanted by `EWCE::load_rdata()` anyway.
* Remove unnecessary parameter descriptions in *MultiEWCE* vignette.
* Change `results_dir` --> `save_dir` throughout.
* `load_example_CTD` --> `load_example_ctd` for function naming consistency.
* Add *README.Rmd*. 

## Bug fixes

* Update to use `EWCE`>1.0.0
* Make more functions internal:
  - `is_not_analysed` 
  - `get_valid_gene_lists` 
* Make paths with `file.path` instead of hard-coded "/"
* Write files to `tempdir()` instead of current working directory.
* Remove unused packages:
  - `cowplot`
  - `ewceData`
* Fix all unit tests (not set up in correct format originally).
