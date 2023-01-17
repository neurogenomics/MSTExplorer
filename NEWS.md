# MultiEWCE 0.1.2

## New features

* New functions:
  - `load_example_results`
  - `ontology_plot`
  - various support functions
* Add `ontologyPlot` as new *Import*. 

## Bug fixes

* `load_example_ctd`
  - Pass `file` to `piggyback::pb_download`

# MultiEWCE 0.1.1

## New features

* Added a `NEWS.md` file to track changes to the package.
* Added `rworkflows`
* `ewce_para`
  - Harmonize arguments with `EWCE`
  - Set default `list_names` arg
  - Output named list of saved files instead of `TRUE`.
  - Add messages that appear in parallel. Make EWCE message silent.
  - Allow users to set seed.
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
