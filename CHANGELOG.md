# nf-core/pathogensurveillance: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.2.0 - Current

## Added

- Added option to use BWA3 for read alignment, increasing speed. This is controlled by the `--aligner` parameter used to switch between BWA and BWA3. 'bwa' must be used with non-AVX2-capable CPUs.

## Changed

- Removed `TRIM_AND_SKETCH` step. Samples that fail assembly are no longer sketched from raw reads and are instead excluded from sketch-based comparisons and downstream analyses that rely on them. Warnings for failed assemblies remain in the report output.
- Assembly RAM is now calculated more accurately to avoid retrying to increase RAM
- Increase number of reference assemblies downloaded by default
- Added `scaledRam` helper function and `--ram_scale`/`--max_ram` parameters to control RAM allocation across processes, analogous to the existing CPU scaling functionality.
- Changed defaults for `--max_cpus` from `0` to `64` and `--max_ram` from `0` to `1000`. Previously, `0` meant "no limit"; now the parameters enforce a hard ceiling and `0` is invalid.
- When there are fewer unique strains than `n_ref_strains`, the pipeline now makes up the difference with more examples of the same strain. This ensure that multiple representatives of each species-level taxon are downloaded when strain-level info is not present.
- Reference assemblies are now downloaded in parallel from NCBI.

### `Fixed`

- Modified Flye module to detect "No disjointigs were assembled" error and exit gracefully instead of retrying.
- Fixed errors when using "excluded" NCBI accessions
- Handle samples that return no sendsketch results.
- Prevent PIRATE and BUSCO phylogeny analyses from running on references associated with samples whose assemblies failed.
- Prevent `SOURMASH_COMPARE` from running on report groups that have no successful assemblies.
- Fixed intermittent PIRATE 'input file name collision' error
- Reference NCBI assemblies excluded by the user are no longer downloaded.
- When multiple assembly versions exist for the same organism (e.g. GCA_000230695.1 vs GCA_000230695.3, or GCA vs GCF), only the best version is retained. RefSeq (GCF) is preferred over GenBank (GCA), and higher version numbers are preferred.

## 1.1.0 - 2026-01-30

### `Added`

- `--max_parallel_downloads` parameter to control how many downloads can occur in parallel. This can be used in HPC/cloud contexts to raise the default maximum that is process-specific and generally below 10 to avoid exceeding API limits when running locally. The `cloud` profile was added to provide environment-specific defaults for this parameter and any other relevant ones in the future.
- `--cpu_scale` parameter to scale the number of CPUs used by multithreaded processes. This allows users to adjust CPU usage across all multithreaded processes
  proportionally (e.g., setting to 0.5 will halve CPU usage, setting to 2 will double it).
- `skip_core_phylogeny` parameter to skip the core phylogeny step if not required.

### `Changed`

- `INITIAL_CLASSIFICATION` now includes taxa with representatives within 1% of the highest ANI match for each rank if there are no matches that pass the quality thresholds.
- New and much better tree plotting widget for the main report.
- 4 samples/references are now required to make variant trees to allow for bootstrapping.

### `Fixed`

- Removed redundant `reference_id` column in output of `PICK_ASSEMBLIES`
- Made `FIND_ASSEMBLIES` not error when exit code is 0 (e.g. when "New version of client..." output in standard error).
- Fixed bakta database download error.
- Updated NCBI tools use to download data, fixing error with downloads.

### `Dependencies`

- Various modules updates.

## 1.0.0 - 2025-06-27

Initial release of nf-core/pathogensurveillance, created with the [nf-core](https://nf-co.re/) template.
