# nf-core/pathogensurveillance: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.1.0

### `Added`

- `--max_parallel_downloads` parameter to control how many downloads can occur in parallel. This can be used in cluster/cloud contexts to raise the default maximum that is process-specific and generally below 10 to avoid exceeding API limits when running locally. The `cloud` and `cluster` profiles were added to provide environment-specific defaults for this parameter and any other relevant ones in the future.
- `--cpu_scale` parameter to scale the number of CPUs used by multithreaded processes. This allows users to adjust CPU usage across all multithreaded processes
  proportionally (e.g., setting to 0.5 will halve CPU usage, setting to 2 will double it).
- `skip_core_phylogeny` parameter to skip the core phylogeny step if not required.

### `Changed`

- `INITIAL_CLASSIFICATION` now includes taxa with representatives within 1% of the highest ANI match for each rank if there are no matches that pass the quality thresholds.
- New and much better tree plotting widget for the main report.

### `Fixed`

- Removed redundant `reference_id` column in output of `PICK_ASSEMBLIES`
- Made `FIND_ASSEMBLIES` not error when exit code is 0 (e.g. when "New version of client..." output in standard error).
- Fixed bakta database download error.
- Updated NCBI tools use to download data, fixing error with downloads.

### `Dependencies`

- Various modules updates.

## 1.0.0 - 2025-06-27

Initial release of nf-core/pathogensurveillance, created with the [nf-core](https://nf-co.re/) template.
