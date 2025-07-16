# nf-core/pathogensurveillance: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.0.0 - 2025-06-27

Initial release of nf-core/pathogensurveillance, created with the [nf-core](https://nf-co.re/) template.

## 1.1.0

### `Added`

- `--max_parallel_downloads` parameter to control how many downloads can occur in parallel. This can be used in cluster/cloud contexts to raise the default maximum that is process-specific and generally below 10 to avoid exceeding API limits when running locally. The `cloud` and `cluster` profiles were added to provide environment-specific defaults for this parameter and any other relevant ones in the future.

### `Fixed`

### `Dependencies`

### `Deprecated`
