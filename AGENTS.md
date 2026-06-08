# AGENTS.md — nf-core/pathogensurveillance

Guidelines for AI-assisted development of this Nextflow/nf-core pipeline.

## Project Overview

- **Name**: `nf-core/pathogensurveillance`
- **Type**: Nextflow DSL2 pipeline
- **Purpose**: Population genomics pipeline for pathogen identification, variant detection, and biosurveillance
- **Repository**: https://github.com/nf-core/pathogensurveillance

## Architecture

- **Local modules** (custom logic): `modules/local/`
- **nf-core modules** (community-curated): `modules/nf-core/`
- **Subworkflows** (local + nf-core): `subworkflows/`
- **Helper scripts** (R/Perl): `bin/`
- **Config hierarchy**: `nextflow.config` → `conf/base.config` → `conf/modules.config`
- **Test data**: from `https://raw.githubusercontent.com/nf-core/test-datasets/refs/heads/pathogensurveillance/`

## Key Files

| File                                | Purpose                                                                                  |
| ----------------------------------- | ---------------------------------------------------------------------------------------- |
| `main.nf`                           | Entry point. Orchestrates init, main workflow, completion                                |
| `workflows/pathogensurveillance.nf` | Main workflow logic. Channels data between subworkflows                                  |
| `nextflow.config`                   | Global parameters, profiles, manifest, plugins (`nf-schema@2.5.1`)                       |
| `nextflow_schema.json`              | Parameter schema for CLI validation and help generation                                  |
| `conf/base.config`                  | Default resource allocations and helper functions (`scaledCpu`, `scaledRam`)             |
| `conf/modules.config`               | Per-module `ext.args`, `publishDir`, resource overrides, and `assemblyCpu`/`assemblyRam` |
| `assets/schema_input.json`          | Samplesheet validation schema                                                            |
| `.nf-core.yml`                      | Lint overrides and template metadata                                                     |
| `nf-test.config`                    | nf-test configuration (profile `test`, ignores `modules/nf-core/**/tests`)               |
| `tests/default.nf.test`             | Pipeline-level snapshot test                                                             |
| `bin/`                              | R and Perl scripts called by local modules                                               |
| `.editorconfig`                     | Indentation rules (4 spaces `.nf`, 2 spaces `.md/.yml/.yaml/.html/.css/.js`)             |

## Code Style & Conventions

- **Indentation**: 4 spaces for `.nf` files; 2 spaces for `.md`, `.yml`, `.yaml`, `.html`, `.css`, `.js` (per `.editorconfig`). Do not modify `.editorconfig` rules.
- **Trailing whitespace**: Remove it (enforced by pre-commit).
- **End-of-file**: Ensure files end with a newline (enforced by pre-commit).
- **Process labels**: reuse these from `conf/base.config` for new modules
- **Channel naming** (from `CONTRIBUTING.md`):
  - Initial: `ch_output_from_<process>`
  - Intermediate: `ch_<previous>_for_<next>`
- **Parameters**:
  - Define defaults in `nextflow.config` under `params {}`
  - Add to `nextflow_schema.json` (use `nf-core pipelines schema build` when available)
  - Use snake_case
- **Resource scaling helpers**:
  - `scaledCpu(baseCpus, attempt)` — use in `conf/modules.config` or `conf/base.config`
  - `scaledRam(baseRam, attempt, exponent=null)` — use for RAM scaling
  - `assemblyCpu(domain, readSize)` / `assemblyRam(domain, readSize, attempt, exitStatus)` — Flye/SPAdes only
- **Process containers**:
  - Prefer `nf-core/modules` containers (Biocontainers/BioConda)
  - Local modules must define `conda` and `container` directives
  - Use Wave/Seqera containers for `arm64` compatibility
- **Error handling**:
  - Use `errorStrategy = { task.exitStatus in ((130..145) + 104 + 175) ? 'retry' : 'finish' }` as default
  - API-limited processes (NCBI) use `maxForks`, `beforeScript` with exponential backoff, and `ignore`/`retry` strategies
- **Publishing**:
  - Configure `publishDir` in `conf/modules.config` via `withName:` selectors
  - Exclude `versions.yml` from publishing: `saveAs: { filename -> filename.equals('versions.yml') ? null : filename }`
- **Version tracking**:
  - Every process should emit `path "versions.yml"` with tool versions under `"${task.process}"`

## AI-Assisted Development Recommendations

### Make Minimal Changes

- Prefer editing existing files over creating new ones.
- When adding a new process, add it to the workflow, create a module in `modules/local/`, register it in `conf/modules.config`, and add parameters to `nextflow.config` + `nextflow_schema.json`.
- Do not modify `modules/nf-core/` or `subworkflows/nf-core/` directly; these are managed via `nf-core modules` commands.

### Respect the Resource Model

- Always use `scaledCpu`/`scaledRam` or `assemblyCpu`/`assemblyRam` instead of hardcoding resources.
- If adding a new multithreaded process, assign it a `process_*` label or use `scaledCpu` in `conf/modules.config`.
- Consider `maxForks` and `beforeScript` for any process that hits external APIs (e.g., NCBI, SRA).

### Keep API Rate Limits in Mind

- Several processes (e.g., `BBMAP_SENDSKETCH`, `FIND_ASSEMBLIES`, `DOWNLOAD_ASSEMBLIES`, `INITIAL_CLASSIFICATION`, `SRATOOLS_FASTERQDUMP`) interact with NCBI or other public servers.
- When modifying these or adding similar ones, preserve or replicate the `maxForks`, `errorStrategy`, and `beforeScript` patterns used in `conf/modules.config`.
- Respect the `params.max_parallel_downloads` parameter for fork limits.

### Maintain Snapshot Compatibility

- The pipeline has nf-test snapshot tests in `tests/default.nf.test` and `tests/default.nf.test.snap`.
- **Do not run nf-test after every change** — it is very slow. Reserve full test runs for CI or pre-PR validation.
- If you change output file names or directory structures, you will likely need to update the snapshot file. Communicate this clearly.

### Version & Changelog

- Update `CHANGELOG.md` for any user-facing change (follow Keep a Changelog format).
- Make change descriptions short and only mention changes relevant to end-users. Dont describe changes to implementation details.

### Documentation

- Update `docs/usage.md` or `docs/output.md` if you change inputs, outputs, or behavior.
- Update `CITATIONS.md` if adding new tools.

## Linting & Validation

### What to run locally (fast)

- **pre-commit** (available in this environment):
  ```bash
  pre-commit run --all-files
  ```
  This runs:
  - `prettier` (formatting)
  - `trailing-whitespace`
  - `end-of-file-fixer`

### What to run before submitting / when preparing a PR

- **nf-core lint** (requires `nf-core` Python package installed):
  ```bash
  nf-core pipelines lint --dir . --markdown lint_results.md
  ```

  - The project has `.nf-core.yml` with specific overrides (e.g., `files_unchanged`, `nf_test_content`, `pipeline_todos`, `template_strings` are disabled).
  - Do not remove these overrides unless you are intentionally template-syncing.
  - Some lint failures can be auto-fixed: `nf-core pipelines lint --fix --dir .`

### What NOT to run after every change

- **nf-test** and **full pipeline runs** are slow. Avoid running them on every edit. Use them:
  - Before opening a PR
  - After major logic changes
  - Via CI (GitHub Actions) rather than locally if possible

