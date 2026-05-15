# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

PerturbPlan is an R package (current version 0.3.1) for **power analysis and experimental design of CRISPR perturbation single-cell experiments** — both perturb-seq (whole-transcriptome) and TAP-seq (targeted). It combines R and C++ (Rcpp) code for efficient statistical computation.

The package is developed by the Katsevich Lab (Wharton/UPenn) and pairs with a separately deployed **PerturbPlan web app** (<https://katsevich-lab-perturbplan.share.connect.posit.cloud/>). The package's current primary public role is **preparing custom reference data for that web app**; it also exposes more advanced programmatic power-analysis APIs that go "beyond the web app."

> Note: There is **no Shiny app bundled in this package** (no `inst/shiny/`, no `launch_app()`). The web app is a separate deployment.

## Common Development Commands

### Build and Check
```bash
# Build the package
R CMD build .

# Check the package (replace version as needed)
R CMD check perturbplan_0.3.1.tar.gz
```

### Development Workflow
```r
# Load package during development
devtools::load_all()

# Run all tests
devtools::test()

# Run a specific test file
devtools::test(filter = "library_computation")

# Generate documentation from roxygen2 comments
devtools::document()

# Check package without building
devtools::check()

# Rebuild the pkgdown site
pkgdown::build_site()
```

## Architecture

The package supports three workflows.

### 1. Pilot-data preprocessing (`R/pilot_data_preprocessing.R`, `R/pilot_data_help.R`)

Turns raw 10x Cell Ranger output into reference data usable by the web app and power functions.

- `reference_data_preprocessing_10x()` — Step 1: aggregates Cell Ranger outputs across SRR runs into a response matrix, a read-UMI table, and a naive mapping-efficiency estimate. (Cell Ranger **count** only — not Cell Ranger multi.)
- `reference_data_processing()` — Step 2: produces `baseline_expression_stats` (negative-binomial expression parameters) and `library_parameters` (read-UMI saturation curve).
- Helpers: `obtain_qc_response_data()`, `obtain_qc_read_umi_table()`, `obtain_mapping_efficiency()`, `obtain_expression_information()`, `library_estimation()`.

### 2. Prospective power analysis / design optimization (`R/power_plan.R`, `R/plan_help.R`)

- `compute_power_plan()` / `compute_power_plan_overall()` / `compute_power_plan_per_grid()` — power for prospective experimental designs.
- `cost_power_computation()` — power analysis with cost minimization across experimental parameters.
- `find_optimal_cost_design()` — binary-search optimization to find designs meeting a power target.
- `obtain_fixed_variable_constraining_cost()` — cost-constrained design helper.
- `extract_fc_expression_info()` / `extract_expression_info()` — sample fold-change/expression info (importance sampling over a gene list).
- `get_pilot_data_from_package()` — load one of the bundled pilot datasets.

### 3. Retrospective (post-hoc) power analysis (`R/power_posthoc.R`, `R/posthoc_help.R`)

- `compute_power_posthoc()` — power for a completed experiment, given per-gRNA cell counts (`cells_per_grna`) and `discovery_pairs`. Returns individual power per perturbation-gene pair plus expected total discoveries. Uses a score-test statistic and BH/Bonferroni FDP estimation.

### Cross-cutting components

- **C++ performance layer** (`src/`): test-statistic distributions (`compute_distribution_teststat_fixed_es.cpp`, `compute_distribution_teststat_random_es.cpp`), BH multiple-testing cutoff (`BH_cutoff.cpp`), QC (`compute_QC.cpp`), library-size saturation curves (`library_size_curves.cpp`), cell-range identification (`identify_cell_range.cpp`), Monte Carlo overall power (`overall_power.cpp`). Bound via Rcpp; see `R/RcppExports.R` / `src/RcppExports.cpp`. Built with `BH`, `Rcpp`, `RcppEigen`, `PoissonBinomial` (`LinkingTo`).
- **Quality control** (`R/QC_computation.R`): pairwise QC checks (minimum non-zero cell counts).
- **Input validation** (`R/check.R`): one validator per major public function.
- **Variable bindings** (`R/perturbplan.R`): `utils::globalVariables()` declarations to suppress NSE R CMD check notes.

### Data

`data/` ships 8 reference pilot datasets, each a list with `baseline_expression_stats`, `library_parameters`, and `mapping_efficiency`:
`K562_Gasperini`, `K562_10x`, `K562_Ray` (TAP-seq), `A549_Sakellaropoulos`, `THP1_Yao`, `T_CD8_Shifrut`, `iPSC_Tian`, `iPSC_neuron_Tian`. Documented in `R/data.R`; regeneration scripts in `inst/data-raw/`.

## Key Design Decisions

- **C++ for hot paths**: Monte Carlo loops, distribution calculations, and multiple-testing corrections are implemented in C++. Prefer C++ implementations over R loops for computationally intensive operations.
- **Modular workflows**: preprocessing, prospective power, and post-hoc power are independent function families with shared C++ and validation layers.
- **`rSAC_fn_wrapper` library model** (since v0.3.0): `library_estimation()` returns an `rSAC_fn_wrapper` list (`method_used`, `UMI_per_cell_at_saturation`, `reads_norm`, `n_cells`, plus method-specific params) using `preseqR` ZTNB/RFA saturation-curve fitting — this replaced the older `minpack.lm` NLS fit. Internal C++ functions consume this wrapper.

## Testing

The package uses testthat (edition 3) with helper functions in `tests/testthat/helper-*.R` for test-data generation. Tests compare analytical computations against simulations to ensure accuracy. Test coverage is ~70%.

## Vignettes

`vignettes/` contains three articles:
- `preprocess-reference.Rmd` — "Prepare Data For Web App" (the primary, currently navbar-linked article).
- `prospective-power.Rmd` — advanced prospective power workflows (perturb-seq vs TAP-seq, comparing designs).
- `posthoc.Rmd` — retrospective power analysis.

`_pkgdown.yml` currently only links `preprocess-reference` in the navbar; the other two are present but commented out.

## Parameter Naming Convention

**IMPORTANT**: Use `TPM_threshold` (Transcripts Per Million) everywhere — never `tmp_threshold`.

- All function parameters, variable names, and documentation should use `TPM_threshold`.
- When modifying existing functions, preserve existing parameter names exactly.
- Note `sequenced_reads_per_cell` refers to **raw** sequencing reads (before mapping), not mapped reads (renamed from `reads_per_cell`/`raw_reads_per_cell` in v0.2.0).

## Common Mistakes to Avoid

- **`TPM_threshold` vs `tmp_threshold`**: always use `TPM_threshold`.
- **Parameter consistency**: when adding parameters, verify spelling matches existing usage and that function definitions and call sites agree.
- **Preserve existing code**: when modifying existing functions, keep parameter names exactly as they are.

## Git Workflow Requirements

**IMPORTANT**: When commit and push is requested, the **entire repository** should be committed and pushed, not just specific changes.

- **Complete Sync**: after commit and push, there should be **no difference** between the local directory and remote repository.
- **Clean Working Tree**: `git status` should show a clean working tree after pushing.
- **Full Commit**: use `git add .` to stage all changes before committing, unless specifically instructed otherwise.
- **Repository Consistency**: the remote repository should always reflect the complete current state of local development.

Active development happens on the `dev` branch; `main` is the release branch.
