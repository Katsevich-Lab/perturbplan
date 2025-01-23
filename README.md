
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PerturbPlan

<!-- badges: start -->
<!-- badges: end -->

The goal of PerturbPlan is to facilitate experimental design and power
analysis for perturb-seq experiments.

``` r
library(perturbplan)
```

In the example below, we will demonstrate how to use
`compute_power_posthoc()` to compute power after an experiment is
conducted (as opposed to prospectively).

## PerturbPlan inputs

### 1. Dataset summaries

Let us consider a toy perturb-seq dataset with 3 enhancers, each
targeted by two gRNAs, and 4 genes. Suppose there are 1000 total cells
in the experiment:

``` r
num_total_cells <- 1000
```

Based on our data, suppose we have information on the number of cells
that received each gRNA:

``` r
cells_per_grna
#>      grna_id grna_target num_cells
#> 1 enh1_grna1        enh1        92
#> 2 enh2_grna1        enh2        96
#> 3 enh3_grna1        enh3        85
#> 4 enh1_grna2        enh1        91
#> 5 enh2_grna2        enh2       104
#> 6 enh3_grna2        enh3        91
```

Furthermore, we have computed the mean and size parameters for the
baseline expression of each gene:

``` r
baseline_expression_stats
#>   response_id expression_mean expression_size
#> 1       gene1        6.677385        6.157879
#> 2       gene2        1.075534        9.396240
#> 3       gene3        1.970892        3.000133
#> 4       gene4        6.635279        3.321322
```

### Analysis choices

Suppose we tested the effects of all elements on all genes. We can
encode this in a `discovery_pairs` data frame, the same as in `sceptre`:

``` r
discovery_pairs
#>    grna_target response_id
#> 1         enh1       gene1
#> 2         enh2       gene1
#> 3         enh3       gene1
#> 4         enh1       gene2
#> 5         enh2       gene2
#> 6         enh3       gene2
#> 7         enh1       gene3
#> 8         enh2       gene3
#> 9         enh3       gene3
#> 10        enh1       gene4
#> 11        enh2       gene4
#> 12        enh3       gene4
```

We have analyzed the data using the complement control group, left-sided
tests, and default values for the pairwise QC parameters:

``` r
control_group <- "complement"
side <- "both"
n_nonzero_trt_thresh <- 7
n_nonzero_cntrl_thresh <- 7
```

### 3. Power analysis parameters

We want to compute the power for enhancer-gene pair assuming that gRNAs
targeting that enhancer have fold changes on the gene with mean and
standard deviations of 0.85 and 0.13, respectively:

``` r
fold_change_mean <- 0.85
fold_change_sd <- 0.13
```

We want to deem perturbation-gene pairs as significant in our power
analysis if they pass a threshold of 0.005:

``` r
cutoff <- 0.005
```

This number may be obtained as the $p$-value cutoff from the original
analysis.

## Running PerturbPlan

We can now compute the power for each enhancer-gene pair using the
`compute_power_posthoc()` function:

``` r
power_results <- compute_power_posthoc(
  num_total_cells = num_total_cells,
  cells_per_grna = cells_per_grna,
  baseline_expression_stats = baseline_expression_stats,
  discovery_pairs = discovery_pairs,
  control_group = control_group,
  side = side,
  n_nonzero_trt_thresh = n_nonzero_trt_thresh,
  n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh,
  fold_change_mean = fold_change_mean,
  fold_change_sd = fold_change_sd,
  cutoff = cutoff
)
```

## PerturbPlan outputs

The result is a list with fields `expected_num_discoveries` and
`individual_power`:

``` r
names(power_results)
#> [1] "individual_power"         "expected_num_discoveries"
```

The field `expected_num_discoveries` contains the expected number of
discoveries across all enhancer-gene pairs:

``` r
power_results$expected_num_discoveries
#> [1] 5.034662
```

The field `individual_power` is a data frame containing the power for
each enhancer-gene pair:

``` r
power_results$individual_power
#> # A tibble: 12 × 3
#>    grna_target response_id power
#>    <fct>       <chr>       <dbl>
#>  1 enh1        gene1       0.599
#>  2 enh2        gene1       0.616
#>  3 enh3        gene1       0.591
#>  4 enh1        gene2       0.253
#>  5 enh2        gene2       0.273
#>  6 enh3        gene2       0.245
#>  7 enh1        gene3       0.312
#>  8 enh2        gene3       0.334
#>  9 enh3        gene3       0.303
#> 10 enh1        gene4       0.499
#> 11 enh2        gene4       0.520
#> 12 enh3        gene4       0.490
```
