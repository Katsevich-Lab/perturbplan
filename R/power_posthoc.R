# This is a Rscript computing the power function using score test

#' Compute power for each perturbation-gene pair
#'
#' @param discovery_pairs A data frame specifying which element-gene pairs to consider, with columns `grna_target` and `response_id`
#' @param cells_per_grna A data frame specifying how many cells contain each gRNA, with columns `grna_id`, `grna_target`, and `num_cells`
#' @param baseline_expression_stats A data frame specifying the baseline expression statistics for each gene, with columns `response_id`, `expression_mean`, and `expression_size`
#' @param control_group A character string specifying the control group, either "complement" or "nt_cells"
#' @param fold_change_mean A numeric value to use for mean effect size for all element-gene pairs
#' @param fold_change_sd A numeric value to use for standard deviation of effect size for all element-gene pairs
#' @param num_total_cells (Required only if control_group == "complement") A positive integer specifying the total number of cells in the experiment
#' @param cutoff (Optional) A numeric value between 0 and 1 to use as the p-value cutoff
#' @param n_nonzero_trt_thresh (Optional) An integer specifying the sceptre QC parameter of the same name; defaults to 7
#' @param n_nonzero_cntrl_thresh (Optional) An integer specifying the sceptre QC parameter of the same name; defaults to 7
#' @param side (Optional) A character string specifying the side of the test, either "left", "right", or "both"; defaults to "both"
#' @param multiple_testing_method (Optional) A character string specifying the multiple testing correction method to use, either "BH" or "bonferroni"; defaults to "BH"
#' @param multiple_testing_alpha (Optional) A numeric value between 0 and 1 specifying the alpha level for multiple testing correction; defaults to 0.1
#'
#' @return A list with two elements: `individual_power` (a data frame with columns `grna_target`, `response_id`, and `power`) and `expected_num_discoveries` (a numeric value)
#'
#' @examples
#' ## --- Toy perturb-seq dataset setup ---
#' # Total number of cells in the experiment
#' num_total_cells <- 1000L
#'
#' # Number of cells receiving each gRNA
#' # (3 enhancers × 2 gRNAs each + 2 non-targeting gRNAs)
#' cells_per_grna <- data.frame(
#'   grna_id      = c("enh1_grna1","enh2_grna1","enh3_grna1","nt_grna1",
#'                    "enh1_grna2","enh2_grna2","enh3_grna2","nt_grna2"),
#'   grna_target  = c("enh1","enh2","enh3","non-targeting",
#'                    "enh1","enh2","enh3","non-targeting"),
#'   num_cells    = c(93L,113L,112L,104L,84L,104L,107L,105L),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Baseline expression statistics (negative binomial mean and size per gene)
#' baseline_expression_stats <- data.frame(
#'   response_id       = paste0("gene", 1:4),
#'   expression_mean   = c(2.002931, 12.326867, 4.014221, 1.460472),
#'   expression_size   = c(0.2967991, 8.3723191, 2.5988431, 2.6746265),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Discovery pairs: test 3 enhancers against 4 genes
#' discovery_pairs <- within(expand.grid(
#'   grna_target = c("enh1","enh2","enh3"),
#'   response_id = paste0("gene", 1:4),
#'   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
#' ), { grna_target <- as.character(grna_target); response_id <- as.character(response_id) })
#'
#' ## --- Analysis choices ---
#' control_group <- "complement"          # use complement control group
#' side <- "left"                         # left-sided test
#' n_nonzero_trt_thresh   <- 7L            # min. nonzero counts in treatment
#' n_nonzero_cntrl_thresh <- 7L            # min. nonzero counts in control
#'
#' ## --- Power analysis parameters ---
#' fold_change_mean <- 0.85                # expected mean fold change
#' fold_change_sd   <- 0.13                # expected SD of fold change
#' cutoff <- 0.005                         # significance threshold (p-value cutoff)
#'
#' ## --- Run PerturbPlan posthoc power analysis ---
#' power_results <- compute_power_posthoc(
#'   num_total_cells = num_total_cells,
#'   cells_per_grna = cells_per_grna,
#'   baseline_expression_stats = baseline_expression_stats,
#'   discovery_pairs = discovery_pairs,
#'   control_group = control_group,
#'   side = side,
#'   n_nonzero_trt_thresh = n_nonzero_trt_thresh,
#'   n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh,
#'   fold_change_mean = fold_change_mean,
#'   fold_change_sd = fold_change_sd,
#'   cutoff = cutoff
#' )
#'
#' # Inspect outputs
#' names(power_results)                    # available fields
#' power_results$expected_num_discoveries  # expected number of discoveries
#' head(power_results$individual_power)    # power per enhancer-gene pair
#'
#' @export
compute_power_posthoc <- function(
    discovery_pairs,
    cells_per_grna,
    baseline_expression_stats,
    control_group,
    fold_change_mean,
    fold_change_sd,
    num_total_cells = NULL,
    cutoff = NULL,
    n_nonzero_trt_thresh = 7L,
    n_nonzero_cntrl_thresh = 7L,
    side = "both",
    multiple_testing_method = "BH",
    multiple_testing_alpha = 0.1) {

  ############################# perform input checks ###########################
  input_check_posthoc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = control_group,
    fold_change_mean = fold_change_mean,
    fold_change_sd = fold_change_sd,
    num_total_cells = num_total_cells,
    cutoff = cutoff,
    n_nonzero_trt_thresh = n_nonzero_trt_thresh,
    n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh,
    side = side,
    multiple_testing_method = multiple_testing_method,
    multiple_testing_alpha = multiple_testing_alpha
  )

  ############################# create enhancer_gene df ########################
  enhancer_gene <- discovery_pairs |>
    # join grna df
    dplyr::left_join(
      cells_per_grna |>
        dplyr::filter(grna_target != "non-targeting") |>
        dplyr::group_by(grna_target) |>
        dplyr::summarize(num_trt_cells = sum(num_cells),
                         num_trt_cells_sq = sum(num_cells^2)) |>
        dplyr::ungroup(),
      "grna_target",
      relationship = "many-to-one"
    ) |>
    # join gene expression df
    dplyr::left_join(
      baseline_expression_stats, "response_id",
      relationship = "many-to-one"
    )

  ############################# obtain number of control cells #################
  if (control_group == "nt_cells") {
    num_cntrl_cells <- cells_per_grna |>
      dplyr::filter(grna_target == "non-targeting") |>
      dplyr::summarize(sum(num_cells)) |>
      dplyr::pull()
    enhancer_gene <- enhancer_gene |> dplyr::mutate(
      num_cntrl_cells = num_cntrl_cells
      )
  } else { # control_group == "complement"
    enhancer_gene <- enhancer_gene |> dplyr::mutate(
      num_cntrl_cells = num_total_cells - num_trt_cells
    )
  }

  ################## transform the scalar-valued effect size mean/sd ###########
  if(is.numeric(fold_change_mean)){
    # create the effect size matrices
    enhancer_gene <- enhancer_gene |> dplyr::mutate(
      fold_change_mean = fold_change_mean,
      fold_change_sd = fold_change_sd
    )
  }else{
    # join the enhancer_gene df and effect size mean(sd) data frames
    enhancer_gene <- enhancer_gene |>
      dplyr::left_join(
        fold_change_mean,
        c("grna_target", "response_id"),
        relationship = "one-to-one"
      ) |>
      dplyr::left_join(
        fold_change_sd,
        c("grna_target", "response_id"),
        relationship = "one-to-one"
      )
  }

  ########################### prepare for multiple testing #####################
  enhancer_gene <- enhancer_gene |>
    dplyr::group_by(grna_target, response_id) |>
    dplyr::mutate(
      # compute mean and sd of the test statistic for each pair
      test_stat_distribution = compute_distribution_teststat(
        num_trt_cells = num_trt_cells,
        num_cntrl_cells = num_cntrl_cells,
        num_trt_cells_sq = num_trt_cells_sq,
        expression_mean = expression_mean,
        expression_size = expression_size,
        fold_change_mean = fold_change_mean,
        fold_change_sd = fold_change_sd
      ),
      # extract mean and sd from test_stat_distribution
      mean_test_stat = unlist(test_stat_distribution)["mean"],
      sd_test_stat = unlist(test_stat_distribution)["sd"],
      # compute QC probability
      QC_prob = compute_QC(
        fold_change_mean = fold_change_mean,
        expression_mean = expression_mean,
        expression_size = expression_size,
        num_cntrl_cells = num_cntrl_cells,
        num_trt_cells = num_trt_cells,
        n_nonzero_trt_thresh = n_nonzero_trt_thresh,
        n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh)
    ) |>
    dplyr::select(-test_stat_distribution) |>
    dplyr::ungroup()

  ########################### correct multiplicity #############################
  # compute cutoff if it is NULL
  if(is.null(cutoff)){
    cutoff <- enhancer_gene |>
      dplyr::summarize(
        cutoff = adjusted_cutoff(mean_list = mean_test_stat,
                                 sd_list = sd_test_stat,
                                 multiple_testing_alpha = multiple_testing_alpha,
                                 multiple_testing_method = multiple_testing_method,
                                 side = side, QC_prob = QC_prob)
      ) |> dplyr::select(cutoff) |> dplyr::pull()
  }
  # compute the adjusted power
  enhancer_gene <- enhancer_gene |>
    dplyr::mutate(
      cutoff = cutoff,
      power = rejection_computation(mean_list = mean_test_stat,
                                    sd_list = sd_test_stat,
                                    side = side,
                                    cutoff = cutoff) * (1 - QC_prob)
    )

  # store individual power and rejection size as the output
  output <- list(
    individual_power = enhancer_gene |> dplyr::select(grna_target, response_id, power),
    expected_num_discoveries = sum(enhancer_gene$power)
  )

  return(output)
}


#' Compute approximate power of a CRISPR screen
#'
#' This function computes the approximate power of detecting an effect (gene perturbation)
#' in a CRISPR screen given various experimental and sequencing parameters. Internally, it
#' calculates an average library size using provided parameters and estimates baseline
#' expression levels. The function then calls \code{\link{compute_power_posthoc}} to obtain
#' power estimates and the expected number of discoveries.
#'
#' @param recovery_rate A numeric value (between 0 and 1) indicating the fraction
#'   of cells that survive and are captured after library preparation.
#' @param num_total_reads A numeric value specifying the total number of reads
#'   generated by sequencing. This is used to estimate the \code{library_size}.
#' @param mapping_efficiency A numeric value (between 0 and 1) indicating the
#'   fraction of reads that successfully map to the transcriptome.
#' @param cells_per_grna A data frame specifying the number of cells per gRNA,
#'   with columns \code{grna_id}, \code{grna_target}, and \code{num_cells}.
#' @param baseline_relative_expression_stats A data frame specifying the relative
#'   expression levels for each gene, with columns \code{response_id} and
#'   \code{relative_expression}.
#' @param fold_change_mean A numeric value indicating the mean fold change effect size
#'   for all gRNA-gene pairs (or a data frame with \code{grna_target} and \code{response_id}
#'   columns for per-pair values).
#' @param fold_change_sd A numeric value indicating the standard deviation of the fold
#'   change effect size for all gRNA-gene pairs (or a data frame with \code{grna_target}
#'   and \code{response_id} columns for per-pair values).
#' @param num_planned_cells A numeric value indicating the total planned number
#'   of cells before losses in library preparation.
#' @param control_group A character string specifying the control group, either
#'   \code{"complement"} or \code{"nt_cells"}. This is passed to
#'   \code{compute_power_posthoc}.
#' @param UMI_per_cell A numeric value specifying the theoretical saturation level (in UMIs)
#'   for each cell.
#' @param variation A numeric value controlling how overdispersion in UMIs per read
#'   is modeled.
#' @param side (Optional) A character string specifying the side of the test, either
#'   \code{"left"}, \code{"right"}, or \code{"both"}. Defaults to \code{"both"}.
#' @param multiple_testing_method (Optional) A character string specifying the multiple
#'   testing correction method, either \code{"BH"} or \code{"bonferroni"}. Defaults
#'   to \code{"BH"}.
#' @param multiple_testing_alpha (Optional) A numeric value (between 0 and 1) specifying the
#'   alpha level for multiple testing correction. Defaults to 0.1.
#' @param cutoff (Optional) A numeric value between 0 and 1 to use as the p-value cutoff.
#'   If \code{NULL}, the function determines it automatically using the specified
#'   \code{multiple_testing_method} and \code{multiple_testing_alpha}.
#' @param discovery_pairs A data frame specifying which gRNA-gene pairs to consider, with
#'   columns \code{grna_target} and \code{response_id}.
#' @param n_nonzero_trt_thresh (Optional) An integer specifying the sceptre QC parameter
#'   of the same name; defaults to 7.
#' @param n_nonzero_cntrl_thresh (Optional) An integer specifying the sceptre QC parameter
#'   of the same name; defaults to 7.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{individual_power}: A data frame with columns \code{grna_target},
#'         \code{response_id}, and \code{power}, giving the power for each pair.
#'   \item \code{expected_num_discoveries}: A numeric value indicating the expected total
#'         number of discoveries.
#' }
#'
#' @keywords internal
power_function <- function(

    ###################### specify experimental design parameters ##############
    recovery_rate,                                              # Library prep parameters
    num_total_reads, mapping_efficiency,                        # Sequencing parameters

    ######################## specify the power-determining parameters ##########
    cells_per_grna, baseline_relative_expression_stats,
    fold_change_mean, fold_change_sd, num_planned_cells, control_group,
    UMI_per_cell, variation,

    ###################### specify test-related parameters #####################
    side = "both", multiple_testing_method = "BH", multiple_testing_alpha = 0.1,
    cutoff = NULL, discovery_pairs,

    ######################## specify QC-related parameters ################
    n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L
){
  ########## compute library size with other power-determing parameters ########

  # compute the number of total cells (singletons) surviving from library preparation
  num_total_cells <- num_planned_cells * recovery_rate

  # compute the reads per cell
  reads_per_cell <- num_total_reads*mapping_efficiency/num_total_cells

  # compute the averaged library size with read per cell
  avg_library_size <- UMI_per_cell * (1 - exp(-reads_per_cell / UMI_per_cell) *
                                        (1 + variation * reads_per_cell^2 / (2*UMI_per_cell^2)))

  ####### perform power calculation with power-determining parameters ##########

  # compute the baseline_expression
  baseline_expression_stats <- baseline_relative_expression_stats |>
    dplyr::mutate(expression_mean = avg_library_size * relative_expression) |>
    dplyr::select(-relative_expression)

  # compute power using function compute_power
  power_result <- compute_power_posthoc(
    discovery_pairs = discovery_pairs,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    control_group = control_group,
    fold_change_mean = fold_change_mean,
    fold_change_sd = fold_change_sd,
    num_total_cells = num_total_cells,
    cutoff = cutoff,
    n_nonzero_trt_thresh = n_nonzero_trt_thresh,
    n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh,
    side = side,
    multiple_testing_method = multiple_testing_method,
    multiple_testing_alpha = multiple_testing_alpha)

  # return the power_result
  return(power_result)
}

