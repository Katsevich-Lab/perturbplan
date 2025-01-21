# This is a Rscript computing the power function using score test

#' Power function for CRISPR screen experimental design.
#'
#'
#'
#' @param perturb_type Either CRISPRi or CRISPRko; currently can only handle CRISPRi.
#'
#' @param effect_size_mean Mean fold change for each (element, gene) (matrix; row: L genomic elements;
#' column: J genes). Rownames and colnames should be specified. The rowname should look like "enh1"
#' and column name should be ensembl form or like "gene1". Alternatively, a scalar value can just be
#' specified if all enhancer-gene pairs share the same effect size.
#' @param effect_size_sd Sd fold change for each (element, gene). Format should be the same as `effect_size_mean`.
#'
#' @param target_cell_df Perturb cell size for L elements (dataframe). The data frame should include
#' three columns: grna_id (value should be like "gRNAd" where d is an integer), grna_target (value should
#' be like "enhl" where l is an integer) and num_cells (value should be integer). Thus each row of the data
#' frame includes the number of cells for each (gRNA, element) pair.
#' @param control_cell_vec Control cell size for L elements (vector; length L). Vector should be named
#' for grna_target.
#'
#' @param library_size Averaged library size (scalar; positive valued).
#' @param relative_expression Baseline relative expression level for J genes (named vector; of length J).
#' Names of the vector should be different genes and the format should be consistent with the column names of
#' `effect_size_mean` and `effect_size_sd`.
#' @param size_parameter Size parameter for J genes (vector; length J). Names should be same as `relative_expression`.
#'
#' @param side Left, right or both.
#' @param multiple_testing_method Multiplicity correction; default is BH.
#' @param multiple_testing_alpha Significance level of interest (FDR of interest).
#' @param n_nonzero_trt_thresh QC threshold for treatment cell.
#' @param n_nonzero_cntrl_thresh QC threshold for control cell.
#'
#' @param UMI_s Estimate for UMI count per singlet for the underlying type of cell.
#' @param recovery_rate Recovery rate for cells surviving the library preparation.
#' @param doublet_rate Doublet rate for droplet containing more than one cell.
#' @param doublet_factor Ratio of averaged UMI count per doublet to UMI per singlet.
#' @param planned_read Planned total sequencing reads.
#' @param mapping_efficiency Mapping efficiency for sequenced reads.
#'
#' @param intermediate_outcome A logic value indicating if only mean, sd of test statistics and QC probability are desired.
#'
#' @return The output format is the same as [compute_power()].
#' @importFrom dplyr if_else group_by ungroup summarise
#' @importFrom stats setNames
#' @export

power_function <- function(
    ######################## specify the type of perturb-seq experiment ########
    perturb_type = "CRISPRi",                                          # currenly only support CRISPRi

    ######################## specify the power-determining parameters ##########
    target_cell_df, control_cell_vec,
    library_size = NULL, relative_expression, size_parameter,
    effect_size_mean = 0.8, effect_size_sd = 0,

    ######################## specify QC-related parameters ################
    n_nonzero_trt_thresh = 7, n_nonzero_cntrl_thresh = 7,

    ###################### specify experimental design parameters ##############
    UMI_s = NULL,                                                      # cell-type specific parameter
    recovery_rate = NULL, doublet_rate = NULL, doublet_factor = NULL,  # Library prep parameters
    planned_read = NULL, mapping_efficiency = NULL,                    # Sequencing parameters

    ###################### specify test-related parameters #####################
    side, multiple_testing_method = "BH", multiple_testing_alpha = 0.1,

    ###################### output asymptotic mean/sd if needed #################
    intermediate_outcome = FALSE
){

  # either library size or UMI_s has to be provided
  if(all(is.null(library_size), is.null(UMI_s))){
    stop("One of library size or UMI per singlet has to be provided!")
  }

  ########## compute library size with other power-determing parameters ########
  if(is.null(library_size)){

    # compute the total number of cell
    cell_per_element <- target_cell_df |>
      dplyr::group_by(gRNA_target) |>
      dplyr::summarise(cell_per_element = sum(num_cells)) |>
      dplyr::ungroup()
    target_cell_vec <- setNames(cell_per_element$cell_per_element, cell_per_element$element)
    total_cell <- (control_cell_vec + target_cell_vec)[1]

    # compute the reads per cell
    read_c <- read_per_cell(planned_read = planned_read,
                            mapping_efficiency = mapping_efficiency,
                            planned_cell = total_cell,
                            recovery_rate = recovery_rate)

    # compute the averaged library size with read per cell
    library_size <- library_computation(UMI_s =  UMI_s,
                                        read_c = read_c,
                                        doublet_rate = doublet_rate,
                                        doublet_factor = doublet_factor)
  }

  ####### perform power calculation with power-determining parameters ##########

  # compute the baseline_expression
  baseline_expression <- setNames(relative_expression * library_size, names(relative_expression))

  # compute power using function compute_power
  power_result <- compute_power(perturb_type = perturb_type,
                                control_cell_vec = control_cell_vec, target_cell_df = target_cell_df,
                                baseline_expression = baseline_expression, size_parameter = size_parameter,
                                effect_size_mean = effect_size_mean, effect_size_sd = effect_size_sd,
                                n_nonzero_trt_thresh = n_nonzero_trt_thresh, n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh,
                                side = side, multiple_testing_method = multiple_testing_method, multiple_testing_alpha = multiple_testing_alpha,
                                intermediate_outcome = intermediate_outcome)

  # return the power_result
  return(power_result)
}


#' Computing power for each (enhancer, gene) pair.
#' @description
#'
#' This function have two functionalities: (1, intermediate goal) compute the mean
#' and sd of the test statistic; (2, final goal) output the power for each enhancer-gene
#' pair.
#'
#' @inheritParams power_function
#' @param cutoff Cutoff for p-values to reject the test.
#' @param baseline_expression Baseline mean expression in the group of control cells.
#'
#' @return A dataframe or a list depending on `intermediate_outcome`:
#' @section Dataframe: if `intermediate_outcome` is TRUE, only mean, sd of the score tests and QC probability will be outputted.
#' @section List: otherwise, the final power for each element-gene pair will be outputted, together with the estimated discovery set.
#'
#' @export
compute_power <- function(perturb_type = "CRISPRi",
                          control_cell_vec, target_cell_df,
                          baseline_expression, size_parameter,
                          effect_size_mean, effect_size_sd,
                          n_nonzero_trt_thresh = 7, n_nonzero_cntrl_thresh = 7,
                          side = "both", multiple_testing_method = "BH", multiple_testing_alpha = 0.1, cutoff = NULL,
                          intermediate_outcome = FALSE){

  ###################################### check the input #######################
  input_check(perturb_type = perturb_type,
              control_cell_vec = control_cell_vec, target_cell_df = target_cell_df,
              baseline_expression = baseline_expression, size_parameter = size_parameter,
              effect_size_mean = effect_size_mean, effect_size_sd = effect_size_sd,
              n_nonzero_trt_thresh = n_nonzero_trt_thresh, n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh,
              side = side, multiple_testing_method = multiple_testing_method,
              multiple_testing_alpha = multiple_testing_alpha, cutoff = cutoff,
              intermediate_outcome = intermediate_outcome)

  ######## specify default order of gRNA_target and gene sets ##################
  default_gRNA_target_set <- names(control_cell_vec)
  default_gene_set <- names(baseline_expression)

  # obtain Enhancer-Gene pair
  E2G_pair <- as.vector(outer(default_gRNA_target_set, default_gene_set, paste, sep = "_"))

  ################## transform the scalar-valued effect size matrix ############
  if(is.numeric(effect_size_mean)){

    # create the effect size matrices
    effect_size_mean <- matrix(effect_size_mean,
                               nrow = length(default_gRNA_target_set),
                               ncol = length(default_gene_set),
                               dimnames = list(
                                 gRNA_target = default_gRNA_target_set,
                                 gene = default_gene_set
                               ))
    effect_size_sd <- matrix(effect_size_sd,
                             nrow = length(default_gRNA_target_set),
                             ncol = length(default_gene_set),
                             dimnames = list(
                               gRNA_target = default_gRNA_target_set,
                               gene = default_gene_set
                             ))

  }

  ##################### compute number of target cells #########################
  cell_per_element <- target_cell_df |>
    dplyr::group_by(gRNA_target) |>
    dplyr::summarise(per_element = sum(num_cells),
                     sq_per_element = sum(num_cells^2)) |>
    dplyr::ungroup()
  target_cell_vec <- setNames(cell_per_element$per_element, cell_per_element$gRNA_target)
  target_cell_vec_sq <- setNames(cell_per_element$sq_per_element, cell_per_element$gRNA_target)

  ####### change the order of vectors or matrices to default order #############
  target_cell_vec <- target_cell_vec[default_gRNA_target_set]
  target_cell_vec_sq <- target_cell_vec_sq[default_gRNA_target_set]
  baseline_expression <- baseline_expression[default_gene_set]
  size_parameter <- size_parameter[default_gene_set]
  effect_size_mean <- effect_size_mean[default_gRNA_target_set, default_gene_set]
  effect_size_sd <- effect_size_sd[default_gRNA_target_set, default_gene_set]

  # divide the power computation depending on perturb_type
  if(perturb_type == "CRISPRi"){

    # compute asy_mean and asy_sd
    asy_quantity <- distribution_teststat(control_cell_vec = control_cell_vec,
                                          target_cell_vec = target_cell_vec,
                                          target_cell_vec_sq = target_cell_vec_sq,
                                          baseline_expression = baseline_expression,
                                          size_parameter = size_parameter,
                                          effect_size_mean = effect_size_mean,
                                          effect_size_sd = effect_size_sd)

    # obtain asy_mean and asy_sd by vectorizing the obtained matrices
    asy_mean <- as.vector(asy_quantity$mean)
    asy_sd <- as.vector(asy_quantity$sd)

  }else{
    stop("Code is not finished for CRISPR perturbation beyond interference!")
  }

  # next step depends if QC_prob is used or not
  if((n_nonzero_trt_thresh == 0) & (n_nonzero_cntrl_thresh == 0)){

    # add QC probability
    QC_prob <- rep(0, length(E2G_pair))

  }else{

    # compute the QC_prob
    QC_prob <- QC_prob(effect_size_mean = effect_size_mean,
                       baseline_expression = baseline_expression,
                       size_parameter = size_parameter,
                       control_cell_vec = control_cell_vec, target_cell_vec = target_cell_vec,
                       n_nonzero_trt_thresh = n_nonzero_trt_thresh, n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh)

    # convert to a vector
    QC_prob <- as.vector(QC_prob)

  }

  # output depending the intermediate outcome is required or not
  if(intermediate_outcome){

    # create a dataframe with outered names
    result_df <- data.frame(
      pair = E2G_pair,
      asy_mean = asy_mean,
      asy_sd = asy_sd,
      QC_prob = QC_prob
    )

  }else{

    # compute the adjusted power and discovery set with QC probability
    adjusted_power_list <- adjusted_power(mean_list = asy_mean,
                                          sd_list = asy_sd,
                                          multiple_testing_alpha = multiple_testing_alpha,
                                          multiple_testing_method = multiple_testing_method,
                                          side = side,
                                          QC_prob = QC_prob,
                                          cutoff = cutoff)

    # create a dataframe with names
    result_df <- list(
      individual_power = data.frame(
        pair = E2G_pair,
        adjusted_power = adjusted_power_list$adjusted_power
      ),
      num_discovery = adjusted_power_list$discovery_size_estimate
    )
  }

  # return the dataframe
  return(result_df)
}
