# This is a Rscript computing the power function using score test

#' Power function for CRISPR screen experimental design.
#'
#'
#'
#' @param perturb_type Either CRISPRi or CRISPRko; currently can only handle CRISPRi.
#'
#' @param effect_size_mean Mean fold change for each (element, gene) (matrix; row: L genomic elements; column: J genes).
#' @param effect_size_sd Sd fold change for each (element, gene) (matrix; row: L genomic elements; column: J genes).
#' @param control_cell_vec Control cell size for L elements (vector; length L).
#' @param target_cell_mat Perturb cell size for L elements (matrix; row: L genomic elements; column: K gRNA libraries). Rownames should be specified for enhancer targets.
#' @param library_size Averaged library size (scalar; positive valued).
#' @param relative_expression Baseline relative expression level for J genes (named vector; of length J). Names of the vector should be different genes.
#' @param size_parameter Size parameter for J genes (vector; length J).
#'
#' @param sideness Left, right or both.
#' @param correction Multiplicity correction; default is BH.
#' @param sig_level Significance level of interest (FDR of interest).
#'
#' @param QC A logic value; if TRUE, then do quality control.
#' @param n_nonzero_trt QC threshold for treatment cell.
#' @param n_nonzero_ctl QC threshold for control cell.
#'
#' @param UMI_s Estimate for UMI count per singlet for the underlying type of cell.
#' @param recovery_rate Recovery rate for cells surviving the library preparation.
#' @param doublet_rate Doublet rate for droplet containing more than one cell.
#' @param doublet_factor Ratio of averaged UMI count per doublet to UMI per singlet.
#' @param planned_read Planned total sequencing reads.
#' @param mapping_efficiency Mapping efficiency for sequenced reads.
#'
#' @param intermediate_outcome A logic value indicating if only mean, sd of test statistics and QC are desired.
#'
#' @return The output format is the same as [compute_power()].
#' @importFrom dplyr if_else
#' @importFrom stats setNames
#' @export

power_function <- function(
    ######################## specify the type of perturb-seq experiment ########
    perturb_type = "CRISPRi",                                          # currenly only support CRISPRi

    ######################## specify the power-determining parameters ##########
    control_cell_vec, target_cell_mat,
    library_size = NULL, relative_expression, size_parameter,
    effect_size_mean, effect_size_sd = matrix(0, nrow = nrow(effect_size_mean), ncol = ncol(effect_size_mean)),

    ######################## specify sceptre-related parameters ################
    QC = TRUE, n_nonzero_trt = 7, n_nonzero_ctl = 7,

    ###################### specify experimental design parameters ##############
    UMI_s = NULL,                                                      # cell-type specific parameter
    recovery_rate = NULL, doublet_rate = NULL, doublet_factor = NULL,  # Library prep parameters
    planned_read = NULL, mapping_efficiency = NULL,                    # Sequencing parameters

    ###################### specify test-related parameters #####################
    sideness, correction = "BH", sig_level = 0.1,

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
    target_cell_vec <- rowSums(target_cell_mat)
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
                                control_cell_vec = control_cell_vec, target_cell_mat = target_cell_mat,
                                baseline_expression = baseline_expression, size_parameter = size_parameter,
                                effect_size_mean = effect_size_mean, effect_size_sd = effect_size_sd,
                                QC = QC, n_nonzero_trt = n_nonzero_trt, n_nonzero_ctl = n_nonzero_ctl,
                                sideness = sideness, correction = correction, sig_level = sig_level,
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
#' @param baseline_expression Baseline gene expression in the group of control cells.
#'
#' @return A dataframe or a list depending on `intermediate_outcome`:
#' @section Dataframe: if `intermediate_outcome` is TRUE, only mean, sd of the score tests and QC probability will be outputted.
#' @section List: otherwise, the final power for each element-gene pair will be outputted, together with the estimated discovery set.
#'
#' @export
compute_power <- function(perturb_type,
                          control_cell_vec, target_cell_mat,
                          baseline_expression, size_parameter,
                          effect_size_mean, effect_size_sd,
                          QC, n_nonzero_trt, n_nonzero_ctl,
                          sideness, correction = NULL, sig_level = NULL, cutoff = NULL,
                          intermediate_outcome = FALSE){

  # compute the number of total cells targeting certain element
  target_cell_vec <- rowSums(target_cell_mat)

  # compute the mean gene expression appearing in the power formula
  mean_expression <- mean_expression_computation(baseline_expression = baseline_expression,
                                                 effect_size_mean = effect_size_mean,
                                                 num_control = control_cell_vec,
                                                 num_trt = target_cell_vec)

  # divide the power computation depending on perturb_type
  if(perturb_type == "CRISPRi"){

    # compute the different variance quantities
    pooled_mean <- mean_expression$pooled_mat
    trt_mean <- mean_expression$trt_mat
    ctl_mean <- mean_expression$ctl_mat

    # compute asy_mean and asy_sd
    asy_quantity <- distribution_teststat(control_cell_vec, target_cell_mat, size_parameter,
                                          effect_size_mean, effect_size_sd,
                                          pooled_mean, trt_mean, ctl_mean)

    # obtain asy_mean and asy_sd
    asy_mean <- asy_quantity$mean
    asy_sd <- asy_quantity$sd

  }else{
    stop("Code is not finished for CRISPR perturbation beyond interference!")
  }

  # next step depends if QC_prob is used or not
  if(QC){

    # compute the QC_prob
    QC_prob <- QC_prob(effect_size_mean = effect_size_mean,
                       baseline_expression = baseline_expression,
                       size_parameter = size_parameter,
                       num_control = control_cell_vec, num_trt = target_cell_vec,
                       n_nonzero_trt = n_nonzero_trt, n_nonzero_ctl = n_nonzero_ctl)

  }else{

    # add QC probability
    QC_prob <- rep(0, length(asy_mean))

  }

  # extract gRNA and gene names and obtain Enhancer-Gene pair
  gRNA_name <- rownames(target_cell_mat)
  gene_name <- names(baseline_expression)
  E2G_pair <- as.vector(outer(gRNA_name, gene_name, paste, sep = "_"))

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
                                          sig_level = sig_level,
                                          correction = correction,
                                          sideness = sideness,
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
