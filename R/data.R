#' A549 10x Chromium reference data for CRISPR power analysis
#'
#' @description
#' Pre-computed pilot data from A549 lung adenocarcinoma cells generated using 
#' 10x Chromium technology. Contains baseline gene expression parameters and 
#' library size information for power analysis of CRISPR-based perturbation 
#' experiments.
#'
#' @format A list with 2 elements:
#' \describe{
#'   \item{baseline_expression}{List containing:
#'     \itemize{
#'       \item \code{baseline_expression}: Data frame with gene expression data (14,179 genes × 3 columns):
#'         \itemize{
#'           \item \code{response_id}: Character vector of Ensembl gene IDs
#'           \item \code{relative_expression}: Numeric vector of relative expression levels (TPM/1e6 scale)
#'           \item \code{expression_size}: Numeric vector of dispersion parameters (theta)
#'         }
#'       \item \code{expression_dispersion_curve}: Function for predicting dispersion from expression level
#'     }
#'   }
#'   \item{library_parameters}{List containing:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Maximum UMI per cell parameter (42,374)
#'       \item \code{variation}: Variation parameter for PCR bias (0.376)
#'     }
#'   }
#' }
#'
#' @details
#' This dataset was generated from A549 lung adenocarcinoma cells using single-cell 
#' RNA sequencing with 10x Chromium technology. The data has been processed to 
#' extract key statistical parameters needed for power analysis calculations:
#'
#' \itemize{
#'   \item Gene expression levels filtered for genes with TPM >= 1
#'   \item Dispersion parameters estimated using negative binomial models
#'   \item Library size parameters fitted using saturation-magnitude curves
#' }
#'
#' @source
#' Processed from single-cell RNA-seq data of A549 cells using 10x Chromium technology.
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for accessing this data programmatically
#' \code{\link{extract_fc_expression_info}} for using this data in power analysis
#'
#' @examples
#' # Load A549 data
#' data(A549_10x)
#' 
#' # Examine structure
#' str(A549_10x, max.level = 2)
#' 
#' # Access baseline expression data
#' baseline <- A549_10x$baseline_expression$baseline_expression
#' head(baseline)
#' 
#' # Access library parameters
#' lib_params <- A549_10x$library_parameters
#' print(lib_params)
#' 
#' # Use with power analysis functions
#' fc_info <- extract_fc_expression_info(
#'   minimum_fold_change = 0.8,
#'   gRNA_variability = 0.1,
#'   biological_system = "A549"
#' )
"A549_10x"

#' K562 10x Chromium reference data for CRISPR power analysis
#'
#' @description
#' Pre-computed pilot data from K562 chronic myelogenous leukemia cells generated 
#' using 10x Chromium technology. Contains baseline gene expression parameters and 
#' library size information for power analysis of CRISPR-based perturbation 
#' experiments.
#'
#' @format A list with 2 elements:
#' \describe{
#'   \item{baseline_expression}{List containing:
#'     \itemize{
#'       \item \code{baseline_expression}: Data frame with gene expression data
#'       \item \code{expression_dispersion_curve}: Function for predicting dispersion from expression level
#'     }
#'   }
#'   \item{library_parameters}{List containing:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Maximum UMI per cell parameter
#'       \item \code{variation}: Variation parameter for PCR bias
#'     }
#'   }
#' }
#'
#' @details
#' This dataset was generated from K562 chronic myelogenous leukemia cells using 
#' single-cell RNA sequencing with 10x Chromium technology. K562 is a widely used 
#' cell line in CRISPR screening experiments.
#'
#' @source
#' Processed from single-cell RNA-seq data of K562 cells using 10x Chromium technology.
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for accessing this data programmatically
#' \code{\link{extract_fc_expression_info}} for using this data in power analysis
#'
#' @examples
#' # Load K562 data
#' data(K562_10x)
#' 
#' # Use with power analysis functions
#' fc_info <- extract_fc_expression_info(
#'   minimum_fold_change = 0.8,
#'   gRNA_variability = 0.1,
#'   biological_system = "K562"
#' )
"K562_10x"

#' THP-1 10x Chromium reference data for CRISPR power analysis
#'
#' @description
#' Pre-computed pilot data from THP-1 monocytic leukemia cells generated using 
#' 10x Chromium technology. Contains baseline gene expression parameters and 
#' library size information for power analysis of CRISPR-based perturbation 
#' experiments.
#'
#' @format A list with 2 elements:
#' \describe{
#'   \item{baseline_expression}{List containing:
#'     \itemize{
#'       \item \code{baseline_expression}: Data frame with gene expression data
#'       \item \code{expression_dispersion_curve}: Function for predicting dispersion from expression level
#'     }
#'   }
#'   \item{library_parameters}{List containing:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Maximum UMI per cell parameter
#'       \item \code{variation}: Variation parameter for PCR bias
#'     }
#'   }
#' }
#'
#' @details
#' This dataset was generated from THP-1 monocytic leukemia cells using single-cell 
#' RNA sequencing with 10x Chromium technology. THP-1 cells are commonly used to 
#' study monocyte and macrophage biology.
#'
#' @source
#' Processed from single-cell RNA-seq data of THP-1 cells using 10x Chromium technology.
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for accessing this data programmatically
#' \code{\link{extract_fc_expression_info}} for using this data in power analysis
#'
#' @examples
#' # Load THP-1 data
#' data(THP_1_10x)
#' 
#' # Use with power analysis functions
#' fc_info <- extract_fc_expression_info(
#'   minimum_fold_change = 0.8,
#'   gRNA_variability = 0.1,
#'   biological_system = "THP-1"
#' )
"THP_1_10x"

#' T_CD8 10x Chromium reference data for CRISPR power analysis
#'
#' @description
#' Pre-computed pilot data from primary CD8+ T cells generated using 10x Chromium 
#' technology. Contains baseline gene expression parameters and library size 
#' information for power analysis of CRISPR-based perturbation experiments.
#'
#' @format A list with 2 elements:
#' \describe{
#'   \item{baseline_expression}{List containing:
#'     \itemize{
#'       \item \code{baseline_expression}: Data frame with gene expression data
#'       \item \code{expression_dispersion_curve}: Function for predicting dispersion from expression level
#'     }
#'   }
#'   \item{library_parameters}{List containing:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Maximum UMI per cell parameter
#'       \item \code{variation}: Variation parameter for PCR bias
#'     }
#'   }
#' }
#'
#' @details
#' This dataset was generated from primary CD8+ T cells using single-cell RNA 
#' sequencing with 10x Chromium technology. CD8+ T cells are important effector 
#' cells in adaptive immunity and cancer immunotherapy.
#'
#' @source
#' Processed from single-cell RNA-seq data of primary CD8+ T cells using 10x Chromium technology.
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for accessing this data programmatically
#' \code{\link{extract_fc_expression_info}} for using this data in power analysis
#'
#' @examples
#' # Load T_CD8 data
#' data(T_CD8_10x)
#' 
#' # Use with power analysis functions
#' fc_info <- extract_fc_expression_info(
#'   minimum_fold_change = 0.8,
#'   gRNA_variability = 0.1,
#'   biological_system = "T_CD8"
#' )
"T_CD8_10x"

#' iPSC 10x Chromium reference data for CRISPR power analysis
#'
#' @description
#' Pre-computed pilot data from induced pluripotent stem cells (iPSCs) generated 
#' using 10x Chromium technology. Contains baseline gene expression parameters and 
#' library size information for power analysis of CRISPR-based perturbation 
#' experiments.
#'
#' @format A list with 2 elements:
#' \describe{
#'   \item{baseline_expression}{List containing:
#'     \itemize{
#'       \item \code{baseline_expression}: Data frame with gene expression data
#'       \item \code{expression_dispersion_curve}: Function for predicting dispersion from expression level
#'     }
#'   }
#'   \item{library_parameters}{List containing:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Maximum UMI per cell parameter
#'       \item \code{variation}: Variation parameter for PCR bias
#'     }
#'   }
#' }
#'
#' @details
#' This dataset was generated from induced pluripotent stem cells (iPSCs) using 
#' single-cell RNA sequencing with 10x Chromium technology. iPSCs are valuable 
#' for studying development and disease modeling.
#'
#' @source
#' Processed from single-cell RNA-seq data of iPSCs using 10x Chromium technology.
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for accessing this data programmatically
#' \code{\link{extract_fc_expression_info}} for using this data in power analysis
#'
#' @examples
#' # Load iPSC data
#' data(iPSC_10x)
#' 
#' # Use with power analysis functions
#' fc_info <- extract_fc_expression_info(
#'   minimum_fold_change = 0.8,
#'   gRNA_variability = 0.1,
#'   biological_system = "iPSC"
#' )
"iPSC_10x"

#' iPSC-derived neuron 10x Chromium reference data for CRISPR power analysis
#'
#' @description
#' Pre-computed pilot data from iPSC-derived neurons generated using 10x Chromium 
#' technology. Contains baseline gene expression parameters and library size 
#' information for power analysis of CRISPR-based perturbation experiments.
#'
#' @format A list with 2 elements:
#' \describe{
#'   \item{baseline_expression}{List containing:
#'     \itemize{
#'       \item \code{baseline_expression}: Data frame with gene expression data
#'       \item \code{expression_dispersion_curve}: Function for predicting dispersion from expression level
#'     }
#'   }
#'   \item{library_parameters}{List containing:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Maximum UMI per cell parameter from S-M curve
#'       \item \code{variation}: Variation parameter characterizing PCR amplification bias
#'     }
#'   }
#' }
#'
#' @source
#' Processed from single-cell RNA-seq data of iPSC-derived neurons using 10x Chromium technology.
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for accessing this data programmatically
#' \code{\link{extract_fc_expression_info}} for using this data in power analysis
#'
#' @examples
#' # Load iPSC-derived neuron data
#' data(iPSC_neuron_10x)
#' 
#' # Examine baseline expression data
#' baseline_expr <- iPSC_neuron_10x$baseline_expression$baseline_expression
#' print(head(baseline_expr))
#' 
#' # Check library parameters
#' lib_params <- iPSC_neuron_10x$library_parameters
#' print(lib_params)
#' 
#' # Use with power analysis functions
#' fc_info <- extract_fc_expression_info(
#'   minimum_fold_change = 0.8,
#'   gRNA_variability = 0.1,
#'   biological_system = "iPSC_neuron"
#' )
"iPSC_neuron_10x"

#' Reference expression datasets metadata
#'
#' @description
#' Metadata table describing the available reference expression datasets in the 
#' perturbplan package. This table maps biological systems to their corresponding 
#' data processing configurations and functions.
#'
#' @format A data frame with 6 rows and 4 columns:
#' \describe{
#'   \item{cell_type}{Character. The biological system name (K562, THP-1, T_CD8, A549, iPSC, iPSC_neuron)}
#'   \item{platform}{Character. The experimental platform used (all "10x" for 10x Chromium)}
#'   \item{config_name}{Character. Configuration variable name for data source paths}
#'   \item{process_function}{Character. Name of the processing function for each dataset}
#' }
#'
#' @details
#' This metadata table is used internally by \code{\link{get_pilot_data_from_package}} 
#' to map biological system names to their corresponding data objects and processing 
#' functions. Each row represents one supported biological system:
#'
#' \itemize{
#'   \item \strong{K562}: Chronic myelogenous leukemia cells
#'   \item \strong{THP-1}: Monocytic leukemia cells  
#'   \item \strong{T_CD8}: Primary CD8+ T cells
#'   \item \strong{A549}: Lung adenocarcinoma cells
#'   \item \strong{iPSC}: Induced pluripotent stem cells
#' }
#'
#' @source
#' Internal metadata for perturbplan package data management.
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for using this metadata to access pilot data
#'
#' @examples
#' # Load metadata
#' data(reference_expression_datasets)
#' print(reference_expression_datasets)
#' 
#' # View available biological systems
#' available_systems <- reference_expression_datasets$cell_type
#' print(available_systems)
"reference_expression_datasets"