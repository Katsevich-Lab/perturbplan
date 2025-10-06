#' A549 Sakellaropoulos Chromium reference data for CRISPR power analysis
#'
#' @description
#' Pre-computed pilot data from A549 lung adenocarcinoma cells generated using
#' 10x Chromium technology. Contains baseline gene expression parameters and
#' library size information for power analysis of CRISPR-based perturbation
#' experiments.
#'
#' @format A list with 3 elements:
#' \describe{
#'   \item{baseline_expression_stats}{Data frame with gene expression data (20,406 genes × 3 columns):
#'     \itemize{
#'       \item \code{response_id}: Character vector of Ensembl gene IDs
#'       \item \code{relative_expression}: Numeric vector of relative expression levels (TPM/1e6 scale)
#'       \item \code{expression_size}: Numeric vector of dispersion parameters (theta)
#'     }
#'   }
#'   \item{library_parameters}{List containing:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Maximum UMI per cell parameter (42,377)
#'       \item \code{variation}: Variation parameter for PCR bias (0.376)
#'     }
#'   }
#'   \item{mapping_efficiency}{Numeric. Mapping efficiency value (0.794)}
#' }
#'
#' @details
#' This dataset was generated from A549 lung adenocarcinoma cells using single-cell
#' RNA sequencing with 10x Chromium technology.
#'
#' \strong{Cells Used in Relative Expression Estimate:} All cells in high-moi condition
#'
#' The data has been processed to extract key statistical parameters needed for
#' power analysis calculations.
#'
#' @source
#' \strong{Paper:} MethNet: a robust approach to identify regulatory hubs and their distal targets
#' from cancer methylomes
#'
#' \strong{Author and Year:} Sakellaropoulos et al., 2024
#'
#' \strong{Journal:} Nature Communications
#'
#' \strong{Accession:} GSE236304
#'
#' \strong{PMID:} 37577603
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for accessing this data programmatically
#'
#' @examples
#' data(A549_Sakellaropoulos)
#' str(A549_Sakellaropoulos)
"A549_Sakellaropoulos"

#' K562 Gasperini reference data for CRISPR power analysis
#'
#' @description
#' Pre-computed pilot data from K562 chronic myelogenous leukemia cells generated
#' using 10x Chromium technology. Contains baseline gene expression parameters and
#' library size information for power analysis of CRISPR-based perturbation
#' experiments.
#'
#' @format A list with 3 elements:
#' \describe{
#'   \item{baseline_expression_stats}{Data frame with gene expression data (19,942 genes × 3 columns):
#'     \itemize{
#'       \item \code{response_id}: Character vector of Ensembl gene IDs
#'       \item \code{relative_expression}: Numeric vector of relative expression levels (TPM/1e6 scale)
#'       \item \code{expression_size}: Numeric vector of dispersion parameters (theta)
#'     }
#'   }
#'   \item{library_parameters}{List containing:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Maximum UMI per cell parameter (59,163)
#'       \item \code{variation}: Variation parameter for PCR bias (0.397)
#'     }
#'   }
#'   \item{mapping_efficiency}{Numeric. Mapping efficiency value (0.711)}
#' }
#'
#' @details
#' This dataset was generated from K562 chronic myelogenous leukemia cells using
#' single-cell RNA sequencing with 10x Chromium technology.
#'
#' \strong{Cells Used in Relative Expression Estimate:} All cells in high-moi condition
#'
#' @source
#' \strong{Paper:} A Genome-wide Framework for Mapping Gene Regulation via Cellular Genetic Screens
#'
#' \strong{Author and Year:} Gasperini et al., 2019
#'
#' \strong{Journal:} Cell
#'
#' \strong{Accession:} GSE120861
#'
#' \strong{PMID:} 30612741
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for accessing this data programmatically
#'
#' @examples
#' data(K562_Gasperini)
#' str(K562_Gasperini)
"K562_Gasperini"

#' K562 10x Genomics reference data for CRISPR power analysis
#'
#' @description
#' Pre-computed pilot data from K562 chronic myelogenous leukemia cells generated
#' using 10x Chromium technology. This is a reference dataset from 10x Genomics.
#' Contains baseline gene expression parameters and library size information for
#' power analysis of CRISPR-based perturbation experiments.
#'
#' @format A list with 3 elements:
#' \describe{
#'   \item{baseline_expression_stats}{Data frame with gene expression data (16,395 genes × 3 columns):
#'     \itemize{
#'       \item \code{response_id}: Character vector of Ensembl gene IDs
#'       \item \code{relative_expression}: Numeric vector of relative expression levels (TPM/1e6 scale)
#'       \item \code{expression_size}: Numeric vector of dispersion parameters (theta)
#'     }
#'   }
#'   \item{library_parameters}{List containing:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Maximum UMI per cell parameter (61,081)
#'       \item \code{variation}: Variation parameter for PCR bias (0.421)
#'     }
#'   }
#'   \item{mapping_efficiency}{Numeric. Mapping efficiency value (0.801)}
#' }
#'
#' @details
#' This dataset was generated from K562 chronic myelogenous leukemia cells using
#' single-cell RNA sequencing with 10x Chromium technology.
#'
#' \strong{Cells Used in Relative Expression Estimate:} All cells in high-moi condition
#'
#' @source
#' \strong{Paper:} 10X Genomics dataset: 10k K562 cells
#'
#' \strong{Author and Year:} 10x Genomics (example data for K562), 2018
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for accessing this data programmatically
#'
#' @examples
#' data(K562_10x)
#' str(K562_10x)
"K562_10x"

#' THP-1 Yao reference data for CRISPR power analysis
#'
#' @description
#' Pre-computed pilot data from THP-1 monocytic leukemia cells generated using
#' 10x Chromium technology. Contains baseline gene expression parameters and
#' library size information for power analysis of CRISPR-based perturbation
#' experiments.
#'
#' @format A list with 3 elements:
#' \describe{
#'   \item{baseline_expression_stats}{Data frame with gene expression data (16,373 genes × 3 columns):
#'     \itemize{
#'       \item \code{response_id}: Character vector of Ensembl gene IDs
#'       \item \code{relative_expression}: Numeric vector of relative expression levels (TPM/1e6 scale)
#'       \item \code{expression_size}: Numeric vector of dispersion parameters (theta)
#'     }
#'   }
#'   \item{library_parameters}{List containing:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Maximum UMI per cell parameter (77,799)
#'       \item \code{variation}: Variation parameter for PCR bias (0.354)
#'     }
#'   }
#'   \item{mapping_efficiency}{Numeric. Mapping efficiency value (0.677)}
#' }
#'
#' @details
#' This dataset was generated from THP-1 monocytic leukemia cells using single-cell
#' RNA sequencing with 10x Chromium technology.
#'
#' \strong{Cells Used in Relative Expression Estimate:} NT (non-targeting) cells in low-moi condition
#'
#' @source
#' \strong{Paper:} Scalable genetic screening for regulatory circuits using compressed
#' Perturb-seq
#'
#' \strong{Author and Year:} Yao et al., 2023
#'
#' \strong{Journal:} Nature Biotechnology
#'
#' \strong{Accession:} GSE221321
#'
#' \strong{PMID:} 36747806
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for accessing this data programmatically
#'
#' @examples
#' data(THP1_Yao)
#' str(THP1_Yao)
"THP1_Yao"

#' T_CD8 Shifrut reference data for CRISPR power analysis
#'
#' @description
#' Pre-computed pilot data from primary CD8+ T cells generated using 10x Chromium
#' technology. Contains baseline gene expression parameters and library size
#' information for power analysis of CRISPR-based perturbation experiments.
#'
#' @format A list with 3 elements:
#' \describe{
#'   \item{baseline_expression_stats}{Data frame with gene expression data (20,935 genes × 3 columns):
#'     \itemize{
#'       \item \code{response_id}: Character vector of Ensembl gene IDs
#'       \item \code{relative_expression}: Numeric vector of relative expression levels (TPM/1e6 scale)
#'       \item \code{expression_size}: Numeric vector of dispersion parameters (theta)
#'     }
#'   }
#'   \item{library_parameters}{List containing:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Maximum UMI per cell parameter (8,801)
#'       \item \code{variation}: Variation parameter for PCR bias (0.297)
#'     }
#'   }
#'   \item{mapping_efficiency}{Numeric. Mapping efficiency value (0.679)}
#' }
#'
#' @details
#' This dataset was generated from primary CD8+ T cells using single-cell RNA
#' sequencing with 10x Chromium technology.
#'
#' \strong{Cells Used in Relative Expression Estimate:} All cells in high-moi condition
#'
#' @source
#' \strong{Paper:} Genome-wide CRISPR Screens in Primary Human T Cells Reveal Key
#' Regulators of Immune Function
#'
#' \strong{Author and Year:} Shifrut et al., 2018
#'
#' \strong{Journal:} Cell
#'
#' \strong{Accession:} GSE119450
#'
#' \strong{PMID:} 30449619
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for accessing this data programmatically
#'
#' @examples
#' data(T_CD8_Shifrut)
#' str(T_CD8_Shifrut)
"T_CD8_Shifrut"

#' iPSC Tian reference data for CRISPR power analysis
#'
#' @description
#' Pre-computed pilot data from induced pluripotent stem cells (iPSCs) generated
#' using 10x Chromium technology. Contains baseline gene expression parameters and
#' library size information for power analysis of CRISPR-based perturbation
#' experiments.
#'
#' @format A list with 3 elements:
#' \describe{
#'   \item{baseline_expression_stats}{Data frame with gene expression data (21,207 genes × 3 columns):
#'     \itemize{
#'       \item \code{response_id}: Character vector of Ensembl gene IDs
#'       \item \code{relative_expression}: Numeric vector of relative expression levels (TPM/1e6 scale)
#'       \item \code{expression_size}: Numeric vector of dispersion parameters (theta)
#'     }
#'   }
#'   \item{library_parameters}{List containing:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Maximum UMI per cell parameter (39,079)
#'       \item \code{variation}: Variation parameter for PCR bias (0.405)
#'     }
#'   }
#'   \item{mapping_efficiency}{Numeric. Mapping efficiency value (0.704)}
#' }
#'
#' @details
#' This dataset was generated from induced pluripotent stem cells (iPSCs) using
#' single-cell RNA sequencing with 10x Chromium technology.
#'
#' \strong{Cells Used in Relative Expression Estimate:} All cells in high-moi condition
#'
#' @source
#' \strong{Paper:} CRISPR Interference-Based Platform for Multimodal Genetic Screens
#' in Human iPSC-Derived Neurons
#'
#' \strong{Author and Year:} Tian et al., 2019
#'
#' \strong{Journal:} Neuron
#'
#' \strong{Accession:} GSE124703
#'
#' \strong{PMID:} 31422865
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for accessing this data programmatically
#'
#' @examples
#' data(iPSC_Tian)
#' str(iPSC_Tian)
"iPSC_Tian"

#' iPSC-derived neuron Tian reference data for CRISPR power analysis
#'
#' @description
#' Pre-computed pilot data from iPSC-derived neurons generated using 10x Chromium
#' technology. Contains baseline gene expression parameters and library size
#' information for power analysis of CRISPR-based perturbation experiments.
#'
#' @format A list with 3 elements:
#' \describe{
#'   \item{baseline_expression_stats}{Data frame with gene expression data (23,882 genes × 3 columns):
#'     \itemize{
#'       \item \code{response_id}: Character vector of Ensembl gene IDs
#'       \item \code{relative_expression}: Numeric vector of relative expression levels (TPM/1e6 scale)
#'       \item \code{expression_size}: Numeric vector of dispersion parameters (theta)
#'     }
#'   }
#'   \item{library_parameters}{List containing:
#'     \itemize{
#'       \item \code{UMI_per_cell}: Maximum UMI per cell parameter (32,285)
#'       \item \code{variation}: Variation parameter for PCR bias (0.443)
#'     }
#'   }
#'   \item{mapping_efficiency}{Numeric. Mapping efficiency value (0.614)}
#' }
#'
#' @details
#' This dataset was generated from iPSC-derived neurons using single-cell RNA
#' sequencing with 10x Chromium technology.
#'
#' \strong{Cells Used in Relative Expression Estimate:} All cells in high-moi condition
#'
#' @source
#' \strong{Paper:} CRISPR Interference-Based Platform for Multimodal Genetic Screens
#' in Human iPSC-Derived Neurons
#'
#' \strong{Author and Year:} Tian et al., 2019
#'
#' \strong{Journal:} Neuron
#'
#' \strong{Accession:} GSE124703
#'
#' \strong{PMID:} 31422865
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for accessing this data programmatically
#'
#' @examples
#' data(iPSC_neuron_Tian)
#' str(iPSC_neuron_Tian)
"iPSC_neuron_Tian"

#' Reference expression datasets metadata
#'
#' @description
#' Metadata table describing the available reference expression datasets in the
#' perturbplan package. This table maps biological systems to their corresponding
#' data processing configurations and functions.
#'
#' @format A data frame with 6 rows and 4 columns:
#' \describe{
#'   \item{cell_type}{Character. The biological system name}
#'   \item{platform}{Character. The experimental platform used (all "10x" for 10x Chromium)}
#'   \item{config_name}{Character. Configuration variable name for data source paths}
#'   \item{process_function}{Character. Name of the processing function for each dataset}
#' }
#'
#' @details
#' This metadata table is used internally by \code{\link{get_pilot_data_from_package}}
#' to map biological system names to their corresponding data objects and processing
#' functions.
#'
#' @source
#' Internal metadata for perturbplan package data management.
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for using this metadata to access pilot data
#'
#' @keywords internal
"reference_expression_datasets"
