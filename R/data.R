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
#'   \item{library_parameters}{List containing saturation curve parameters estimated using preseqR.
#'     Data-adaptively selected method: ZTNB. Contains:
#'     \itemize{
#'       \item \code{method_used}: "ZTNB" - Zero-truncated negative binomial method
#'       \item \code{L}: 1.19e9 - Total expected distinct UMIs at saturation
#'       \item \code{size}: 1.64 - ZTNB shape parameter
#'       \item \code{mu}: 0.244 - ZTNB mean parameter
#'       \item \code{reads_norm}: 11,759 - Normalization constant (reads per cell)
#'       \item \code{n_cells}: 24,744 - Number of cells in pilot data
#'       \item \code{UMI_per_cell_at_saturation}: 48,179 - Maximum UMI per cell at infinite sequencing depth
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
#'   \item{library_parameters}{List containing saturation curve parameters estimated using preseqR.
#'     Data-adaptively selected method: ZTNB. Contains:
#'     \itemize{
#'       \item \code{method_used}: "ZTNB" - Zero-truncated negative binomial method
#'       \item \code{L}: 4.68e8 - Total expected distinct UMIs at saturation
#'       \item \code{size}: 2.14 - ZTNB shape parameter
#'       \item \code{mu}: 0.395 - ZTNB mean parameter
#'       \item \code{reads_norm}: 23,835 - Normalization constant (reads per cell)
#'       \item \code{n_cells}: 7,749 - Number of cells in pilot data
#'       \item \code{UMI_per_cell_at_saturation}: 60,414 - Maximum UMI per cell at infinite sequencing depth
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
#'   \item{library_parameters}{List containing saturation curve parameters estimated using preseqR.
#'     Data-adaptively selected method: ZTNB. Contains:
#'     \itemize{
#'       \item \code{method_used}: "ZTNB" - Zero-truncated negative binomial method
#'       \item \code{L}: 4.86e8 - Total expected distinct UMIs at saturation
#'       \item \code{size}: 1.67 - ZTNB shape parameter
#'       \item \code{mu}: 0.451 - ZTNB mean parameter
#'       \item \code{reads_norm}: 28,745 - Normalization constant (reads per cell)
#'       \item \code{n_cells}: 7,632 - Number of cells in pilot data
#'       \item \code{UMI_per_cell_at_saturation}: 63,710 - Maximum UMI per cell at infinite sequencing depth
#'     }
#'   }
#'   \item{mapping_efficiency}{Numeric. Mapping efficiency value (0.700)}
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
#'   \item{baseline_expression_stats}{Data frame with gene expression data (16,366 genes × 3 columns):
#'     \itemize{
#'       \item \code{response_id}: Character vector of Ensembl gene IDs
#'       \item \code{relative_expression}: Numeric vector of relative expression levels (TPM/1e6 scale)
#'       \item \code{expression_size}: Numeric vector of dispersion parameters (theta)
#'     }
#'   }
#'   \item{library_parameters}{List containing saturation curve parameters estimated using preseqR.
#'     Data-adaptively selected method: ZTNB. Contains:
#'     \itemize{
#'       \item \code{method_used}: "ZTNB" - Zero-truncated negative binomial method
#'       \item \code{L}: 4.96e8 - Total expected distinct UMIs at saturation
#'       \item \code{size}: 1.79 - ZTNB shape parameter
#'       \item \code{mu}: 0.124 - ZTNB mean parameter
#'       \item \code{reads_norm}: 10,930 - Normalization constant (reads per cell)
#'       \item \code{n_cells}: 5,608 - Number of cells in pilot data
#'       \item \code{UMI_per_cell_at_saturation}: 88,390 - Maximum UMI per cell at infinite sequencing depth
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
#'   \item{library_parameters}{List containing saturation curve parameters estimated using preseqR.
#'     Data-adaptively selected method: ZTNB. Contains:
#'     \itemize{
#'       \item \code{method_used}: "ZTNB" - Zero-truncated negative binomial method
#'       \item \code{L}: 1.11e8 - Total expected distinct UMIs at saturation
#'       \item \code{size}: 2.69 - ZTNB shape parameter
#'       \item \code{mu}: 2.17 - ZTNB mean parameter
#'       \item \code{reads_norm}: 19,435 - Normalization constant (reads per cell)
#'       \item \code{n_cells}: 12,340 - Number of cells in pilot data
#'       \item \code{UMI_per_cell_at_saturation}: 8,957 - Maximum UMI per cell at infinite sequencing depth
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
#'   \item{library_parameters}{List containing saturation curve parameters estimated using preseqR.
#'     Data-adaptively selected method: ZTNB. Contains:
#'     \itemize{
#'       \item \code{method_used}: "ZTNB" - Zero-truncated negative binomial method
#'       \item \code{L}: 5.16e8 - Total expected distinct UMIs at saturation
#'       \item \code{size}: 2.28 - ZTNB shape parameter
#'       \item \code{mu}: 0.55 - ZTNB mean parameter
#'       \item \code{reads_norm}: 21,301 - Normalization constant (reads per cell)
#'       \item \code{n_cells}: 13,321 - Number of cells in pilot data
#'       \item \code{UMI_per_cell_at_saturation}: 38,699 - Maximum UMI per cell at infinite sequencing depth
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
#'   \item{library_parameters}{List containing saturation curve parameters estimated using preseqR.
#'     Data-adaptively selected method: ZTNB. Contains:
#'     \itemize{
#'       \item \code{method_used}: "ZTNB" - Zero-truncated negative binomial method
#'       \item \code{L}: 3.42e8 - Total expected distinct UMIs at saturation
#'       \item \code{size}: 2.55 - ZTNB shape parameter
#'       \item \code{mu}: 0.779 - ZTNB mean parameter
#'       \item \code{reads_norm}: 23,512 - Normalization constant (reads per cell)
#'       \item \code{n_cells}: 11,350 - Number of cells in pilot data
#'       \item \code{UMI_per_cell_at_saturation}: 30,173 - Maximum UMI per cell at infinite sequencing depth
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

#' K562 Ray TAP-seq reference data for CRISPR power analysis
#'
#' @description
#' Pre-computed pilot data from K562 chronic myelogenous leukemia cells generated
#' using TAP-seq (targeted perturb-seq) with 10x Chromium technology. Contains
#' baseline gene expression parameters and library size information for power
#' analysis of CRISPR-based perturbation experiments. TAP-seq uses targeted
#' sequencing to profile a focused gene panel, providing cost-effective power
#' analysis for experiments targeting specific pathways or gene sets.
#'
#' @format A list with 3 elements:
#' \describe{
#'   \item{baseline_expression_stats}{Data frame with gene expression data (303 genes × 3 columns):
#'     \itemize{
#'       \item \code{response_id}: Character vector of Ensembl gene IDs
#'       \item \code{relative_expression}: Numeric vector of relative expression levels (TPM/1e6 scale)
#'       \item \code{expression_size}: Numeric vector of dispersion parameters (theta)
#'     }
#'   }
#'   \item{library_parameters}{List containing saturation curve parameters estimated using preseqR.
#'     Data-adaptively selected method: RFA (Rational Function Approximation). Contains:
#'     \itemize{
#'       \item \code{method_used}: "RFA" - Rational function approximation method
#'       \item \code{valid_estimator}: TRUE - RFA estimator is valid
#'       \item \code{coefs_real}: Real parts of RFA coefficients
#'       \item \code{coefs_imag}: Imaginary parts of RFA coefficients
#'       \item \code{poles_real}: Real parts of RFA poles
#'       \item \code{poles_imag}: Imaginary parts of RFA poles
#'       \item \code{reads_norm}: 16,790 - Normalization constant (reads per cell)
#'       \item \code{n_cells}: 8,278 - Number of cells in pilot data
#'       \item \code{UMI_per_cell_at_saturation}: 3,684 - Maximum UMI per cell at infinite sequencing depth
#'     }
#'   }
#'   \item{mapping_efficiency}{Numeric. Mapping efficiency value (0.371)}
#' }
#'
#' @details
#' This dataset was generated using DC TAP-seq (Direct-Capture Targeted Perturb-seq),
#' an enhanced version of targeted perturb-seq that integrates CRISPR-based perturbations
#' with direct-capture single-cell RNA sequencing. By capturing guide RNAs alongside
#' targeted gene transcripts within the same sequencing reaction, DC TAP-seq enables
#' high-throughput, unbiased mapping of distal regulatory element–gene interactions
#' with improved sensitivity and reduced technical noise. This approach allows
#' simultaneous measurement of perturbation identity and gene expression in thousands
#' of single cells, facilitating large-scale functional dissection of noncoding regions
#' at single-cell resolution.
#'
#' \strong{Cells Used in Relative Expression Estimate:} All cells in high-moi condition
#'
#' @source
#' \strong{Paper:} An unbiased survey of distal element-gene regulatory interactions with direct-capture targeted Perturb-seq
#'
#' \strong{Author and Year:} Ray et al., 2025
#'
#' \strong{Journal:} (Publication details pending)
#'
#' \strong{Accession:} GSE303901
#'
#' \strong{PMID:} 41000760
#'
#' @seealso
#' \code{\link{get_pilot_data_from_package}} for accessing this data programmatically
#'
#' @examples
#' data(K562_Ray)
#' str(K562_Ray)
"K562_Ray"

#' Reference expression datasets metadata
#'
#' @description
#' Metadata table describing the available reference expression datasets in the
#' perturbplan package. This table maps biological systems to their corresponding
#' data processing configurations and functions.
#'
#' @format A data frame with 8 rows and 5 columns:
#' \describe{
#'   \item{dataset_name}{Character. The unique identifier for each dataset}
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
