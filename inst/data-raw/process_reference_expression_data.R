library(perturbplan)
library(tibble)   # tribble()
library(purrr)    # pwalk()
library(usethis)  # use_data()

overwrite <- FALSE   # set TRUE to rebuild even if .rda exists
source("data-raw/cellwise_reference_processing.R")
# 1. Define and save your spec ------------------------------------------
reference_expression_datasets <- tribble(
  ~ dataset_name, ~cell_type,                       ~platform, ~config_name, ~process_function,
  "K562_Gasperini", "K562",                            "10x",      "LOCAL_GASPERINI_2019_SRA_DATA_DIR", "process_k562_gasperini",
  "K562_10x", "K562",                            "10x",      "LOCAL_10X_2018_DATA_DIR", "process_k562_10x",
  "THP1_Yao", "THP-1",                           "10x",      "LOCAL_YAO_2023_DATA_DIR", "process_thp1_10x",
  "T_CD8_Shifrut", "T_CD8",                         "10x",      "LOCAL_SHIFRUT_2018_DATA_DIR", "process_t_cd8_10x",
  "A549_Sakellaropoulos", "A549",                           "10x",      "LOCAL_SAKELLAROPOULOS_2024_RAW_DATA_DIR", "process_a549_10x",
  "iPSC_Tian", "iPSC",                           "10x",      "LOCAL_TIAN_2019_RAW_DATA_DIR_1", "process_ipsc_10x",
  "iPSC_neuron_Tian", "iPSC_neuron",                    "10x",      "LOCAL_TIAN_2019_RAW_DATA_DIR_2", "process_ipsc_neuron_10x"
)

rda_file <- file.path("data", "reference_expression_datasets.rda")
if (!file.exists(rda_file)) {
  usethis::use_data(reference_expression_datasets, overwrite = overwrite)
} else {
  message("reference_expression_datasets.rda exists; skipping.")
}

# 2. Processing loop ----------------------------------------------------
pwalk(
  reference_expression_datasets,
  function(dataset_name, cell_type, platform, config_name, process_function) {

    rda_file <- file.path("data", paste0(dataset_name, ".rda"))

    if (!overwrite && file.exists(rda_file)) {
      message(dataset_name, " exists; skipping.")
      return(invisible(NULL))
    }

    # Locate path
    path <- .get_config_path(config_name)

    # Dynamically call the processing function
    obj <- get(process_function)(path)

    # Save object
    assign(package_name, obj, envir = .GlobalEnv)
    save(list = dataset_name,
         file = file.path("data", paste0(dataset_name, ".rda")),
         envir = .GlobalEnv, compress = "xz")
    # do.call(usethis::use_data, list(list = package_name, overwrite = TRUE, envir = .GlobalEnv))
  }
)
