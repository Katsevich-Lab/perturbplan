library(perturbplan)
library(tibble)   # tribble()
library(purrr)    # pwalk()
library(usethis)  # use_data()

overwrite <- FALSE   # set TRUE to rebuild even if .rda exists
source("data-raw/cellwise_reference_processing.R")
# 1. Define and save your spec ------------------------------------------
reference_expression_datasets <- tribble(
  ~cell_type,                       ~platform, ~config_name, ~process_function,
  "K562",                            "10x",      "LOCAL_GASPERINI_2019_SRA_DATA_DIR", "process_k562_10x",
  "THP-1",                           "10x",      "LOCAL_YAO_2023_DATA_DIR", "process_thp1_10x",
  "T_CD8",                         "10x",      "LOCAL_SHIFRUT_2018_DATA_DIR", "process_t_cd8_10x",
  "A549",                           "10x",      "LOCAL_SAKELLAROPOULOS_2024_RAW_DATA_DIR", "process_a549_10x",
  "iPSC",                           "10x",      "LOCAL_TIAN_2019_RAW_DATA_DIR_1", "process_ipsc_10x",
  "iPSC_neuron",                    "10x",      "LOCAL_TIAN_2019_RAW_DATA_DIR_2", "process_ipsc_neuron_10x"
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
  function(cell_type, platform, config_name, process_function) {

    # Create a valid object name for the dataset
    package_name <- paste0(cell_type, "_", platform) |>
      stringr::str_replace_all("[^A-Za-z0-9_]", "_")

    rda_file <- file.path("data", paste0(package_name, ".rda"))

    if (!overwrite && file.exists(rda_file)) {
      message(package_name, " exists; skipping.")
      return(invisible(NULL))
    }

    # Locate path
    path <- .get_config_path(config_name)

    # Dynamically call the processing function
    obj <- get(process_function)(path)

    # Save object
    assign(package_name, obj, envir = .GlobalEnv)
    save(list = package_name,
         file = file.path("data", paste0(package_name, ".rda")),
         envir = .GlobalEnv)
    # do.call(usethis::use_data, list(list = package_name, overwrite = TRUE, envir = .GlobalEnv))
  }
)
