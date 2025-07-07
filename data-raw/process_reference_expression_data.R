library(perturbplan)
library(tibble)   # tribble()
library(purrr)    # pwalk()
library(usethis)  # use_data()

overwrite <- FALSE   # set TRUE to rebuild even if .rda exists

# 1. Define and save your spec ------------------------------------------
reference_expression_datasets <- tribble(
  ~cell_type,                       ~platform, ~config_name, ~process_function,
  "K562",                            "10x",      "LOCAL_10X_2022_DATA_DIR", "process_k562_10x",
  # "MCF 7 (treated with IFN-gamma)", "Parse",    "MCF_IFN_2021", "process_mcf7_parse",
  "THP-1",                           "10x",    "LOCAL_YAO_2023_DATA_DIR", "process_thp1_10x",
)

# Save the spec into data/reference_expression_datasets.rda
usethis::use_data(reference_expression_datasets, overwrite = overwrite)

# 2. Processing loop ----------------------------------------------------
pwalk(
  reference_expression_datasets,
  function(cell_type, platform, config_name, overwrite) {
    
    # Build the sanitized object name
    package_name <- paste0(cell_type, "_", platform) |>
      stringr::str_replace_all("[^A-Za-z0-9_]", "_")
    
    rda_file <- file.path("data", paste0(package_name, ".rda"))
    
    if (!overwrite && file.exists(rda_file)) {
      message(package_name, " exists; skipping.")
      return(invisible(NULL))
    }
    
    path <- .get_config_path(config_name)
    
    obj <- if (tolower(platform) == "10x") {
      reference_data_preprocessing_10x(path)
    } else {
      ds <- get_reference_datasets(path)
      reference_data_preprocessing(ds$response_matrix, ds$read_umi_table)
    }
    
    assign(package_name, obj, envir = .GlobalEnv)
    usethis::use_data(list = package_name, overwrite = overwrite)
  },
  overwrite = overwrite
)