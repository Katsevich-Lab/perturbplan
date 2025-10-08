# make_cellranger_tiny.R
# Build a tiny, self-consistent subset from a full Cell Ranger run:
# (1) outs/filtered_feature_bc_matrix/{matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz}
# (2) outs/molecule_info.h5  (minimal datasets: count, umi, barcode_idx, feature_idx, gem_group, barcodes, features/id)
# (3) outs/filtered_feature_bc_matrix.h5  (only matrix/barcodes)
# (4) outs/metrics_summary.csv  (all columns from a template; only "Number of Reads" filled with generated value)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(R.utils)
  library(rhdf5)
})

# ========= EDIT THESE PATHS =========
SRC_DIR       <- file.path(.get_config_path("LOCAL_GASPERINI_2019_SRA_DATA_DIR") ,"processed", "at_scale", "run1","SRR7967494")
DEST_DIR      <- "inst/extdata/cellranger_tiny" # destination tiny dataset root
TEMPLATE_CSV  <- file.path(SRC_DIR,"outs","metrics_summary.csv")
# ====================================

# ---- ensure destination dirs ----
dir.create(file.path(DEST_DIR, "outs", "filtered_feature_bc_matrix"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DEST_DIR, "outs"), showWarnings = FALSE)

# =========================
# (1) filtered_feature_bc_matrix: read -> tiny subset -> write
# =========================
src_matdir <- file.path(SRC_DIR, "outs", "filtered_feature_bc_matrix")
stopifnot(dir.exists(src_matdir))

mtx_path <- file.path(src_matdir, "matrix.mtx.gz")
fea_path <- file.path(src_matdir, "features.tsv.gz")
bar_path <- file.path(src_matdir, "barcodes.tsv.gz")

X_full   <- as(Matrix::readMM(mtx_path), "dgCMatrix")       # genes x cells
features <- data.table::fread(fea_path, header = FALSE)     # V1=id, V2=name, V3=feature_type
barcodes <- data.table::fread(bar_path, header = FALSE)[[1]]# robust gz read

stopifnot(nrow(X_full) == nrow(features), ncol(X_full) == length(barcodes))

# set target sizes
G <- min(5L, nrow(X_full))   # number of rows (genes) to keep
C <- min(8L, ncol(X_full))   # number of cols (cells) to keep

# choose rows that have at least one non-zero
row_has <- Matrix::rowSums(X_full) > 0
if (!any(row_has)) stop("No non-zero rows in input matrix.")
keep_genes <- which(row_has)[seq_len(min(G, sum(row_has)))]

# for each kept row, pick ONE non-zero column
#    - must be unique across rows (no duplicates allowed)
keep_cols <- integer(0)
for (r in keep_genes) {
  nz <- which(X_full[r, ] != 0)
  # exclude already used columns
  pick <- setdiff(nz, keep_cols)
  if (!length(pick)) {
    stop("Row ", r, " cannot be assigned a unique non-zero column.")
  }
  keep_cols <- c(keep_cols, pick[1L])
}

# if we still have fewer than C columns,
#    add extra non-zero columns (not used yet, no duplicates)
if (length(keep_cols) < C) {
  cand <- setdiff(which(Matrix::colSums(X_full[keep_genes, , drop=FALSE] != 0) > 0), keep_cols)
  if (length(cand)) {
    hits <- Matrix::colSums(X_full[keep_genes, cand, drop=FALSE] != 0)
    add  <- cand[order(-as.numeric(hits))]
    add  <- head(add, C - length(keep_cols))
    keep_cols <- c(keep_cols, add)
  }
}

X_small        <- X_full[keep_genes, keep_cols, drop = FALSE]
features_small <- features[keep_genes, , drop = FALSE]
barcodes_small <- barcodes[keep_cols]

X_small <- Matrix::drop0(X_small)              # drop explicit 0s
X_small <- methods::as(X_small, "dgCMatrix")   # force into dgCMatrix
stopifnot(sum(X_small != 0) > 0)  # ensure not all-zero

# write matrix.mtx.gz
mm_tmp <- tempfile(fileext = ".mtx")
Matrix::writeMM(X_small, file = mm_tmp)
mm_gz  <- file.path(DEST_DIR, "outs", "filtered_feature_bc_matrix", "matrix.mtx.gz")
R.utils::gzip(mm_tmp, destname = mm_gz, overwrite = TRUE)

# write features.tsv.gz  (no header)
fea_tmp <- tempfile(fileext = ".tsv")
data.table::fwrite(features_small, fea_tmp, sep = "\t", col.names = FALSE)
fea_gz  <- file.path(DEST_DIR, "outs", "filtered_feature_bc_matrix", "features.tsv.gz")
R.utils::gzip(fea_tmp, destname = fea_gz, overwrite = TRUE)

# write barcodes.tsv.gz  (single column, matches matrix columns)
bar_tmp <- tempfile(fileext = ".tsv")
writeLines(barcodes_small, con = bar_tmp)
bar_gz  <- file.path(DEST_DIR, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
R.utils::gzip(bar_tmp, destname = bar_gz, overwrite = TRUE)

# ============================================================
# (2) Make tiny molecule_info.h5 — filter rows by cell_id & gene_id
#      cell_id = paste(barcodes[barcode_idx+1], gem_group, sep="-")
#      gene_id = features$id[feature_idx+1]
#      Keep rows where cell_id ∈ barcodes_small  AND  gene_id ∈ features_small[[1]]
# ============================================================
mol_src <- file.path(SRC_DIR, "outs", "molecule_info.h5")
stopifnot(file.exists(mol_src))
mol_dst <- file.path(DEST_DIR, "outs", "molecule_info.h5")

# Read needed datasets from full molecule_info.h5
count        <- rhdf5::h5read(mol_src, "count")
umi          <- rhdf5::h5read(mol_src, "umi")
mol_barcodes <- rhdf5::h5read(mol_src, "barcodes")          # bare barcodes (no "-1")
barcode_idx  <- rhdf5::h5read(mol_src, "barcode_idx")       # 0-based
gem_group    <- rhdf5::h5read(mol_src, "gem_group")         # integer
feature_idx  <- rhdf5::h5read(mol_src, "feature_idx")       # 0-based
feat_ids     <- rhdf5::h5read(mol_src, "features/id")       # character vector

# Construct full-length cell_id & gene_id like your function does
cell_id_full <- paste(mol_barcodes[barcode_idx + 1L], gem_group, sep = "-")
gene_id_full <- feat_ids[feature_idx + 1L]

# The tiny ffm barcodes should match the TSV barcodes_small (which include "-1")
qc_cells_target <- unique(barcodes_small)          # EXPECT strings like "AAAC...-1"
keep_gene_ids   <- unique(features_small[[1]])     # first column is gene_id

# Keep only molecules where both conditions hold
keep_rows <- which(cell_id_full %in% qc_cells_target & gene_id_full %in% keep_gene_ids)

# Build tiny vectors; reindex barcode_idx/feature_idx into compact 0-based space
# (a) compact *barcode* space based on barcodes_small (strip "-1" to map to mol_barcodes)
bare_from_tsv <- sub("-[0-9]+$", "", barcodes_small)
bar_keep0     <- match(bare_from_tsv, mol_barcodes) - 1L     # 0-based indices in full space
stopifnot(all(!is.na(bar_keep0)))

# (b) compact *feature* space based on features_small[[1]]
feat_keep0    <- match(keep_gene_ids, feat_ids) - 1L         # 0-based indices in full space
stopifnot(all(!is.na(feat_keep0)))

# Maps from full-space indices to compact indices (0-based)
bar_compact_map  <- setNames(seq_along(bar_keep0) - 1L, bar_keep0)   # names are full 0-based
feat_compact_map <- setNames(seq_along(feat_keep0) - 1L, feat_keep0)

# Compute compact indices for kept rows
bar_idx_full_for_rows  <- barcode_idx[keep_rows]
feat_idx_full_for_rows <- feature_idx[keep_rows]
bar_idx_compact  <- as.integer(bar_compact_map[as.character(bar_idx_full_for_rows)])
feat_idx_compact <- as.integer(feat_compact_map[as.character(feat_idx_full_for_rows)])
stopifnot(all(!is.na(bar_idx_compact)), all(!is.na(feat_idx_compact)))

# Minimal datasets for tiny molecule_info.h5
count_s       <- as.integer(c(2,1,1,2,rep(1,length(keep_rows)-4)))
umi_s         <- as.integer(umi[keep_rows])
barcode_idx_s <- bar_idx_compact
feature_idx_s <- feat_idx_compact
gem_group_s   <- as.integer(gsub("^.*-", "", cell_id_full[keep_rows]))  # derive from kept cell_ids
barcodes_s    <- mol_barcodes[bar_keep0 + 1L]            # bare barcodes in compact order
features_id_s <- feat_ids[feat_keep0 + 1L]               # gene IDs in compact order

# Helper: create compressed dataset & write
create_and_write <- function(file, name, data, chunk_len = 1024L, level = 9L) {
  if (is.matrix(data)) {
    dims  <- dim(data)
    chunk <- pmax(1L, pmin(dims, c(3L, min(1024L, dims[2]))))
    rhdf5::h5createDataset(file, name, dims = dims,
                           storage.mode = if (is.integer(data)) "integer" else "character",
                           chunk = chunk, level = level)
  } else {
    len   <- length(data)
    chunk <- min(chunk_len, max(1L, len))
    rhdf5::h5createDataset(file, name, dims = len,
                           storage.mode = if (is.integer(data)) "integer" else "character",
                           chunk = chunk, level = level)
  }
  rhdf5::h5write(data, file, name)
}

# Write tiny molecule_info.h5
if (file.exists(mol_dst)) file.remove(mol_dst)
rhdf5::h5createFile(mol_dst)
create_and_write(mol_dst, "count",       count_s)
create_and_write(mol_dst, "umi",         umi_s)
create_and_write(mol_dst, "barcode_idx", barcode_idx_s)   # compact 0-based
create_and_write(mol_dst, "feature_idx", feature_idx_s)   # compact 0-based
create_and_write(mol_dst, "gem_group",   gem_group_s)
create_and_write(mol_dst, "barcodes",    barcodes_s)      # bare barcodes in compact order
rhdf5::h5createGroup(mol_dst, "features")
create_and_write(mol_dst, "features/id", features_id_s)
rhdf5::H5close()

# ============================================================
# (3) filtered_feature_bc_matrix.h5 — only "matrix/barcodes"
#     Must match the TSV barcodes used for qc filtering in your function.
# ============================================================
ffm_dst <- file.path(DEST_DIR, "outs", "filtered_feature_bc_matrix.h5")
if (file.exists(ffm_dst)) file.remove(ffm_dst)
rhdf5::h5createFile(ffm_dst)
rhdf5::h5createGroup(ffm_dst, "matrix")
# Use exactly the TSV barcodes (already include the "-gem_group" suffix)
create_and_write(ffm_dst, "matrix/barcodes", unique(barcodes_small))
rhdf5::H5close()

# ============================================================
# (4) metrics_summary.csv — clone template columns; fill "Number of Reads"
#     Use the sum of counts in the *kept molecule rows* (post-filter).
# ============================================================
number_of_reads <- sum(as.numeric(count_s), na.rm = TRUE)

if (!file.exists(TEMPLATE_CSV)) {
  stop("Template CSV not found: ", TEMPLATE_CSV)
}
tpl <- read.csv(TEMPLATE_CSV, check.names = FALSE)
row_list <- as.list(rep(NA, length(tpl))); names(row_list) <- names(tpl)
df_out <- as.data.frame(row_list, check.names = FALSE, stringsAsFactors = FALSE)

col_nr <- "Number of Reads"
if (!(col_nr %in% names(df_out))) {
  stop("Column 'Number of Reads' not found in template CSV.")
}
df_out[[col_nr]] <- number_of_reads

out_csv <- file.path(DEST_DIR, "outs", "metrics_summary.csv")
write.csv(df_out, out_csv, row.names = FALSE, quote = TRUE)

cat("Tiny dataset ready at:\n", DEST_DIR, "\n",
    "- outs/filtered_feature_bc_matrix/{matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz}\n",
    "- outs/molecule_info.h5 (rows filtered by cell_id & gene_id)\n",
    "- outs/filtered_feature_bc_matrix.h5 (matrix/barcodes)\n",
    "- outs/metrics_summary.csv (Number of Reads = ", format(number_of_reads, scientific = FALSE), ")\n", sep = "")
