

# ===============================================================
# File: matrix_to_seurat_rds.R
# File Created: Th Sep 2023
# Author: Zhang Bin
# -----
# Last Modified: Fri Sep 08 2023
# Modified By: Zhang Bin
# ===============================================================

rm(list=ls())
gc()

library(Seurat)
library(tidyverse)

args <- commandArgs(T)
if (length(args) < 4) {
    stop("USAGE: Rscript matrix_to_seurat_rds.R sample_name count_matrix_path.tsv barcode_path.tsv out_file.rds [meta_data_path.csv]\n")
}
sample <- args[1]
count_matrix_path <- args[2]
barcode_path <- args[3]
out <- args[4]
meta_data_path <- if (length(args) >= 5) args[5] else NA

usage <- str_c(
    "USAGE:\n",
    "\tRscript matrix_to_seurat_rds.R sample_name count_matrix_path.tsv barcode_path.tsv out_file.rds [meta_data_path.csv]",
    "\n",
    "\tmeta_data_path: meta.data is formatted as CSV, which should have a cell name column named 'CN' that is the same as the cell name in the barcode file."
)

# check
if (!file.exists(count_matrix_path))
    stop(str_c("ERROR: cannot find count matrix: ", count_matrix_path, "\n\n", usage))
if (!file.exists(barcode_path))
    stop(str_c("ERROR: cannot find barcode: ", barcode_path, "\n\n", usage))
if (!is.na(meta_data_path) && !file.exists(meta_data_path))
    stop(str_c("ERROR: cannot find meta data: ", meta_data_path, "\n\n", usage))

message(str_c(
    "NOTE: transform count matrix to Seurat obj!\n",
    "\tcount matrix:\t", count_matrix_path,"\n",
    "\tbarcode:\t", barcode_path,"\n",
    "\tmeta.data:\t", ifelse(is.na(meta_data_path), "", meta_data_path), "\n"
))


# CN: CellName, CB: CellBarcode, CI:CellId
barcode <- read_tsv(barcode_path, col_names = c("CN", "CB")) %>%
    mutate(
        CI = str_extract(CN, "\\d+$"),
        CI = str_c(sample, "_", CB, "-", CI)
    )
CB2CI <- setNames(barcode$CI, barcode$CB)

# meta.data
if (!is.na(meta_data_path)) {
    meta.data <- read_csv(meta_data_path) 
    if(colnames(meta.data)[1] %in% c("barcode")) {
        colnames(meta.data)[1] <- "CB"
    }
    meta.data <- barcode %>%
        left_join(meta.data) %>%
        column_to_rownames("CI")
} else {
    meta.data <- barcode %>%
        column_to_rownames("CI")
}

# matrix
count_matrix <- read_tsv(count_matrix_path) %>%
    column_to_rownames("GeneID") %>%
    as.matrix()
colnames(count_matrix) <- unname(CB2CI[colnames(count_matrix)])
count_matrix[is.na(count_matrix)] <- 0

# create seurat obj
rds <- CreateSeuratObject(
        counts = count_matrix,
        project = sample,
        meta.data = meta.data
    )
rds
rds %>%
    saveRDS(out)

message(str_c(
    "NOTE:\n",
    "\toutput:\t", out
))

