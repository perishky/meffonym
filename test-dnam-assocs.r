#!/usr/bin/env R

if (!require("data.table", quietly = TRUE)) install.packages("data.table", quietly = TRUE)
if (!require("remotes", quietly = TRUE)) install.packages("remotes", quietly = TRUE)
if (!require("ewaff", quietly = TRUE)) remotes::install_github("perishky/ewaff", quiet = TRUE)
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
data.dir <- args[1]
files <- list(
  Methylation_matrix = file.path(data.dir, "Methylation_matrix.csv"),
  Illumina_matrix = file.path(data.dir, "Illumina_matrix.csv"),
  dnam_score = file.path(data.dir, "DNA_Methylation_Scores.csv"),
  camda_score = file.path(data.dir, "camda_matrix.csv"),
  phenotype = file.path(data.dir,"phenotype.csv"),
  ECC = file.path(data.dir,"estimate_cell_counts.csv"),
  reads = file.path(data.dir,"cutadapt_filtered_reads_plot.txt"),
  removed_samples = file.path(data.dir,"removed_samples.txt"))

pheno <- data.frame(fread(files$phenotype))

expected_columns <- c("Sex", "Country.x", "Age_Rerc", "Cancer_status", "BMI_C", "Smoke_status", "Alc_Lifetime", "batch")
cat("Expected columns are: ", paste(expected_columns, collapse=", "), "\n")
missing <- setdiff(expected_columns, colnames(pheno))
extra <- setdiff(colnames(pheno), expected_columns)
message({
  if (length(missing)) message("Warning: Missing columns in phenotype file: ", paste(missing, collapse=", "))
  if (length(extra)) message("Note: Extra columns in phenotype file: ", paste(extra, collapse=", "))
  if (!length(missing) && !length(extra)) message("All expected columns are present in the phenotype file. Proceed please.")
})

reads<-fread(file.path(files$reads))
reads$Sample <- paste0("X", gsub("-", ".", reads$Sample))
colnames(reads)[colnames(reads) == "Reads passing filters"] <- "reads"
reads<-reads %>% mutate(Sample = gsub("_R[12]", "", Sample)) %>% distinct(Sample, .keep_all = TRUE) 
reads <- reads %>% mutate(Sample = gsub(".*_(.+)", "\\1", Sample)) 

pheno <- merge(pheno,reads,by="Sample")
rownames(pheno) <- pheno$Sample
pheno$Sample <- NULL

ECC <- fread(file.path(files$ECC))
ECC <- data.frame(ECC)
cell_types <- ECC$V1
ECC$V1 <- NULL
ECC <- t(ECC)
colnames(ECC) <- cell_types
rownames(ECC) <- gsub(".*_(S\\d+)", "\\1", rownames(ECC))
pheno <- cbind(pheno, ECC[pheno$Sample,])

removed_samples <- readr::read_lines(files$removed_samples)
pheno <- pheno[!pheno$Sample %in% removed_samples,]

models <- read.csv(text = "
var,model
Alc_Lifetime,methylation~Alc_Lifetime+Smoke_status+BMI_C+Age_Rerc+Sex+reads
Sex,methylation~Sex+reads
Country.x,methylation~Country.x
Age_Rerc,methylation~Age_Rerc+reads
Cancer_status,methylation~Cancer_status+Age_Rerc+reads
")
models$model <- paste(models$model, "+CD4T+CD8T+NK+Mono+Bcell+Granulocytes")

process_ewaff <- function(pheno, meth, output_folder, summary_folder = NULL, manifest = NULL, models) {
  common_samples <- intersect(rownames(pheno), colnames(meth))
  stopifnot(length(common_samples) > 1)
  pheno <- pheno[common_samples,]
  meth <- meth[,common_samples]
  results <- list()
  for (i in seq_len(nrow(models))) {
    formula <- as.formula(models$model[i])
    var <- models$var[i]
    sites.ret <- ewaff.sites(
      formula,
      variable.of.interest = var,
      methylation = as.matrix(meth),
      data = pheno,
      method = "glm")
    if (!is.null(manifest)) {
      sum.ret <- ewaff.summary(sites.ret, manifest$chr, manifest$start, meth_matrix)
      ewaff.report(sum.ret, output.file = paste0(output_folder, "/report_", var, ".html"))
      results[[var]]$sum.ret <- sum.ret
    } 
    table_file <- if (!is.null(summary_folder)) {
      dir.create(summary_folder, showWarnings = FALSE, recursive = TRUE)
      paste0(summary_folder, "/summary_statistics_", var, ".csv")
    } else {
      paste0(output_folder, "/sites_ret_", var, ".csv")
    } 
    fwrite(cbind("rownames" = rownames(sites.ret$table), sites.ret$table), file = table_file, row.names = FALSE)
    results[[var]]$sites.ret <- sites.ret
  } 
  return(results)
}

prepare_data <- function(file, id_col=NULL, manifest_cols = NULL) {
  data <- data.frame(fread(file))
  if (is.null(manifest_cols)) 
    manifest <- data.frame(chr="score", start=1:nrow(data), end=1:nrow(data))
  else 
    manifest <- data[manifest_cols]
  if (is.null(id_col))
    rownames(data) <- paste(manifest$chr, manifest$start, sep=":")
  else
    rownames(data) <- data[[id_col]]
  return(list(data = data, manifest = manifest))}

datasets <- list(
  Methylation_matrix = list(file = files$Methylation_matrix, manifest_cols = c("chr", "start", "end")),
  Illumina_matrix = list(file = files$Illumina_matrix, id_col = "CpGs", manifest_cols = c("chr", "start", "end")),
  dnam_score = list(file = files$dnam_score, id_col = "V1"),
  camda_score = list(file = files$camda_score, manifest_cols = c("chr", "start", "end")))

results <- lapply(names(datasets), function(name) {
  dataset <- datasets[[name]]
  output_folder <- paste0("output/", name)
  dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
  summary_folder <- paste0(output_folder, "/summary_statistics_", name)
  data_prepared <- prepare_data(dataset$file, dataset$id_col, dataset$manifest_cols)
  result <- process_ewaff(
    pheno=pheno,
    meth = data_prepared$data,
    output_folder = output_folder,
    summary_folder = summary_folder,
    manifest = data_prepared$manifest,
    models = models)
  save_file <- paste0(output_folder, "/results_", name, ".rda")
  save(list = c("result"), file = save_file) 
  return(result)
})
names(results) <- names(datasets)
