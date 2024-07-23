# MIT License
#
# Copyright (c) 2023 Eric BAREKE (eric.bareke@mcgill.ca), Emma Carlson (emma.carlson@mail.mcgill.ca), and Majewski Lab (McGill University)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# How to run it ?
# Rscript DGEAfromRes.r -f /path/to/your/input/files -t 2 -p 1

# List of required packages
libraries <- c(
  "clusterProfiler",
  "enrichplot",
  "ggplot2",
  "DESeq2",
  "EnhancedVolcano",
  "grDevices",
  "RColorBrewer",
  "pheatmap",
  "ggrepel",
  "dplyr",
  "argparse",
  "here",
  "org.Hs.eg.db",
  "pathfindR"
)

# Function to check, install, and load packages
check_install_load <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
#    # Install the package if not installed
#    if (package %in% BiocManager::installed()) {
#      BiocManager::install(package)
#    } else {
#      install.packages(package)
#    }
  }

  # Load the package
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

# Apply the function to each package in the list
lapply(libraries, check_install_load)


# Function to parse command line arguments
parse_args <- function() {
  parser <- ArgumentParser()
  parser$add_argument("-f", "--filePath", type = "character", help = "Path to the input files")
  parser$add_argument("-t", "--fc_cutoff", type = "numeric", help = "Fold change cutoff", default = 1.5)
  parser$add_argument("-p", "--adjusted_p_value_cutoff", type = "numeric", help = "Adjusted p-value cutoff", default = 0.05)

  args <- parser$parse_args()
  return(args)
}

# Parse command line arguments
opt <- parse_args()


############# Function Definitions #############

# Function for saving plots in multiple formats
save_plot <- function(plot, filename) {
  png_filename <- paste0(filename, ".png")
#  pdf_filename <- paste0(filename, ".pdf")

  # Save as PNG
  ggsave(png_filename, plot = plot, path = outputPath, dpi = 600, device = 'png')

  # Save as PDF
#  ggsave(pdf_filename, plot = plot, path = outputPath, dpi = 600, device = 'pdf')
}

# Perform inner joins
merge_and_rename <- function(data, subset_data, suffix) {
  merged_data <- merge(counts(dds, normalized = FALSE), subset_data, by = "row.names", all = FALSE)
  colnames(merged_data)[1] <- "GeneSymbol"
  colnames(merged_data)[-1] <- paste0(colnames(merged_data)[-1], suffix)
  return(merged_data)
}

# Write results to Excel-compatible files with raw counts and Gene Symbol
write_results <- function(data, filename) {
  write.csv(data, file = here::here(outputPath, filename), row.names = FALSE, quote = FALSE)
}

# Function to perform GO enrichment analysis
perform_GO_enrichment <- function(gene_list, direction, ont_type, output_path, adjusted_p_value_cutoff = 0.1) {
  GO_results <- enrichGO(
    gene = gene_list,
    universe = background,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ont_type,
    pAdjustMethod = "fdr",
    minGSSize = 10,
    maxGSSize = 2000,
    pvalueCutoff = adjusted_p_value_cutoff,
    qvalueCutoff = adjusted_p_value_cutoff,
    readable = TRUE
  )

  # Plot and save charts...
  save_enrichment_plots(GO_results, direction, "GO", ont_type, output_path)
}

perform_KEGG_enrichment <- function(gene_list, direction, output_path, adjusted_p_value_cutoff = 0.1) {
  KEGG_results <- enrichKEGG(
    gene = gene_list,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = adjusted_p_value_cutoff,
    pAdjustMethod = "BH",
    universe = background,
    minGSSize = 10,
    maxGSSize = 2000,
    qvalueCutoff = adjusted_p_value_cutoff,
    use_internal_data = FALSE
  )

  # Plot and save charts...
  save_enrichment_plots(KEGG_results, direction, "KEGG", NULL, output_path)
}

# Save enrichment plots
save_enrichment_plots <- function(results, direction, analysis_type, ont_type, output_path) {
  bar_plot <- barplot(results, showCategory = 15)
  dot_plot <- dotplot(results, showCategory = 15)
  # Save plots and data in different formats with unique names
  finame <- paste0("leadingEdge_", analysis_type, ifelse(!is.null(direction), paste0("_", direction), ""), ifelse(!is.null(ont_type), paste0("_", ont_type), ""), "_results.csv")
  write_results(results,finame)
  save_plot(bar_plot, paste0("barPlot_", analysis_type, ifelse(!is.null(direction), paste0("_", direction), ""), ifelse(!is.null(ont_type), paste0("_", ont_type), "")))
  save_plot(dot_plot, paste0("dotPlot_", analysis_type, ifelse(!is.null(direction), paste0("_", direction), ""), ifelse(!is.null(ont_type), paste0("_", ont_type), "")))

}

# Assign values to variables
filePath <- opt$filePath
FC_cutoff <- as.numeric(opt$fc_cutoff)
adjusted_p_value_cutoff <- as.numeric(opt$adjusted_p_value_cutoff)

# Set input and output paths
inputPath <- file.path(filePath, "inputFiles")
outputPath <- file.path(filePath, "outputFiles")

# Create inputFiles folder if it doesn't exist
if (!dir.exists(inputPath)) {
  dir.create(inputPath, recursive = TRUE)
} else {
  cat("Input directory already exists. Skipping creation.\n")
}

# Create outputFiles folder if it doesn't exist
if (!dir.exists(outputPath)) {
  dir.create(outputPath, recursive = TRUE)
} else {
  cat("Output directory already exists. Skipping creation.\n")
}


# Copy count matrix and sample information files to inputFiles folder
lessBindingFiles <- list.files(pattern = "*less_input.t.*")
moreBindingFiles <- list.files(pattern = "*more_input.t.*")

if (length(lessBindingFiles) > 0) {
  # Specify the destination folder for count matrix files
  less_binding_matrix <- file.path(inputPath, lessBindingFiles)

  # Copy count matrix files with overwrite
  file.copy(lessBindingFiles, less_binding_matrix, overwrite = TRUE)
}

if (length(moreBindingFiles) > 0) {
  # Specify the destination folder for sample information files
  more_binding_matrix <- file.path(inputPath, moreBindingFiles)

  # Copy sample information files with overwrite
  file.copy(moreBindingFiles, more_binding_matrix, overwrite = TRUE)
}

lessBindingFile <- list.files(path = inputPath, pattern = ".*less_input.t.*")
moreBindingFile <- list.files(path = inputPath, pattern = ".*more_input.t.*")

if (length(lessBindingFile) > 0) {
  # Specify the destination folder for count matrix files
  less_binding_matrix <- file.path(inputPath, lessBindingFiles)

  # Copy count matrix files with overwrite
  file.copy(lessBindingFile, less_binding_matrix, overwrite = TRUE)
}

if (length(moreBindingFile) > 0) {
  # Specify the destination folder for sample information files
  more_binding_matrix <- file.path(inputPath, moreBindingFiles)

  # Copy sample information files with overwrite
  file.copy(moreBindingFile, less_binding_matrix, overwrite = TRUE)
}


lessBinding <- read.table(file.path(inputPath, lessBindingFiles), header = FALSE)
moreBinding <- read.table(file.path(inputPath, moreBindingFiles), header = FALSE)
backgroundL <- read.table("/Users/ebareke/Desktop/Katie_ATACSeq-CUTnTAG_V4_Feb-02-2024/Quantification/Background.txt", header = FALSE)

colnames(lessBinding) <- c("GeneSymbol")
colnames(moreBinding) <- c("GeneSymbol")

colnames(backgroundL) <- c("GeneSymbol")

# Define input parameters as needed...

universe <- backgroundL$GeneSymbol
background <- mapIds(org.Hs.eg.db, keys = universe, column = "ENTREZID", keytype = "SYMBOL")

background <- background[!is.na(background)]

listDN <- lessBinding$GeneSymbol
listDN <- mapIds(org.Hs.eg.db, keys = listDN, column = "ENTREZID", keytype = "SYMBOL")
listDN <- listDN[!is.na(names(listDN))]

# Main analysis for downregulated genes
perform_GO_enrichment(listDN, "Less", "BP", outputPath, adjusted_p_value_cutoff)
perform_GO_enrichment(listDN, "Less", "MF", outputPath, adjusted_p_value_cutoff)
perform_GO_enrichment(listDN, "Less", "CC", outputPath, adjusted_p_value_cutoff)
perform_KEGG_enrichment(listDN, "Less", outputPath, adjusted_p_value_cutoff)

# Define other input parameters as needed

listUP <- moreBinding$GeneSymbol
listUP <- mapIds(org.Hs.eg.db, keys = listUP, column = "ENTREZID", keytype = "SYMBOL")
listUP <- listUP[!is.na(names(listUP))]

# Main analysis for upregulated genes
perform_GO_enrichment(listUP, "More", "BP", outputPath, adjusted_p_value_cutoff)
perform_GO_enrichment(listUP, "More", "MF", outputPath, adjusted_p_value_cutoff)
perform_GO_enrichment(listUP, "More", "CC", outputPath, adjusted_p_value_cutoff)
perform_KEGG_enrichment(listUP, "More", outputPath, adjusted_p_value_cutoff)
