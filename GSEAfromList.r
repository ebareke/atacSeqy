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
# Rscript GSEAfromList.r -f /path/to/your/input/files -p 1

# List of required packages
libraries <- c(
  "DESeq2",
  "org.Hs.eg.db",
  "tidyverse",
  "fgsea",
  "ggplot2",
  "ggrepel",
  "reshape2",
  "argparse",
  "here",
  "data.table",
  "grDevices"
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
  parser$add_argument("-p", "--adjusted_p_value_cutoff", type = "numeric", help = "Adjusted p-value cutoff", default = 0.1)
  args <- parser$parse_args()
  return(args)
}

# Parse command line arguments
opt <- parse_args()


# Function Definitions

# Function for saving plots in multiple formats
save_plot <- function(plot, filename) {
  png_filename <- paste0(filename, ".png")
  #pdf_filename <- paste0(filename, ".pdf")

  # Save as PNG
  ggsave(png_filename, plot = plot, path = outputPath, width = 24, height = 16, dpi = 600, device = 'png')

  # Save as PDF
  #ggsave(pdf_filename, plot = plot, path = outputPath,  width = 16, height = 12, dpi = 600, device = 'pdf')
}

# Write results to Excel-compatible files with raw counts and Gene Symbol
write_results <- function(data, filename) {
  write.table(data, file = here::here(outputPath, filename), sep = "\t", row.names = FALSE, quote = FALSE)
}

save_results <- function(data, filename) {
  fwrite(data, file = here::here(outputPath, filename), sep = "\t")
}

#### Perform GSEA
GSEA = function(gene_list, GMT_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(fgsea)

  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGMT = fgsea::gmtPathways(GMT_file)

  fgRes <- fgsea::fgseaMultilevel(pathways = myGMT,
                           stats = gene_list,
                           minSize=20, ## minimum gene set size
                           maxSize=2000) %>%
                  as.data.frame() %>%
                  dplyr::filter(padj < !!pval) %>%
                  arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))

  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGMT,
                                      stats = gene_list)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))

  fgRes$Enrichment = ifelse(fgRes$NES > 0, "More Binding", "Less Binding")
  filtRes = rbind(head(fgRes, n = 20),
                  tail(fgRes, n = 20 ))

  total_up = sum(fgRes$Enrichment == "More Binding")
  total_down = sum(fgRes$Enrichment == "Less Binding")
  header = paste0("Top 40 Significant Enrichments (FDR = 1) - [Total: More Bound = ", total_up, ", Less Bound = ", total_down, "]")

  colos = setNames(c("firebrick2", "dodgerblue2"),
                 c("More Binding", "Less Binding"))

 g1= ggplot(filtRes, aes(reorder(pathway, NES), NES, size = 1.2)) +
  geom_point( aes(fill = Enrichment, size = size), shape=21) +
  scale_fill_manual(values = colos ) +
  scale_size_continuous(range = c(2,10)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x="Enrichments", y="Normalized Enrichment Score",
       title=header) +
        theme_bw() +
    theme(axis.text.x=element_text(size=rel(2.5))) +
    theme(axis.text.y=element_text(size=rel(2.5))) +
    theme(axis.title.x=element_text(size=rel(2.5))) +
    theme(axis.title.y=element_text(size=rel(2.5))) + 
    theme(legend.position="none", plot.title = element_text(hjust=0.5,size = 16, face = "bold"), legend.title=element_text(size=14, face = "bold"))



  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}

# Assign values to variables
filePath <- opt$filePath
adjusted_p_value_cutoff <- as.numeric(opt$adjusted_p_value_cutoff)

# Set input and output paths
inputPath <- file.path(filePath, "inputFiles")
outputPath <- file.path(filePath, "gseaFiles")

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


# Copy DESeq result file to inputFiles folder
DESeq2ResFiles <- list.files(pattern = "*sea_input.t.*")

if (length(DESeq2ResFiles) > 0) {
  # Specify the destination folder for count matrix files
  deseq2_result_matrix <- file.path(inputPath, DESeq2ResFiles)

  # Copy count matrix files with overwrite
  file.copy(DESeq2ResFiles, deseq2_result_matrix, overwrite = TRUE)
}

DESeq2ResFile <- list.files(path = inputPath, pattern = ".*sea_input.t.*")

if (length(DESeq2ResFile) > 0) {
  # Specify the destination folder for count matrix files
  deseq2_result_matrix <- file.path(inputPath, DESeq2ResFiles)

  # Copy count matrix files with overwrite
  file.copy(DESeq2ResFile, deseq2_result_matrix, overwrite = TRUE)
}


DESeq2Res <- read.table(file.path(inputPath, DESeq2ResFiles), header = FALSE)

colnames(DESeq2Res) <- c("GeneSymbol","log2FoldChange")

gene_list <- DESeq2Res$log2FoldChange
names(gene_list) = DESeq2Res$GeneSymbol

gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

# List of GMT files
gmt_files <- c(
  "c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt",
  "c3.tft.v2023.2.Hs.symbols.gmt",
  "c5.go.bp.v2023.2.Hs.symbols.gmt",
  "c5.go.mf.v2023.2.Hs.symbols.gmt",
  "c6.all.v2023.2.Hs.symbols.gmt",
  "c7.all.v2023.2.Hs.symbols.gmt",
  "h.all.v2023.2.Hs.symbols.gmt"
)

# Loop over GMT files
for (gmt_file in gmt_files) {
  # Construct the pathway object
  pathways <- paste("/Users/ebareke/Desktop/Katie_ATACSeq-CUTnTAG_V4_Feb-02-2024/Quantification/signatures/", gmt_file, sep = "")

  # Perform fgsea analysis
  fgsea_res = GSEA(gene_list, pathways, pval = adjusted_p_value_cutoff)

  # Extract GSEA result
  res <- fgsea_res$Results

  # Save the results in EXCEL-compatible format
   save_results(res, paste(gsub("\\.gmt$", "", gmt_file), "_gsea_results.txt", sep = ""))

  # Create a plot
  plot <- fgsea_res$Plot

  # Save the plot in both formats
  save_plot(plot, paste(gsub("\\.gmt$", "", gmt_file), "_gseaEnrichmentPlot", sep = ""))
}

# Delete the Rplots.pdf file
# file.remove("Rplots.pdf")
