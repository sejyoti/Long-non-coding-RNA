
setwd("/home/sejyoti/Downloads/Gastic cancer cell Line")
# Load required libraries
library(DESeq2)

# Read data from CSV file
file_path <- "/home/sejyoti/Downloads/Gastic cancer cell Line/Thirdset.xlsx"
data <- read.csv(file_path, header = TRUE)

# Remove duplicate rows based on Gene ID
unique_data <- data[!duplicated(data$Gene.ID), ]

# Continue with data preparation
count_matrix <- as.matrix(unique_data[, -1])  # Remove the Gene.ID column
rownames(count_matrix) <- unique_data$Gene.ID

# Create a sample table with conditions
sample_table <- data.frame(
  condition = c("AGS_rep1", "AGS_rep2", "KATO-III_rep1", "KATO-III_rep2"),
  replicate = c("rep1", "rep2", "rep1", "rep2")
)
rownames(sample_table) <- colnames(count_matrix)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_table, design = ~condition)

# Perform normalization (try different methods)
dds <- DESeq(dds, betaPrior = TRUE)

# Conduct differential expression analysis
results <- results(dds)

# Remove rows with NA values in padj or log2FoldChange
non_na_results <- results[!is.na(results$padj) & !is.na(results$log2FoldChange), ]

# ... Previous code ...

# Adjust filtering criteria to capture more genes
significant_results <- non_na_results[non_na_results$padj < 0.5 & abs(non_na_results$log2FoldChange) > 0.02, ]

# Separate upregulated and downregulated genes
upregulated_results <- significant_results[significant_results$log2FoldChange > 0, ]
downregulated_results <- significant_results[significant_results$log2FoldChange < 0, ]

# ... Rest of the code ...

# Save results to CSV files
write.csv(upregulated_results, file = upregulated_file_path, row.names = TRUE)
write.csv(downregulated_results, file = downregulated_file_path, row.names = TRUE)

# Load the 'fgsea' package if not already loaded
if (!requireNamespace("fgsea", quietly = TRUE))
  install.packages("fgsea", dependencies = TRUE)

library(fgsea)

# Create a vector of gene symbols for fgsea
gene_symbols <- rownames(upregulated_results)

# Example gene sets (replace with your own gene sets)
gene_sets <- list(
  upregulated_genes = upregulated_results$GeneID,
  downregulated_genes = downregulated_results$GeneID
)

# Create a named vector of log2FoldChange values using gene IDs as names
stats_named <- setNames(upregulated_results$log2FoldChange, upregulated_results$GeneID)

# Run GSEA
gsea_results <- fgsea(pathways = gene_sets, stats = stats_named, nPerm = 1000)
# Check for non-matching gene IDs between stats_named and gene_sets
missing_gene_ids <- setdiff(names(stats_named), unique(unlist(gene_sets)))
if (length(missing_gene_ids) > 0) {
  cat("Gene IDs missing in gene_sets:", paste(missing_gene_ids, collapse = ", "), "\n")
}




# Print the GSEA results
print(gsea_results)

