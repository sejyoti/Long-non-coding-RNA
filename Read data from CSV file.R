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
  condition = c("HFE145", "HFE145", "AGS", "AGS"),
  replicate = c("rep1", "rep2", "rep1", "rep2")
)
rownames(sample_table) <- colnames(count_matrix)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_table, design = ~condition)

# Perform normalization
dds <- DESeq(dds)

# Conduct differential expression analysis
results <- results(dds)

# Filter significant differentially expressed genes
significant_results <- results[results$padj < 0.5 & abs(results$log2FoldChange) > 1, ]


# Check for missing values in the results data frame
missing_values <- apply(is.na(results), 2, any)
missing_columns <- colnames(results)[missing_values]

if (length(missing_columns) > 0) {
  print(paste("Missing values found in columns:", paste(missing_columns, collapse = ", ")))
  # Handle missing values if necessary
} else {
  # Proceed with filtering
  significant_results <- results[results$padj < 0.05 & abs(results$log2FoldChange) > 1, ]
}

# Remove rows with missing values
results_clean <- results[complete.cases(results), ]

# Proceed with filtering
significant_results <- results_clean[results_clean$padj < 0.05 & abs(results_clean$log2FoldChange) > 1, ]

# Separate upregulated and downregulated genes
upregulated_results <- significant_results[significant_results$log2FoldChange > 0, ]
downregulated_results <- significant_results[significant_results$log2FoldChange < 0, ]


# Provide a complete file path including the file name and extension
upregulated_file_path <- "/home/sejyoti/Downloads/Gastic cancer cell Line/upregulated_results.csv"
downregulated_file_path <- "/home/sejyoti/Downloads/Gastic cancer cell Line/downregulated_results.csv"

# Save results to CSV files
write.csv(upregulated_results, file = upregulated_file_path, row.names = TRUE)
write.csv(downregulated_results, file = downregulated_file_path, row.names = TRUE)




# install.packages("ggplot2")
library(ggplot2)
# Convert DESeqResults to data frame
volcano_data <- as.data.frame(significant_results)

# Create a new column for color based on fold change
volcano_data$color <- ifelse(volcano_data$log2FoldChange > 1, "Upregulated",
                             ifelse(volcano_data$log2FoldChange < -1, "Downregulated", "Not Significant"))

# Create the volcano plot using the new data frame
volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = c("red", "blue", "black")) +
  labs(x = "log2 Fold Change", y = "-log10(Adjusted p-value)", title = "Volcano Plot")

print(volcano_plot)



install.packages("pheatmap")

# Extract expression data for significant genes
significant_gene_ids <- rownames(significant_results)
significant_expression <- count_matrix[significant_gene_ids, ]

# Load pheatmap library
library(pheatmap)


# Create a heatmap of expression data
heatmap_expression <- pheatmap(significant_expression, 
                               color = colorRampPalette(c("blue", "white", "red"))(50),
                               scale = "row",  # Scale rows to have similar variance
                               main = "Heatmap of Differentially Expressed Genes",
                               fontsize_row = 8,  # Adjust row label font size
                               fontsize_col = 8   # Adjust column label font size
)


# Display the heatmap
heatmap_expression



# Assuming you have the significant_results data frame and count_matrix available

# Load necessary libraries
library(ggplot2)

# Create a new data frame containing expression levels and conditions for significant genes
significant_expression_data <- data.frame(
  Gene_ID = rownames(significant_results),
  log2FoldChange = significant_results$log2FoldChange,
  pvalue = significant_results$pvalue,
  condition = sample_table$condition[rownames(significant_results)]
)

# Create a scatter plot
scatter_plot <- ggplot(significant_expression_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = condition), size = 3, alpha = 0.7) +
  scale_color_manual(values = c("HFE145" = "blue", "AGS" = "red")) +
  labs(x = "log2 Fold Change",
       y = "-log10(p-value)",
       title = "Scatter Plot of Differentially Expressed Genes",
       color = "Condition") +
  theme_minimal()

# Display the scatter plot
print(scatter_plot)

genes to test 
