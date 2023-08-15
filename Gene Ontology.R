setwd("/home/sejyoti/Downloads/Gastic cancer cell Line/hjc 1/HFE145_GO_G1-HFE145_G2M")

# Load required libraries
library(ggplot2)
library(viridis)

# Read the data
enrichment_results <- read.csv("enrichment(KEGG).csv")
fdr_palette <- viridisLite::viridis(100, option = "D")

# Create the plot
plot <- ggplot(enrichment_results, aes(x = -log10(FDR), y = reorder(Pathway, -log10(FDR)), size = nGenes, fill = FDR)) +
  geom_point(pch = 21, alpha = 0.7) +
  scale_fill_viridis(trans = "log", option = "D", direction = -1) +
  scale_size_continuous(range = c(1, 7)) +
  labs(title = "KEGG Enrichment Analysis",
       x = "-log10(FDR)", y = "Pathway", size = "nGenes", fill = "FDR") +
  theme_bw()  # Use theme_bw() for a white background

# Save the plot as a PDF
ggsave("enrichment_bubble_plot.pdf", plot, width = 10, height = 20, device = "pdf")

