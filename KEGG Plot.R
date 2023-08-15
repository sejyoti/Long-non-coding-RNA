
# Load the data from the CSV file
data <- read.csv('/home/sejyoti/Downloads/Gastic cancer cell Line/hjc 1/AGS-KATOIII/enrichment(kegg).csv')

# Load required packages
library(ggplot2)

# Load the data from the CSV file
data <- read.csv('/home/sejyoti/Downloads/Gastic cancer cell Line/hjc 1/AGS-KATOIII/enrichment(kegg).csv')

# Create the enrichment plot
enrichment_plot <- ggplot(data, aes(x = -log10(FDR), y = reorder(`Pathway`, -log10(FDR)))) +
  geom_point(aes(size = `nGenes`, color = `group`), shape = 16, alpha = 0.7) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white")) +
  labs(
    x = "-log₁₀(FDR)",
    y = "Pathway",
    color = "Group",
    title = "KEGG Pathway Enrichment",
    subtitle = "Upregulated and Downregulated Pathways",
    caption = "Data source: Your Source"
  ) +
  geom_text(aes(label = paste("nGenes:", `nGenes`)), hjust = -0.2, size = 3, color = "black") +
  geom_text(aes(label = paste("FoldEnriched:", round(`Fold enriched`, 2))), hjust = -0.2, vjust = -0.2, size = 3, color = "black") +
  annotate("text", x = 0, y = -1, label = "Red: Upregulated\nBlue: Downregulated", hjust = 0, vjust = 0, size = 3, color = "black")

# Create a new graphics device
dev.new()
# Save the plot as an image file
ggsave("enrichment_plot_csv.png", enrichment_plot, width = 10, height = 6, dpi = 300)

