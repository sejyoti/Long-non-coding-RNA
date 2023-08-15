if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")


library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(AnnotationDbi)


# Load upregulated and downregulated gene lists
upregulated <- read.csv("/home/sejyoti/Downloads/Gastic cancer cell Line/upregulated_results.csv")
downregulated <- read.csv("/home/sejyoti/Downloads/Gastic cancer cell Line/downregulated_results.csv")

# Load your background dataset and reshape it to long format
background_data <- read.csv("/home/sejyoti/Downloads/Gastic cancer cell Line/GSE134150_all.counts.csv")
background_long <- pivot_longer(background_data, cols = -GeneID, names_to = "Sample", values_to = "Expression")

# Trim whitespace from gene identifiers in background_long
background_long$GeneID <- trimws(background_long$GeneID)

# Perform GSEA for upregulated genes
enrich_result_up <- enricher(gene = upregulated$GeneID,
                             TERM2GENE = background_long,
                             universe = unique(background_long$GeneID),
                             pvalueCutoff = 0.05,
                             pAdjustMethod = "BH")

# Perform GSEA for downregulated genes
enrich_result_down <- enricher(gene = downregulated$GeneID,
                               TERM2GENE = background_long,
                               universe = unique(background_long$GeneID),
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH")

# Extract lncRNA, protein-coding RNA, and mRNA gene sets from upregulated enrichment results
lncRNA_up <- enrich_result_up$geneSets[grepl("lncRNA", enrich_result_up$Description)]
protein_coding_up <- enrich_result_up$geneSets[grepl("protein_coding", enrich_result_up$Description)]
mRNA_up <- enrich_result_up$geneSets[grepl("mRNA", enrich_result_up$Description)]

# Extract lncRNA, protein-coding RNA, and mRNA gene sets from downregulated enrichment results
lncRNA_down <- enrich_result_down$geneSets[grepl("lncRNA", enrich_result_down$Description)]
protein_coding_down <- enrich_result_down$geneSets[grepl("protein_coding", enrich_result_down$Description)]
mRNA_down <- enrich_result_down$geneSets[grepl("mRNA", enrich_result_down$Description)]

# Print names of the extracted gene sets
print("Upregulated Gene Sets:")
print(paste("lncRNA:", names(lncRNA_up)))
print(paste("Protein-Coding:", names(protein_coding_up)))
print(paste("mRNA:", names(mRNA_up)))

print("Downregulated Gene Sets:")
print(paste("lncRNA:", names(lncRNA_down)))
print(paste("Protein-Coding:", names(protein_coding_down)))
print(paste("mRNA:", names(mRNA_down)))

