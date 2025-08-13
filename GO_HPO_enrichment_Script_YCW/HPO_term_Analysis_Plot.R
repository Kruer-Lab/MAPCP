# Step 1: Install and load required packages
if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
if (!requireNamespace("HGNChelper", quietly = TRUE)) install.packages("HGNChelper")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")

library(viridis)
library(ggplot2)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(HGNChelper)
library(org.Hs.eg.db)
library(magrittr)

# Step 2: Prepare helper function
prepare_genes_for_enrich <- function(gene_symbols) {
  gene_symbols_unique <- unique(gene_symbols)
  cat("Original input genes:", length(gene_symbols), "\n")
  cat("Unique genes (after deduplication):", length(gene_symbols_unique), "\n")
  
  validated <- HGNChelper::checkGeneSymbols(gene_symbols_unique)
  corrected <- validated$Suggested.Symbol[validated$Approved]
  cat("Approved HGNC symbols:", length(corrected), "\n\n")
  
  return(corrected)
}

# Step 3: Read gene lists
module_genes <- read.csv("Module_Genes_List/midnightblue_gene_symbols_filtered.csv", 
                         stringsAsFactors = FALSE)$Gene_Symbol
module_genes_clean <- prepare_genes_for_enrich(module_genes)

background_genes <- read.csv("Module_Genes_List/grey_gene_symbols_filtered.csv", 
                             stringsAsFactors = FALSE)$Gene_Symbol
background_genes_clean <- prepare_genes_for_enrich(background_genes)

universe_fixed <- union(module_genes_clean, background_genes_clean)

# Step 4: Load HPO gene-term mappings
hpo_data <- read.delim("genes_to_phenotype.txt", sep = "\t", header = TRUE, quote = "", fill = TRUE, stringsAsFactors = FALSE)
colnames(hpo_data) <- c("ncbi_gene_id", "gene_symbol", "hpo_id", "hpo_name", "frequency", "disease_id")

TERM2GENE <- dplyr::select(hpo_data, hpo_id, gene_symbol) %>% distinct()
TERM2NAME <- dplyr::select(hpo_data, hpo_id, hpo_name) %>% distinct()

# Optional: Checking step-How many genes in HPO TERM gene
cat("Numbers of Overlapped Genes in HPO Term Gene Sets:", sum(module_genes_clean %in% TERM2GENE$gene_symbol), "\n\n")

# Step 5: Run enrichment analysis
enrich_result <- enricher(
  
  gene = module_genes_clean,
  universe = universe_fixed,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# Step 6: Prepare and save result
enrich_df <- as.data.frame(enrich_result)
enrich_df$GeneRatioNum <- sapply(strsplit(enrich_df$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
enrich_df$BgRatioNum <- sapply(strsplit(enrich_df$BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
enrich_df$enrichment <- as.numeric(enrich_df$GeneRatioNum) / as.numeric(enrich_df$BgRatioNum)
enrich_df$neg_log10_fdr <- -log10(as.numeric(enrich_df$qvalue))

write.csv(enrich_df, file = "HPO_enrichment_results.csv", row.names = FALSE)

# Step 7: Plot top 20 terms
top_terms <- enrich_df %>% arrange(qvalue) %>% slice_head(n = 20)
top_terms <- top_terms %>% mutate(Description = str_trunc(Description, width = 75, side = "right"))

p <- ggplot(top_terms, aes(x = enrichment, y = reorder(Description, neg_log10_fdr))) +
  geom_point(aes(size = Count, color = neg_log10_fdr)) +
  scale_color_viridis_c(name = "-log10(FDR q-value)") +
  scale_size_continuous(range = c(3, 8), name = "Gene Count") +
  labs(
    x = "Fold Enrichment (Count / Gene Ratio)",
    y = "HPO Term",
    title = "HPO Term Enrichment Analysis"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  )

ggsave("HPO_Enrichment_Dotplot.pdf", plot = p, width = 10, height = 8)
