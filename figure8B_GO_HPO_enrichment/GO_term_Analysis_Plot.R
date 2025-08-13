# Step 1: Install and load required packages
library(viridis)
library(ggplot2)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(GOSemSim)
library(org.Hs.eg.db)
library(HGNChelper)
library(magrittr)

###############
# Set Options #
###############
# Install first if needed: install.packages("optparse")
# library(optparse)

# option_list <- list(
#   make_option(c("-i", "--input"), type = "character", help = "Input target gene list", metavar = "FILE"),
#   make_option(c("-b", "--background"), type = "character", help = "Background gene list", metavar = "FILE")
# )

# opt <- parse_args(OptionParser(option_list = option_list))

# input_file <- opt$input
# background_file <- opt$background

# cat("Target gene file:", input_file, "\n")
# cat("Background gene file:", background_file, "\n")

#############
# Functions #
#############
prepare_genes_for_enrichGO <- function(gene_symbols, OrgDb) {
  # Remove duplicates before validation
  gene_symbols_unique <- unique(gene_symbols)
  cat("Original input genes:", length(gene_symbols), "\n")
  cat("Unique genes (after deduplication):", length(gene_symbols_unique), "\n")
  cat("Removed duplicates:", length(gene_symbols) - length(gene_symbols_unique), "\n\n")
  
  # Step 1: Validate symbols
  validated <- HGNChelper::checkGeneSymbols(gene_symbols_unique)
  
  # Step 2: Keep only approved or corrected symbols
  corrected <- validated$Suggested.Symbol[validated$Approved]
  cat("Approved HGNC symbols:", length(corrected), "\n\n")
  
  # Step 3: Convert to ENSEMBL (or ENTREZID)
  gene_df <- clusterProfiler::bitr(corrected,
                                   fromType = "SYMBOL",
                                   toType = "ENSEMBL",  # Change to "ENTREZID" if preferred
                                   OrgDb = OrgDb)
  
  # Step 4: Deduplicate SYMBOLs to 1:1 mapping
  gene_df_unique <- gene_df %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  cat("Final unique mappings (SYMBOL to ENSEMBL):", nrow(gene_df_unique), "\n")
  
  return(gene_df_unique)
}

##############
# Read files #
##############
# Step 2: Read the module and background gene lists
module_genes <- read.csv("Module_Genes_List/midnightblue_gene_symbols_filtered.csv", 
                         stringsAsFactors = FALSE)$Gene_Symbol
module_genes_df <- prepare_genes_for_enrichGO(module_genes, OrgDb = org.Hs.eg.db)


background_genes <- read.csv("Module_Genes_List/grey_gene_symbols_filtered.csv", 
                             stringsAsFactors = FALSE)$Gene_Symbol
background_genes_df <- prepare_genes_for_enrichGO(background_genes, OrgDb = org.Hs.eg.db)

universe_fixed <- union(module_genes_df$ENSEMBL, background_genes_df$ENSEMBL) # enrichGO assume the gene were the target out from the universe list

####################
# GO Term Analysis #
####################
# Step 3: Run GO enrichment analysis with enrichGO
enrich_result <- enrichGO(
  gene = module_genes_df$ENSEMBL,           # Foreground genes (module genes)
  universe = universe_fixed,   # Background genes
  OrgDb = org.Hs.eg.db,          # Organism database (human)
  ont = "BP",                    # Ontology: Biological Process (BP), Molecular Function (MF), Cellular Component (CC)
  pvalueCutoff = 0.05,           # P-value cutoff
  qvalueCutoff = 0.05,           # FDR q-value cutoff
  pAdjustMethod = "BH",          # FDR adjustment method (Benjamini-Hochberg)
  keyType = "ENSEMBL",            # Gene identifier type (assuming gene symbols)
  readable = TRUE                # Map gene IDs to symbols for readability
)

# Debug: Check the enrichResult object
print("enrichResult object:")
print(head(enrich_result@result))
print("Number of enriched terms:")
print(nrow(enrich_result@result))

# If no enriched terms are found, stop and investigate
if (nrow(enrich_result@result) == 0) {
  stop("No enriched GO terms found. Check your gene lists or adjust the pvalueCutoff/qvalueCutoff.")
}

# Save result to csv
# Convert enrich_result to a data frame
enrich_result_df <- as.data.frame(enrich_result)

# Convert GeneRatio and BgRatio from strings (e.g., "5/300") to numeric and then get enrichment
enrich_result_df <- enrich_result_df %>%
  dplyr::mutate(
    gene_ratio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    bg_ratio = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    enrichment = gene_ratio / bg_ratio
  )

# Save to CSV
write.csv(enrich_result_df, file = "enrichGO_with_enrichment.csv", row.names = FALSE)

#########################
# Simplify the GO Terms #
#########################
# Step 4: Simplify GO terms using clusterProfiler's simplify function
# Prepare semantic similarity data
sem_data <- tryCatch(
  {
    GOSemSim::godata("org.Hs.eg.db", ont = "BP")
  },
  error = function(e) {
    message("Failed to create semData with GOSemSim::godata: ", e$message)
    message("Proceeding with simplify without semData (will attempt to fetch GO data automatically).")
    return(NULL)
  }
)

# Debug: Check if semData was created successfully
if (!is.null(sem_data)) {
  print("semData created successfully.")
  go_terms_in_data <- enrich_result@result$ID
  terms_found <- go_terms_in_data %in% names(sem_data@IC)
  print("GO terms found in semData:")
  print(table(terms_found))
} else {
  print("semData is NULL. Simplify will attempt to fetch GO data automatically.")
}

# Run simplify with a low cutoff
set.seed(123)
simplified_go <- tryCatch(
  {
    simplify(
      enrich_result,
      cutoff = 0.3,
      by = "pvalue",
      select_fun = min,
      measure = "Wang",
      semData = sem_data
    )
  },
  error = function(e) {
    message("Simplify failed: ", e$message)
    message("Proceeding without simplification.")
    return(enrich_result)
  }
)

# Debug: Check the simplified results
print("Simplified GO terms:")
print(head(simplified_go@result))
print("Number of simplified terms:")
print(nrow(simplified_go@result))

# Save result to csv
# Convert simplified result to a data frame
simplified_df <- as.data.frame(simplified_go)

# Add gene_ratio, bg_ratio, and enrichment columns
simplified_df <- simplified_df %>%
  dplyr::mutate(
    gene_ratio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    bg_ratio = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    enrichment = gene_ratio / bg_ratio
  )

# Save to CSV
write.csv(simplified_df, file = "simplified_GO_enrichment.csv", row.names = FALSE)

###################
# Plot-Simplified #
###################
# Step 5: Prepare data for plotting
# Add -log10(FDR q-value) for coloring
simplified_data <- simplified_df %>%
  mutate(neg_log10_fdr = -log10(qvalue))

# Debug: Check the simplified data
print("Simplified data:")
print(head(simplified_data))
print("Number of rows in simplified_data:")
print(nrow(simplified_data))

# If simplified_data is empty, use the original results
if (nrow(simplified_data) == 0) {
  message("No terms remain after simplification. Using original terms.")
  simplified_data <- as.data.frame(enrich_result@result) %>%
    mutate(neg_log10_fdr = -log10(qvalue))
}

# Step 6: Select the top 20 terms by FDR q-value
simplified_data_top <- simplified_data %>%
  arrange(qvalue) %>%
  slice_head(n = 20)

# Debug: Check the top 20 terms
print("Top 20 terms:")
print(head(simplified_data_top))
print("Number of rows in simplified_data_top:")
print(nrow(simplified_data_top))

# If simplified_data_top is empty, stop and investigate
if (nrow(simplified_data_top) == 0) {
  stop("No terms remain after selecting top 20. Check your data.")
}

# Shorten long descriptions for better readability
# Shorten long descriptions for better readability
simplified_data_top <- simplified_data_top %>%
  mutate(Description = str_trunc(Description, width = 75, side = "right"))

# Step 7: Create the dot plot
p <- ggplot(simplified_data_top, aes(x = enrichment, y = reorder(Description, neg_log10_fdr))) +
  geom_point(aes(size = Count, color = neg_log10_fdr)) +
  scale_color_viridis_c(name = "-log10(FDR q-value)") +
  scale_size_continuous(range = c(3, 8), name = "Gene Count") +
  labs(
    x = "Fold Enrichment (Count / Gene Ratio)",
    y = "GO Term",
    title = "GO Term Enrichment Analysis (Simplified)"
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

# Save the plot as a PDF
ggsave("GO_Enrichment_Dotplot_Simplified.pdf", plot = p, width = 10, height = 8)


#############
# Plot-Full #
#############
# Add -log10(qvalue) for coloring
enrich_result_df <- enrich_result_df %>%
  mutate(neg_log10_fdr = -log10(qvalue))

# Select top 20 by qvalue
enrich_result_df_top <- enrich_result_df %>%
  arrange(qvalue) %>%
  slice_head(n = 20) %>%
  mutate(Description = str_trunc(Description, width = 75, side = "right"))

# Create the dot plot
p_full <- ggplot(enrich_result_df_top, aes(x = enrichment, y = reorder(Description, neg_log10_fdr))) +
  geom_point(aes(size = Count, color = neg_log10_fdr)) +
  scale_color_viridis_c(name = "-log10(FDR q-value)") +
  scale_size_continuous(range = c(3, 8), name = "Gene Count") +
  labs(
    x = "Fold Enrichment (GeneRatio / BgRatio)",
    y = "GO Term",
    title = "GO Term Enrichment Analysis"
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

# Save the plot
ggsave("GO_Enrichment_Dotplot_Full.pdf", plot = p_full, width = 10, height = 8)

