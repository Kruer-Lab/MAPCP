#source("/gpfs/ysm/project/kahle/sp592/Walker2019/Walker_Enrichment_Plots.R", echo = T)

library(tidyverse)
library(magrittr)
library(useful)
library(org.Hs.eg.db)
library(HGNChelper)
library(ggpubr)
library(cowplot)
#library(miSciTools)
library(egg)
library(ggpubr)
options(stringsAsFactors=FALSE)

setwd(".")

rm(list=ls())

load("MAP_CP_Walker2019_Results.RData")

# Function for module plot
p.value.plots.module <- function(gene.set, disease.label) {
  # Verify column names in dis.res
  if (!"p" %in% colnames(dis.res)) {
    if ("p_value" %in% colnames(dis.res)) {
      dis.res <- dis.res %>% rename(p = p_value)
    } else if ("P" %in% colnames(dis.res)) {
      dis.res <- dis.res %>% rename(p = P)
    } else {
      stop("Column 'p' not found in dis.res. Available columns: ", paste(colnames(dis.res), collapse = ", "))
    }
  }
  if (!"p.transform" %in% colnames(dis.res)) {
    if ("p_transform" %in% colnames(dis.res)) {
      dis.res <- dis.res %>% rename(p.transform = p_transform)
    } else {
      stop("Column 'p.transform' not found in dis.res. Available columns: ", paste(colnames(dis.res), collapse = ", "))
    }
  }
  
  # Get unique modules with p.transform >= 2, or all modules if none found
  unique_module <- if (length(unique(dis.res %>% filter(p.transform >= 2) %>% distinct(module) %>% pull(module))) == 0) {
    unique(dis.res$module)
  } else {
    unique(dis.res %>% filter(p.transform >= 2) %>% distinct(module) %>% pull(module))
  }
  # print("Using unique_module:")
  # print(unique_module)
  
  # Subset the data
  dis.res.sub <- dis.res %>%
    filter(module %in% unique_module) %>%
    filter(disease %in% names(gene.set)) %>%
    mutate(disease = factor(disease, levels = names(gene.set)))
  
  # Validate disease matching
  if (!any(dis.res$disease %in% names(gene.set))) {
    stop("No diseases in dis.res$disease match names(gene.set). Available diseases: ", 
         paste(unique(dis.res$disease), collapse = ", "), 
         ". Expected: ", paste(names(gene.set), collapse = ", "))
  }
  
  # Debug: Check columns and rows
  # print("Columns in dis.res.sub:")
  # print(colnames(dis.res.sub))
  # print("Number of rows in dis.res.sub:")
  # print(nrow(dis.res.sub))
  # print("First few rows of dis.res.sub:")
  # print(head(dis.res.sub))
  
  # Precompute the label, size, fontface, and color
  dis.res.sub <- dis.res.sub %>%
    mutate(
      # Label: Show significant values with *, non-significant without
      label = ifelse(
        .data$p < 0.05 / length(unique(dis.res$module)),
        paste(round(.data$p.transform, digits = 1), "*", sep = ""),
        round(.data$p.transform, digits = 1)
      ),
      # Text size: 3 for significant, scaled for non-significant
      text_size = ifelse(
        .data$p < 0.05 / length(unique(dis.res$module)),
        3,
        scales::rescale(.data$p.transform, to = c(1, 2.5))
      ),
      # Text fontface: bold for significant, italic for non-significant
      text_fontface = ifelse(
        .data$p < 0.05 / length(unique(dis.res$module)),
        "bold",
        "italic"
      ),
      # Text color: Gradient from grey to black based on p.transform
      text_color = ifelse(
        .data$p < 0.05 / length(unique(dis.res$module)),
        "black",
        colorRampPalette(c("grey70", "black"))(100)[
          as.integer(scales::rescale(.data$p.transform, to = c(1, 100)))
        ])
    )
  
  # Create the plot
  ggplot(
    dis.res.sub, 
    aes(x = module, y = disease)
  ) +
    geom_tile(aes(fill = p.transform), color = "white") +
    scale_fill_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 1), low = "aliceblue", high = "darkseagreen4") +
    scale_y_discrete(labels = gsub("  ", "\n", gene.set)) +
    theme_classic() + 
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1, colour = "black"),
      axis.text.y = element_text(angle = 0, size = 8, colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(angle = 90, size = 10, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
      legend.position = "top", 
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    geom_text(
      aes(
        label = label,
        color = text_color  # Map the precomputed color
      ),
      size = dis.res.sub$text_size,  # Use precomputed size (outside aes)
      fontface = dis.res.sub$text_fontface  # Use precomputed fontface (outside aes)
    ) +
    scale_colour_identity() +  # Use the precomputed colors directly
    ylab(disease.label) +
    NULL 
}

## Function for cell type plot
p.value.plots.celltype <- function(gene.set, disease.label, x.axis.labels) {
  # Subset the data
  hydro.cell.type <- hydro.cell.type %>% filter(disease %in% names(gene.set))
  hydro.cell.type$disease <- hydro.cell.type$disease %>% factor(levels = names(gene.set))
  
  # Define cell type order and labels
  order.cells <- c("vRG", "oRG", "PgS", "PgG2M", "IP", "ExN", "ExM", "ExM-U", "ExDp1", "ExDp2", "InMGE", "InCGE", "OPC", "End", "Per", "Mic")
  renamed.cells <- c("Ventricular Radial Glia (vRG)",
                     "Outer Radial Glia (oRG)",
                     "Cycling Progenitor: S phase (PgS)",
                     "Cycling Progenitor: G2/M phase (PgG2M)",
                     "Intermediate Progenitor (IP)",
                     "Excitatory Neuron: Migrating (ExN)",
                     "Excitatory Neuron: Maturing (ExM)",
                     "Excitatory Neuron: Maturing - Upper Enriched Layer (ExM-U)",
                     "Excitatory Neuron: Deep Layer 1 (ExDp1)",
                     "Excitatory Neuron: Deep Layer 2 (ExDp2)",
                     "Interneuron: MGE (InMGE)",
                     "Interneuron: CGE (InCGE)",
                     "Oligodendrocyte Precursor (OPC)",
                     "Endothelial (End)",
                     "Pericyte (Per)",
                     "Microglia (Mic)"
  )
  
  hydro.cell.type$cell.type <- hydro.cell.type$cell.type %>% factor(levels = order.cells)
  
  # Debug: Check columns and rows
  # print("Columns in hydro.cell.type:")
  # print(colnames(hydro.cell.type))
  # print("Number of rows in hydro.cell.type:")
  # print(nrow(hydro.cell.type))
  # print("First few rows of hydro.cell.type:")
  # print(head(hydro.cell.type))
  
  # Precompute the label, size, fontface, and color
  hydro.cell.type <- hydro.cell.type %>%
    mutate(
      # Label: Show significant values with *, non-significant without
      label = ifelse(
        .data$p.ind < 0.05 / length(unique(.data$cell.type)),
        paste(round(.data$p.transform, digits = 1), "*", sep = ""),
        round(.data$p.transform, digits = 1)
      ),
      # Text size: 3 for significant, scaled for non-significant
      text_size = ifelse(
        .data$p.ind < 0.05 / length(unique(.data$cell.type)),
        3,
        scales::rescale(.data$p.transform, to = c(1, 2.5))
      ),
      # Text fontface: bold for significant, italic for non-significant
      text_fontface = ifelse(
        .data$p.ind < 0.05 / length(unique(.data$cell.type)),
        "bold",
        "italic"
      ),
      # Text color: Gradient from grey to black based on p.transform
      text_color = ifelse(
        .data$p.ind < 0.05 / length(unique(.data$cell.type)),
        "black",
        colorRampPalette(c("grey70", "black"))(100)[
        as.integer(scales::rescale(.data$p.transform, to = c(1, 100)))
      ])
    )
  
  # Create the plot
  ggplot(
    hydro.cell.type, 
    aes(x = cell.type, y = disease)
  ) +
    geom_tile(aes(fill = p.transform), color = "white") +
    scale_fill_continuous(limits=c(0, 20), breaks=seq(0, 20, by=5), low = "aliceblue", high = "rosybrown4") +
    scale_x_discrete(labels = renamed.cells) +
    scale_y_discrete(labels = gsub("  ", "\n", gene.set)) +
    theme_classic() + 
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1, colour = "black"), 
      axis.text.y = element_text(angle = 0, size = 8, colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(angle = 90, size = 10, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
      legend.position = "top", 
      legend.title = element_text(size = 12, face = "bold"), 
      legend.text = element_text(size = 10)
    ) +
    geom_text(
      aes(
        label = label,
        color = text_color  # Map the precomputed color
      ),
      size = hydro.cell.type$text_size,  # Use precomputed size (outside aes)
      fontface = hydro.cell.type$text_fontface  # Use precomputed fontface (outside aes)
    ) +
    scale_colour_identity() +  # Use the precomputed colors directly
    ylab(disease.label) +
    NULL
}

## Plotting
# Custom display names for each panel
label_map <- c(
  "PanelApp_Baylor_Invitae_Intersection_Risk" = "Clinical CP Panels",
  "PMID32989326_Risk" = "Published CP Risk Genes",
  "MAP_CP_Known_or_Novel" = "MAP CP (Clinician-curated)",
  "MAP_CP_Risk" = "MAP CP (Constraint-based)"
)

# Prepare CONTROL gene sets
control_gene_sets <- setNames(
  paste0(label_map[control_panels]),
  # paste0(label_map[control_panels], " (n=", lengths(gene_list[control_panels]), ")"),
  control_panels
)

# Prepare CASE gene sets
case_gene_sets <- setNames(
  paste0(label_map[case_panels]),
  # paste0(label_map[case_panels], " (n=", lengths(gene_list[case_panels]), ")"),
  case_panels
)

###############
# Module Plot #
###############
## Disease Control Panel
control_module_plots <- list()  # Initialize list to store ggplot objects

for (i in seq_along(control_gene_sets)) {
  panel_name <- names(control_gene_sets)[i]
  gene_set <- control_gene_sets[i]
  
  # Generate the plot (this returns a ggplot object)
  # Generate the plot
  plot <- gene_set %>%
    p.value.plots.module(disease.label = "") +     # "Known Diseases\nGene"
    scale_x_discrete(labels = NULL)
  
  # Add title and colorbar only for the first panel
  if (i == 1) {
    plot <- plot +
      ggtitle("Module Enrichment") +
      guides(fill = guide_colourbar(title = expression('-log'[10]*'p'))) +
      NULL
  } else {
    plot <- plot +
      guides(fill = F) +
      NULL
  }
  
  # ✅ Correct: store the ggplot object directly
  control_module_plots[[i]] <- plot  
}

## Case Panel
case_module_plot <-
  case_gene_sets %>% p.value.plots.module(disease.label = "") +           # "Cerebral\nPalsy"
  guides(fill=F) +
  NULL

## Optionally delete the existing file manually
if (file.exists("Module_Enrichment_Plots.pdf")) {
  file.remove("Module_Enrichment_Plots.pdf")
}

## Save as pdf
pdf(
  file="Module_Enrichment_Plots.pdf",
  height = 5,
  width = 5,
  onefile = F
)
# Combine control plots and case plot dynamically
plots_to_arrange <- c(control_module_plots, list(case_module_plot))
# do.call() handles unpacking the list into separate arguments properly.
do.call(egg::ggarrange, c(plots_to_arrange, list(
  nrow = length(control_module_plots) + 1,
  heights = c(rep(1, length(control_module_plots)), 3)
)))
dev.off()

#################
# Cellytpe Plot #
#################
## Disease Control Panel
control_celltype_plots <- list()  # Initialize list to store ggplot objects

for (i in seq_along(control_gene_sets)) {
  panel_name <- names(control_gene_sets)[i]
  gene_set <- control_gene_sets[i]
  
  # Generate the plot (this returns a ggplot object)
  # Generate the plot
  plot <- gene_set %>%
    p.value.plots.celltype(disease.label = "") +          # "Known Diseases\nGene"
    scale_x_discrete(labels = NULL)
  
  # Add title and colorbar only for the first panel
  if (i == 1) {
    plot <- plot +
      ggtitle("Module Enrichment") +
      guides(fill = guide_colourbar(title = expression('-log'[10]*'p'))) +
      NULL
  } else {
    plot <- plot +
      guides(fill = F) +
      NULL
  }
  
  # ✅ Correct: store the ggplot object directly
  control_celltype_plots[[i]] <- plot  
}

## Case Panel
case_celltype_plot <-
  case_gene_sets %>% p.value.plots.celltype(disease.label = "") +           # "Cerebral\nPalsy"
  guides(fill=F) +
  NULL

## Optionally delete the existing file manually
if (file.exists("Celltype_Enrichment_Plots.pdf")) {
  file.remove("Celltype_Enrichment_Plots.pdf")
}

## Save as pdf
pdf(
  file="Celltype_Enrichment_Plots.pdf",
  height = 8,
  width = 10,
  onefile = F
)
# Combine control plots and case plot dynamically
plots_to_arrange <- c(control_celltype_plots, list(case_celltype_plot))
# do.call() handles unpacking the list into separate arguments properly.
do.call(egg::ggarrange, c(plots_to_arrange, list(
  nrow = length(control_celltype_plots) + 1,
  heights = c(rep(1, length(control_celltype_plots)), 3)
)))
dev.off()

