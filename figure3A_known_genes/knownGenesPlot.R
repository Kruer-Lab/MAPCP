library(tidyverse)
library(RColorBrewer)

denovo <- read.csv("MAPCP_Variants_Table_03.01.2025_FINAL_denovo.tsv", sep='\t')
recessive <- read.csv("MAPCP_Variants_Table_03.01.2025_FINAL_recessive.tsv", sep='\t')

# Only keep known top candidates
denovo <- denovo[denovo$Top.Candidate=="Yes",]
recessive <- recessive[recessive$Top.Candidate=="Yes",]
recessive <- recessive[grepl("Known", recessive$Gene.Category),]

# Sanity check for duplicate samples
denovo$ID %in% recessive$ID
recessive$ID %in% denovo$ID

# Merge tables
denovo <- data.frame(sapply(denovo, as.character))
merged <- bind_rows(denovo, recessive)
merged$PredominantMotorDisorder <- str_trim(merged$PredominantMotorDisorder)

# Bar plot of genes by inheritance mode
ggplot(merged, aes(x=Type)) +
  geom_bar() +
  ggtitle("Known disease genes by inheritance mode") +
  theme_minimal()

# Stacked bar plot by inheritance mode and predominant motor disorder
ggplot(merged, aes(fill=PredominantMotorDisorder, x=Type)) +
  geom_bar(position="stack", stat="count") +
  theme_minimal()

# Clean up
merged$Type <- factor(merged$Type, levels = c("Xlink", "CompHet", "De novo", "Hom"))
merged$PredominantMotorDisorder <- factor(merged$PredominantMotorDisorder, 
                                          levels = c("Not available",
                                                     "Dyskinetic-chorea, athetosis, and/or ballism",
                                                     "Ataxic", "Dyskinetic-Dystonia",
                                                     "Hypotonic", "Spastic"))

# Colorblind friendly palette
color_palette <- c("#D55E00", "#0072B2", "#F0E442", "#009E73", "#56B4E9", "#E69F00")

ggplot(merged, aes(fill=PredominantMotorDisorder, x=Type)) +
  geom_bar(position="stack", stat="count") +
  theme_minimal() +
  ggtitle("Known disease genes by mode of inheritance and predominant movement disorder of proband") +
  xlab("Inheritance mode") +
  ylab("Number of cases") +
  labs(fill="Predominant movement disorder") +
  scale_fill_manual(
    values=color_palette,
    guide = guide_legend(reverse=TRUE)
  ) +
  scale_x_discrete(labels =c("Xlink"="X-Linked", 
                             "CompHet"= "Compound\nHeterozygous",
                             "De novo" = "De novo",
                             "Hom" = "Homozygous\nRecessive")) +
  scale_y_continuous(
    limits=c(0, 90),
    breaks=c(seq(0, 90, by=15))
  )

# Save as png 800x600