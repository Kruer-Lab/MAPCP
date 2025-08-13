library(readxl)
library(tidyverse)

denovo <- read.csv("MAPCP_Variants_Table_03.01.2025_FINAL_denovo.tsv", sep='\t')
recessive <- read.csv("MAPCP_Variants_Table_03.01.2025_FINAL_recessive.tsv", sep='\t')

denovo <- denovo[!(is.na(denovo$ID) | denovo$ID == ""),]
recessive <- recessive[!(is.na(recessive$ID) | recessive$ID == ""),]

genemap2 <- read.csv("genemap2_01272025.txt", sep='\t', skip=3)
genemap2 <- genemap2[-(18626:18701),]

# Recessive variants
mergedRecessive <- merge(recessive, genemap2, by.x="Gene", by.y="Approved.Gene.Symbol", all.x=TRUE)
write.table(mergedRecessive, "mergedVariantsPhenotypesRecessive.tsv", sep='\t', row.names=FALSE)

# De novo variants
mergedDeNovo <- merge(denovo, genemap2, by.x="Gene", by.y="Approved.Gene.Symbol", all.x=TRUE)
write.table(mergedDeNovo, "mergedVariantsPhenotypesDeNovo.tsv", sep='\t', row.names=FALSE)


#### After manual review for motor phenotype ####
motorRecessive <- read_excel("multilocus.xlsx", sheet="recessive",
                                   skip=0, col_names=TRUE)

# Keep homozygous recessive and motor phenotypes
motorRecessive <- motorRecessive[motorRecessive$Type=="Hom" & motorRecessive$`Motor Component`=="Y",]

# Remove duplicate genes within a family
motorRecessive <- motorRecessive[!duplicated(motorRecessive[1:2]),]

# Get number of times multiple genes have appeared
occurrenceBreakdown <- as.data.frame(table(table(motorRecessive$ID)))

ggplot(occurrenceBreakdown, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  ggtitle("Occurrence of recessive movement genes in families") +
  xlab("Number of recessive movement genes in families") +
  ylab("Number of families")
# Save as 500x500

# Analyze by ROH similarity

# Keep only samples with >1 genes
numOccurrences <- data.frame(table(motorRecessive$ID))
numOccurrences[numOccurrences$Freq > 1,]
dropSingletons <- motorRecessive[motorRecessive$ID %in%
                                   numOccurrences$Var1[numOccurrences$Freq > 1],]

write.table(dropSingletons, "dropSingletons.tsv", sep='\t', row.names=FALSE)

# After checking for overlapping ROHs
rohOverlap <- duplicated(dropSingletons[c(2,19)])
sameROH <- dropSingletons[dropSingletons$ID %in% dropSingletons[rohOverlap, 2]$ID &
                            dropSingletons$ROH.Size.of.Variant.in.Proband..Mb. %in%
                            dropSingletons[rohOverlap, 19]$ROH.Size.of.Variant.in.Proband..Mb.,]

occurrenceBreakdown$Overlap <- c(0,4,2,0,2)
occurrenceBreakdown$NonOverlap <- occurrenceBreakdown$Freq - occurrenceBreakdown$Overlap

pivot <- occurrenceBreakdown %>% select(Var1, Overlap, NonOverlap) %>%
  pivot_longer(!Var1)

# Stacked bar
pivot$name <- factor(pivot$name, levels = c("Overlap", "NonOverlap"))

ggplot(pivot, aes(x=Var1, y=value, fill=name)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  ggtitle("Number of recessive movement disorder genes by family") +
  xlab("Number of genes") +
  ylab("Number of families") +
  labs(fill="ROH Region") +
  scale_fill_manual(labels=c("Overlap" = "2 variants in same ROH region",
                             "NonOverlap" = "Different ROH regions"),
    values = c("#E69F00", "#56B4E9"),
    guide = guide_legend(reverse=FALSE)
  )

# Save as 500x500