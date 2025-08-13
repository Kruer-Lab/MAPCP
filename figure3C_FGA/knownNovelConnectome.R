# Functional Genomic Alignment (FGA) code from 
# https://lab.rockefeller.edu/casanova/assets/file/FGA_plot.R.txt
library(ape)
library(biomaRt)
library(HGNChelper)

m <- as.matrix(read.table("distance_matrix.txt", head=T, row.names=1))
arbol <- nj(as.dist(m))
arbol$edge.length <- arbol$edge.length + abs(min(arbol$edge.length))+1

# Load in geneDx panel, known genes, novel candidates
geneList <- read.csv("gene_list.txt", header=FALSE)
geneList$Color = c(rep("Black", 1184), rep("Red", 107), rep("Green", 27))
geneList$TipLabel = geneList$V1
geneList$TipLabel[0:1184] <- ""
geneList$Size = c(rep(0, 1184), rep(1, 107), rep(1, 27))

# Match gene names with those in our paper
geneList[geneList$TipLabel=="SPG20", "TipLabel"] <- "SPART"
geneList[geneList$TipLabel=="ADCK3", "TipLabel"] <- "COQ8A"
geneList[geneList$TipLabel=="EPRS1", "TipLabel"] <- "EPRS"
geneList[geneList$TipLabel=="CPSF3L", "TipLabel"] <- "INTS11"
geneList[geneList$TipLabel=="RARS", "TipLabel"] <- "RARS1"
geneList[geneList$TipLabel=="H3F3A", "TipLabel"] <- "H3-3A"
geneList[geneList$TipLabel=="ASNA1", "TipLabel"] <- "GET3"

# Remove duplicate genes keeping color for called genes over geneDx genes
deduplicated <- geneList[!duplicated(geneList$V1, fromLast=TRUE),]

# Reorder colors to match tree
orderedColors <- deduplicated[match(arbol$tip.label, deduplicated$V1),]
write.table(orderedColors$V1, "orderedColors.tsv", sep='\t')

plot(arbol,edge.width=1, tip.color=orderedColors$Color, 
     type="fan", show.tip.label=FALSE)

arbol$tip.label <- orderedColors$TipLabel

pdf('FGA_tree.pdf',width=10, height=10,pointsize=0.1)
plot(arbol,edge.width=1, tip.color=orderedColors$Color, 
     type="fan", show.tip.label=TRUE)

dev.off()

# Add GO terms to tree

# Check if annotations are still good
checkedSymbols <- checkGeneSymbols(orderedColors$V1)
fixedSymbols <- merge(orderedColors, checkedSymbols, by.x="V1", by.y="x", sort=FALSE)
fixedSymbols[fixedSymbols$Suggested.Symbol=="B3GNT2 /// B4GAT1", "Suggested.Symbol"] <- "B4GAT1"
fixedSymbols[fixedSymbols$Suggested.Symbol=="EPRS1 /// QARS1", "Suggested.Symbol"] <- "QARS1"

# Get annotations
ensembl <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart = "genes")
annotation <- getBM(filters = 'hgnc_symbol', 
                    attributes=c('hgnc_symbol', 'entrezgene_id'),
                    values = fixedSymbols$Suggested.Symbol,
                    mart = ensembl)
annotatedData <- merge(fixedSymbols , annotation, by.x="Suggested.Symbol", by.y="hgnc_symbol", all.x=TRUE)
write.table(annotatedData, "annotatedData.tsv", sep='\t', row.names = FALSE)

# Add node labels
testTree <- arbol
testTree$node.label <- 1:testTree$Nnode
pdf('FGA_treeNewNodes.pdf',width=10, height=10,pointsize=0.1)
plot(testTree,edge.width=1, tip.color=orderedColors$Color, 
     type="fan",show.tip.label=TRUE, show.node.label=TRUE)
dev.off()

pdf('FGA_treeNewNodesDefaultNoNodes.pdf',width=10, height=10,pointsize=0.1)
plot(testTree,edge.width=1, tip.color=orderedColors$Color, 
     show.tip.label=TRUE, show.node.label=FALSE)
dev.off()

pdf('FGA_treeNewNodesDefault.pdf',width=10, height=10,pointsize=0.1)
plot(testTree,edge.width=1, tip.color=orderedColors$Color, 
     show.tip.label=TRUE, show.node.label=TRUE)
dev.off()

# Get all subtrees using corrected names for DAVID input
testTree$tip.label <- fixedSymbols$Suggested.Symbol
allSubtrees <- subtrees(testTree)

node182 <- allSubtrees[182][[1]]$tip.label
node203 <- allSubtrees[203][[1]]$tip.label
node232 <- allSubtrees[232][[1]]$tip.label
node241 <- allSubtrees[241][[1]]$tip.label
node269 <- allSubtrees[269][[1]]$tip.label
node202 <- allSubtrees[202][[1]]$tip.label
node183 <- allSubtrees[183][[1]]$tip.label
node10 <- allSubtrees[10][[1]]$tip.label
node125 <- allSubtrees[125][[1]]$tip.label
node126 <- allSubtrees[126][[1]]$tip.label
node6 <- allSubtrees[6][[1]]$tip.label


# Get significant by DAVID
node6Enrichment <- read.csv("davidOutput/node6David.txt", sep='\t')
node6Enrichment <- node6Enrichment[node6Enrichment$Benjamini<0.05,]
node6Enrichment$GO<- sapply(strsplit(node6Enrichment$Term, split="~"), "[[", 1)
write.table(node6Enrichment, "davidOutput/signif/node6Enrichment.tsv", sep='\t', row.names=FALSE)

node10Enrichment <- read.csv("davidOutput/node10David.txt", sep='\t')
node10Enrichment <- node10Enrichment[node10Enrichment$Benjamini<0.05,]
node10Enrichment$GO<- sapply(strsplit(node10Enrichment$Term, split="~"), "[[", 1)
write.table(node10Enrichment, "davidOutput/signif/node10Enrichment.tsv", sep='\t', row.names=FALSE)

node125Enrichment <- read.csv("davidOutput/node125David.txt", sep='\t')
node125Enrichment <- node125Enrichment[node125Enrichment$Benjamini<0.05,]
node125Enrichment$GO<- sapply(strsplit(node125Enrichment$Term, split="~"), "[[", 1)
write.table(node125Enrichment, "davidOutput/signif/node125Enrichment.tsv", sep='\t', row.names=FALSE)

node126Enrichment <- read.csv("davidOutput/node126David.txt", sep='\t')
node126Enrichment <- node126Enrichment[node126Enrichment$Benjamini<0.05,]
node126Enrichment$GO<- sapply(strsplit(node126Enrichment$Term, split="~"), "[[", 1)
write.table(node126Enrichment, "davidOutput/signif/node126Enrichment.tsv", sep='\t', row.names=FALSE)

node182Enrichment <- read.csv("davidOutput/node182David.txt", sep='\t')
node182Enrichment <- node182Enrichment[node182Enrichment$Benjamini<0.05,]
node182Enrichment$GO<- sapply(strsplit(node182Enrichment$Term, split="~"), "[[", 1)
write.table(node182Enrichment, "davidOutput/signif/node182Enrichment.tsv", sep='\t', row.names=FALSE)

node183Enrichment <- read.csv("davidOutput/node183David.txt", sep='\t')
node183Enrichment <- node183Enrichment[node183Enrichment$Benjamini<0.05,]
node183Enrichment$GO<- sapply(strsplit(node183Enrichment$Term, split="~"), "[[", 1)
write.table(node183Enrichment, "davidOutput/signif/node183Enrichment.tsv", sep='\t', row.names=FALSE)

node202Enrichment <- read.csv("davidOutput/node202David.txt", sep='\t')
node202Enrichment <- node202Enrichment[node202Enrichment$Benjamini<0.05,]
node202Enrichment$GO<- sapply(strsplit(node202Enrichment$Term, split="~"), "[[", 1)
write.table(node202Enrichment, "davidOutput/signif/node202Enrichment.tsv", sep='\t', row.names=FALSE)

node203Enrichment <- read.csv("davidOutput/node203David.txt", sep='\t')
node203Enrichment <- node203Enrichment[node203Enrichment$Benjamini<0.05,]
node203Enrichment$GO<- sapply(strsplit(node203Enrichment$Term, split="~"), "[[", 1)
write.table(node203Enrichment, "davidOutput/signif/node203Enrichment.tsv", sep='\t', row.names=FALSE)

node232Enrichment <- read.csv("davidOutput/node232David.txt", sep='\t')
node232Enrichment <- node232Enrichment[node232Enrichment$Benjamini<0.05,]
node232Enrichment$GO<- sapply(strsplit(node232Enrichment$Term, split="~"), "[[", 1)
write.table(node232Enrichment, "davidOutput/signif/node232Enrichment.tsv", sep='\t', row.names=FALSE)

node241Enrichment <- read.csv("davidOutput/node241David.txt", sep='\t')
node241Enrichment <- node241Enrichment[node241Enrichment$Benjamini<0.05,]
node241Enrichment$GO<- sapply(strsplit(node241Enrichment$Term, split="~"), "[[", 1)
write.table(node241Enrichment, "davidOutput/signif/node241Enrichment.tsv", sep='\t', row.names=FALSE)

node269Enrichment <- read.csv("davidOutput/node269David.txt", sep='\t')
node269Enrichment <- node269Enrichment[node269Enrichment$Benjamini<0.05,]
node269Enrichment$GO<- sapply(strsplit(node269Enrichment$Term, split="~"), "[[", 1)
write.table(node269Enrichment, "davidOutput/signif/node269Enrichment.tsv", sep='\t', row.names=FALSE)
