library(readr)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)

# Note to Peter: Datasets for plotting made available in Supplementary Tables (not numbered yet). Filenames used here are original table names

## For Plotting GO Term Enrichments from DAVID Overrepresentation Analysis ##
GOTermsDf=read.csv("Data/mapcp/mapcp_GO_Pathway_analyses/DAVID_TopRecessiveCausativeGenes_GOTerms_04222025.csv") # GO Terms Data from DAVID will go here
GOTermsDf_sorted = GOTermsDf[order(GOTermsDf$PValue),]

GOTermsDf_sorted$GeneRatio = GOTermsDf_sorted$Count/GOTermsDf_sorted$List.Total # Calculating gene ratio
# Getting top 25 terms sorted by GeneRatio
GOTermsDf_Top25 = GOTermsDf_sorted[1:25,] # Getting top 25 terms sorted by GeneRatio
GOTermsDf_Top25 = GOTermsDf_Top25[order(GOTermsDf_Top25$GeneRatio, decreasing = FALSE),]
GOTermsDf_Top25$TermOrder = c(1:25)

GOTermsDf_Top25$Term = factor(GOTermsDf_Top25$Term, levels = GOTermsDf_Top25$Term)
  
# Plot figure
plt=ggplot(GOTermsDf_Top25, aes(x=`GeneRatio`, y=`Term`, color=`PValue`, size=`Count`))
plt=plt+geom_point()
plt=plt+theme(axis.text.y = element_text(color = "black", size = 12))
plt=plt+scale_color_gradient(low = "blue", high = "red")+theme_classic() #Add color-gradient
ggsave("GO_Terms_Enrichment.pdf", plot = plt) # # save plot
############################################################################


## For Plotting KEGG Pathway Term Enrichments from DAVID Overrepresentation Analysis ##
KEGGTermsDf=read.csv("Data/mapcp/mapcp_GO_Pathway_analyses/DAVID_TopRecessiveCausativeGenes_KEGGTerms_04222025.csv") # KEGG Pathway Terms Data from DAVID will go here
KEGGTermsDf_sorted = KEGGTermsDf[order(KEGGTermsDf$PValue),]

KEGGTermsDf_sorted$GeneRatio = KEGGTermsDf_sorted$Count/KEGGTermsDf_sorted$List.Total # Calculating gene ratio
# Getting top 25 terms sorted by GeneRatio
KEGGTermsDf_Top25 = KEGGTermsDf_sorted[1:25,] 
KEGGTermsDf_Top25 = KEGGTermsDf_Top25[order(KEGGTermsDf_Top25$GeneRatio, decreasing = FALSE),]
KEGGTermsDf_Top25$TermOrder = c(1:25)

# For only getting KEGG Terms
KEGGTermsDf_Ex = KEGGTermsDf_sorted[KEGGTermsDf_sorted$Category=='KEGG_PATHWAY',]
KEGGTermsDf_Top25 = KEGGTermsDf_Ex[order(KEGGTermsDf_Ex$GeneRatio, decreasing = FALSE),]

KEGGTermsDf_Top25$Term = factor(KEGGTermsDf_Top25$Term, levels = KEGGTermsDf_Top25$Term)

# Plot figure
plt1=ggplot(KEGGTermsDf_Top25, aes(x=`GeneRatio`, y=`Term`, color=`PValue`, size=`Count`))
plt1=plt1+geom_point()
plt1=plt1+theme(axis.text.y = element_text(color = "black", size = 12))
plt1=plt1+scale_color_gradient(low = "blue", high = "red")+theme_classic() #Add color-gradient
ggsave("KEGG_Pathway_Terms_Enrichment.pdf", plot = plt1) # save plot
###################################################################
