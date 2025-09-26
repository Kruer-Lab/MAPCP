library(readr)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readxl)


## For Plotting GO Term Enrichments from DAVID Overrepresentation Analysis ##
GOTermsDf=read_excel("/path/to/directory/SupplementaryTables_Revised.xlsx", sheet = "Supplementary Table 15", skip=1, col_names = TRUE) # Change /path/to/directory to where the supplemental file is
GOTermsDf_sorted = GOTermsDf[order(GOTermsDf$PValue),]

#Getting top 25 terms from p-value sorted data
GOTermsDf_Top25 = GOTermsDf_sorted[1:25,] 
#Sorting GO terms based on fold enrichment
GOTermsDf_Top25 = GOTermsDf_Top25[order(GOTermsDf_Top25$Fold.Enrichment, decreasing = FALSE),]
GOTermsDf_Top25$TermOrder = c(1:25)

GOTermsDf_Top25$Term = factor(GOTermsDf_Top25$Term, levels = GOTermsDf_Top25$Term)

# Plot figure
#plt=ggplot(GOTermsDf_Top25, aes(x=`GeneRatio`, y=`Term`, color=`PValue`, size=`Count`))
plt=ggplot(GOTermsDf_Top25, aes(x=`Fold.Enrichment`, y=`Term`, color=`PValue`, size=`Count`))
plt=plt+geom_point()
plt=plt+theme(axis.text.y = element_text(color = "black", size = 12))
plt=plt+scale_color_gradient(low = "turquoise4", high = "salmon")+theme_classic() #Add color-gradient
ggsave("GO_Terms_Enrichment.pdf", plot = plt) # # save plot
############################################################################


## For Plotting KEGG Pathway Term Enrichments from DAVID Overrepresentation Analysis ##
KEGGTermsDf=read_excel("/path/to/directory/SupplementaryTables_Revised.xlsx", sheet = "Supplementary Table 16", skip=1, col_names = TRUE)  # Change /path/to/directory to where the supplemental file is
KEGGTermsDf_sorted = KEGGTermsDf[order(KEGGTermsDf$PValue),]

# Getting top 25 terms from p-value sorted data
KEGGTermsDf_Top25 = KEGGTermsDf_sorted[1:25,] 
#Sorting KEGG terms based on fold enrichment
KEGGTermsDf_Top25 = KEGGTermsDf_Top25[order(KEGGTermsDf_Top25$Fold.Enrichment, decreasing = FALSE),]
KEGGTermsDf_Top25$TermOrder = c(1:25)

# For only getting KEGG Terms
KEGGTermsDf_Ex = KEGGTermsDf_sorted[KEGGTermsDf_sorted$Category=='KEGG_PATHWAY',]
#KEGGTermsDf_Top25 = KEGGTermsDf_Ex[order(KEGGTermsDf_Ex$GeneRatio, decreasing = FALSE),]
KEGGTermsDf_Top25 = KEGGTermsDf_Ex[order(KEGGTermsDf_Ex$Fold.Enrichment, decreasing = FALSE),]

KEGGTermsDf_Top25$Term = factor(KEGGTermsDf_Top25$Term, levels = KEGGTermsDf_Top25$Term)

# Plot figure
#plt1=ggplot(KEGGTermsDf_Top25, aes(x=`GeneRatio`, y=`Term`, color=`PValue`, size=`Count`))
plt1=ggplot(KEGGTermsDf_Top25, aes(x=`Fold.Enrichment`, y=`Term`, color=`PValue`, size=`Count`))
plt1=plt1+geom_point()
plt1=plt1+theme(axis.text.y = element_text(color = "black", size = 12))
plt1=plt1+scale_color_gradient(low = "turquoise4", high = "salmon")+theme_classic() #Add color-gradient
ggsave("KEGG_Pathway_Terms_Enrichment.pdf", plot = plt1) # save plot
###################################################################
