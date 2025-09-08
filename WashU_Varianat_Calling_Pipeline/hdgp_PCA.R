#!/usr/bin/env Rscript

##################################
# Purpose: To handle the PCA coordinate file output from Laser/Trace 
# Author: Sheng Chih Jin <shengchih.jin@yale.edu>
# Modifier: Yung-Chun Wang <yung-chun@wustl.edu> 
# Language: R
# Version: 2
# Comment: Require Study.ProPC.coord (Study) and HGDP_938.RefPC.coord (Reference)
# Last Modified Date: 10-10-2022
###################################

## R-Command Line info
### https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/

library("optparse")

ref_coord = "/storage1/fs1/jin810/Active/testing/yung-chun/database/reference/LASER/HGDP_938.RefPC.coord"

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="full directory for the input file [study.coord].", metavar="character"),
  make_option(c("-r", "--ref"), type="character", default=ref_coord, 
              help="output file prefix name [default=HGDP_938.RefPC.coord].", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}

## Create output folder to store results
out_dir = dirname(opt$input)

### Test name
input_name = basename(opt$input)
test_name = strsplit(input_name, split="[.]")[[1]][1]

folder_name = paste(test_name, "PCA_results", sep ="-", collapse = NULL)
results = paste(out_dir, folder_name, sep = "/", collapse = NULL)
dir.create(results)

## Read reference coordinate
Ref = read.table(file=opt$ref, header=T, stringsAsFactors=FALSE)

## Read sample coordinate
mydata = read.table(file=opt$input, header=TRUE, sep="\t", stringsAsFactors=FALSE)

## Change data type
Ref$PC1 <- as.numeric(as.character(Ref$PC1))
Ref$PC2 <- as.numeric(as.character(Ref$PC2))
Ref$PC3 <- as.numeric(as.character(Ref$PC3))

## Create each category group
pca.afr = Ref[which(Ref[[1]] %in% c("BantuSouthAfrica","BantuKenya","BiakaPygmy","Mandenka","MbutiPygmy","San","Yoruba")),]
pca.ceu = Ref[which(Ref[[1]] %in% c("Adygei","Basque","French","Italian","Orcadian","Russian","Sardinian","Tuscan")),]
pca.middleeast = Ref[which(Ref[[1]] %in% c("Bedouin","Druze","Mozabite","Palestinian")),]
pca.southasia = Ref[which(Ref[[1]] %in% c("Balochi","Brahui","Burusho","Hazara","Kalash","Makrani","Pathan","Sindhi","Uygur")),]
pca.eastasia = Ref[which(Ref[[1]] %in% c("Cambodian","Dai","Daur","Han","Han-NChina","Hezhen","Japanese","Lahu","Miao","Mongola","Naxi","Oroqen","She","Tu","Tujia","Xibo","Yakut","Yi")),]
pca.oceania = Ref[which(Ref[[1]] %in% c("Melanesian","Papuan")),]
pca.america = Ref[which(Ref[[1]] %in% c("Colombian","Karitiana","Maya","Pima","Surui")),]

## Set range of PC1, PC2, PC3
small.cex <- 0.7
rangePC1 <- c(min(c(Ref$PC1,mydata$PC1)), max(c(Ref$PC1,mydata$PC1)));
rangePC2 <- c(min(c(Ref$PC2,mydata$PC2)), max(c(Ref$PC2,mydata$PC2)));
rangePC3 <- c(min(c(Ref$PC3,mydata$PC3)), max(c(Ref$PC3,mydata$PC3)));

######### Cases at the front PC1 vs PC2
## Call the pdf command to start the plot
pdf(paste(results,"LASER_PCA_Analysis_GMFK3_b38.pdf", sep = "/", collapse = NULL), width=7, height=7.9)
par(mar=c(5,5,4,2))

## PCreate the plot
# Africans
plot(pca.afr$PC1, pca.afr$PC2, xlim= rangePC1, ylim=rangePC2, xlab="PC1", ylab="PC2", cex=small.cex, cex.lab=2, cex.axis=1.7, cex.main=2.5, main="PCA of Cases",asp=1,col="orange", pch=1)
# Europeans
points(pca.ceu$PC1, pca.ceu$PC2, col="blue", cex=small.cex, pch=1)
# Middle easterns
points(pca.middleeast$PC1, pca.middleeast$PC2, col="yellow", cex=small.cex, pch=1)
# South asians
points(pca.southasia$PC1,pca.southasia$PC2, col="red", cex=small.cex, pch=1)
# East asians
points(pca.eastasia$PC1,pca.eastasia$PC2, col="pink", cex=small.cex, pch=1)
# Oceania
points(pca.oceania$PC1, pca.oceania$PC2, col="light green", cex=small.cex, pch=1)
# America
points(pca.america$PC1, pca.america$PC2, col="purple", cex=small.cex, pch=1)
# Cases
points(mydata$PC1, mydata$PC2, col="black", cex=small.cex, pch=1)

legend("bottomleft",bty='n',  c("Cases", "European", "African", "American","Easet Asian","South Asian","Middle Eastern","Oceania"), 
		col=c("black", "blue", "orange", "purple","pink","red","yellow","light green"), pch=1, cex=0.6, ncol=1)

## Run dev.off() to create the file!
dev.off()

######### Cases at the background PC1 vs PC2
pdf(paste(results,"LASER_PCA_Analysis_GMKF3_PC1_PC2.pdf", sep = "/", collapse = NULL), width=7, height=7.9)
par(mar=c(5,5,4,2))

# Cases
plot(mydata$PC1, mydata$PC2, xlim= rangePC1, ylim=rangePC2, xlab="PC1", ylab="PC2", cex=small.cex, cex.lab=2, cex.axis=1.7, cex.main=2.5, main="PCA of Cases",asp=1,col="black", pch=1)
# Europeans
points(pca.ceu$PC1, pca.ceu$PC2, col="blue", cex=small.cex, pch=1)
# Middle easterns
points(pca.middleeast$PC1, pca.middleeast$PC2, col="yellow", cex=small.cex, pch=1)
# South asians
points(pca.southasia$PC1,pca.southasia$PC2, col="red", cex=small.cex, pch=1)
# East asians
points(pca.eastasia$PC1,pca.eastasia$PC2, col="pink", cex=small.cex, pch=1)
# Oceania
points(pca.oceania$PC1, pca.oceania$PC2, col="light green", cex=small.cex, pch=1)
# America
points(pca.america$PC1, pca.america$PC2, col="purple", cex=small.cex, pch=1)
# Africans
points(pca.afr$PC1, pca.afr$PC2, col="orange",cex=small.cex, pch=1)

legend("bottomleft",bty='n',  c("Cases", "European", "African", "American","Easet Asian","South Asian","Middle Eastern","Oceania"), 
       col=c("black", "blue", "orange", "purple","pink","red","yellow","light green"), pch=1, cex=0.6, ncol=1)

dev.off()

######### Cases at the background PC1 vs PC3
pdf(paste(results,"LASER_PCA_Analysis_GMFK3_PC1_PC3.pdf", sep = "/", collapse = NULL), width=7, height=7.9)
par(mar=c(5,5,4,2))

# Cases
plot(mydata$PC1, mydata$PC3, xlim= rangePC1, ylim=rangePC3, xlab="PC1", ylab="PC3", cex=small.cex, cex.lab=2, cex.axis=1.7, cex.main=2.5, main="PCA of Cases",asp=1,col="black", pch=1)
# Europeans
points(pca.ceu$PC1, pca.ceu$PC3, col="blue", cex=small.cex, pch=1)
# Middle easterns
points(pca.middleeast$PC1, pca.middleeast$PC3, col="yellow", cex=small.cex, pch=1)
# South asians
points(pca.southasia$PC1,pca.southasia$PC3, col="red", cex=small.cex, pch=1)
# East asians
points(pca.eastasia$PC1,pca.eastasia$PC3, col="pink", cex=small.cex, pch=1)
# Oceania
points(pca.oceania$PC1, pca.oceania$PC3, col="light green", cex=small.cex, pch=1)
# America
points(pca.america$PC1, pca.america$PC3, col="purple", cex=small.cex, pch=1)
# Africans
points(pca.afr$PC1, pca.afr$PC3, col="orange",cex=small.cex, pch=1)

legend("bottomleft",bty='n',  c("Cases", "European", "African", "American","Easet Asian","South Asian","Middle Eastern","Oceania"), 
       col=c("black", "blue", "orange", "purple","pink","red","yellow","light green"), pch=1, cex=0.6, ncol=1)

dev.off()

###### Write the SampleID and its correpsonding ethinicity into txt file
## Classificatio
# African
index = which(mydata$PC2 < mydata$PC1-280)
ID = mydata[index,]$indivID
Eth = rep("African",length(ID))
Dataframe.af = data.frame(ID=ID, Ethinicity=Eth)
file = paste(results,"African.txt", sep = "/", collapse = NULL)
write.table(Dataframe.af, file=file, col.names=F, row.names=F, sep="\t", quote=F)

# South Asians
index = which(mydata$PC2 >= -50 & mydata$PC2 < 150 & mydata$PC2 < -mydata$PC1 + 190)
ID = mydata[index,]$indivID
Eth = rep("South_Asian",length(ID))
Dataframe.sa = data.frame(ID=ID, Ethinicity=Eth)
file = paste(results,"South_Asian.txt", sep = "/", collapse = NULL)
write.table(Dataframe.sa, file=file, col.names=F, row.names=F, sep="\t", quote=F)

# Middle East
index = which(mydata$PC2 >= mydata$PC1 -280 & mydata$PC2 < mydata$PC1 + 125 & mydata$PC2 >= -mydata$PC1 + 190)
index.middleeast = index
ID = mydata[index,]$indivID
Eth = rep("Middle_East",length(ID))
Dataframe.me = data.frame(ID=ID, Ethinicity=Eth)
file = paste(results,"Middle_East.txt", sep = "/", collapse = NULL)
write.table(Dataframe.me, file=file, col.names=F, row.names=F, sep="\t", quote=F)

# European
index = which(mydata$PC2 >= 150)
index <- index[!index %in% index.middleeast]
ID = mydata[index,]$indivID
Eth = rep("European",length(ID))
Dataframe.e = data.frame(ID=ID, Ethinicity=Eth)
file = paste(results,"European.txt", sep = "/", collapse = NULL)
write.table(Dataframe.e, file=file, col.names=F, row.names=F, sep="\t", quote=F)

# East asians
index = which(mydata$PC3 > mydata$PC1 + 180 & mydata$PC2 < -50) 
ID = mydata[index,]$indivID
Eth = rep("East_Asian",length(ID))
Dataframe.ea = data.frame(ID=ID, Ethinicity=Eth)
file = paste(results,"East_Asian.txt", sep = "/", collapse = NULL)
write.table(Dataframe.ea, file=file, col.names=F, row.names=F, sep="\t", quote=F)

# Americans
index = which(mydata$PC3 < -220 & mydata$PC2 < -50)
ID = mydata[index,]$indivID
Eth = rep("American",length(ID))
Dataframe.a = data.frame(ID=ID, Ethinicity=Eth)
file = paste(results,"American.txt", sep = "/", collapse = NULL)
write.table(Dataframe.a, file=file, col.names=F, row.names=F, sep="\t", quote=F)

# Merge dataframe
Dataframe = rbind(Dataframe.af, Dataframe.sa, Dataframe.me, Dataframe.e, Dataframe.ea, Dataframe.a)
file = paste(results,"All_Ethnicity.txt", sep = "/", collapse = NULL)
write.table(Dataframe, file=file, col.names=F, row.names=F, sep="\t", quote=F)
