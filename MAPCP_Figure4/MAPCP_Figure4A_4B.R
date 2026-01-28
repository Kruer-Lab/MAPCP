library(stringr)
library(readr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readxl)

## Figure 4A ##
rohData <- read_excel("/path/to/directory/SupplementaryTables_11.05.25.xlsx", sheet = "Supplementary Table 18", skip=1, col_names = TRUE) #change /path/to/directory to where the supplementary tables are

# Calculating quartile values for each cohort
gccpQuartile <- quantile(rohData$Tot_roh_proband_Mb[rohData$Cohort=='internal_control_cohort'], probs = c(0,0.25,0.5,0.75,1))
consanguineousNoQuartile <- quantile(rohData$Tot_roh_proband_Mb[rohData$Cohort=='consanguinity_no'], probs = c(0,0.25,0.5,0.75,1))
consanguineousYesQuartile <- quantile(rohData$Tot_roh_proband_Mb[rohData$Cohort=='consanguinity_yes'], probs = c(0,0.25,0.5,0.75,1))

# Boxplot for each cohort
plt1 <- ggplot(rohData, aes(Cohort, Tot_roh_proband_Mb))
plt1 <- plt1+geom_boxplot(outlier.shape = NA, position = "dodge2", aes(colour = Cohort),show.legend = FALSE) 
plt1 <- plt1+geom_jitter(width = 0.2, alpha = 0.2)
plt1 <- plt1+theme_classic()+xlab(" ")+ylab("Total RoH in Proband (Mb)")
plt1 <- plt1+scale_x_discrete(labels = c("Denied \n Consaguinity", "Accepted \n Consanguinity", "Control"))
plt1 <- plt1+ggtitle("Distribution of Total RoH Sizes in probands w.r.t Reported Consanguinity")
plt1 <- plt1+geom_label(x=rohData$Cohort[56], y=-5, label="n=26")
plt1 <- plt1+geom_label(x=rohData$Cohort[95], y=-5, label="n=158")
plt1 <- plt1+geom_label(x=rohData$Cohort[1], y=-5, label="n=55")

# plotting p-values
myComp <- list(c("consanguinity_no","consanguinity_yes"), c("consanguinity_yes","internal_control_cohort"), c("consanguinity_no","internal_control_cohort"))
plt1 <- plt1+stat_compare_means(comparisons = myComp, method="wilcox.test")
plt1
ggsave("Total_ROH_Size_vs_Reported_Consanguinity.pdf", plot = plt1)
############


## Figure 4B ##
ConsangunityDf = read_excel("/path/to/directory/SupplementaryTables_11.05.25.xlsx", sheet = "Supplementary Table 19", skip=1, col_names = TRUE) #change /path/to/directory to where the supplementary tables are

plt1 <- ggplot(ConsangunityDf, aes(x = Expected_Percentage_of_ROH_burden, y = Observed_Percentage_of_ROH_burden))
plt1 <- plt1+geom_point()+geom_smooth(aes(Expected_Percentage_of_ROH_burden,Observed_Percentage_of_ROH_burden), method="lm")
plt1 <- plt1+theme_classic()
plt1

# Performing paired t-test
t_test_result <- t.test(ConsangunityDf$Observed_Percentage_of_ROH_burden, ConsangunityDf$Expected_Percentage_of_ROH_burden, paired = TRUE)
print(t_test_result)

pval <- signif(t_test_result$p.value, 4)
tval <- signif(t_test_result$statistic, 4)
degfree <- t_test_result$parameter

plt1 <- plt1+annotate("text", x = 1.5, y = max(ConsangunityDf$Observed_Percentage_of_ROH_burden), label = paste0("t-test p value = ",pval), size = 5)
plt1 <- plt1+annotate("text", x = 1.5, y = max(ConsangunityDf$Observed_Percentage_of_ROH_burden)+1, label = paste0("t-test t value (",degfree,") = ",tval), size = 5)
plt1
ggsave("Total_ROH_Size_Observed_vs_Expected.pdf", plot = plt1)
###############


## Figure 4B ##
ConsangunityDf = read_csv("Data/mapcp/Tables/MAPCP_ROH_KnownPedigreeFams_Obs_vs_Exp_ROC_08202025.csv")

plt1 <- ggplot(ConsangunityDf, aes(x = Expected, y = Observed))
plt1 <- plt1+geom_point()+geom_smooth(aes(Expected,Observed), method="lm")
plt1 <- plt1+theme_classic()
plt1

# Performing paired t-test
t_test_result = t.test(ConsangunityDf$Observed, ConsangunityDf$Expected, paired = TRUE)
print(t_test_result)

pval <- signif(t_test_result$p.value, 4)
tval <- signif(t_test_result$statistic, 4)
degfree <- t_test_result$parameter

plt1 <- plt1+annotate("text", x = 1.5, y = max(ConsangunityDf$Observed), label = paste0("t-test p value = ",pval), size = 5)
plt1 <- plt1+annotate("text", x = 1.5, y = max(ConsangunityDf$Observed)+1, label = paste0("t-test t value (",degfree,") = ",tval), size = 5)
ggsave("Total_ROH_Size_Observed_vs_Expected.pdf", plot = plt1)

###############

