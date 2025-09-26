args=commandArgs(trailingOnly=TRUE) # importing arguments from bash script
len=length(args)


if (len < 6) {
    stop("Error: minimum required files not provided")
} else {
    work_folder <- args[1]
    in_folder <- args[2]
    filename_proband <- args[3]
    unaffectedFam <- args[4]
    affectedFam <- args[5]
    familyID <- args[6]
}


argsDf <- data.frame(variables=c("work_folder", "in_folder", "filename_proband", "unaffectedFam", "affectedFam", "familyID"), Values=c(args[1],args[2],args[3],args[4],args[5],args[6]))


####### Check for empty entries #########
## check entry of proband members
checkProband <- function(filename_proband) {
if (missing(filename_proband) || filename_proband == "" || grepl("^\\s*$", filename_proband)) {
    stop("Error: variable proband_rohFile cannot be empty. If none available, enter FALSE")
  }    
}

## check entry of unaffected members
checkUnaffected <- function(unaffectedFam) {
if (missing(unaffectedFam) || unaffectedFam == "" || grepl("^\\s*$", unaffectedFam)) {
    stop("Error: variable unaffectedFam_rohFile cannot be empty. If none available, enter FALSE")
  }    
}

## check entry of affected members
checkAffected <- function(affectedFam) {
if (missing(affectedFam) || affectedFam == "" || grepl("^\\s*$", affectedFam)) {
    stop("Error: variable affectedFam_rohFile cannot be empty. If none available, enter FALSE")
  }    
}



####### Display provided sample information before proceeding ############
cat(paste0("The sample/family ID is: ",familyID), "\n")
cat("\n")
cat(paste0("The proband is: ",filename_proband), "\n")
cat("\n")
if (unaffectedFam!="FALSE") {
    unaffectedFam_rohFiles <- unlist(strsplit(unaffectedFam, " "))
    cat("Unaffected Family ROH Files:\n", paste(unaffectedFam_rohFiles, collapse = "\n"), "\n")
    cat("\n")
} else {
    cat("No unaffected family members provided", "\n")
    cat("\n")
}

if (affectedFam!="FALSE") {
    affectedFam_rohFiles <- unlist(strsplit(affectedFam, " "))
    cat("Affected Family ROH Files:\n", paste(affectedFam_rohFiles, collapse = "\n"), "\n")
    cat("\n")
} else {
    cat("No affected family members provided", "\n")
    cat("\n")
}


########### Get date and time stamps #############
date <- format(Sys.Date()) 
time <- format(Sys.time(), "%H_%M_%S")


########## Install and load packages from CRAN to local folder if package not located in .libPaths #########
install_if_missing <- function(packages) {
    newPackages <- packages[!(packages %in% installed.packages(lib.loc=.libPaths())[,"Package"])]
    
    if (length(newPackages)) {
        newPackages1 <- packages[!(packages %in% installed.packages(lib.loc=paste0(work_folder,"/R_libs/"))[,"Package"])]
        if (length(newPackages1)) {
            install.packages(newPackages1, paste0(work_folder,"/R_libs/"), repos = "https://cloud.r-project.org", dependencies = TRUE)
            library(newPackages1, lib=paste0(work_folder,"/R_libs"), character.only = TRUE)
            cat(paste(newPackages1, "package installed and loaded locally"), "\n")
        } else {
            library(packages, lib=paste0(work_folder,"/R_libs"), character.only = TRUE)
            cat(paste(packages, "package loaded locally"),"\n")
        }

    } else {
        installedPackagePath <- unique(installed.packages()[packages,"LibPath"])
        library(packages, lib=installedPackagePath, character.only = TRUE)
        cat(paste(packages, "package loaded from system"), "\n")
    }
}

install_if_missing("tidyverse")
install_if_missing("NCmisc")
install_if_missing("ggrepel")
install_if_missing("labeling")
#####################################################################


######### Create output directory for each run #########
dirPathName <- paste0(work_folder,"/outputs/",familyID,"_",date,"_",time)
if (dir.exists(dirPathName)==FALSE) {
    dir.create(dirPathName)
} else {
    stop("Error: Directory already exists")
}

############ Reading proband data ############
proband <- read_tsv(paste0(in_folder,"/",filename_proband),comment="##")
#proband=read_tsv(filename_proband,comment="##")

colnames(proband)[colnames(proband) == "#Chr"] = "Chr"
probandNewDf <- as.data.frame(proband)




############### Performing segregation of roh regions ########################
if (unaffectedFam=="FALSE" && affectedFam=="FALSE") {
    cat("Only proband provided. No segregation performed", "\n")
    ### Save updated data as a csv file 
    probandNewDf$`Total ROH Size` = rep(sum(probandNewDf$Size), times=nrow(probandNewDf))
    new_file <- paste(dirPathName,"/",familyID,"_probandOnly.csv",sep="")
    #write.csv(probandNewDf, file=new_file, row.names = FALSE) 

} else if (unaffectedFam=="FALSE" & affectedFam!="FALSE") { # If only affected family members (siblings) are included
    cat("Performing segregation. Only commonly shared regions will be included.", "\n")
    for (i in affectedFam_rohFiles) {
        # A nice progress bar
        max=100; for (cc in 1:max) { loop.tracker(cc,max); wait(0.004,"s") }

        finDf <- data.frame(Chr = character(), Begin = numeric(), End = numeric(), Size = numeric(), Percentage_homozygosity = numeric())  # Initialize dataframe
        affFam <- as.data.frame(read_tsv(paste0(in_folder,"/",i), comment="##")) # Read roh data for affected members
        colnames(affFam)[colnames(affFam) == "#Chr"] = "Chr"
        for (j in 1:nrow(probandNewDf)) {
            for (k in 1:nrow(affFam)) {
                if (probandNewDf$Chr[j]==affFam$Chr[k] & probandNewDf$Begin[j]<=affFam$End[k] & probandNewDf$End[j]>=affFam$Begin[k]) {
                    beginVal <- max(probandNewDf$Begin[j], affFam$Begin[k])
                    endVal <- min(probandNewDf$End[j], affFam$End[k])
                    rohSize <- endVal-beginVal
                    tempDf <- data.frame(Chr=probandNewDf$Chr[j],Begin=beginVal,End=endVal,Size=rohSize/1000000,Percentage_homozygosity=probandNewDf$Percentage_homozygosity[j])
                    finDf <- rbind(finDf, tempDf)
                }
            }
        }
        probandNewDf <- finDf # Update proband roh regions after each segregation
    }
    ### Save updated data as a csv file 
    probandNewDf$`Total ROH Size`  <-  rep(sum(probandNewDf$Size), times=nrow(probandNewDf))
    new_file <- paste(dirPathName,"/",familyID,"_proband&affected.csv",sep="")
    write.csv(probandNewDf, file=new_file, row.names = FALSE) 


} else if (unaffectedFam!="FALSE" & affectedFam=="FALSE") { # If only unaffected family members (parents/siblings) are included
    cat("Performing segregation. All commonly shared regions will be excluded.", "\n")
    for (i in unaffectedFam_rohFiles) {
        # A nice progress bar
        max=100; for (cc in 1:max) { loop.tracker(cc,max); wait(0.004,"s") }

        finDf <- data.frame(Chr = character(), Begin = numeric(), End = numeric(), Size = numeric(), Percentage_homozygosity = numeric())  # Initialize dataframe
        unaffFam <- as.data.frame(read_tsv(paste0(in_folder,"/",i), comment="##")) # Read roh data for unaffected members
        colnames(unaffFam)[colnames(unaffFam) == "#Chr"] = "Chr"
        for (j in 1:nrow(probandNewDf)) {
            # Function to check overlap
            overlap <- any(sapply(1:nrow(unaffFam),function(k) {
                begin2 <- unaffFam$Begin[k]
                end2 <- unaffFam$End[k]
                chr <- unaffFam$Chr[k]
                (probandNewDf$Begin[j]<=end2 & probandNewDf$End[j]>=begin2 & probandNewDf$Chr[j]==chr)
            }))
            overlap_indices <- sapply(1:nrow(unaffFam), function(k){
                if (probandNewDf$Begin[j]<=unaffFam$End[k] & probandNewDf$End[j]>=unaffFam$Begin[k] & probandNewDf$Chr[j]==unaffFam$Chr[k]) {
                    return(k)
                } else {
                    return(NA)
                }
            })

            ### Check for overlap
            if (!overlap) {
                tempDf <- data.frame(Chr=probandNewDf$Chr[j],Begin=probandNewDf$Begin[j],End=probandNewDf$End[j],Size=(probandNewDf$End[j]-probandNewDf$Begin[j])/1000000,Percentage_homozygosity=probandNewDf$Percentage_homozygosity[j])
                finDf <- rbind(finDf, tempDf)
                #print(tempDf)
            } else {
                ind <- overlap_indices[!is.na(overlap_indices)]
                if (probandNewDf$Begin[j] < unaffFam$Begin[ind]) { # Check for non-overlapping on the left side
                    tempDf <- data.frame(Chr=probandNewDf$Chr[j],Begin=probandNewDf$Begin[j],End=min(probandNewDf$End[j],unaffFam$Begin[ind]-1),Size=(probandNewDf$End[j]-probandNewDf$Begin[j])/1000000,Percentage_homozygosity=probandNewDf$Percentage_homozygosity[j])
                }
                if (probandNewDf$End[j] > unaffFam$End[ind]) {  # Check for non-overlapping on the right side
                    tempDf <- data.frame(Chr=probandNewDf$Chr[j],Begin=max(probandNewDf$Begin[j],unaffFam$End[ind]+1),End=probandNewDf$End[j],Size=(probandNewDf$End[j]-probandNewDf$Begin[j])/1000000,Percentage_homozygosity=probandNewDf$Percentage_homozygosity[j])
                }
                finDf <- rbind(finDf, tempDf)
                #tempDf=data.frame(Chr=probandNewDf$Chr[j],Begin=probandNewDf$Begin[j],End=probandNewDf$End[j],Size=(probandNewDf$End[j]-probandNewDf$Begin[j])/1000000,Percentage_homozygosity=probandNewDf$Percentage_homozygosity[j])
            }
        }
        probandNewDf <- finDf # Update proband roh regions after each segregation
    }
    ### Save updated data as a csv file 
    probandNewDf$`Total ROH Size`  <-  rep(sum(probandNewDf$Size), times=nrow(probandNewDf))
    new_file <- paste(dirPathName,"/",familyID,"_proband&unaffected.csv",sep="")
    write.csv(probandNewDf, file=new_file, row.names = FALSE) 


} else if (unaffectedFam!="FALSE" & affectedFam!="FALSE") { # If both affected and unaffected family members (parents/siblings) are included 
    cat("Performing segregation. Only commonly shared regions will be included from affected and all shared regions from unaffected will be excluded.", "\n")
    # Perform segregation for affected family members first
    for (i in affectedFam_rohFiles) {
        # A nice progress bar
        max=100; for (cc in 1:max) { loop.tracker(cc,max); wait(0.004,"s") }

        finDf <- data.frame(Chr = character(), Begin = numeric(), End = numeric(), Size = numeric(), Percentage_homozygosity = numeric())  # Initialize dataframe
        affFam <- as.data.frame(read_tsv(paste0(in_folder,"/",i), comment="##")) # Read roh data for affected members
        colnames(affFam)[colnames(affFam) == "#Chr"] = "Chr"
        for (j in 1:nrow(probandNewDf)) {
            for (k in 1:nrow(affFam)) {
                if (probandNewDf$Chr[j]==affFam$Chr[k] & probandNewDf$Begin[j]<=affFam$End[k] & probandNewDf$End[j]>=affFam$Begin[k]) {
                    beginVal <- max(probandNewDf$Begin[j], affFam$Begin[k])
                    endVal <- min(probandNewDf$End[j], affFam$End[k])
                    rohSize <- endVal-beginVal
                    tempDf <- data.frame(Chr=probandNewDf$Chr[j],Begin=beginVal,End=endVal,Size=rohSize/1000000,Percentage_homozygosity=probandNewDf$Percentage_homozygosity[j])
                    finDf <- rbind(finDf, tempDf)
                }
            }
        }
        probandNewDf <- finDf # Update proband roh regions after each segregation
    }

    # Perform segregation of unaffected family members after segregation of affected members
    for (i in unaffectedFam_rohFiles) {
        # A nice progress bar
        max=100; for (cc in 1:max) { loop.tracker(cc,max); wait(0.004,"s") }

        finDf <- data.frame(Chr = character(), Begin = numeric(), End = numeric(), Size = numeric(), Percentage_homozygosity = numeric()) # Initialize dataframe
        unaffFam <- as.data.frame(read_tsv(paste0(in_folder,"/",i), comment="##")) # Read roh data for unaffected members
        colnames(unaffFam)[colnames(unaffFam) == "#Chr"] = "Chr"
        for (j in 1:nrow(probandNewDf)) {
            # Function to check overlap
            overlap <- any(sapply(1:nrow(unaffFam),function(k) {
                begin2 <- unaffFam$Begin[k]
                end2 <- unaffFam$End[k]
                chr <- unaffFam$Chr[k]
                (probandNewDf$Begin[j]<=end2 & probandNewDf$End[j]>=begin2 & probandNewDf$Chr[j]==chr)
            }))
            overlap_indices <- sapply(1:nrow(unaffFam), function(k){
                if (probandNewDf$Begin[j]<=unaffFam$End[k] & probandNewDf$End[j]>=unaffFam$Begin[k] & probandNewDf$Chr[j]==unaffFam$Chr[k]) {
                    return(k)
                } else {
                    return(NA)
                }
            })

            ### Check for overlap
            if (!overlap) {
                tempDf <- data.frame(Chr=probandNewDf$Chr[j],Begin=probandNewDf$Begin[j],End=probandNewDf$End[j],Size=(probandNewDf$End[j]-probandNewDf$Begin[j])/1000000,Percentage_homozygosity=probandNewDf$Percentage_homozygosity[j])
                finDf <- rbind(finDf, tempDf)
                #print(tempDf)
            } else {
                ind <- overlap_indices[!is.na(overlap_indices)]
                if (probandNewDf$Begin[j] < unaffFam$Begin[ind]) { # Check for non-overlapping on the left side
                    tempDf <- data.frame(Chr=probandNewDf$Chr[j],Begin=probandNewDf$Begin[j],End=min(probandNewDf$End[j],unaffFam$Begin[ind]-1),Size=(probandNewDf$End[j]-probandNewDf$Begin[j])/1000000,Percentage_homozygosity=probandNewDf$Percentage_homozygosity[j])
                }
                if (probandNewDf$End[j] > unaffFam$End[ind]) {  # Check for non-overlapping on the right side
                    tempDf <- data.frame(Chr=probandNewDf$Chr[j],Begin=max(probandNewDf$Begin[j],unaffFam$End[ind]+1),End=probandNewDf$End[j],Size=(probandNewDf$End[j]-probandNewDf$Begin[j])/1000000,Percentage_homozygosity=probandNewDf$Percentage_homozygosity[j])
                }
                finDf <- rbind(finDf, tempDf)
                #tempDf=data.frame(Chr=probandNewDf$Chr[j],Begin=probandNewDf$Begin[j],End=probandNewDf$End[j],Size=(probandNewDf$End[j]-probandNewDf$Begin[j])/1000000,Percentage_homozygosity=probandNewDf$Percentage_homozygosity[j])
            }
        }
        probandNewDf <- finDf # Update proband roh regions after each segregation
    }
    ### Save updated data as a csv file 
    probandNewDf$`Total ROH Size` <- rep(sum(probandNewDf$Size), times=nrow(probandNewDf))
    new_file <- paste(dirPathName,"/",familyID,"_affected&unaffected.csv",sep="")
    write.csv(probandNewDf, file=new_file, row.names = FALSE)
}

#print(probandNewDf)
##############################################################################################################################