#!/bin/bash

echo


############### Required variables ##################
# Working dir (path to fastmap folder with databases, scripts, and outputs. Default: current directory)
workDir=$(pwd)/Segregation_Analysis

# Full path to input dir with autozygosity information
inputDir=${workDir}/testFiles

# Proband roh info file (Only one file, required, currently only tsvs are applicable)
probandFile="testProbandData.tsv"

# Roh info files for unaffected family members, space separated, currently only tsvs are applicable (if none, indicate FALSE)
unaffectedFiles="testUnaffectedMember1Data.tsv testUnaffectedMember2Data.tsv"

# Roh info files for affected family members, space separated, currently only tsvs are applicable (if none, indicate FALSE)
affectedFiles="testAffectedMemberData.tsv"

# Csv file with all variants information. Please use full path
variantFilepath=${workDir}/testFiles/TestVariantData.csv

# Family/Sample Identifier used for identifying variants from variant list and for naming outputs
sampleID="Test"


# Parse command-line arguments for passing values to individual variables 
while [[ $# -gt 0 ]]; do
  case $1 in
    -workDir|--workDir)
      workDir=$2
      shift 2
      ;;
    -inputDir|--inputDir)
      inputDir=$2
      shift 2
      ;;
    -probandFile|--probandFile)
      probandFile=$2
      shift 2
      ;;
    -unaffectedFiles|--unaffectedFiles)
      unaffectedFiles=$2
      shift 2
      ;;
    -affectedFiles|--affectedFiles)
      affectedFiles=$2
      shift 2
      ;;
    -sampleID|--sampleID)
      sampleID=$2
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done


## Script format
Rscript ${workDir}/ROHSegregationAnalysis.R ${workDir} ${inputDir} ${probandFile} "${unaffectedFiles}" "${affectedFiles}" ${sampleID} 2>output.log