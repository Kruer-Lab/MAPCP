#!/bin/bash

#SBATCH -o /path/to/directory/slurm/out/slurm.%j.out #output file (Set your own sbatch output folder here)
#SBATCH -e /path/to/directory/slurm/error/slurm.%j.err #error file (Set your own sbatch error folder here)
#SBATCH -t 0-04:00:00 #d-hh:mm:ss
#SBATCH -N 1 #Number of nodes
##SBATCH --exclusive #for exlcusive node usage

# This code can perform autozygosity analysis on multiple trio vcfs using AutoMap (PMID: 33483490)
# Before using AutoMap, please clone the AutoMap code repository from Github using the following command
# git clone https://github.com/mquinodo/AutoMap
# for multiple vcf files, use --vcflist instead of --vcf, input must be a txt file containing the names of all VCF files
# for output including chromosome X, use --chrX option
# Following modules must be loaded on the hpc environment before autozygosity analysis using AutoMap:
# bcftools v1.15.1, bedtools v2.30.0, perl v5.26.0, gatk/4.2.5.0, r v4.1.0
# If running locally on a Linux system, ensure to install the aforementioned tools beforehand

source /path/to/directory/Autozygosity/AutozygositySettings.sh # Location of Autozygosity folder

echo
# ************************* Required Modules *********************************
module load bcftools-1.14-gcc-11.2.0
module load bedtools2-2.30.0-gcc-11.2.0
module load perl-5.26.2-gcc-12.1.0
module load gatk-4.2.6.1-gcc-11.2.0
module load r-4.2.2-gcc-11.2.0

# ************************* Optional Arguments for AutoMap ************************
dp=8 # minimal variant depth (DP), default = 8
binom=0.000001 # minimal p-value (binomial), default = 0.000001
pclow=0.25 # min alternative reads ratio for het variants (percaltlow), default = 0.25
pchigh=0.75 # max alternative reads ratio for het variants (percalthigh), default = 0.75
win=7 # sliding window size (window), default = 7
winthres=5 # window threshold (windowthres), default = 5
minsz=1.0 # min ROH size detected in Mb (minsize), default = 1.0
minvr=25 # min # of variant detected in ROH (minvar), default = 25
minpc=88 # min % of homozygous variant detected in ROH
mgap=10 # max gap between 2 variants in ROH in Mb (maxgap), default = 10
ext=1.0 # max extension at both ROH boundaries (extend), default = 1.0

# Excel file path with all variants
filepath_xl=/data/kruerlab/pbisarad/fastmap/VarList/All_variants_mapcp.xlsx

# calculating ROH of each individual id in the combined trio vcf file
for i in ${path_name[0]};do
vcf_file=${i}
vcf_folder=$(basename -a ${vcf_file} | sed 's/.trio.raw.vcf//')
mkdir -p "${output_dir}/${vcf_folder}/"
echo
echo "******************Autozygosity analysis of multi sample vcf using AutoMap**********************" 
# Note: Change directory to location of AutoMap here
bash /path/to/directory/AutoMap/AutoMap_v1.2.sh --vcf $vcf_file --multivcf --genome hg38 --out $output_dir/$vcf_folder/ \
--DP $dp --binomial $binom --percaltlow $pclow --percalthigh $pchigh --window $win --windowthres $winthres --minsize $minsz --minvar $minvr --minperc $minpc --maxgap $mgap --extend $ext
echo "DONE"
echo "**********************end of sample analysis******************************"
wait
done