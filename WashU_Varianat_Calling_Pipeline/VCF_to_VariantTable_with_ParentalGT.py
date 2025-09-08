#!/usr/bin/python

###################################
# Purpose: Filter annovar-annotated vcf file, extract mutation carriers, and get rid of 0/0 and ./.
# Author: Sheng Chih Jin <shengchih.jin@yale.edu>
# Modifier: Yung-Chun Wang <yung-chun@wustl.edu>
# Language: Python3
# Version: 4
# Comment: Add sibling information after parent's information. Pedigree file should be `FamID, Proband, Mother, Father, Sibling1, Sibling2...` Remove the singleton output.
# Last Modified Date: 02-18-2023
###################################

##Setup the important stuff
#Import necessary packages
import sys, os, numbers, linecache, gzip, shutil

if len(sys.argv) != 3:
	sys.stderr.write("Usage: python VCF_to_VariantTable_with_ParentalGT_v3.py <vcf file> <ped file>\n")
	sys.stderr.write("Usage: python VCF_to_VariantTable_with_ParentalGT_v3.py exome_calls_pass_step2_normalized_anno.hg19_multianno_0.001.vcf 29Trios.ped \n")
	sys.exit(-1)

cwd = os.getcwd()

######################################-Function-##########################################
# Test Wheather the input is gzip file, True=gzip file, False=un gzip file
def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

# Convert DP or GQ into 0 if '.' is found
def convertQual(DP):
	if DP == '.':
		return 0
	else:
		return float(DP)

# Get GT information
def getGT(sample,data,store):
    if sample in store:
        sample_info = data[store.index(sample)]
        sample_GT = sample_info.split(':')[0]
    else:
        sample_GT = "NA"
    return sample_GT

# Deal with different length of information stored in format (ex: GT:RGQ:DP:GQ:PID:PS:SB:PGT:AD:PL)
def makeinfo(sample,data,store,format):
    if sample in store:
        sample_info = data[store.index(sample)]
        sample_GT = sample_info.split(':')[0]
        if len(sample_info.split(':')) < len(format):
            sample_AD = 'NA'
            sample_DP = 'NA'
            sample_GQ = 'NA'
        else:
            if 'AD' in format:
                sample_AD = sample_info.split(':')[format.index('AD')]
            else:
                sample_AD = 'NA'

            if 'DP' in format:
                sample_DP = sample_info.split(':')[format.index('DP')]
            else:
                sample_DP = 'NA'

            if 'GQ' in format:
                sample_GQ = sample_info.split(':')[format.index('GQ')]
            else:
                sample_GQ = 'NA'
    else:
        sample_GT = 'NA'
        sample_AD = 'NA'
        sample_DP = 'NA'
        sample_GQ = 'NA'

    # if sample_DP == '0' or sample_DP == '.' or sample_DP == 'NA' or sample_AD == 'NA':
    #     sample_VAF = 'NA'
    # else:
    #     sample_VAF = str(float(sample_AD.split(',')[1])/float(sample_DP))
    try:
        if (
            sample_AD != 'NA' and ',' in sample_AD and
            sample_DP not in ['0', '.', 'NA']
        ):
            ref, alt = sample_AD.split(',')
            sample_VAF = str(float(alt) / float(sample_DP))
        else:
            sample_VAF = 'NA'
    except (IndexError, ValueError, ZeroDivisionError):
        sample_VAF = 'NA'
    
    out.write('\t'.join([sample,sample_GT,sample_AD,sample_VAF,sample_DP,sample_GQ])+'\t')
######################################-Function-##########################################

###################################-Prepare steps-########################################
# get the pedigree list and generate a full list of sample name that is not singeleton
ped_list = []
fam_list = []

pedfile = sys.argv[2]
with open(pedfile,'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            else:
                data = line.strip().split('\t')
                ped_list.append(data[1:])
                fam_list.extend(data[1:])
fam_list=sorted(set(fam_list))[1:]

# Calculate the maxium number of members in a family
for family in ped_list:
    num=len(family)

# Test if the input vcf file is gzip or not, if so it will ungzip it
vcf_f = sys.argv[1]
vcf_nf = os.path.splitext(vcf_f)[0] #-----file name without gz
# Unzip the gz file
if is_gz_file(vcf_f):
	with gzip.open(vcf_f, 'rb') as f_in:
    		with open(vcf_nf, 'wb') as f_out:
        		shutil.copyfileobj(f_in, f_out)
	vcf_file = vcf_nf
else:
	vcf_file = vcf_f	

# get the full sample list from the VCF
with open(vcf_file, 'r') as file:
	for line in file:
		if line[0] == '#':
			if '#CHROM' in line:
				data = line.strip().split('\t')
				list = data[9:]	

# get the singleton list by excluding samples from ped (trios)
singleton_list = []
for ID in list:
        if ID not in fam_list:
                singleton_list.append(ID)

# output file name
out = open(vcf_file[:-4]+'_variantTable.txt','w')
###################################-Prepare steps-########################################

header = []
header_flag = 'F'
order = []

with open(vcf_file, 'r') as file:
    for line in file:
        data = line.strip().split('\t')
        if line.startswith('#CHROM'):
            store = data   # store saves each header
        elif line.startswith('##'):
            continue
        else:
            #chr = data[0]
            #pos = data[1] 
            info = data[7].split(';')  # info stores the info of each line
            Info = {}  # Info is a dictionary; key is the metric; value is the variant info
            for item in info:
                if len(item.split('=')) == 2:
                    Info[item.split('=')[0]] = item.split('=')[1]
                    order.append(item.split('=')[0])

            # print out header
            if header_flag == 'F':
                for key in order:
                    header.append(key)
                
            # Write header into txt file
                n=1
                for i in range(num):
                    if i == 0:
                        out.write('ProbandID\tProband_GT\tProband_AD\tProband_VAF\tProband_DP\tProband_GQ\t')
                    elif i == 1:
                        out.write('MotherID\tMother_GT\tMother_AD\tMother_VAF\tMother_DP\tMother_GQ\t')
                    elif i == 2: 
                        out.write('FatherID\tFather_GT\tFather_AD\tFather_VAF\tFather_DP\tFather_GQ\t')
                    else:
                        Sibling = 'Sibling' + str(n)
                        out.write(Sibling+'_ID'+'\t'+Sibling+'_GT'+'\t'+Sibling+'_AD'+'\t'+Sibling+'_VAF'+'\t'+Sibling+'_DP'+'\t'+Sibling+'_GQ'+'\t')
                        n+=1
                
                out.write('\t'.join(store[0:7])+'\t'+'\t'.join(header)+'\n')
                
                header_flag = 'T'
                
            format = data[store.index('FORMAT')].split(':')
            for family in ped_list:
                for nm in range(len(family)):
                    if nm == 0:
                        proband = family[nm]
                        if proband != ".":
                            proband_GT = getGT(proband,data,store)
                        else:
                            proband_GT = "NA"
                    elif nm == 1:
                        mother = family[nm]
                        if mother != ".":
                            mother_GT = getGT(mother,data,store)
                        else:
                            mother_GT = "NA"
                    elif nm == 2:
                        father = family[nm]
                        if father != ".":
                            father_GT = getGT(father,data,store)
                        else:
                            father_GT = "NA"
                    
                if proband not in store:
                    continue
                
                if proband_GT == "NA":
                    continue
                else:
                    if (mother_GT != "NA" and father_GT != "NA"):
                        if (proband_GT == '0/0' or proband_GT == '0|0' or proband_GT == './.') and (father_GT == '0/0' or father_GT == '0|0' or father_GT == './.') and (mother_GT == '0/0' or mother_GT == '0|0' or mother_GT == './.'):
                            continue
                        else:
                            makeinfo(proband,data,store,format)
                            makeinfo(mother,data,store,format)
                            makeinfo(father,data,store,format)  
        
                    elif (mother_GT == "NA" and father_GT == "NA"):
                        if (proband_GT == '0/0' or proband_GT == '0|0' or proband_GT == './.'):
                            continue
                        else:
                            makeinfo(proband,data,store,format)

                            mother_info=['NA' for i in range(6)]
                            out.write('\t'.join(mother_info)+'\t')

                            father_info=['NA' for i in range(6)]
                            out.write('\t'.join(father_info)+'\t')

                    elif (mother_GT == "NA" and father_GT != "NA"):
                        if (proband_GT == '0/0' or proband_GT == '0|0' or proband_GT == './.') and (father_GT == '0/0' or father_GT == '0|0' or father_GT == './.'):
                            continue
                        else:
                            makeinfo(proband,data,store,format)

                            mother_info=['NA' for i in range(6)]
                            out.write('\t'.join(mother_info)+'\t')

                            makeinfo(father,data,store,format)
                    
                    elif (mother_GT != "NA" and father_GT == "NA"):
                        if (proband_GT == '0/0' or proband_GT == '0|0' or proband_GT == './.') and (mother_GT == '0/0' or mother_GT == '0|0' or mother_GT == './.'):
                            continue
                        else:
                            makeinfo(proband,data,store,format)
                            makeinfo(mother,data,store,format)

                            father_info=['NA' for i in range(6)]
                            out.write('\t'.join(father_info)+'\t')                          
               
                for i in range(len(family)):
                    if i > 2:
                        sibling = family[i]
                        if sibling == ".":
                            sibling_info=['NA' for i in range(6)]
                            out.write('\t'.join(sibling_info)+'\t')
                        else:
                            makeinfo(sibling,data,store,format)
               
                #print out line contents
                out.write('\t'.join(data[0:7]))
                
                for i in header:
                    if i not in Info:
                        out.write('\t'+'NA')
                    else:
                        out.write('\t'+Info[i])
                out.write('\n')
                
            # for each singleton proband
            #for proband in singleton_list:
            #    proband_GT = getGT(proband,data,store)

            #    if (proband_GT == '0/0' or proband_GT == '0|0' or proband_GT == './.'):
            #        continue
                
                #makeinfo(proband,data,store,format)
                
                # Add NA info for parents and siblings
            #    for i in range(num):
            #        if i > 0:
            #            na_info=['NA' for i in range(6)]
                        #out.write('\t'.join(na_info)+'\t')

                #out.write('\t'.join(data[0:7]))
            #    for i in header:
            #        if i not in Info:
                        #out.write('\t'+'NA')
            #        else:
                        #out.write('\t'+Info[i])
                #out.write('\n')


out.close()

# Remove ungzip file
if is_gz_file(vcf_f):
    os.remove(vcf_nf)
