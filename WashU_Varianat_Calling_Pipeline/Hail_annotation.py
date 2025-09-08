#!/usr/bin/python
###################################
# Purpose: Hail annotation steps for the variant Tables (txt format)
# Author: Yung-Chun Wang <yung-chun@wustl.edu>
# Modifier: 
# Language: Python
# Version: 2
# Comment: This script for running annotation for variant Table.txt
# Created Date: 04-06-2023
# Last Modified Date: 09-03-2025
# Updated: for MAP_CP
###################################

###########
# MODULES #
###########
########################################
import os
import sys
import subprocess
import time
from datetime import datetime
#import pyspark
import argparse
import hail as hl

## Add my personal module
sys.path.append('/home/yung-chun/2_python_scripts/my_hail_module')
import my_hail_function_general as mygl
import my_hail_function_variantTable_annotation as myva
import my_hail_table_list as myhl
########################################

###########
# Date #
###########
########################################
today = datetime.now().strftime("%Y-%m-%d")
########################################

###########
# Jobname #
###########
########################################
job_name = "annotation"
########################################

#############
# ARGUMENTS #
#############
########################################
### ARGV begin:
parser = argparse.ArgumentParser(description='Annotation for the variant Table, and then save to hail table and txt.', prog="HAIL_variantTable_annotation")

### contents
parser.add_argument(
    "-i", "--input_file"
    , dest='file_fullpath'
    , required=True
    , type=mygl.hail_file_exists
    , help="Variant table path:\n" 
            "This program take input file in txt format.\n"
            "If for 'step2' and 'step3', the input should be hail table.")
parser.add_argument(
    "-o", "--output"
    , dest='output_fullpath'
    , required=False
    , default=None
    , type=mygl.is_directory_check
    , help="Output path: Provide the absolute path for folder to store the output file.")
parser.add_argument(
    "-m", "--mode"
    , dest='processing_mode'
    , required=True
    , default='all'
    , nargs='+'
    , choices=['all',
               'step1',
               'step2',
               'step3',
               'cadd_null',
               'plotread',
               'blat_anno',
               'deepvariant_anno',
               'recurrence']
    , help="Choose the mode of annotation steps:\n"
            "- 'all' performs all three annotation steps (VEP, damage_protein, and variant_types).\n"
            "- 'step1' performs the first annotation step (VEP).\n"
            "- 'step2' performs the second annotation step (damage_protein) and requires a CADD web file to be provided.\n"
            "- 'step3' performs the third annotation step (variant_types) and requires a phenotype file to be provided.\n"
            "- 'cadd_null' performs the cadd_null out tsv for uploading to cadd web, make sure the input file is hail_table after step1.\n"
            "- 'plotread' generates the plotread format txt file (ID, ProbandID, #CHROM, POS, Proband_GT).")
			
## Select the varaints only show "Pass" in either my_denovo, deepvaraint, blat
parser.add_argument(
    "-s", "--select"
    , dest='select_mode'
    , required=False
    , default=None
    , nargs='+'
    , choices=['denovo',
               'blat',
               'deepvariant']
    , help="Choose the mode of annotation steps:\n"
            "- 'denovo' only keep the denovo positive variants.\n"
            "- 'blat' only keep vairants pass blat.\n"
            "- 'deepvariant' only keep vairants pass deepvariant.")
parser.add_argument(
    "-c", "--cadd_web_file"
    , dest='cadd_web_file'
    , required=False
    , help="Provide the web cadd file for processing damage protein characterization (step2).")
parser.add_argument(
    "-p", "--phenotype_file"
    , dest='phenotype_file'
    , required=False
    , help="Provide the phenotype file for cleaning (step1), and calling the variants types (step3).")
## denovo_file
parser.add_argument(
    "-d", "--denovo_file"
    , dest='denovo_file'
    , required=False
    , help="Provide the denovo calling file from hl.denovo (step3).")
## upd_file
parser.add_argument(
    "-upd", "--upd_file"
    , dest='upd_file'
    , required=False
    , help="Provide the upd calling file (step3).")
## unique_FamID
parser.add_argument(
    "-uf", "--unique_FamID"
    , dest='unique_FamID'
    , required=False
    , help="Provide the unique_FamID file for filtering the unique family ID in joint variant table (step1).")
## exclude_samples
parser.add_argument(
    "-es", "--exclude_samples"
    , dest='exclude_samples'
    , required=False
    , help="Provide the exclude_samples file for exluding the sample from joint variant table (step1).")
## true_set_file (positive control)
parser.add_argument(
    "-ts", "--true_set_file"
    , dest='true_set_file'
    , required=False
    , help="Provide the true_set_file for annotating the positive vairaints (step3).")
## project_key to find the correspnding files
parser.add_argument(
    "-k", "--project_key"
    , dest='project_key'
    , required=False
    , help="Provide the project_key to search for the phenotype, cadd_web, and denovo file in hail_list.")
## deepvariant parsing file
parser.add_argument(
    "-dp", "--deepvariant_parse_file"
    , dest='deep_file'
    , required=False
    , help="Provide the deep_file for annotating the positive vairaints (step4).")
## blat parsing file
parser.add_argument(
    "-bl", "--blat_parse_file"
    , dest='blat_file'
    , required=False
    , help="Provide the blat_file for annotating the positive vairaints (step4).")
## blat manually correction file
parser.add_argument(
    "-blc", "--blat_correct_file"
    , dest='blat_correction_file'
    , required=False
    , help="Provide the blat_correction_file for annotating (step4).")
## Provide the columns name for the output
parser.add_argument(
    "-col", "--columns"
    , dest='select_columns'
    , required=False
    , default=None
    , nargs="+"
    , help="Specific columns to select for the output.")
## Version
parser.add_argument(
    "-v", "--version"
    , version='%(prog)s 1.0'
    , action='version'
    , help="Print out program version.")

args = parser.parse_args()

##################
# CHECKING STEPS #
###################
########################################
## Mode selection checking:
mode = set(args.processing_mode)
print(f"Your mode selection is {mode}")

valid_modes = [
    {'all'},
    {'step1'},
    {'step2'},
    {'step3'},
    {'step1', 'step2'},
    {'step2', 'step3'},
    {'step1', 'step2', 'step3'},
    {'cadd_null'},
    {'plotread'},
    {'blat_anno'},
    {'deepvariant_anno'},
    {'deepvariant_anno', 'plotread'},
    {'blat_anno', 'plotread'},
    {'blat_anno', 'deepvariant_anno'},
    {'blat_anno', 'deepvariant_anno', 'plotread'},
    {'recurrence'}
]
valid_mode_options = [sorted(option) for option in valid_modes]

if not any(mode == set(m) for m in valid_modes):
    print(f"Error: Invalid mode provided by arugments: {mode}.\n"
          f"Please choose from {valid_mode_options}")
    exit(1)

# Check if args.select_mode is not None before processing
if args.select_mode:
    select_mode = set(args.select_mode)
    print(f"Your select mode selection is {select_mode}")

    valid_select_modes = [
        {'denovo'},
        {'blat'},
        {'deepvariant'},
        {'denovo', 'blat'},
        {'denovo', 'deepvariant'},
        {'blat', 'deepvariant'},
        {'denovo', 'blat', 'deepvariant'}
    ]
    valid_select_mode_options = [sorted(option) for option in valid_select_modes]

    if not any(select_mode == set(m) for m in valid_select_modes):
        print(f"Error: Invalid select mode provided by arguments: {select_mode}.\n"
              f"Please choose from {valid_select_mode_options}")
        exit(1)
else:
    print("No select mode provided. Use default or specific behavior based on the rest of your program.")

## Check the file requirement for each step
project_name = args.project_key

cadd_web_file = None
phenotype_file = None
denovo_file = None
upd_file = None

#######
# WES #
#######
if 'all' in mode:
    # Check if cadd_web_file and phenotype_file was provided in module or by arguments
    if not project_name and (not args.cadd_web_file or not args.phenotype_file):
        print('Error: cadd_web (-c), phenotype (-p) and hail denovo (-d) files must be provided by arguments or specify the project ' 
              'name that can be used to search in hail_table_list for the completed annotation.')
        exit(1)
    elif not project_name:
        cadd_web_file = args.cadd_web_file
        phenotype_file = args.phenotype_file
        # denovo_file = args.denovo_file
    else:
        cadd_web_file = myhl.project_annotaion_related_file_paths[project_name]['cadd_web_file']
        phenotype_file = myhl.project_annotaion_related_file_paths[project_name]['phenotype_file']
        # denovo_file = myhl.project_annotaion_related_file_paths[project_name]['denovo_file']

if 'step1' in mode:
    # Check if cadd_web_file was provided in module or by arguments
    if not project_name and not args.phenotype_file:
        print('Error: phenotype (-p) file must be provided by arguments or specify the project ' 
              'name that can be used to search in hail_table_list for processing step1.')
        exit(1)
    elif not project_name:
        phenotype_file = args.phenotype_file
    else:
        phenotype_file = myhl.project_annotaion_related_file_paths[project_name]['phenotype_file']

if 'step2' in mode:
    # Check if cadd_web_file was provided in module or by arguments
    if not project_name and not args.cadd_web_file:
        print('Error: cadd_web (-c) file must be provided by arguments or specify the project ' 
              'name that can be used to search in hail_table_list for processing step2.')
        exit(1)
    elif not project_name:
        cadd_web_file = args.cadd_web_file
    else:
        cadd_web_file = myhl.project_annotaion_related_file_paths[project_name]['cadd_web_file']
 
# if 'step3' in mode:
#     # Check if denovo_file was provided in module or by arguments
#     if not project_name and not args.denovo_file:
#         print('Error: hail denovo (-d) file must be provided by arguments or specify the project ' 
#               'name that can be used to search in hail_table_list for processing step3.')
#         exit(1)
#     elif not project_name:
#         denovo_file = args.denovo_file
#     else:
#         denovo_file = myhl.project_annotaion_related_file_paths[project_name]['denovo_file']

#################
# Deep and Blat #
#################
if 'deepvariant_anno' in mode:
    # Check if cadd_web_file was provided in module or by arguments
    if not project_name and not args.deep_file:
        print('Error: deep parsing (-dp) file must be provided by arguments or specify the project ' 
              'name that can be used to search in hail_table_list for processing step4.')
        exit(1)
    elif not project_name:
        deep_file = args.deep_file
    else:
        deep_file = myhl.project_annotaion_related_file_paths[project_name]['deep_file']
 
if 'blat_anno' in mode:
    # Check if cadd_web_file was provided in module or by arguments
    if not project_name and not args.blat_file:
        print('Error: blat parsing (-bl) file must be provided by arguments or specify the project ' 
              'name that can be used to search in hail_table_list for processing step4.')
        exit(1)
    elif not project_name:
        blat_file = args.blat_file
    else:
        blat_file = myhl.project_annotaion_related_file_paths[project_name]['blat_file']

########################################

###############
# INPUT/OUPUT #
###############
########################################
## Input
input_file = args.file_fullpath
unique_FamID = args.unique_FamID
exclude_samples = args.exclude_samples
true_set_file = args.true_set_file
select_columns = args.select_columns


blat_correction_file = args.blat_correction_file

## Change mode for all to each steps
if mode == {'all'}:
    mode = {'step1', 'step2', 'step3'}

## Final Output
if len(mode) > 1:
    mode_name = '_'.join(sorted(mode))
else:
    mode_name = ''.join(mode)
name_pattern = "annotation_" + mode_name

output_path = args.output_fullpath
final_output_file = mygl.create_output_name(input_file, name_pattern, '.ht', output_path)

########################################

##################
# PRINT SETTINGS #
##################
########################################
print('Your input file is: ', input_file)
print('Cadd: ', cadd_web_file)
print('Phenotype: ', phenotype_file)
print('Denovo: ', denovo_file)
# print('unique_FamID: ', unique_FamID)
# print('exclude_samples: ', exclude_samples)
# print('true_set: ', true_set_file)
print('Your final output stored here: ',final_output_file)
########################################

###################
# HAIL INITIATION #
###################
########################################
job_name = job_name + "_" + mode_name
mygl.hail_initialize(job_name)
########################################

#######
# WES #
#######
## Step1: input should be variantable .txt format
if 'step1' in mode:
    name_pattern = "annotation_step1"
    output_step1 = mygl.create_output_name(input_file, name_pattern, '.ht', output_path)
    output_file = output_step1
    print("Output for step1: ", output_step1)
    
    if unique_FamID is not None and exclude_samples is not None: 
        myva.step1_annotation(input_file, output_file, phenotype_file, unique_FamID, exclude_samples)
    else:
        myva.step1_annotation(input_file, output_file, phenotype_file)

if 'step2' in mode:
    if 'step1' in mode:
        input_file = output_step1
        name_pattern = "step2"
    else:
        name_pattern = "annotation_step2"

    output_step2 = mygl.create_output_name(input_file, name_pattern, '.ht', output_path)
    output_file = output_step2
    print("Output for step2: ", output_step2)
    
    myva.step2_annotation(input_file, output_file, cadd_web_file)

if 'step3' in mode:
    if 'step2' in mode:
        input_file = output_step2
        name_pattern = "step3"
    else:
        name_pattern = "annotation_step3"

    output_step3 = mygl.create_output_name(input_file, name_pattern, '.ht', output_path)
    output_file = output_step3
    print("Output for step3: ", output_step3)

    myva.step3_annotation(input_file, output_file, denovo_file, upd_file, true_set_file)

###################################
# BLAT and DEEPVARIANT annotation #
###################################
if 'deepvariant_anno' in mode:
    print("Input file for deepvariant annotation is: ", input_file)

    name_pattern = "anno_deepvariant"
    
    output_deep = mygl.create_output_name(input_file, name_pattern, '.ht', output_path)
    output_file = output_deep
    print("Output for deepvariant annotation: ", output_deep)
    
    if not mygl.file_exists(output_file):
        myva.step4_annotation_deep_blat(input_file, output_file, deep_file=deep_file)
    else:
        print(f"File {output_file} already exists. Skipping step.")
    
    ## reset input, output for the next step
    input_file, output_file = mygl.process_input_output(input_file=input_file, output_file=output_file)

if 'blat_anno' in mode:
    print("Input file for blat annotation is: ", input_file)

    if not blat_correction_file:
        name_pattern = "anno_blat"
    else:
        name_pattern = "anno_blat_correct"

    output_blat = mygl.create_output_name(input_file, name_pattern, '.ht', output_path)
    output_file = output_blat
    print("Output for blat annotation: ", output_file)
    
    if not mygl.file_exists(output_file):
        if not blat_correction_file:
            myva.step4_annotation_deep_blat(input_file, output_file, blat_file=blat_file)
        else:
            myva.step4_annotation_deep_blat(input_file, output_file, blat_file=blat_file, blat_correction_file=blat_correction_file)
    else:
        print(f"File {output_file} already exists. Skipping step.")
    
    ## reset input, output for the next step
    input_file, output_file = mygl.process_input_output(input_file=input_file, output_file=output_file)

########
# ELSE #
########
if 'cadd_null' in mode:
    ht = hl.read_table(input_file)
    myva.cadd_null_out_tsv(ht, input_file)

if 'plotread' in mode:
    ## If combine with other steps, like blat,deepvariant annotation, it will use the output file as input_file
    print("Input file for plotread file preparation step is: ", input_file)
    
    ht = hl.read_table(input_file)

    if select_columns is None:
        select_all=True
    else:
        select_all=False

    if 'denovo' in select_mode:
        ht = ht.filter(hl.is_missing(ht.my_denovo), keep=False)
    
    if 'blat' in select_mode:
        ht = ht.filter((ht.Blat == "Pass"), keep=True)

    if 'deepvariant' in select_mode:
        ht = ht.filter((ht.Deepvariant == "Pass"), keep=True)

    output_file =  mygl.create_output_name(input_file, pattern="forplotread",extension=".txt")
    
    if not mygl.file_exists(output_file):
        mygl.out_for_downstream(ht, out_file=output_file, select_all=select_all, cols=select_columns)
    else:
        print(f"File {output_file} already exists. Skipping step.")

    ## reset input, output for the next step
    input_file, output_file = mygl.process_input_output(input_file=input_file, output_file=output_file)

if 'recurrence' in mode:
    name_pattern = "recurrence"
    output_file = mygl.create_output_name(input_file, name_pattern, '.ht', output_path)

    # Add recurrent information
    ht = recur_byID(ht)
    ht = recur_bygene(ht)
    ht = recur_proteindmg_bygene(ht)

    ## Save output file
    mygl.out_save(ht, out_file=output_file)

    ## reset input, output for the next step
    input_file, output_file = mygl.process_input_output(input_file=input_file, output_file=output_file)