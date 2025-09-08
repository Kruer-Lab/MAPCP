#!/usr/bin/python
###################################
# Purpose: Functions for the Hail annotation steps for WES_CONSANG
# Author: Yung-Chun Wang <yung-chun@wustl.edu>
# Modifier: 
# Language: Python
# Version: 1
# Comment: This funciton only process hail table read from vairant Table
# Created Date: 03-22-2023
# Last Modified Date: 04-04-2023
###################################

import sys
import os
import csv
import subprocess
import re
import tempfile
import pandas as pd

import hail as hl
from hail.typecheck import typecheck, nullable, typecheck_method
#from hail.typecheck import (typecheck, nullable, anytype, enumeration, tupleof,func_spec, oneof, arg_check, args_check, anyfunc,sequenceof)


## Add my personal module
sys.path.append('/home/yung-chun/2_python_scripts/my_hail_module')
import my_hail_table_list as myhl
import my_hail_function_general as mygl

#################
# Genotype List #
#################
non_lst = ['./.','.|.']
non_lst_hl = hl.literal(non_lst)

ref_lst = ['0/0','0|0']
ref_lst_hl = hl.literal(ref_lst)

het_lst = ['1/0','0/1','1|0','0|1']
het_lst_hl = hl.literal(het_lst)

hom_lst = ['1/1','1|1']
hom_lst_hl = hl.literal(hom_lst)

#############
# FUNCTIONS #
#############
#######################################
## Clean up the variant table
def remove_malformed_rows(input_path):
    """
    Reads a TSV file and removes rows that do not match the header's field count.
    Returns the path to a cleaned temporary file.
    """
    with open(input_path, "r") as f:
        lines = f.readlines()

    header = lines[0].strip().split("\t")
    expected_n_fields = len(header)

    # Create a temp file to hold cleaned content
    tmpfile = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.tsv')
    
    # Write header
    tmpfile.write("\t".join(header) + "\n")
    
    # Write only valid rows
    for line in lines[1:]:
        if len(line.strip().split("\t")) == expected_n_fields:
            tmpfile.write(line)

    tmpfile.close()
    return tmpfile.name

## Read txt file and key by choosen columnID
@typecheck(file=str,
           keys=str)
def read_txt2ht_simple(file, keys):
    """
    Import the txt file and then key by the choosen columnID

    :param file: /path/to/variant/txt_file
    :return: hail table.
    """
    _ht = hl.import_table(file,delimiter='\t').key_by(keys)

    return _ht

## Import variant Table to hail table
@typecheck(file=str)
def read_variantTable(file):
    """
    Import the variant Tabel into hail table with setting the Proband_AD to str and POS to int

    :param file: /path/to/variant/Table
    :return: hail table.
    """
    rg38 = hl.get_reference('GRCh38')
    
    types={
            "Proband_AD":hl.tstr, 'POS':hl.tint, 'END':hl.tint
          }
    _ht = hl.import_table(file,delimiter='\t',types=types)

    _ht = _ht.annotate(
                    locus = hl.locus(_ht['#CHROM'], _ht.POS, reference_genome=rg38),
                    alleles=hl.array([_ht.REF, _ht.ALT])
    )
    
    _ht = _ht.key_by('locus','alleles')

    return _ht

def read_variant_table_step_blat_deep(file):
    if isinstance(file, str) and mygl.file_exists(file) and file.endswith('.ht'):
        _ht = hl.read_table(file).key_by('ProbandID','locus')
    elif isinstance(file, str) and mygl.file_exists(file) and file.endswith('.txt'):
        _ht = read_ht_from_txt(file)
    
    return _ht

## Key the hail table read from variant Table by locus and alleles
@typecheck(ht=hl.Table)
def key_by_alleles_locus_ht(ht):
    """
    Modify the hail table read from the variant Table by generating the locus and allels, and then key by [locus, alleles].

    :param ht: hail table
    :return: keyed hail table
    """
    rg38 = hl.get_reference('GRCh38')

    _ht = ht.annotate(
                locus = hl.locus(ht['#CHROM'], ht.POS, reference_genome=rg38),
                alleles = hl.array([ht.REF, ht.ALT]),
                Proband_DP_int = hl.if_else(
                                            (ht.Proband_AD != ".") ,
                                            hl.array([hl.int(ht.Proband_AD.split(",")[0]), hl.int(ht.Proband_AD.split(",")[1])]),
                                            hl.array([hl.missing(hl.tint),hl.missing(hl.tint)])
                                            )
                )

    _ht = _ht.key_by('locus','alleles')
    
    return _ht

## Transform the bed/txt annotation file into hail table
@typecheck(input_file=str,
           output_path=str,
           key_by=str,
           columns=str,
           values=str)
def transform_annotation_bed_ht(input_file: str,
                                output_path: str = '/storage1/fs1/jin810/Active/References/2023_old_references/hail/Ref_HailFormat',
                                key_by: str = "interval",
                                columns: str = 'info', 
                                values: str = 'in_bed'):
    """
    Transform the bed or txt annotation file to the hail table. The file should either followed these two settings. 1: 4 column- chr, start, end, information or 2. 5 column- index, chr, start, end, information.

    :param input_file: tab-delimited files
    :param output_path: path to save hail table
    :key_by: interval or locus
    :param columns: Comma-separated string of column names.
    :param values: Comma-separated string of values corresponding to each column.
    :return: hail table
    """
    
    #############
    # FUNCTIONS #
    #############
    
    ## Chekc if it is gzip
    def is_gz_file(filepath):
        with open(filepath, 'rb') as test_f:
            return test_f.read(2) == b'\x1f\x8b'

    ## Read the file without line start with "#"
    def filter_file_content(file_path):
        with open(file_path, 'r') as f:
            content = '\n'.join([line.strip() for line in f if not line.startswith('#')])
        
        return content	

    ## Get column number function
    def get_col_num(file):
        with open(file) as f:
            # Read the first line
            first_line = f.readline().strip()

            # Split the line based on one or more spaces or a tab
            first_row = re.split(r'\s+|\t', first_line)
            num_cols = len(first_row)

        return num_cols

    ## Determine the types and colnames function
    def get_type_header(file):
        types = {}  # Default empty dictionary
        header_name = {}  # Default empty dictionary

        num_cols = get_col_num(file)

        if num_cols == 5:
            types = {'f0':hl.tint, 'f1':hl.tstr, 'f2':hl.tint, 'f3':hl.tint, 'f4':hl.tstr }
            header_name = {'f0':'index', 'f1':'contig', 'f2':'start', 'f3':'end', 'f4':'info'}
        elif num_cols == 4:
            types = {'f0':hl.tstr, 'f1':hl.tint, 'f2':hl.tint, 'f3':hl.tstr }
            header_name = {'f0':'contig', 'f1':'start', 'f2':'end', 'f3':'info'}
        elif num_cols == 3:
            types = {'f0':hl.tstr, 'f1':hl.tint, 'f2':hl.tint}
            header_name = {'f0':'contig', 'f1':'start', 'f2':'end'}
        else:
            raise ValueError(f"Unexpected number of columns {num_cols}. Expected 3, 4, or 5 columns.")

        return types, header_name

    ## Add column for those bed file only have chr star end
    def add_columns_to_table(ht: hl.Table, columns: str, values: str) -> hl.Table:
        """
        Add new columns with specified names and values to a Hail table.

        :param ht: The Hail table to modify.
        :param columns: Comma-separated string of column names.
        :param values: Comma-separated string of values corresponding to each column.
        :return: The modified Hail table.
        """
        column_list = columns.split(',')
        value_list = values.split(',')

        if len(column_list) != len(value_list):
            raise ValueError("Number of columns must match number of values.")

        for column, value in zip(column_list, value_list):
            ht = ht.annotate(**{column: hl.literal(value)})

        return ht

    #########
    # START #
    #########

    ## Tracking the decompress
    decompressed = False

    ## Check if the input file is gzipped either by extension or by content
    if is_gz_file(input_file):

        print("Input file is gzipped. Decompressing to standard output...")

        ### Decompress the gzip file to a temporary file
        temp_file_decomp = input_file + ".decompressed"

        with open(temp_file_decomp, 'w') as f_out:
            subprocess.run(['gzip', '-dc', input_file], stdout=f_out)

        input_file = temp_file_decomp

        ## Change decompress status flag
        decompressed = True

    ## Remove the linw with "#"
    content = filter_file_content(input_file)
    temp_file_removehash = input_file + ".remove_hash"

    with open(temp_file_removehash, 'w') as f_out:
        f_out.write(content)

    input_file = temp_file_removehash

    ## Read reference
    rg38 = hl.get_reference('GRCh38')
    
    ## get output name
    name = os.path.basename(input_file)
    out_hl = mygl.create_output_name(name, pattern="hail", extension=".ht", output_path=output_path)

    ## Import file
    types, header_name = get_type_header(input_file)
    ht = hl.import_table(input_file, delimiter="\t", types=types, no_header=True, force_bgz=True).rename(header_name)

    ## To handle the column number is 3, add info column
    if get_col_num(input_file) == 3:
        ht = add_columns_to_table(ht, columns, values)
    
    ## Add essential information to the hail table
    ht = ht.annotate(pos = hl.int((ht.start + ht.end)/2))
    
    ht = ht.annotate(
                     locus = hl.locus(ht.contig, ht.pos, reference_genome=rg38),
                     interval = hl.interval(
                                           hl.locus(ht.contig, ht.start, reference_genome=rg38),
                                           hl.locus(ht.contig, ht.end, reference_genome=rg38)
                                           )
                  )
    ht = ht.key_by(key_by)
    
    ht = ht.select('contig', 'start', 'end', 'pos', 'locus', 'info')
    
    ## Save output
    mygl.out_save(ht, out_hl)

    ## Remove the temp file
    if decompressed:
        os.remove(temp_file_decomp)

    os.remove(temp_file_removehash)

## Read file into ht and then rekey to alleles and locus
@typecheck(file=str)
def read_variantTable_ht_rekey(file):
    """
    Import the variant Tabel into hail table. It will set the following info to a specific data types. And then generate the locus and alleles for setting as key.
    - POS: int
    - Proband/Mother/Father_DP: int
    - Proband/Mother/Father_VAF: float
    - Proband_GQ: int
    
    :param file: /path/to/variant/Table
    :return: hail table
    """
    rg38 = hl.get_reference('GRCh38')

    # Clean input file first
    cleaned_file = remove_malformed_rows(file)

    types = {
    "cadd": hl.tfloat,
    "vep_consequence": hl.tarray(hl.tstr),
    "gnomad_genomes_AF": hl.tfloat,
    "bravo8_AF": hl.tfloat,
    "Proband_AD": hl.tstr,
    "POS": hl.tint,
    "caddv16_phred": hl.tfloat,
    "Mother_AD": hl.tstr,
    "Father_AD": hl.tstr,
    "MQ": hl.tfloat,
    "Proband_DP_int": hl.tarray(hl.tint),
    "Mother_DP_int": hl.tarray(hl.tint),
    "Father_DP_int": hl.tarray(hl.tint),
    "Proband_VAF": hl.tfloat,
    "Mother_VAF": hl.tfloat,
    "Father_VAF": hl.tfloat
    }

    _ht = hl.import_table(cleaned_file,delimiter='\t',types=types)

    # Set checking step for proband_DP_int (if the table was from hail denvo after modify_hldenovo_for_outplotread, it should have *_DP_int)
    ## Get the schema of the table
    schema = _ht.row
    ## Check if the fields already exist
    has_Proband_DP_int = 'Proband_DP_int' in schema
    has_Mother_DP_int = 'Mother_DP_int' in schema
    has_Father_DP_int = 'Father_DP_int' in schema
    has_Proband_AD = 'Proband_AD' in schema
    has_Mother_AD = 'Mother_AD' in schema
    has_Father_AD = 'Father_AD' in schema
    has_Proband_DP = 'Proband_DP' in schema
    has_Mother_DP = 'Mother_DP' in schema
    has_Father_DP = 'Father_DP' in schema
    has_Proband_VAF = 'Proband_VAF' in schema
    has_Mother_VAF = 'Mother_VAF' in schema
    has_Father_VAF = 'Father_VAF' in schema

    # Remove rows where any AD field contains "."
    _ht = _ht.filter(
    ~(_ht.Proband_AD.contains(".") |
      _ht.Mother_AD.contains(".") |
      _ht.Father_AD.contains(".")),
    keep=True
    )
    
    # Proband_DP_int
    if has_Proband_AD and not has_Proband_DP_int:
        _ht = _ht.annotate(Proband_DP_int=hl.if_else(
            _ht.Proband_AD != ".",
            hl.array([hl.int(_ht.Proband_AD.split(",")[0]), hl.int(_ht.Proband_AD.split(",")[1])]),
            hl.array([hl.missing(hl.tint), hl.missing(hl.tint)])
    ))
    
    # Mother_DP_int
    if has_Mother_AD and not has_Mother_DP_int:
        _ht = _ht.annotate(Mother_DP_int=hl.if_else(
            _ht.Mother_AD != ".",
            hl.array([hl.int(_ht.Mother_AD.split(",")[0]), hl.int(_ht.Mother_AD.split(",")[1])]),
            hl.array([hl.missing(hl.tint), hl.missing(hl.tint)])
        ))

    # Father_DP_int
    if has_Father_AD and not has_Father_DP_int:
        _ht = _ht.annotate(Father_DP_int=hl.if_else(
            _ht.Father_AD != ".",
            hl.array([hl.int(_ht.Father_AD.split(",")[0]), hl.int(_ht.Father_AD.split(",")[1])]),
            hl.array([hl.missing(hl.tint), hl.missing(hl.tint)])
        ))
    
    # Proband_VAF
    if has_Proband_AD and has_Proband_DP and not has_Proband_VAF:
        _ht = _ht.annotate(Proband_VAF=hl.if_else(
            (_ht.Proband_AD != ".") & (_ht.Proband_DP != "."),
            hl.int(_ht.Proband_AD.split(",")[1]) / hl.int(_ht.Proband_DP),
            hl.missing(hl.tfloat)
        ))
    
    # Mother_VAF
    if has_Mother_AD and has_Mother_DP and not has_Mother_VAF:
        _ht = _ht.annotate(Mother_VAF=hl.if_else(
            (_ht.Mother_AD != ".") & (_ht.Mother_DP != "."),
            hl.int(_ht.Mother_AD.split(",")[1]) / hl.int(_ht.Mother_DP),
            hl.missing(hl.tfloat)
        ))
    
    # Father_VAF
    if has_Father_AD and has_Father_DP and not has_Father_VAF:
        _ht = _ht.annotate(Father_VAF=hl.if_else(
            (_ht.Father_AD != ".") & (_ht.Father_DP != "."),
            hl.int(_ht.Father_AD.split(",")[1]) / hl.int(_ht.Father_DP),
            hl.missing(hl.tfloat)
        ))
    
    # Deal with NaN value in the VAF and DP_int
    _ht = _ht.annotate(
    Proband_VAF=hl.if_else(hl.is_nan(_ht.Proband_VAF), hl.missing(hl.tfloat), _ht.Proband_VAF),
    Mother_VAF=hl.if_else(hl.is_nan(_ht.Mother_VAF), hl.missing(hl.tfloat), _ht.Mother_VAF),
    Father_VAF=hl.if_else(hl.is_nan(_ht.Father_VAF), hl.missing(hl.tfloat), _ht.Father_VAF),
    
    Proband_DP_int=hl.if_else(
        hl.any(lambda x: hl.is_nan(x), _ht.Proband_DP_int),
        hl.array([hl.missing(hl.tint), hl.missing(hl.tint)]),
        _ht.Proband_DP_int
    ),
    Mother_DP_int=hl.if_else(
        hl.any(lambda x: hl.is_nan(x), _ht.Mother_DP_int),
        hl.array([hl.missing(hl.tint), hl.missing(hl.tint)]),
        _ht.Mother_DP_int
    ),
    Father_DP_int=hl.if_else(
        hl.any(lambda x: hl.is_nan(x), _ht.Father_DP_int),
        hl.array([hl.missing(hl.tint), hl.missing(hl.tint)]),
        _ht.Father_DP_int
    )
    )


    # Perform the remaining annotations unconditionally
    _ht = _ht.annotate(
                    locus = hl.locus(_ht['#CHROM'], _ht.POS, reference_genome=rg38),
                    alleles = hl.array([_ht.REF, _ht.ALT]),
                    Proband_GQ = hl.if_else(
                                                (_ht.Proband_GQ != "."),
                                                hl.int(_ht.Proband_GQ),
                                                hl.missing(hl.tint)
                                            ),
                    Mother_GQ = hl.if_else(
                                                (_ht.Mother_GQ != "."),
                                                hl.int(_ht.Mother_GQ),
                                                hl.missing(hl.tint)
                                            ),
                    Father_GQ = hl.if_else(
                                                (_ht.Father_GQ != "."),
                                                hl.int(_ht.Father_GQ),
                                                hl.missing(hl.tint)
                                            )
                    )

    _ht = _ht.key_by('locus', 'alleles')

    return _ht

## Read upd_bed file into ht and then rekey to interval
@typecheck(file=str)
def read_ht_locus_from_triomix_upd_bed(file):
    """
    Read the UPD bed file (#CHROM,Start,End,ProbandID,Contamination,UPD_side(5MB)).

    :param file: /path/to/file
    :return: hail table
    """    
    rg38 = hl.get_reference('GRCh38')
    
    types={
           'Start':hl.tint, 'End':hl.tint    
          }

    _ht = hl.import_table(file,delimiter='\t',types=types)
  
    _ht = _ht.annotate(
                      interval = hl.interval(
                                       hl.locus(_ht['#CHROM'], _ht.Start, reference_genome=rg38),
                                       hl.locus(_ht['#CHROM'], _ht.End, reference_genome=rg38)
                                       )
                      )
    
    _ht = _ht.key_by('SampleID')

    return _ht

## To filter the hail table by providing FamID list, also exclude the samples from the list provided
@typecheck(ht=hl.Table,
           unique_FamID=str,
           exclude_samples=str)
def select_FamID(ht, unique_FamID, exclude_samples=None):
    """
    To filter the variant Table based on the FamID provided in the file. This will include all the samples that contains the FamID that is specified in the unique_FamID file. However, it will exclude the samples that was specified in the exclude_samples file.

    :param ht: hail table
    :param unique_FamID: Files containg the Family ID. ex: CP-F0001.
    :param exclude_samples: Files containing the Sample ID. ex: CP-F0001-003A.
    :return: hail table
    """
    # Step 1: Get the family list
    if isinstance(unique_FamID, list):
        Fam_lst = unique_FamID
    else:
        with open(unique_FamID, 'r') as file:
            Fam_lst = [line.strip() for line in file if not line.startswith('#')]
    
    # Step 2: Get the sample list based on the Family list
    # Step 2.1: Get unique sample list from variantTable
    UniqueSample_from_ht = ht.aggregate(hl.agg.collect_as_set(ht.ProbandID))

    # Step 2.2: Get sample list based on the Family list
    Sample_lst = [s for s in UniqueSample_from_ht if any(pattern in s for pattern in Fam_lst)]
    
    # Step2.3 Exclude samples if specified
    if exclude_samples is not None:
        if isinstance(exclude_samples, str):
            with open(exclude_samples, 'r') as file:
                exclude_lst = [line.strip() for line in file if not line.startswith('#')]
        else:
            exclude_lst = exclude_samples
        Sample_lst = [s for s in Sample_lst if s not in exclude_lst]
    
    # Step 2.4: Transform the sample list to a hail set
    Sample_set = hl.literal(set(Sample_lst))

    # Step 3: Subset the variantTable
    ht = ht.filter(Sample_set.contains(ht.ProbandID))

    return ht

## To remove some vairants that has no ALT and REF
### Specifal case for some variant that don't have Alt 
@typecheck(ht=hl.Table)
def remove_no_ref_alt_gt(ht):
    """
    Remove the variants that didn't have info for ALT and REF, and it was assinged as "*"
    Remove Proband_GT is ./. or .|.

    :param ht: hail table
    :return: hail table
    """
    ht = ht.filter( (ht.ALT == "*") | (ht.REF == "*") , keep=False)
    
    # Filter out if any GT column contains "."
    _ht = ht.filter(
        ht.Proband_GT.contains(".") |
        ht.Mother_GT.contains(".") |
        ht.Father_GT.contains("."),
        keep=False
    )
    
    return _ht

## Annotate the hail table with gene_name, bravo, gnomad, metasvm, cadd, feature typ.
#@mygl.set_ht_key_and_reset
@typecheck(ht=hl.Table)
def ht_annotate_gene_bravo_gnomad_meta_cadd_feature_type(ht):
    """
    Annotate the hail table with the following information:
    - SNP/INDEL
    - Gene name (gencode_gene_id_name_ht)
    - MAF from bravo_v8 (combined_reference_v2_ht)
    - MAF from gnomad_genome_v3 (combined_reference_v2_ht)
    - MetaSVM (combined_reference_v2_ht)
    - cadd (cadd_v16_ht)
    - clinvar (clinvar_ht)
    - feature_type (gencode_gene_id_name_ht)
    - low_complexity_region ()
    
    :param ht: hail table
    :return: hail table

    Notes: 
          the table for annotation should be named as same as the indicated inside ().
    """
    gencode_gene_id_name_ht = hl.read_table(myhl.gencode_gene_id_name_ht_path)
    combined_reference_v2_ht = hl.read_table(myhl.combined_reference_data_v2_ht_path)
    cadd_v16_ht = hl.read_table(myhl.cadd_v16_ht_path)
    clinvar_ht = hl.read_table(myhl.clinvar_ht_path)
    low_complex_ht = hl.read_table(myhl.low_complex_path)

    ht = ht.annotate(
                    SNP_INDEL = hl.if_else(
                                      hl.is_snp(ht.REF,ht.ALT) == True
                                      , "SNP"
                                      , "INDEL"
                                      ),
                    gene_name = hl.if_else(
                                        hl.is_defined(gencode_gene_id_name_ht[ht.locus].geneName),
                                        gencode_gene_id_name_ht[ht.locus].geneName,
                                        hl.missing(hl.tstr)
                                        )
                    )
    ht = ht.annotate(
                    bravo8_AF = hl.if_else(
                                        hl.is_defined(combined_reference_v2_ht[ht.key].topmed.AF),
                                        combined_reference_v2_ht[ht.key].topmed.AF,
                                        0
                                        ),
                    gnomad_genomes_AF = hl.if_else(
                                        hl.is_defined(combined_reference_v2_ht[ht.key].gnomad_genomes.AF),
                                        combined_reference_v2_ht[ht.key].gnomad_genomes.AF,
                                        0
                                        ),
                    MetaSVM = hl.if_else(
                                        hl.is_defined(combined_reference_v2_ht[ht.key].dbnsfp.MetaSVM_pred),
                                        combined_reference_v2_ht[ht.key].dbnsfp.MetaSVM_pred,
                                        hl.missing(hl.tstr)
                                        )
                    )
    ht = ht.annotate(
                    cadd = hl.if_else(
                                        hl.is_defined(cadd_v16_ht[ht.key].PHRED),
                                        cadd_v16_ht[ht.key].PHRED,
                                        hl.missing(hl.tfloat)
                                        )
                    )
    ht = ht.annotate(
                    clinvar = hl.if_else(
                                        hl.is_defined(clinvar_ht[ht.key].info.CLNSIG),
                                        clinvar_ht[ht.key].info.CLNSIG,
                                        hl.array([hl.missing(hl.tstr)])
                                        )
                    )
    ht = ht.annotate(
                    feature_type = hl.set(hl.if_else(
                                        hl.is_defined(gencode_gene_id_name_ht[ht.locus].featureType),
                                        gencode_gene_id_name_ht[ht.locus].featureType,
                                        hl.array([hl.missing(hl.tstr)])
                                        ))
                    )
    ht = ht.annotate(
                    low_complexity = hl.if_else(
                                        hl.is_defined(low_complex_ht[ht.locus]),
                                        "low_complexity",
                                        hl.missing(hl.tstr)
                                        )
                    )   
    
    return ht

## Annotate the hail table with omim information by gene_name
@typecheck(ht=hl.Table)
def ht_annotate_omim(ht):
    """
    Annotate the hail table with the following information:
    - omim_AR (omim_ar_ht)
    - omim (omim_ht)

    :param ht: hail table
    :return: hail table

    Notes: 
          the table for annotation should be named as same as the indicated inside ().
    """
    omim_ht = hl.read_table(myhl.omim_path)
    omim_ar_ht = hl.read_table(myhl.omim_ar_path)

    _ht = ht.annotate(
                    omim_AR = hl.if_else(
                                        hl.is_defined(omim_ar_ht[ht.gene_name].info),
                                        omim_ar_ht[ht.gene_name].info,
                                        hl.missing(hl.tstr)
                                        ),
                    omim = hl.if_else(
                                    hl.is_defined(omim_ht[ht.gene_name].info),
                                    omim_ht[ht.gene_name].info,
                                    hl.missing(hl.tstr)
                                    )
    )

    return _ht

## Additional annotation steps for cerebral palsy (CP) project
@typecheck(ht=hl.Table)
def ht_annotate_CP_related_DB(ht):
    """
    Annotate the hail table with the following information:
    - HSP CP related gene (hspgene_ht)
    - GeneDx CP related gene (genedx_cp_ht)

    :param ht: hail table
    :return: hail table

    Notes: 
          the table for annotation should be named as same as the indicated inside ().
    """
    hspgene_ht = hl.read_table(myhl.hspgene_path)
    genedx_cp_ht = hl.read_table(myhl.genedx_cp_path)

    _ht = ht.annotate(
                    HSP_gene = hl.if_else(
                                hl.is_defined(hspgene_ht[ht.locus].geneName),
                                hspgene_ht[ht.locus].geneName,
                                hl.missing(hl.tstr)
                                ),
                    GeneDx_CP = hl.if_else(
                                hl.is_defined(genedx_cp_ht[ht.locus].geneName),
                                genedx_cp_ht[ht.locus].geneName,
                                hl.missing(hl.tstr)
                                )
                        )
    
    return _ht

## Remove the variants without gene name and chrY showed in Female sample 
@typecheck(ht=hl.Table)
def remove_no_genename(ht):
    """
    Remove the vairants as the following information:
    - Variants has no gene name

    :param ht: hail table
    :return: hail table

    Notes:
          Make sure the variant table has been annotated with gene name (gencode_gene_id_name_ht).
    """
    ht = ht.filter(hl.is_missing(ht['gene_name']), keep=False)
    
    return ht


@typecheck(ht=hl.Table)
def remove_femalechrY(ht):
    """
    Remove the vairants as the following information:
    - Variants located at ChrY in Female samples

    :param ht: hail table
    :return: hail table

    Notes:
          Make sure the variant table has been annotated with gender information (Pheno_ht).
    """    
    #ht = ht.filter(hl.is_defined(ht['gene_name']), keep=True)
    ht = ht.filter(
                   ((ht.gender_mf == "Female") & (ht['#CHROM'] == "chrY")), keep=False
                  )
    
    return ht

######################
# VEP protein damage #
######################
## Annotate vep consequence information
@typecheck(ht=hl.Table)
def ht_annotation_vep_consequence(ht):
    """
    Annotate vep consequence information.
    - vep_consequence (vep104)

    :param ht: hail table
    :return: hail table
    """
    vep104 = myhl.vep104

    ## Store the hail table key information and set it to key_by alleles and locus (required by vep)
    key_info = ht.key.keys()
    ht = ht.key_by('locus','alleles')

    ## Annotation with VEP 104
    _ht = hl.vep(ht, vep104, csq=True, tolerate_parse_error=True)

    ## Get varaints category information (eq. start-loss, stop-gain, etc) from vep_consequence
    _ht = _ht.annotate(
                    vep_consequence = hl.set(hl.flatmap(lambda x: hl.if_else(x.contains("|"), x.split("\|")[1].split("\&"), hl.missing(hl.tarray(hl.tstr))), 
                                                                      _ht.vep))
                        )
    _ht = _ht.drop(_ht.vep)

    ## re-key to the original setting
    _ht = _ht.key_by(*key_info)
    
    return _ht

## Determine the protein damage variants based on the information of vep consequence, cadd, and metasvm 
@typecheck(ht=hl.Table)
def ht_annotate_proteindmg(ht):
    """
    Annotate the protein damage class information to the variant table

    :param ht: hail table
    :return: hail table

    Notes:
          The variant table should have vep_consequence, CADD, and MetaSVM information.
    """
    lof_consequences = hl.set(['splice_acceptor_variant', 'splice_donor_variant', 'frameshift_variant', 'stop_lost', 'stop_gained', 'start_lost'])
    Dmis_consequences = hl.set(['missense_variant', 'protein_altering_variant'])
    Syn_consequences = hl.set(['synonymous_variant'])

    _ht = ht.annotate(
                    protein_damage = hl.case()
                    .when(
                         hl.any(hl.map(lambda x: lof_consequences.contains(x), ht.vep_consequence)), "Lof"
                         )
                    .when(
                         (hl.any(hl.map(lambda x: Dmis_consequences.contains(x), ht.vep_consequence)) &
                           ((ht.caddv16_phred < 20) & (ht.MetaSVM != 'D'))), "Mis"
                         )
                    .when(
                         (hl.any(hl.map(lambda x: Dmis_consequences.contains(x), ht.vep_consequence)) & 
                           ((ht.caddv16_phred >= 20) | (ht.MetaSVM == 'D'))), "DMis"
                         )
                    .when(
                         hl.any(hl.map(lambda x: Syn_consequences.contains(x), ht.vep_consequence)), "Syn"
                         )
                    .default(
                            hl.missing(hl.tstr)
                            )
                    )
    return _ht

############
# CADD WEB #
############
## Output the variants has no cadd value for cadd web annotation
@typecheck(ht=hl.Table,
                  file=str)
def cadd_null_out_tsv(ht,file):
    """
    Save the variants without cadd-annotated information into txt file for cadd web annotation. 

    :param ht: hail table
    :param file: /path/to/variantTable
    :return: hail table

    Notes:
          - file was used to create the output name. ('file'_cadd_null.tsv)
          - maxium file size uploaded to cadd web is 2MB. (Example code in bash: split -b 2m $file -d "${file::-4}"_ or split -n 2 $file -d "${file::-4}"_ ) 
    """
    ht_noCADD = ht.filter(hl.is_missing(ht.cadd)).key_by()
    ht_noCADD_out = ht_noCADD.select("#CHROM","POS","ID","REF","ALT")
    
    ## out name
    #noCADD_outdir = file[:-4]+"_cadd_null.tsv"
    noCADD_outdir = mygl.create_output_name(file, pattern='cadd_null', extension='.tsv')

    ht_noCADD_out.export(noCADD_outdir,delimiter='\t')

    print("Output for cadd_null tsv file is at: ", noCADD_outdir)

## Read cadd_web.tsv file into hail table for annotation usage
@typecheck(cadd_file=str)
def read_webcadd(cadd_file):
    """
    Read in the tsv file generated from cadd web service. And then key by alleles, locus. 

    :param file: /path/to/cadd_web.tsv.bgz
    :return: hail table
    
    Notes:
          - Need to block gzip the cadd_web.tsv file. Example: bgzip -@4 $f.
    """
    rg38 = hl.get_reference('GRCh38')
    
    types = {'f0':hl.tstr, 'f1':hl.tint, 
               'f2':hl.tstr, 'f3':hl.tstr, 
               'f4':hl.tfloat, 'f5':hl.tfloat}
    header_name = {'f0':'contig', 'f1':'pos', 
                   'f2':'ref', 'f3':'alt', 
                   'f4':'rawscore', 'f5':'phred'}

    _cadd_web_ht = hl.import_table(cadd_file, delimiter='\t', 
                              comment=('#'), no_header=True, 
                              types=types, force_bgz=True).rename(header_name)

    #### Add locus and alleles information, and then key_by(locus,alleles)

    _cadd_web_ht = _cadd_web_ht.annotate(hg38_contig = hl.format('chr%s', _cadd_web_ht.contig))
    _cadd_web_ht = _cadd_web_ht.annotate(
                                      locus = hl.locus(_cadd_web_ht.hg38_contig,
                                                       _cadd_web_ht.pos, reference_genome=rg38),
                                      alleles = hl.array([_cadd_web_ht.ref, _cadd_web_ht.alt])
                                      )

    _cadd_web_ht = _cadd_web_ht.key_by('locus','alleles')
    
    return _cadd_web_ht

#####################
# Enhancer/Promoter #
#####################
@typecheck(ht=hl.Table)
def ht_annotate_vista_ccre_gse(ht):
    """
    Annotate the enhancer/promoter information to the variant table

    :param ht: hail table
    :return: hail table

    Notes:
          The variant table should have locus information.
    """
    ## VISTA
    vista_ht = hl.read_table(myhl.vista_path)
    ## ccREs
    ccre_ht = hl.read_table(myhl.ccre_path)
    ## GSE149268 (front lob ocr/pre)
    gse_ocr_ht = hl.read_table(myhl.gse_ocr_path)
    gse_pre_ht = hl.read_table(myhl.gse_pre_path)

    ## Annotation
    _ht = ht.annotate(
                VISTA = hl.if_else(
                                    hl.is_defined(vista_ht[ht.locus].info)
                                    , vista_ht[ht.locus].info
                                    , hl.missing(hl.tstr)
                                    ),
                CCRE = hl.if_else(
                                    hl.is_defined(ccre_ht[ht.locus].info)
                                    , ccre_ht[ht.locus].info
                                    , hl.missing(hl.tstr)
                                    ),
                GSE_OCR = hl.if_else(
                                    hl.is_defined(gse_ocr_ht[ht.locus].info)
                                    , gse_ocr_ht[ht.locus].info
                                    , hl.missing(hl.tstr)
                                    ),
                GSE_PRE = hl.if_else(
                                    hl.is_defined(gse_pre_ht[ht.locus].info)
                                    , gse_pre_ht[ht.locus].info
                                    , hl.missing(hl.tstr)
                                    )
                )

    return _ht

#################
# For Comp_side #
#################
## Step1: Get coulmn information with ID_GENE
@typecheck(ht=hl.Table)
def get_ID_gene(ht): 
    """
    Annotate the table with ID_GENE information for the comp_het vairants identification.

    :param ht: hail table
    :return: hail table 
    """
    _ht = ht.annotate(
                     ProbandID_GENE = ht.ProbandID + "_" + ht.gene_name
                    )
    return _ht

## Add count proband_gene column for those variants annotated with Comp_side(Comp_het)
@typecheck(ht=hl.Table)
def ht_annotate_ID_gene_count_comp_side(ht):
    """
    Get the count of how many variants with same ID and Gene name was annotated as Comp_side.

    :param ht: hail table
    :return: hail table

    Notes:
          This function was used in this function 'ht_annotate_comp_side'. 
    """
    ht_2 = ht.filter(~hl.is_missing(ht.ProbandID_GENE))
    dict_r_ID_GENE_comp_side = hl.dict(ht_2.aggregate(hl.agg.filter(ht_2['comp_side'].contains("Comp_side"), hl.agg.counter(ht_2['ProbandID_GENE']))))
    _ht = ht.annotate(
                    count_ID_GENE_comp_side = hl.if_else( 
                                                      dict_r_ID_GENE_comp_side.contains(ht['ProbandID_GENE']),
                                                      dict_r_ID_GENE_comp_side[ht['ProbandID_GENE']],
                                                      hl.missing(hl.tint)
                                                        )
                    )
    ht_2 = None
    dict_r_ID_GENE_comp_side = None
    
    return _ht

## Step2: Annotate the variants that fall into the comp_side criteria
@typecheck(ht=hl.Table)
def ht_annotate_comp_side(ht):
    """
    Annotate the variants that fall into the comp_side criteria:
    - MAF <= 0.001 both in bravo_v8 and gnomad_genome_v3.1
    - Proband_DP >= 8
    - Protein_damage_class = "Lof","DMis"
    - Remove ChrX and ChrY from Male samples (not included in this function, please see 'ht_annotate_my_comphet')

    :param ht: hail table
    :return: hail table 

    Notes:
          The comp_side will be diveded into comp_side_fa inherited from father or comp_side_ma inherited from mother.
    """
    _ht = ht.annotate(
                    comp_side = hl.case()
                    .when(
                         (
                             het_lst_hl.contains(ht.Proband_GT)
                         )
                       & (
                             (ref_lst_hl.contains(ht.Mother_GT)) & (het_lst_hl.contains(ht.Father_GT))
                         )
                       & (
                             (ht['bravo8_AF'] <= 0.001) & (ht['gnomad_genomes_AF'] <= 0.001)
                         )
                       & (
                             ht['Proband_DP_int'][0] + ht['Proband_DP_int'][1] >= 8
                         )
                       & (
                             (ht.protein_damage == "Lof") | (ht.protein_damage == "DMis")
                         )
                        ,"Comp_side_fa" 
                         )
                    .when(
                         (
                             het_lst_hl.contains(ht.Proband_GT)
                         )
                       & (
                            (ref_lst_hl.contains(ht.Father_GT)) & (het_lst_hl.contains(ht.Mother_GT))
                         )
                       & (
                             (ht['bravo8_AF'] <= 0.001) & (ht['gnomad_genomes_AF'] <= 0.001)
                         )
                       & (
                             ht['Proband_DP_int'][0] + ht['Proband_DP_int'][1] >= 8
                         )
                       & (
                             (ht.protein_damage == "Lof") | (ht.protein_damage == "DMis")
                         )
                        ,"Comp_side_ma" 
                         )
                    .default(
                            hl.missing(hl.tstr)
                            )
                    ) 

    _ht = ht_annotate_ID_gene_count_comp_side(_ht)
    
    return _ht

## Get the filtered hail table containing the trans-comp_het
@typecheck(ht=hl.Table)
def get_trans_comp_het_ht(ht):
    """
    Create a hail table that only contains the trans comp_het variants: have at least 2 comp_side variants in same gene within an individual, and these comp_side variants should at least coming from each parents.

    :param ht: hail table
    :return: hail table 

    Notes:
          - True: chr1_11111_A_G    S101_GENE    comp_side_fa; chr1_11432_C_T    S101_GENE    comp_side_ma
          - False: chr3_54816_T_C    S201_GENE    comp_side_fa; chr3_56418_A_T    S201_GENE    comp_side_fa 
    """
    _ht = ht.filter(~hl.is_missing(ht.comp_side))
    
    _ht = _ht.group_by(_ht.ProbandID_GENE).aggregate(trans_comp_side=hl.agg.collect_as_set(_ht.comp_side))
    
    _ht = _ht.filter(hl.set(['Comp_side_fa', 'Comp_side_ma']).is_subset(_ht.trans_comp_side))
    
    return _ht

#################
# Variant Types #
#################
## Annotate the variants that fall into the hom_alt criteria
@typecheck(ht=hl.Table)
def ht_annotate_my_homalt(ht):
    """
    Annotate the variants that fall into the hom_alt criteria:
    - MAF <= 0.001 both in bravo_v8 and gnomad_genome_v3.1
    - Proband_DP >= 8
    - Protein_damage_class = "Lof","DMis"
    - Remove ChrX and ChrY from Male samples

    :param ht: hail table
    :return: hail table 
    """
    _ht = ht.annotate(
            my_homalt = hl.case()
                    .when(
                         (
                             hom_lst_hl.contains(ht.Proband_GT)
                         )
                       & (
                             (ht['bravo8_AF'] <= 0.001) & (ht['gnomad_genomes_AF'] <= 0.001)
                         )
                       & (
                             ht['Proband_DP_int'][0] + ht['Proband_DP_int'][1] >= 8
                         )
                       & (
                             (ht.protein_damage == "Lof") | (ht.protein_damage == "DMis")
                         )
                      & ~(
                             (ht.gender_mf == "Male") & ((ht['#CHROM'] == "chrX") | (ht['#CHROM'] == "chrY"))
                         )
                          , "Hom_alt"
                         )
                    .default(
                            hl.missing(hl.tstr)
                            )
                    )
    
    return _ht

## Annotate the variants that fall into the comp_het criteria
@typecheck(ht=hl.Table,
           ht_trans_comp_het=hl.Table)
def ht_annotate_my_comphet(ht,ht_trans_comp_het):
    """
    Annotate the variants that fall into the comp_het criteria:
    - Identified in trans_comp_side hail table
    - Remove ChrX and ChrY from Male samples

    :param ht: hail table
    :return: hail table 

    Notes:
          More information was describe in 'ht_annotate_comp_side'.
    """
    _ht = ht.annotate(
            my_comphet = hl.case()
                    .when(
                         (
                             hl.is_defined(ht_trans_comp_het[ht.ProbandID_GENE])
                         ) 
                       & (
                             ht['comp_side'].contains("Comp_side")
                         )
                      & ~(
                             (ht.gender_mf == "Male") & ((ht['#CHROM'] == "chrX") | (ht['#CHROM'] == "chrY"))
                         )
                         , "Comp_het"
                         )
                    .default(
                            hl.missing(hl.tstr)
                            )
                    )
    
    return _ht

## Annotate the variants that fall into the denovo(heterozygous genotype) criteria
@typecheck(ht=hl.Table)
def ht_annotate_my_denovo(ht):
    """
    Annotate the variants that fall into the denovo(GT-het) criteria:
    - Proband_GT is heterozygous in alternative allel while Parent_GT are homozygous in reference allel
    - MAF <= 0.0004 both in bravo_v8 and gnomad_genome_v3.1
    - Proband DP >= 10 and Proband AD_alt >= 5
    - Proband AP_alt >= 10 and Proband_VAF >= 0.2
    - Proband AD_alt < 10 and Proband_VAF >= 0.28
    - Parents DP >= 10 and Parents VAF < 0.035
    - Protein_damage_class = "Lof","DMis","Mis","Syn"

    :param ht: hail table
    :return: hail table 
    """

    _ht = ht.annotate(
            my_denovo = hl.case()                    
                    .when(
                         (
                             het_lst_hl.contains(ht.Proband_GT)
                         )
                       & (
                             (ref_lst_hl.contains(ht.Mother_GT)) & (ref_lst_hl.contains(ht.Father_GT)) 
                         )
                       & (
                             ht.denovo_filtration == "Pass"
                         )
                         , "denovo_het"
                         )
                    .default(
                            hl.missing(hl.tstr)
                            )
                    )
    
    return _ht

## Annotate the variants that fall into the hemizygous (both chrX and chrY) criteria
@typecheck(ht=hl.Table)
# def ht_annotate_my_hemi_XY(ht):
#     """
#     Annotate the variants that fall into the hemizygous criteria:
#     - ChrX and ChrY in Male samples
#     - MAF <= 0.00005 both in bravo_v8 and gnomad_genome_v3.1
#     - Proband DP >= 8
#     - Proband GQ >= 20
#     - MQ >= 40
#     - Proband_VAF >= 0.95
#     - Protein_damage_class = "Lof","DMis"

#     :param ht: hail table
#     :return: hail table
#     """
#     _ht = ht.annotate(
#             my_hemi_XY = hl.case()              
#                     .when(
#                          (
#                              hom_lst_hl.contains(ht.Proband_GT)
#                          )
#                        & (
#                              (ht.gender_mf == "Male") & ((ht['#CHROM'] == "chrX") | (ht['#CHROM'] == "chrY"))
#                          )
#                        & (
#                              (ht['bravo8_AF'] <= 0.0005) & (ht['gnomad_genomes_AF'] <= 0.0005)
#                          )
#                        & (
#                              ht['Proband_DP_int'][0] + ht['Proband_DP_int'][1] >= 8
#                          )
#                        & (
#                              ht['Proband_GQ'] >= 20
#                          )
#                        & (
#                              ht['MQ'] >= 40
#                          )
#                        & (
#                              ht['Proband_VAF'] >= 0.95
#                          )
#                        & (
#                              (ht.protein_damage == "Lof") | (ht.protein_damage == "DMis")
#                          )
#                          , "hemizygous"
#                          )
#                     .default(
#                             hl.missing(hl.tstr)
#                             )
#                     )
    
#     return _ht

def ht_annotate_my_hemi_XY(ht):
    """
    Annotate hemizygous variants on chrX or chrY in males with:
    - Proband GT homozygous
    - rare in population
    - high DP, GQ, VAF
    - optional MQ ≥ 40 (if available)
    - protein_damage in Lof or DMis
    """
    has_MQ = 'MQ' in ht.row

    hemi_criteria = (
        (hom_lst_hl.contains(ht.Proband_GT)) &
        (ht.gender_mf == "Male") &
        ((ht['#CHROM'] == "chrX") | (ht['#CHROM'] == "chrY")) &
        (ht['bravo8_AF'] <= 0.0005) &
        (ht['gnomad_genomes_AF'] <= 0.0005) &
        ((ht['Proband_DP_int'][0] + ht['Proband_DP_int'][1]) >= 8) &
        (ht['Proband_GQ'] >= 20) &
        (ht['Proband_VAF'] >= 0.95) &
        ((ht.protein_damage == "Lof") | (ht.protein_damage == "DMis"))
    )

    if has_MQ:
        hemi_criteria = hemi_criteria & (ht['MQ'] >= 40)

    _ht = ht.annotate(
        my_hemi_XY = hl.if_else(hemi_criteria, "hemizygous", hl.missing(hl.tstr))
    )

    return _ht

##########
# Denovo #
##########
## Annotate the variants with the hard filtration criteria: wes denovo
@typecheck(ht=hl.Table)
def ht_annotate_denovo_filtration(ht):
    """
     Annotate the variants that fall into the following criteria:
    - MAF <= 0.0004 both in bravo_v8 and gnomad_genome_v3.1
    - Proband DP >= 10 and Proband AD_alt >= 5
    - Proband AP_alt >= 10 and Proband_VAF >= 0.2
    - Proband AD_alt < 10 and Proband_VAF >= 0.28
    - Parents DP >= 10 and Parents VAF < 0.035
    - Protein_damage_class = "Lof","DMis","Mis","Syn"

    :param ht: hail table
    :return: hail table
    """
    _ht = ht.annotate(
            denovo_filtration = hl.case()
                    .when(
                         (
                             (ht['bravo8_AF'] <= 0.0004) & (ht['gnomad_genomes_AF'] <= 0.0004)
                         )
                       & (
                             (ht['Proband_DP_int'][1] >= 5) & (ht['Proband_DP_int'][0] + ht['Proband_DP_int'][1] >= 10)
                         )
                       & (
                             ((ht['Proband_VAF'] >= 0.2) & (ht['Proband_DP_int'][1] >= 10)) | ((ht['Proband_VAF'] >= 0.28) & (ht['Proband_DP_int'][1] < 10))
                         )
                       & (
                             (ht['Mother_DP_int'][0] >= 10) & (ht['Father_DP_int'][0] >= 10)
                         )
                       & (
                             (ht['Mother_VAF'] < 0.035) & (ht['Father_VAF'] < 0.035)
                         )
                       & (
                             (ht.protein_damage == "Lof") | (ht.protein_damage == "DMis") | (ht.protein_damage == "Mis") | (ht.protein_damage == "Syn")
                         )
                         , "Pass"
                         )
                 .default(
                             hl.missing(hl.tstr)
                         )
                 )
    
    return _ht

## Process triodenovo variant Table
### For setting the key of the triodenovo_ht
@typecheck(triodenovo_ht=hl.Table)
def key_by_alleles_locus_id_ht(triodenovo_ht):
    """
    Setting the key of the triodenovo_ht: 'locus','alleles','ProbandID'

    :param triodenovo_ht: hail table
    :return: hail table
    """
    rg38 = hl.get_reference('GRCh38')

    _ht = triodenovo_ht.annotate(
                locus = hl.locus(triodenovo_ht['#CHROM'], triodenovo_ht.POS, reference_genome=rg38),
                alleles = hl.array([triodenovo_ht.REF, triodenovo_ht.ALT])
                )

    _ht = _ht.key_by('locus','alleles','ProbandID')
    
    return _ht 

### Annotate the variants that was found in trio_denovo
@typecheck(ht=hl.Table,
           triodenovo_ht=hl.Table)
def ht_annotate_triodenovo(ht, triodenovo_ht):
    """
    Annotate the variants that was found in the trio_denovo hail table

    :param ht: hail table
    :param triodenovo_ht: denovo hail table generated from trio_denovo
    :return: hail table
    """
    ht = ht.key_by('locus', 'alleles', 'ProbandID')
    
    _ht = ht.annotate(
            trio_denovo = hl.if_else(
                                    (hl.is_defined(triodenovo_ht[ht.key]) & (ht.denovo_filtration == "Pass"))
                                  , "Pass"
                                  , hl.missing(hl.tstr)
                                  )
                 )
    
    _ht = _ht.key_by('locus', 'alleles')
    
    return _ht

####################
# Positive Control #
####################
@typecheck(ht=hl.Table,
           ht_pc=hl.Table)
def annotate_positive_control(ht,ht_pc):
    """
    Annotate the variants that was found in the positive control hail table

    :param ht: hail table
    :param ht_pc: hail table containing the positive control variants 
    :return: hail table
    """
    ht = ht.key_by('locus', 'alleles', 'ProbandID')
    
    _ht = ht.annotate(
                     positive_control = hl.if_else(
                                         hl.is_defined(ht_pc[ht.key].Inheritance)
                                       , ht_pc[ht.key].Inheritance
                                       , hl.missing(hl.tstr)
                     )
    
    )
    
    _ht = _ht.key_by('locus', 'alleles')
    
    return _ht

##########
# OTHERS #
##########
## Select variant types (exome)
@typecheck(ht=hl.Table)
def ht_select_VariantTypes(ht):
    """
    Only keep the variants that has been identified as a certian type of variants (hom_alt, denovo, comp_het, and hemizygous)

    :param ht: hail table
    :return: hail table
    """
    _ht = ht.filter( (hl.is_missing(ht.my_homalt)) 
                   & (hl.is_missing(ht.my_comphet))
                   & (hl.is_missing(ht.my_denovo))
                   & (hl.is_missing(ht.my_hemi_XY))
                #    & (hl.is_missing(ht.hail_denovo))
                #    & (hl.is_missing(ht.trio_denovo))
                #    & (hl.is_missing(ht.positive_control))
                   , keep=False)
    
    return _ht

#######################
# PROJECTS ANNOTATION #
#######################
class annotation_project:
    def __init__(self, phenotype_path=None, cadd_web_path=None, denovo_path=None, upd_path=None):
        ## Attributes
        self.phenotype_path = phenotype_path
        self.cadd_web_path = cadd_web_path
        self.denovo_path = denovo_path
        self.upd_path= upd_path

        ## read phenotype file
        if phenotype_path is not None:
            self.phenotype_ht = read_txt2ht_simple(self.phenotype_path,'ProbandID')
		
        ## read cadd_web file
        if cadd_web_path is not None: 
            self.cadd_web_ht = read_webcadd(self.cadd_web_path)

        ## read denovo file
        if denovo_path is not None:
            self.denovo_ht = hl.read_table(self.denovo_path)

        ## read upd_bed_merge file
        if upd_path is not None:
            self.upd_ht = read_ht_locus_from_triomix_upd_bed(self.upd_path)

    ## Additional annotation steps for phenotype information
    @typecheck_method(ht=hl.Table)
    def annotation_pheno(self, ht):
        """
        Annotate the hail table with the following information:
        - gender (Pheno_ht)
        - race (Pheno_ht)
        - disease_info (Pheno_ht)

        :param ht: hail table
        :return: hail table

        Notes: 
            the table for annotation should be named as same as the indicated inside ().
        """
        phenotype_ht = self.phenotype_ht

        _ht = ht.annotate(
            gender_mf = hl.if_else(
                hl.is_defined(phenotype_ht[ht.ProbandID].Sex),
                phenotype_ht[ht.ProbandID].Sex,
                hl.missing(hl.tstr)
            ),
            race_mf = hl.if_else(
                hl.is_defined(phenotype_ht[ht.ProbandID].Ethnicity),
                phenotype_ht[ht.ProbandID].Ethnicity,
                hl.missing(hl.tstr)
            ),
            disease_info_mf = hl.if_else(
                hl.is_defined(phenotype_ht[ht.ProbandID].Disease_info),
                phenotype_ht[ht.ProbandID].Disease_info,
                hl.missing(hl.tstr)
            ),
            # CP_mf = hl.if_else(
            #     hl.is_defined(phenotype_ht[ht.ProbandID].CP),
            #     phenotype_ht[ht.ProbandID].CP,
            #     hl.missing(hl.tstr)
            # ),
            # Consanguinity_mf = hl.if_else(
            #     hl.is_defined(phenotype_ht[ht.ProbandID].Consanguinity),
            #     phenotype_ht[ht.ProbandID].Consanguinity,
            #     hl.missing(hl.tstr)
            # )
        )

        return _ht

    ## Annotate cadd web information
    @typecheck_method(ht=hl.Table)
    def ht_annotate_caddweb(self, ht):
        """
        Annotate the CADD web information to those variants has null cadd value.

        :param ht: hail table
        :return: hail table

        Notes:
            - Phred score will be annotated. 
            - Name of caddweb hail table should be 'cadd_v16_ht'.
        """
        cadd_v16_ht = hl.read_table(myhl.cadd_v16_ht_path)
        cadd_web_ht = self.cadd_web_ht

        _ht = ht.annotate(
                    caddv16_phred = hl.if_else(
                                            hl.is_defined(cadd_v16_ht[ht.key].PHRED),
                                            cadd_v16_ht[ht.key].PHRED,
                                            hl.if_else(
                                                        hl.is_defined(cadd_web_ht[ht.key].phred),
                                                        cadd_web_ht[ht.key].phred,
                                                        hl.missing(hl.tfloat)
                                                        )
                                            )
                    ) 
        
        return _ht

    ## Annotate the varaints that was found in the hl.denovo hail table
    @typecheck_method(ht=hl.Table)
    def ht_annotate_hail_denovo(self, ht):
        """
        Annotate the variants that was found in the hl.denovo hail table

        :param ht: hail table
        :param denovo_ht: denovo hail table generated from hl.denovo
        :return: hail table
        """
        denovo_ht = self.denovo_ht

        ht = ht.key_by('locus', 'alleles', 'ProbandID')

        _ht = ht.annotate(
                        Denovo_hail_confidence = hl.if_else(
                                                hl.is_defined(denovo_ht[ht.key].confidence)
                                                , denovo_ht[ht.key].confidence
                                                , hl.missing(hl.tstr)
                                                )
                      , Denovo_hail_p = hl.if_else(
                                                hl.is_defined(denovo_ht[ht.key].p_de_novo)
                                                , denovo_ht[ht.key].p_de_novo
                                                , hl.missing(hl.tfloat)
                                                )
                      , hail_denovo = hl.if_else(
                                                (hl.is_defined(denovo_ht[ht.key]) & (ht.denovo_filtration == "Pass"))
                                                , "Pass"
                                                , hl.missing(hl.tstr)
                                                )

                        )

        _ht = _ht.key_by('locus', 'alleles')
    
        return _ht

    @typecheck_method(ht=hl.Table)
    def ht_annotate_my_upd(self, ht):
        """
        Annotate the variants that fall into the upd region:
        - Identified in upd bed hail table

        :param ht: hail table
        :return: hail table 
        """
        upd_ht = self.upd_ht
        
        ## Get upd groupd by SampleID table with a set of intervals
        upd_ht_group = upd_ht.group_by(upd_ht.SampleID).aggregate(intervals=hl.agg.collect_as_set(upd_ht.interval))
        
        ## Annotate the set of intervals back to hail table
        ht = ht.annotate(
                        my_upd_interval = hl.if_else(
                                                hl.is_defined(upd_ht_group[ht.ProbandID].intervals),
                                                upd_ht_group[ht.ProbandID].intervals,
                                                hl.set([hl.interval(hl.locus(hl.missing(hl.tstr), hl.missing(hl.tint)),
                                                            hl.locus(hl.missing(hl.tstr), hl.missing(hl.tint)))])
                                                )
                    )
        
        _ht = ht.annotate(
                        my_upd = hl.if_else(
                                            hl.any(hl.map(lambda x: x.contains(ht.locus), ht.my_upd_interval)),
                                            "Pass",
                                            hl.missing(hl.tstr)
                                            ),
                        contamination = hl.if_else(
                                            hl.is_defined(upd_ht[ht.ProbandID].Contamination),
                                            upd_ht[ht.ProbandID].Contamination,
                                            hl.missing(hl.tstr)
                                            ),
                        my_upd_5MB = hl.if_else(
                                            hl.any(hl.map(lambda x: x.contains(ht.locus), ht.my_upd_interval)),
                                            upd_ht[ht.ProbandID].gt_5Mb,
                                            hl.missing(hl.tstr)
                                            )
                        )

        return _ht

## Process with no caddweb info
@typecheck_method(ht=hl.Table)
def ht_annotate_caddv16_only(self, ht):
    """
    Annotate CADD PHRED scores from the local cadd_v16_ht only.
    If not available in cadd_v16_ht, result is missing.

    :param ht: Hail Table
    :return: Hail Table with caddv16_phred column
    """
    _ht = ht.annotate(
            caddv16_phred = ht.cadd
        )

    return _ht


#####################
# DENOVO HAIL TABLE #
#####################
# def concat_hail_table_from_list(file_list):
#     ht_union = hl.Table
#     for i in file_list:
#         file = i
#         base = os.path.basename(file)
#         print(base + ' is running')
        
#         ht = hl.read_table(file)
#         #print("current set rows: ", ht.count())
#         ht_union = ht_union.union(ht)
#         #print("total row: ", ht_union.count())

#     #out_file =
#     mygl.out_save(ht_union, out_file)


####################
# STEPS ANNOTATION #
####################
## Step1: output the cadd_null variants for web annotation
@mygl.calculate_running_time
@typecheck(input_file=str,
           output_file=nullable(str),
           phenotype_file=nullable(str),
           unique_FamID=nullable(str),
           exclude_samples=nullable(str))
def step1_annotation(input_file, output_file=None, phenotype_file=None, unique_FamID=None, exclude_samples=None):
    """
    Annotate the variant table txt file for the step1 annotation, including:
    - remove variants with missing information in Proband_GT, REF, and ALT.
    - annotate gene_name, MAF (bravoV8, gnomadV3), metaSVM, OMIM, and feature type.
    - remove the variants without gene_name annotation, and the chrY variants in female samples.
    - annotate with CP-related dataset (optional.)
    - annotate with vep consequence information.
    - save the output into hail table.
    - save the cadd_null variants to .tsv file for web cadd annotation.
    """
    rg38 = hl.get_reference('GRCh38')

    ## Set the output file if it is not provided
    if output_file is None:
        output_file = mygl.create_output_name(input_file,"step1",".ht")
    
    ## Read-in the txt file to hail table
    if isinstance(input_file, str) and mygl.file_exists(input_file) and input_file.endswith('.ht'):
        print("read the hail input table.")
        ht = hl.read_table(input_file)
        ht = ht.key_by()
        ht = ht.annotate( locus = hl.locus(ht.locus.contig, ht.locus.position, reference_genome=rg38))
        ht = ht.key_by('locus', 'alleles')
    elif isinstance(input_file, str) and mygl.file_exists(input_file) and input_file.endswith('.txt'):
        print("read the input txt file.")
        ht = read_variantTable_ht_rekey(input_file)
    else:
        raise ValueError("Input_file must be the path/file of Hail Table or variant table txt file.")

    ## Add ALT, and REF if doesn't exist
    if not hasattr(ht.row, 'ALT'):
        ht = ht.annotate( ALT = ht.alleles[1] )
    if not hasattr(ht.row, 'REF'):
        ht = ht.annotate( REF = ht.alleles[0] )

    ## Get Proband_GT from GT column (for entries_ht derived from mt)
    if not hasattr(ht.row, 'Proband_GT'):
        ht = ht.annotate( Proband_GT = hl.str(ht.GT) )

    ## Temporary: Select the samples want to include in the analysis
    if unique_FamID is not None:
        ht = select_FamID(ht, unique_FamID, exclude_samples)
    
    ## remove none information
    ht = remove_no_ref_alt_gt(ht)
    
    ## annotate phenotype
    if phenotype_file is not None:
        p1 = annotation_project(phenotype_path=phenotype_file)
        ht = p1.annotation_pheno(ht)
    
    ## Save the temperory output
    output_file_temp1 = mygl.create_output_name(input_file,"step1_temp1",".ht")
    if not mygl.file_exists(output_file_temp1):
        mygl.out_save(ht, out_file=output_file_temp1)
    else:
        print("step1_temp1 file exists, proceed to population info annotation")
    
    ht = hl.read_table(output_file_temp1)

    ## annotate essential informations
    ht = ht_annotate_gene_bravo_gnomad_meta_cadd_feature_type(ht)
    ht = ht_annotate_omim(ht)
    
    ## remove none information for gene_name and chrY variants in female samples
    if phenotype_file is not None:
        ht = remove_femalechrY(ht)
    
    ## Temporary: annotate CP-related dataset
    ht = ht_annotate_CP_related_DB(ht)

    ## Save the temperory output
    output_file_temp2 = mygl.create_output_name(input_file,"step1_temp2",".ht")
    if not mygl.file_exists(output_file_temp2):
        mygl.out_save(ht, out_file=output_file_temp2)
    else:
        print("step1_temp2 file exists, proceed to vep annotation")
    
    ht = hl.read_table(output_file_temp2)

    ## annotate the vep consequence information
    ht = ht_annotation_vep_consequence(ht)

    ## Save the step1 output into hail table
    mygl.out_save(ht, out_file=output_file)

    ## Remove temp files
    if mygl.file_exists(output_file_temp1):
        mygl.remove_item(output_file_temp1)

    if mygl.file_exists(output_file_temp2):
        mygl.remove_item(output_file_temp2)

    ## Read-in output_file
    ht = hl.read_table(output_file)

    ## Save the cad_null variants to tsv
    cadd_null_out_tsv(ht, output_file)

    ## Print completion
    print("Step1 is completed.")

## Step2: annotate the protein damage information
@mygl.calculate_running_time
@typecheck(input_file=str,
           output_file=nullable(str),
           cadd_web_file=nullable(str))
def step2_annotation(input_file, output_file=None, cadd_web_file=None):
    """
    Annotate the variant table for the step2 annotation, including:
    - annotate the cadd value from website annotation information.
    - annotate protein damage information.
    - save the output into hail table.
    """
    ## Read the step1 output hail table
    ht = hl.read_table(input_file)
    
    ## Set the output file if it is not provided
    if output_file is None:
        output_file = mygl.create_output_name(input_file,"step2",".ht")

    ## Annotate cadd_web information
    if cadd_web_file is not None:
        p2 = annotation_project(cadd_web_path=cadd_web_file)
        ht = p2.ht_annotate_caddweb(ht)
    else:
        ht = ht_annotate_caddv16_only(ht)

    ## Annotate protein damage information
    ht = ht_annotate_proteindmg(ht)
    
    ## Save the step2 output into hail table
    mygl.out_save(ht, out_file=output_file)

    ## Print completion
    print("Step2 is completed.")

## Step3: annotate the variant types (hom_alt, comphet, denovo, hemizygous) (WES)
@mygl.calculate_running_time
@typecheck(input_file=str,
           output_file=nullable(str),
           denovo_file=nullable(str),
           upd_file=nullable(str),
           true_set_file=nullable(str))
def step3_annotation(input_file, output_file=None, denovo_file=None, upd_file=None, true_set_file=None):
    """
    Annotate the variant table for the step3 annotation, including:
    - annotate the different variant types (denovo, comphet, hom_alt, hemi_xy)
    - save the output into hail table and plotread.txt file
    """
    ## Read the step1 output hail table
    ht = hl.read_table(input_file)

    ## Set the output file if it is not provided
    if output_file is None:
        output_file = mygl.create_output_name(input_file,"step3",".ht")

    ## For comp-het:
    ### 1st step: generate new column with ID and gene togather
    ht = get_ID_gene(ht)
    ### 2nd step: annotate comp_side information and get trans_comp_het hail table
    ht = ht_annotate_comp_side(ht)
    ht_trans_comp_het = get_trans_comp_het_ht(ht)
    ## Annotate denovo, comphet, hom_alt, hemizygous
    ht = ht_annotate_my_homalt(ht)
    ht = ht_annotate_my_comphet(ht, ht_trans_comp_het)
    ht = ht_annotate_my_hemi_XY(ht)
    ## Annotate denovo filtration criteria results
    ht = ht_annotate_denovo_filtration(ht)
    ht = ht_annotate_my_denovo(ht)

    ## Annotate hail denovo
    p3 = annotation_project(denovo_path=denovo_file, upd_path=upd_file)
    if denovo_file is not None:
        ht = p3.ht_annotate_hail_denovo(ht)
    ## Annotate upd information
    if upd_file is not None:
        ht = p3.ht_annotate_my_upd(ht)

    ## Annotate positive control
    if true_set_file is not None:
        ht_pc = hl.read_table(true_set_file)
        ht = annotate_positive_control(ht,ht_pc)
    
    ## Further select
    ht = ht_select_VariantTypes(ht)

    ## Save the step3 output into hail table
    mygl.out_save(ht, out_file=output_file)

    ## Read-in output_file
    ht = hl.read_table(output_file)

    ## Save out for plotread
    out_plotread =  mygl.create_output_name(output_file, pattern="forplotread",extension=".txt")
    mygl.out_for_downstream(ht,out_plotread,select_all=True)

    ## Print completion
    print("Step3 is completed.")

#####################################
# STEPS BLAT/DEEPVARIANT ANNOTATION #
#####################################
#############
# FUNCTIONS #
#############
## Read Blat/Deep txt file
@typecheck(file=str)
def transform_txt_ht(file):
    """
    Read the txt file into hail table and then key by ProbandID and locus.

    :param file: /path/to/file
    :return: hail table
    """
    rg38 = hl.get_reference('GRCh38')
    _ht = hl.import_table(file, delimiter="\t")
    _ht = _ht.annotate(
                      ProbandID=_ht['#Variants_Chr_Pos'].split(":")[0],
                      locus = hl.locus(_ht['#Variants_Chr_Pos'].split(":")[1], hl.int(_ht['#Variants_Chr_Pos'].split(":")[2]), reference_genome=rg38)
                      ).key_by('ProbandID','locus')
    
    return _ht

@typecheck(file=str)
def read_ht_from_txt(file):
    """
    Read the Blat correction file (ProbandID,#CHROM,POS,Blat).

    :param file: /path/to/file
    :return: hail table
    """
    rg38 = hl.get_reference('GRCh38')
    types={
           'POS':hl.tint
          }

    _ht = hl.import_table(file,delimiter='\t',types=types)

    ## Set the locus, and key_by locus
    _ht = _ht.annotate(
                      locus=hl.locus(_ht['#CHROM'], _ht.POS, reference_genome=rg38)
                      ).key_by('ProbandID','locus')
    
    return _ht

@typecheck(Blat_ht=hl.Table,
           Blat_corr_ht=nullable(hl.Table))
def blat_correction(Blat_ht, Blat_corr_ht=None):
    """
    Correct the No_value in Blat_validation file from the corrected blat file (manually check on UCSC blat)

    :param Blat_ht: hail table
    :return: hail table
    """

    if Blat_corr_ht is not None:
        _ht = Blat_ht.annotate(
                            new_info = hl.if_else(hl.is_defined(Blat_corr_ht[Blat_ht.key].Blat),
                                                    Blat_corr_ht[Blat_ht.key].Blat,
                                                    Blat_ht['Pass_or_Failed']
                                                )
                            )
    else:
        _ht = Blat_ht.annotate(
                            new_info = Blat_ht['Pass_or_Failed']                                
                            )
    
    return _ht

## Annotation blat, deep info
@typecheck(ht=hl.Table,
           Blat_ht=nullable(hl.Table),
           Deep_ht=nullable(hl.Table))
def annotation_blat_deep(ht, Blat_ht=None, Deep_ht=None):
    """
    Annotate the blat and deepvariant information back to hail table.

    :param ht: hail table
    :param Blat_ht: hail table (optional)
    :param Deep_ht: hail table (optional)
    :return: hail table
    """
    
    if Blat_ht is not None:
        ht = ht.annotate(
            Blat=hl.if_else(
                hl.is_defined(Blat_ht[ht.key].new_info),
                Blat_ht[ht.key].new_info,
                hl.missing(hl.tstr)
            )
        )

    if Deep_ht is not None:
        ht = ht.annotate(
            Deepvariant=hl.if_else(
                hl.is_defined(Deep_ht[ht.key].Pass_or_Failed),
                Deep_ht[ht.key].Pass_or_Failed,
                hl.missing(hl.tstr)
            )
        )
    
    return ht

### Hail table filteration function:
@typecheck(ht=hl.Table,
           Blat_ht=nullable(hl.Table),
           Deep_ht=nullable(hl.Table))
def keep_pass_blat_deep(ht, Blat_ht=None, Deep_ht=None):

    if Blat_ht is not None:
        ht = ht.filter((ht.Blat == "Pass"), keep=True)
    
    if Deep_ht is not None:
        ht = ht.filter((ht.Deepvariant == "Pass"), keep=True)
    
    return ht

@typecheck(ht=hl.Table)
def keep_proteindmg(ht):
    _ht =  ht.filter((ht.protein_damage == "Lof") | (ht.protein_damage == "DMis"), keep=True)
    return _ht

### Get recurrent information
#### Only calculate the infromation from Blat and Deep pass variants
@typecheck(ht=hl.Table,
           Blat_ht=nullable(hl.Table),
           Deep_ht=nullable(hl.Table))
def recur_byID(ht, Blat_ht=None, Deep_ht=None):
    ht_2 = keep_pass_blat_deep(ht, Blat_ht, Deep_ht)
    dict_r_ID = hl.dict(ht_2.aggregate(hl.agg.counter(ht_2['ID'])))
    _ht = ht.annotate(
                    recurrent_ID = hl.if_else(
                                             dict_r_ID.contains(ht['ID']),
                                             dict_r_ID[ht['ID']],
                                             1
                                             )
                    )
    ht_2 = None
    dict_r_ID = None

    return _ht

@typecheck(ht=hl.Table,
           Blat_ht=nullable(hl.Table),
           Deep_ht=nullable(hl.Table)) 
def recur_bygene(ht, Blat_ht=None, Deep_ht=None):
    ht_2 = keep_pass_blat_deep(ht, Blat_ht, Deep_ht)
    ht_2 = remove_no_genename(ht_2)
    dict_r_GENE = hl.dict(ht_2.aggregate(hl.agg.counter(ht_2['gene_name'])))
    _ht = ht.annotate(
                    recurrent_GENE = hl.if_else(
                                                dict_r_GENE.contains(ht['gene_name']),
                                                dict_r_GENE[ht['gene_name']],
                                                hl.missing(hl.tint)
                                                )
                    )
    ht_2 = None
    dict_r_GENE = None

    return _ht

@typecheck(ht=hl.Table,
           Blat_ht=nullable(hl.Table),
           Deep_ht=nullable(hl.Table))
def recur_proteindmg_bygene(ht, Blat_ht=None, Deep_ht=None):
    ht_2 = keep_pass_blat_deep(ht, Blat_ht, Deep_ht)                     
    ht_2 =  keep_proteindmg(ht_2)
    ht_2 = remove_no_genename(ht_2)
    dict_r_GENE_dmg = hl.dict(ht_2.aggregate(hl.agg.counter(ht_2['gene_name'])))
    _ht = ht.annotate(
                    recurrent_GENE_dmg = hl.if_else(
                                                    dict_r_GENE_dmg.contains(ht['gene_name']),
                                                    dict_r_GENE_dmg[ht['gene_name']],
                                                    hl.missing(hl.tint)
                                                   )
                    )
    ht_2 = None
    dict_r_GENE_dmg = None

    return _ht

from hail.typecheck import typecheck, nullable

#########
# STEPS #
#########
## Step4 for blat_deep: annotate the blat and deep information
#@mygl.calculate_running_time
@typecheck(input_file=str,
           output_file=nullable(str),
           deep_file=nullable(str),
           blat_file=nullable(str),
           blat_correction_file=nullable(str),
           test_mode=str)
def step4_annotation_deep_blat(input_file, output_file=None, deep_file=None, blat_file=None, blat_correction_file=None, test_mode="no"):
    """
    Annotate the variant table for the step4 annotation, including:
    - annotate the blat result.
    - annotate the deepvariant result.

    :param input_file: hail table path
    :param output_file: hail table output path (optional)
    :param deep_file: text file path for deepvariant (optional)
    :param blat_file: text file path for blat (optional)
    :param blat_correction_file: text file path for blat corrections (optional)
    :return: hail table
    """

    # Error Handling
    if not deep_file and not blat_file:
        raise ValueError("Either deep_file or blat_file should be provided.")
    
    # Test mode and output_file check
    if test_mode.lower() == "no" and not output_file:
        raise ValueError("Output_file should be provided if not in test mode.")

    # Initialize the Hail Tables
    Blat_ht = None
    Deep_ht = None
    Blat_corr_ht = None
    
    # Read-in Blat and Deepvariant Table
    if blat_correction_file:
        Blat_corr_ht = read_ht_from_txt(blat_correction_file)
    
    if blat_file:
        Blat_ht = transform_txt_ht(blat_file)        
        Blat_ht = blat_correction(Blat_ht, Blat_corr_ht)
    
    if deep_file:
        Deep_ht = transform_txt_ht(deep_file)
    
    # Read-in Variant Table
    ht = read_variant_table_step_blat_deep(input_file)
    
    # Annotate the main table
    ht = annotation_blat_deep(ht, Blat_ht=Blat_ht, Deep_ht=Deep_ht)
    
    # Save the output based on test_mode
    if test_mode.lower() == "no":
        mygl.out_save(ht, out_file=output_file)
        print("Step4_blat_deep is completed.")
    elif test_mode.lower() == "yes":
        return ht
    else:
        raise ValueError("test_mode should be either 'no' or 'yes'")
    
    return None