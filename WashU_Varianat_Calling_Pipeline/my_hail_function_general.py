#!/usr/bin/python
###################################
# Purpose: Functions for the Hail general usage
# Author: Yung-Chun Wang <yung-chun@wustl.edu>
# Modifier: 
# Language: Python
# Version: 1
# Comment: 
# Created Date: 02-02-2023
# Last Modified Date: 04-04-2023
###################################

import os
import sys
import glob
import argparse
import subprocess
import pyspark
from pyspark.sql import SparkSession
import time
import shutil
import multiprocessing
from datetime import datetime
import psutil
import pandas as pd

import hail as hl
import hail.expr as expr
import hail.expr.types as t
import hail.expr.functions as functions
import hail.genetics as genetics
import hail.ir as ir
from hail.typecheck import typecheck, typecheck_method, nullable, anytype, enumeration, tupleof, func_spec, oneof, arg_check, args_check, anyfunc, sequenceof
from typing import Union
from hail.table import Table
from hail.matrixtable import MatrixTable

## Add my personal module
sys.path.append('/home/yung-chun/2_python_scripts/my_hail_module')
import my_hail_table_list as myhl

#############
# DECORATOR #
#############
#######################################
## Decorator: Calculate running time function
def calculate_running_time(func):
    """
    Calculate the running time for a function

    :param fun: function
    :return: elapse time (hour/minute/second)
    """
    def wrapper(*args, **kwargs):
        """
        wrapper the function of calculate_running_time
        """
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        hours, rem = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(rem, 60)
        print(f"Running time for function {func.__name__}: {int(hours)}h {int(minutes)}m {seconds:.2f}s")
        return result
    return wrapper

## Decorator: Reset the hail table key and then set back
### ///It seems to be useless since it will make it running longer than expected
def set_ht_key_and_reset(func):
    """
    Set the hail key to 'locus', and 'alleles'.
    And then reset the key back to original setting after the fucntion.
    
    :param fun: function
    :return: hail table
    """
    def wrapper(ht):
        """
        wrapper the function of reset table key
        """
         # Get the current key information
        key_info = ht.key.keys()
        
        # Change the key to 'locus' and 'alleles'
        ht = ht.key_by('locus', 'alleles')
        
        # Call the decorated function
        modified_ht = func(ht)
        
        # Revert the key to the original setting
        modified_ht = modified_ht.key_by(*key_info)
        
        return modified_ht   
    return wrapper

#######################################

#############
# FUNCTIONS #
#############
#######################################
## Hail initiate
def hail_initialize(jobname: str = None, mem: int = 32):
    """
    To initialize the hail environment
    """
    today = datetime.now().strftime("%Y-%m-%d")
    print("Program Starting Date: {}".format(today))

    HAIL_HOME = subprocess.getoutput("pip3 show hail | grep Location | awk -F' ' '{print $2 \"/hail\"}'")
    os.environ["HAIL_DIR"] = f"{HAIL_HOME}/backend"
    # os.environ["HAIL_DIR"] = "/opt/conda/lib/python3.7/site-packages/hail/backend"
    # os.environ["JAVA_HOME"] = "/usr/lib/jvm/java-11-openjdk-amd64"
    os.environ["JAVA_HOME"] = "/usr/lib/jvm/java-8-openjdk-amd64/jre"

    threads = int(os.environ['LSB_MAX_NUM_PROCESSORS'])
    print("Run with {} CPUs.".format(threads))

    # Ensure hail-all-spark.jar exists
    hail_jars = f"{HAIL_HOME}/backend/hail-all-spark.jar"
    if not os.path.exists(hail_jars):
        raise FileNotFoundError(f"Cannot find hail-all-spark.jar at {hail_jars}")
    print(f"Found hail-all-spark.jar at: {hail_jars}")

    # Convert 'mem' to integer if it's not already
    mem = int(mem)

    # Calculate memory settings
    spare_memory = 6 #int(mem * 0.25) + 1
    driver_memory = int(mem * 0.25)  # Use 25% of total memory for driver
    executor_memory = mem - driver_memory - spare_memory  # Subtract driver memory and overhead for executors

    # Spark Configuration
    conf = pyspark.SparkConf().setAll([
        ('spark.master', 'local[{}]'.format(threads)),
        ('spark.app.name', 'Hail'),
        ('spark.jars', f'file://{hail_jars}'),
        ('spark.driver.extraClassPath', f'file://{hail_jars}'),
        ('spark.executor.extraClassPath', f'file://{hail_jars}'),   #'./hail-all-spark.jar'
        ('spark.serializer', 'org.apache.spark.serializer.KryoSerializer'),
        ('spark.kryo.registrator', 'is.hail.kryo.HailKryoRegistrator'),
        ### https://discuss.hail.is/t/turning-run-combiner-performance-for-hail-local-mode/2318
        ('spark.driver.memory', '{}g'.format(driver_memory)),
        ('spark.executor.memory', '{}g'.format(executor_memory)),
        ('spark.driver.maxResultSize', '{}g'.format(driver_memory)),
        ### For IOException: No space left on device
        # https://spark.apache.org/docs/2.3.0/configuration.html
        ('spark.local.dir', '/storage1/fs1/jin810/Active/testing/yung-chun/hail_tmp'),
        ### Error: KryoException: Buffer overflow. Available: 0, required: 15069596
        ('spark.kryoserializer.buffer', '64m'),
        ('spark.kryoserializer.buffer.max', '1024m'),
        ('spark.driver.memoryOverhead', '512m'),
        ('spark.executor.memoryOverhead', '1g'),
        ### To output some essential info for memory leakage (Garbage collection logging)
        ('spark.executor.extraJavaOptions', '-XX:+UseParallelGC -XX:+PrintGCDetails -XX:+PrintGCTimeStamps ' +
         '-XX:+PrintGCDateStamps -XX:+UseGCLogFileRotation -XX:NumberOfGCLogFiles=10 -XX:GCLogFileSize=100M ' +
         '-XX:ConcGCThreads={}'.format(threads)),
    ])

    # Initialize Spark Session
    spark = SparkSession.builder.config(conf=conf).getOrCreate()
    sc = spark.sparkContext

    # Log file setup
    log_base_path = "/storage1/fs1/jin810/Active/testing/yung-chun/Hail_logs"
    log_file_name = "hail_vds_new_combiner_{}.log".format(jobname if jobname else today)
    logfile = os.path.join(log_base_path, log_file_name)
    logfile = create_unique_filename(logfile)

    # Hail Initialization
    hl.init(default_reference='GRCh38', sc=sc, log=logfile)
    # hl.init(master='local[16]') -> to use all the CPU

    return spark, sc

## Get the num of cores for repartition the hail table and matrix table
def repartition_based_on_cores(data: hl.MatrixTable or hl.Table, partitions: int = None) -> hl.MatrixTable or hl.Table:
    """
    Repartition a Hail Table or MatrixTable. 
    If partitions is not provided, it uses twice the number of CPU cores.

    :param data: (hl.MatrixTable or hl.Table): Hail object to repartition.
    :param partitions: (int, optional): Number of partitions. Default to None.

    :return: hl.MatrixTable or hl.Table: Repartitioned Hail object.
    """
    
    if not partitions:
        num_cores = multiprocessing.cpu_count()
        partitions = num_cores * 2

    if isinstance(data, (hl.Table, hl.MatrixTable)):
        return data.repartition(partitions)
    else:
        raise ValueError("Provided data must be either a Hail Table or MatrixTable.")

## Check whether a file exists and return boolean
### From https://broadinstitute.github.io/gnomad_methods/_modules/gnomad/utils/file_utils.html#check_file_exists_raise_error
def file_exists(fname: str) -> bool:
    """
    Check whether a file exists.

    Supports either local or Google cloud (gs://) paths.
    If the file is a Hail file (.ht, .mt, .bm, .parquet, .he, and .vds extensions), it
    checks that _SUCCESS is present.

    :param fname: File name.
    :return: Whether the file exists.
    """
    fext = os.path.splitext(fname)[1]
    if fext in {".ht", ".mt", ".bm", ".parquet", ".he"}:
        paths = [f"{fname}/_SUCCESS"]
    elif fext == ".vds":
        paths = [f"{fname}/reference_data/_SUCCESS", f"{fname}/variant_data/_SUCCESS"]
    else:
        paths = [fname]

    if fname.startswith("gs://"):
        exists_func = hl.hadoop_exists
    else:
        exists_func = os.path.isfile

    exists = all([exists_func(p) for p in paths])

    return exists

## Check existence of files and return the path
def hail_file_exists(fname: str):
    # https://github.com/broadinstitute/gnomad_methods/blob/382fc2c7976d58cc8983cc4c9f0df5d8d5f9fae3/gnomad/utils/file_utils.py#L103
    """
    Check the input filename and see if it is valid and exists

    :param fname: File name.
    :return: input path.
    """
    path = str(fname)
    if file_exists(path):
        return path
    else:
        msg = "IT's NOT a valid Input format"
        raise argparse.ArgumentTypeError(msg)

## Check the the provided input is an directory
def is_directory_check(path: str):
    if not os.path.isdir(path):
        msg = "The provided input is not a directory or not exist."
        raise argparse.ArgumentTypeError(msg)
    else:
        return path

## Re-assign the input and output files for the step by step process in hail
def process_input_output(input_file, output_file=None):
    """
    Processes the input and output file names. 
    
    If an output_file is provided, its value will be assigned to the input_file, 
    and the output_file variable will be cleared (set to None). 
    If no output_file is provided, the function simply returns the original input_file 
    and a None value for output_file.
    
    Parameters:
    - input_file (str): The name of the input file.
    - output_file (str, optional): The name of the output file. Defaults to None.
    
    Returns:
    - tuple: A tuple containing two values: 
             1. The updated input file name (or the original name if no output_file was provided).
             2. None, since the output_file variable is cleared.
    """

    # If output_file is provided, assign its value to input_file
    if output_file is not None:
        input_file = output_file
        output_file = None  # Clear the output_file variable

    return input_file, output_file

## Generate a list of files from the given directory and pattern
def get_file_list_from_dir_pattern(path: str, pattern: str=None):
    """
    Generate a list of the file path from the directory and pattern provided

    :param path: File path.
    :param pattern: Pattern of the file. Default: *
    :return: list of file paths.
    """
    if pattern is None:
        pattern = "*"
    else:
        pattern = "*" + pattern

    file_list = glob.glob(path + '/' + pattern)

    return file_list

## Generate a list of file path from the given file
def get_file_list_from_file(filename: str):
    """
    Generate a list contianing the file paths from a file

    :param filename: /path/to/the/file
    :return: list of file paths.
    """
    file = open(filename, "r")
    data = file.readlines()
    #print("\nData of the File in List:", data)
    
    new_list = []
    for item in data:
        new_list.append(item.replace("\n", ""))
    #print("New List:", new_list)
    
    file.close()

    return new_list

## Check output file name and return the path
def output_filename_check(fname: str):
    """
    Check the output filename and see if the folder exists

    :param fname: File name.
    :return: output path.
    """
    fullname = str(fname)
    folder = os.path.dirname(fullname)
    filename = os.path.basename(fullname)
    if not os.path.isdir(folder):
        msg = "The output folder does NOT EXIST! Please check it manually."
        raise argparse.ArgumentTypeError(msg)
    elif (not filename) or (filename == ""):
        msg = "Your output file name IS BLANK! Please check it manually."
        raise argparse.ArgumentTypeError(msg)
    else:
        return fullname

## Remove the files
def remove_file(file: str):
    """
    Remove the files that provided

    :param file: file path
    
    """
    try:
        os.remove(file)
        print(f"Successfully removed the temporary file: {temp_file_path}")
    except OSError as e:
        print(f"Error deleting the file: {temp_file_path} - {e}")

def remove_item(item_path: str):
    """
    Remove the specified file or directory.
    
    Args:
    - item_path (str): Path to the file or directory to be removed.

    Example usage:
    - remove_item("/path/to/file_or_directory")
    """
    
    # Check if the provided path exists
    if not os.path.exists(item_path):
        print("Specified path does not exist.")
        return
    
    # Remove file or directory based on its type
    if os.path.isfile(item_path):
        os.remove(item_path)
        print(f"File '{item_path}' removed.")
    elif os.path.isdir(item_path):
        shutil.rmtree(item_path)
        print(f"Directory '{item_path}' and its contents removed.")
    else:
        print("Unknown item type. Unable to remove.")

# Load the files based on the different format input
def load_data(filepath: str) -> Union[hl.Table, hl.MatrixTable]:
    """
    Load data from a file into a Hail table or matrix based on file extension.
    
    :param filepath: path for the files, it could be ht, mt, vcf.bgz, vcf.gz or txt.
    :return: hail table or matrix table
    """
    if filepath.endswith('.ht'):
        return hl.read_table(filepath)
    elif filepath.endswith('.mt'):
        return hl.read_matrix_table(filepath)
    elif filepath.endswith(('.vcf.bgz', '.vcf.gz')):
        return hl.import_vcf(filepath, force_bgz=True)
    elif filepath.endswith('.txt'):
        # Handle .txt import here, perhaps with hl.import_table()
        return hl.import_table(filepath)
    else:
        raise ValueError(f"Unsupported file type: {filepath}")

## Get the list of unique values form the selected column
def count_unique_columnID(data,name=None,axis="col"):
    """
    Count the unique values from the selected column.

    :param data: could be hail table or matrix table
    :param name: name of the column to count unique values for (default: "ProbandID" for tables, "s" for matrices)
    :param axis: axis along which to aggregate (default: "row", can be "row" or "col")
    :return: list of unique values, and the count for the unique values.
    """
    if isinstance(data, Table):
        if name is None:
            name = "ProbandID"
        col_expr = hl.str(data[name])
        sample_lst = data.aggregate(hl.agg.collect_as_set(col_expr))
    elif isinstance(data, MatrixTable):
        if name is None:
            name =  "s"
        col_expr = hl.str(data[name])
            
        if axis == "row":
            sample_lst = data.aggregate_rows(hl.agg.collect_as_set(col_expr))
        elif axis == "col":
            sample_lst = data.aggregate_cols(hl.agg.collect_as_set(col_expr))
        else:
            raise ValueError("axis must be 'row' or 'col'")
    else:
        raise ValueError("Input data must be a Hail Table or MatrixTable")

    num = len(sample_lst)
    
    return sample_lst, hl.eval(num)

## Make directory if it doesn't exist
def make_dir_if_not_exist(_dir: str):
    """
    Create directory if it doen't exist
    
    :param _dir: must be a string
    """
    if not isinstance(_dir, str):
        raise ValueError("Directory must be a string")
      
    # Check whether the specified path exists or not
    _is_exist = os.path.exists(_dir)

    # Create a new directory because it does not exist
    if not _is_exist:
        os.makedirs(_dir)

## Create the unique name for the file
def create_unique_filename(filename: str) -> str:
    """
    Generate unique file name. 
    If the file already exists, it will add suffix.
    
    :param filename: must be a string

    Note:
        Filename: Date_analysis_result.txt
        Newname: Date_analysis_result_001.txt
    """
    if not isinstance(filename, str):
        raise TypeError("Input must be a string.")
        
    base, ext = os.path.splitext(filename)
    i = 0
    while os.path.exists(filename):
        filename = f"{base}_{str(i).zfill(3)}{ext}"
        i += 1
    return filename

## Generate outputname based on the input_file and pattern
def create_output_name(input_file: str, 
                       pattern: str, 
                       extension: str = None, 
                       output_path: str = None,
                       unique: bool = False) -> str:
    """
    Generate the output file name based on the input file name and pattern provided

    :param input_file: /path/to/file
    :param pattern: pattern want to add at the tail of file name beofre extention (ex: out_forplotread, pass_qc, etc)
    :param extention:
    :param output_path:
    :return: output file name

    Note:
         Will create the unique file name if it is already existed.
    """
    basename = os.path.basename(input_file)

    if output_path is not None:
        input_file = os.path.join(output_path, basename)
    else:
        input_file = input_file

    num_ext = len(input_file.split('.'))
    if num_ext > 1:
        # input file has an extension
        base, ext = os.path.splitext(input_file)
        for i in range(num_ext-2):
            base, additional_ext = os.path.splitext(base)
            ext = additional_ext + ext
        if extension is None:
            output = f"{base}_{pattern}{ext}"
        else:
            output = f"{base}_{pattern}{extension}"
    else:
        # input file has no extension
        if extension is None:
            output = f"{input_file}_{pattern}"
        else:
            output = f"{input_file}_{pattern}{extension}"

    # Get unique outputname
    if unique is True:
        output = create_unique_filename(output)

    return output

## Save the output to either txt, ht, or mt
@calculate_running_time
def out_save(dataset, out_file):
    """
    Save the files based on the input and the output extentions

    :param dataset: could be hail table or matrix table
    :param out_file: must have extention {.ht, .txt, .mt, }
    
    Note:
        Hail table can only be saved to 'txt' or 'ht'.
        Matrix table can only be save to 'mt'.
    """
    ext = out_file.rsplit(".", 1)[-1]
    
    if ext not in {'txt', 'bed', 'ht', 'mt'}:
        raise ValueError("Output file must have extension '.txt', '.bed', '.ht', or '.mt'")

    out_dir=os.path.dirname(out_file)
    make_dir_if_not_exist(out_dir)
    out_file = create_unique_filename(out_file)

    if isinstance(dataset, hl.Table):
        print('input is hail table')
        if ext == "txt" or ext == "bed":
            print('save to txt: ', out_file)
            dataset.export(out_file, delimiter='\t')
        elif ext == "ht":
            print('save to ht', out_file)
            dataset.write(out_file)
        else:
            raise ValueError("Output file must have extension '.txt', or '.ht'")
    elif isinstance(dataset, hl.MatrixTable):
        print('input is hail matrix table')
        if ext == "mt":
            print('save to mt', out_file)
            dataset.write(out_file)
        else:
            raise ValueError("Output file must have extension '.mt'")
    else:
        raise ValueError("Input dataset must be a Hail Table or MatrixTable object")

## Modify the colname to be albe use out_plotread
### "ID","ProbandID","#CHROM","POS","Proband_GT"
@typecheck(ht_denovo=hl.Table)
def modify_hldenovo_for_outplotread(ht_denovo):
    """
    Modify the hail table get from hl.denovo enable to be output for plotread

    :param_ht_denovo: hail table (output from hl.denovo)
    """
    # if not isinstance(ht_denovo, hl.Table):
    #     raise TypeError("ht should be hail table")

    ht_denovo = ht_denovo.annotate(
                                  ProbandID =  ht_denovo.proband.s
                                , MotherID = ht_denovo.mother.s
                                , FatherID = ht_denovo.father.s
                                , Proband_GT = hl.str(ht_denovo.proband_entry.GT)
                                , Mother_GT = hl.str(ht_denovo.mother_entry.GT)
                                , Father_GT = hl.str(ht_denovo.father_entry.GT)
                                , Proband_DP_int = ht_denovo.proband_entry.AD
                                , Mother_DP_int = ht_denovo.mother_entry.AD
                                , Father_DP_int = ht_denovo.father_entry.AD
                                , Proband_VAF = (ht_denovo.proband_entry.AD[1] / ht_denovo.proband_entry.DP)
                                , Mother_VAF = (ht_denovo.mother_entry.AD[1] / ht_denovo.mother_entry.DP)
                                , Father_VAF = (ht_denovo.father_entry.AD[1] / ht_denovo.father_entry.DP)
                                , Proband_GQ = ht_denovo.proband_entry.GQ
                                , Mother_GQ = ht_denovo.mother_entry.GQ
                                , Father_GQ = ht_denovo.father_entry.GQ
                                , CHROM = ht_denovo.locus.contig
                                , POS = ht_denovo.locus.position
                                , REF = ht_denovo.alleles[0]
                                , ALT = ht_denovo.alleles[1]
    )

    ht_denovo = ht_denovo.annotate(
                                  ID = hl.str('_').join(hl.array([ht_denovo.CHROM, hl.str(ht_denovo.POS), ht_denovo.alleles[0], ht_denovo.alleles[1]]))
    )
    
    ht_denovo = ht_denovo.rename({'CHROM':'#CHROM'})
    
    return ht_denovo

## out for plotread version 2
@calculate_running_time
@typecheck(ht=hl.Table,
           out_file=nullable(str),  # out_file can be None or a string
           select_all=bool,         # select_all is a boolean
           cols=nullable(oneof(str, tuple)),  # cols can be None, a string, or a tuple
           mode=str) # mode=plotread or deapsea
def out_for_downstream(ht: hl.Table, out_file=None, select_all=False, cols=None, mode='plotread'):
    """
    Save the files for different downstream purpose. 
    - plotread (ID, ProbandID, Chr, Pos, Proband_GT)
    - annovar (Chr, Pos, ID, Ref, Alt, Qual, Filter, AC)
    - deapsea (Chr, Pos, ID, Ref, Alt)

    :param ht: hail table (read from variantTable)
    :param out_file: output files, should be .txt extension
    :param select_all: if True, select all columns from the Hail table
    :param cols: add more column ID want to keep (expected string or tuple)
    :param mode: 'plotread', 'deepsea', 'annovar'. Default='plotread'
    """
    # Step 1: Set the default key to locus and alleles
    ht = ht.key_by('locus', 'alleles')

    # Step 2: Check if 'ID' exists in the table, if not, create it as 'Chr_Pos_Ref_Alt'
    if 'ID' not in ht.row:
        print("'ID' column is missing. Annotating 'ID' as 'Chr_Pos_Ref_Alt'.")
        required_cols = ['#CHROM', 'POS', 'REF', 'ALT']
        missing_cols = [col for col in required_cols if col not in ht.row]
        if missing_cols:
            raise ValueError(f"Missing required columns for ID generation: {missing_cols}")
        ht = ht.annotate(ID=hl.str(ht['#CHROM']) + "_" + hl.str(ht['POS']) + "_" + hl.str(ht['REF']) + "_" + hl.str(ht['ALT']))

    # Step 3: Key columns
    key_cols = ["locus", "alleles"]

    # Step 4: Define default columns for each mode
    if mode == 'deapsea':
        default_cols = ["#CHROM", "POS", "ID", "REF", "ALT"]
    elif mode == 'annovar':
        # Annotate QUAL, FILTER with ".", AC with "1", FORMAT, and Sample
        ht = ht.annotate(QUAL=".", FILTER=".", AC="1", FORMAT="GT:AD:DP:GQ:PL", Sample="0/0:0,0:0:0:0,0,0")
        default_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "AC", "FORMAT", "Sample"]
    else:  # plotread
        default_cols = ["ID", "ProbandID", "#CHROM", "POS", "Proband_GT"]

    # Step 5: Select the default columns
    _ht = ht.select(*default_cols)

    # Step 6: If select_all is True, select all additional columns from ht
    if select_all:
        print("Select all the columns.")
        all_cols = [field for field in ht.row_value]
        additional_cols = [col for col in all_cols if col not in default_cols and col not in key_cols]
        _ht = _ht.key_by(*ht.key).join(ht.select(*additional_cols), how='left')

    # Step 7: If specific columns are provided, join only those columns
    elif cols is not None:
        print(f"Only select these columns: {cols}")
        if isinstance(cols, str):  # Convert string to list if only one column is provided
            cols = [cols]
        _ht = _ht.key_by(*ht.key).join(ht.select(*cols), how='left')

    # Step 8: Handle duplicate removal and re-keying
    schema = _ht.row
    has_ProbandID = 'ProbandID' in schema
    if has_ProbandID:
        print("Re-keying by 'ProbandID' and 'ID'.")
        _ht = _ht.key_by('ProbandID', 'ID')
    else:
        print("Re-keying by 'ID'.")
        _ht = _ht.key_by('ID')

    _ht = _ht.distinct()

    # Step 9: Remove the key from the table and drop key columns
    _ht = _ht.key_by()
    _ht = _ht.drop(*key_cols)

    # Step 10: Save or return the table
    if out_file is not None:
        print(f"Saving table to {out_file}")
        out_save(_ht, out_file)
    else:
        return _ht

## Process the deepsea count file
def processed_deepsea_count_file(file_path: str, output: str = 'df') -> Union[str, pd.DataFrame, hl.Table]:
    """
    Processes a DeepSEA count file, finding the first matching column for each row based on the absolute difference 
    specified in the "seqclass_max_absdiff" column, and outputs the result as a CSV file, pandas DataFrame, or Hail Table.

    :param file_path: Path to the input TSV file.
    :param output: Specifies the output format. Options are 'tsv' for TSV file output, 'df' for returning a pandas DataFrame,
                   and 'ht' for returning a Hail Table. Default is 'df'.
    :return: Depending on the output parameter, returns the path to the output CSV file (if output='tsv'), 
             a pandas DataFrame (if output='df'), or a Hail Table (if output='hail').

    Note:
         The function adds a new column 'matched_class' to the output, which contains the name of the first column 
         where the absolute value matches that of the "seqclass_max_absdiff" column, for each row.
         If the 'tsv' option is selected, the output file name is generated by adding "_processed" before the file 
         extension of the input file name. For unique file generation, ensure the 'unique' flag is handled 
         externally if required.
    
    Usage:
        deepsea_file = 'path_to_your_file.tsv'
        output_format = 'tsv'  # Options: 'tsv', 'df', 'ht'
        result = processed_deepsea_count_file(deepsea_file, output=output_format)

    """
    # Read the TSV file into a DataFrame
    df = pd.read_csv(file_path, sep='\t')

    def find_matching_column(row):
        # Convert to float for comparison, handling potential conversion errors
        try:
            match_value = abs(float(row["seqclass_max_absdiff"]))
        except ValueError:
            return None  # Return None if conversion fails
        
        for col in df.columns[1:]:  # Adjusted to skip specific columns if needed
            if df[col].dtype in ['float64', 'float32']:
                try:
                    if abs(float(row[col])) == match_value:
                        return col
                except ValueError:
                    continue

        return None

    df['matched_class'] = df.apply(find_matching_column, axis=1)

    if output == 'tsv':
        # Generate the output file path by adding "_processed" before the file extension
        output_file_path = create_output_name(input_file=file_path, pattern="processed", unique=False)

        # Write the modified DataFrame to a new CSV file
        print(f"Processed file saved to: {output_file_path}")
        df.to_csv(output_file_path, sep='\t', index=False)
        return None
    elif output == 'ht':
        # Convert the pandas DataFrame to Hail Table
        #hl.init()  # Initialize Hail
        ht = hl.Table.from_pandas(df).key_by('name')
        return ht
    else:
        # Return the pandas DataFrame
        return df

## Filter the hail table or matrix table based on the provided criteria
@calculate_running_time
def data_filter_and_select(data: Union[str, hl.Table, hl.MatrixTable], columns: str = None, values: str = None, bed: str = None, reference_genome: str = 'GRCh38', reverse: bool = False) -> Union[hl.Table, hl.MatrixTable]:
    """
    Filter rows based on specified column-value pairs and/or bed region.

    :param data: Path or Hail table/matrix.
    :param columns: Comma-separated string of column names (optional).
    :param values: Comma-separated string of values corresponding to each column (optional).
    :param bed: Path for the hail table containing bed region (optional).
    :param reference_genome: Reference genome identifier. Default is 'GRCh38'.
    :param reverse: If True, exclude regions specified in the bed file.
    :return: Filtered table or matrix.
    """
    ref_genome = hl.get_reference(reference_genome)
    
    # Load data if filepath provided
    if isinstance(data, str):
        data = load_data(data)

    #################
    # Column-Values #
    #################
    # Filter based on column-value pairs
    if columns and values:
        column_list = columns.split(',')
        value_list = values.split(',')

        if len(column_list) != len(value_list):
            raise ValueError("Number of columns must match number of values.")
        
        conditions = []
        for column, value in zip(column_list, value_list):
            if isinstance(data, hl.Table):
                conditions.append(data[column] == hl.literal(value))
            elif isinstance(data, hl.MatrixTable):
                conditions.append(data.row[column] == hl.literal(value))

        combined_condition = hl.all(conditions)
        if reverse:
            combined_condition = ~combined_condition

        if isinstance(data, hl.Table):
            data = data.filter(combined_condition)
        elif isinstance(data, hl.MatrixTable):
            data = data.filter_rows(combined_condition)

    #######
    # BED #
    #######
    # Add locus column if necessary and it doesn't exist
    if bed and isinstance(data, hl.Table) and 'locus' not in data.row:
        data = data.annotate(locus=hl.locus(data['#CHROM'], hl.int(data.POS), reference_genome=ref_genome))

    # Slice into bed region
    if bed:
        bed_ht = hl.read_table(bed)
        locus_condition = hl.is_defined(bed_ht[data.locus])
        if reverse:
            locus_condition = ~locus_condition

        if isinstance(data, hl.Table):
            data = data.filter(locus_condition)
        elif isinstance(data, hl.MatrixTable):
            data = data.filter_rows(locus_condition)

    return data

def merge_ht_tables(ht1, ht2):
    """
    Merges two Hail Tables with distinct rows and different column sets.
    Aligns the schema of both tables by adding missing columns with missing values,
    then concatenates the rows of both tables.

    Parameters:
    ht1 (hail.Table): The first Hail Table to merge.
    ht2 (hail.Table): The second Hail Table to merge.

    Returns:
    hail.Table: A new Hail Table containing all rows and columns from both input tables.
                Columns unique to one table are filled with missing values in rows from the other table.
    """
    # Get the column fields from both tables
    columns_ht1 = list(ht1.row_value)
    columns_ht2 = list(ht2.row_value)

    # Combine and deduplicate the column names while preserving order
    combined_columns = list(dict.fromkeys(columns_ht1 + columns_ht2))

    # Add missing columns to each table
    for col in combined_columns:
        if col not in columns_ht1:
            ht1 = ht1.annotate(**{col: hl.missing(ht2[col].dtype)})
        if col not in columns_ht2:
            ht2 = ht2.annotate(**{col: hl.missing(ht1[col].dtype)})

    # Select columns in the determined order
    ht1_select = ht1.select(*combined_columns)
    ht2_select = ht2.select(*combined_columns)

    # Now, union the tables
    ht_combined = ht1_select.union(ht2_select)

    return ht_combined

## For annotation 
def annotate_hail_table(
    ht: hl.Table, 
    source_ht: hl.Table, 
    columns: list = 'all', 
    ht_key: list = ['locus', 'alleles'], 
    source_ht_key: list = ['locus', 'alleles']
) -> hl.Table:
    """
    Annotate a Hail table with information from a source Hail table using specified keys.
    
    :param ht: Hail table to annotate
    :param source_ht: Source Hail table containing the annotation information
    :param columns: List of columns to annotate from the source table (default is 'all', which means all columns)
    :param ht_key: List of key columns to use for the ht table (default is ['locus', 'alleles'])
    :param source_ht_key: List of key columns to use for the source_ht table (default is ['locus', 'alleles'])
    :return: Annotated Hail table with key reverted to ['locus', 'alleles']
    """
    
    # Step 1: Validate the keys exist in ht and source_ht
    for key in ht_key:
        if key not in ht.row:
            raise ValueError(f"Key '{key}' is not a field in the target Hail table (ht). Available fields: {list(ht.row)}")
    
    for key in source_ht_key:
        if key not in source_ht.row:
            raise ValueError(f"Key '{key}' is not a field in the source Hail table (source_ht). Available fields: {list(source_ht.row)}")
    
    # Step 2: Re-key ht and source_ht
    ht = ht.key_by(*ht_key)
    source_ht = source_ht.key_by(*source_ht_key)

    # Step 3: If 'all' is specified, use all columns from source table excluding keys
    if columns == 'all':
        columns = [col for col in source_ht.row if col not in source_ht.key]  # Exclude key columns
    
    # Step 4: Validate that the user-specified columns exist in the source table
    missing_columns = [col for col in columns if col not in source_ht.row]
    if missing_columns:
        raise ValueError(f"The following columns do not exist in the source table: {missing_columns}")
    
    # Step 5: Select only the specified columns from source_ht
    source_ht = source_ht.select(*columns)
    
    # Step 6: Annotate the main table (ht) with source_ht using the specified keys
    annotated_ht = ht.annotate(**source_ht[ht.key])
    
    # Step 7: Revert the key for ht back to ['locus', 'alleles'] (only if 'locus' and 'alleles' exist)
    if 'locus' in annotated_ht.row and 'alleles' in annotated_ht.row:
        annotated_ht = annotated_ht.key_by('locus', 'alleles')
    else:
        print(f"Warning: 'locus' and 'alleles' keys were not found in the table after annotation, so the key was not reset.")

    return annotated_ht

########
# TEMP #
########
###################
## Fix the rsid column (matrix table)

def re_annotate_rsid(data):
    if isinstance(data, hl.Table):
        data = data.annotate( ID = hl.str('_').join([data.locus.contig, hl.str(data.locus.position), data.alleles[0], data.alleles[1]]) )
    elif isinstance(data, hl.MatrixTable):
        data = data.annotate_rows( ID = hl.str('_').join([data.locus.contig, hl.str(data.locus.position), data.alleles[0], data.alleles[1]]))
    
    return data