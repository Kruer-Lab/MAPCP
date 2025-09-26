home_dir=$(pwd)
echo "The home directory name is: $home_dir" # Current working directory

input_dir=/path/to/directory # Specify the directory where trio sample folders are
echo "The input directory is: $input_dir"
echo
output_dir=${home_dir}/Outputs
if [ -d "$output_dir" ]; then
    echo "Directory exists: $output_dir"
else
    echo "Directory does not exist. Creating: $output_dir"
    mkdir -p "$output_dir"
fi


echo

# ********** Loop to extract relevant vcf file path from multiple sample folders **********
# Please be mindful that this script was written to parse the naming system that our lab uses. For combined trio/multiplex vcfs they are named and stored as the following example: F0000.trio/F00000.trio.raw.vcf
sample_dirs=${input_dir}/F* # Reads all folders with extension .trio. Change extension here if naming system is different.
sample_names=$(basename -a ${sample_dirs}) # Names of individual samples
path_name=()
for i in ${sample_dirs[@]};do
data_file=$(basename -a ${i})
dir_path=( $input_dir"/"${data_file}"/"${data_file}".trio.raw.vcf" ) # change file type here if naming system is diifferent
path_name=(${path_name[@]} ${dir_path}) # Append individual paths to variable
done
# *********************************** End of loop *****************************************

echo "Total number of sample vcf files: ${#path_name[@]}" # Shows total number of current vcf filepaths
echo
echo ${path_name[@]} > AutoMap_Samples.txt # Saves filepaths of individual vcf files in a text file

echo $home_dir
echo $input_dir
echo $output_dir
echo $path_name
echo $sample_names

# Working, Input, and Output directories
export home_dir
export input_dir
export output_dir
export path_name
export sample_names