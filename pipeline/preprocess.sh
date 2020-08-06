#!/usr/bin/env bash

input_csv=$1
output_dir=$2

git_path=${HOME}/projects/PhageLoader
process_bin="${git_path}/bin/preprocess"

if [ ! -f "$input_csv" ]; then
  echo "Error : input csv not found at path : ${input_csv}";
  exit -1;
fi


if [ -z "$output_dir" ]; then
  echo "Error : Output directory path required";
  exit -1;
fi

if [ -d "$output_dir" ]; then
  echo "Overwriting Output Directory : ${output_dir}";
  rm -rf $output_dir
fi

output_csv=${output_dir}/raw_counts.csv
output_arma=${output_dir}/raw_counts.arma
output_sample_names=${output_dir}/sample_names.txt
output_peptide_names=${output_dir}/peptide_names.txt

# create output directory
echo ""
echo "Creating Output Directory";
mkdir $output_dir;

# create sample name list
echo ""
echo "Creating Sample List";
time head -n 1 $input_csv | cut -d ',' -f 2- | tr "," "\n" > ${output_sample_names}

# create peptide name list
echo ""
echo "Creating Peptide List";
time cut -d ',' -f 1 $input_csv | tail -n +2 > ${output_peptide_names}

# create raw count matrix
echo ""
echo "Formatting Input Matrix";
time tail -n +2 $input_csv | cut -d ',' -f 2- > ${output_csv}

# create arma format matrix
echo ""
echo "Converting CSV -> ARMA"
time $process_bin -i ${output_csv} -o ${output_arma}
