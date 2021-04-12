#!/bin/bash
set -uex

# This script can be used to download project data files
# Usage: $bash get_data.sh HSC RNASeq ENCFF247FEJ

# Script commant line arguments

# Cell Line (one of HSC, CMP, CFUE, Erythroblast)
cell_line=$1
# Experiment type (one of RNASeq or ATACSeq)
exp_type=$2
# Encode id to retrieve dataset
encode_id=$3
# Output Directory to store dataset
output_dir="${4:-"../raw/"}"
# creating file tag based on experiment type
if [[ $exp_type = "RNASeq" ]]
then
	echo Yes
	flag_type="tsv"
elif [[ $exp_type = "ATACSeq" ]]
then
    flag_type="bigBed"
else
    echo "Error: wrong experiment type"
    exit 1
fi
# setting the output directory prefix based on cell line
output_dir_prefix=$output_dir$cell_line
# making the output directory if it does not exist
mkdir -p $output_dir_prefix

# Downloading data using wget
echo $output_dir_prefix
wget https://www.encodeproject.org/files/$encode_id/@@download/$encode_id.$flag_type -P $output_dir_prefix
