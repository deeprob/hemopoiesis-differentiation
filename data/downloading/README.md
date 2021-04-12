# Downloading RNA-Seq and ATAC-Seq data

## ENCODE ID stored in text file: encodeIds.txt
<encodeIds.txt> file: consists of the IDs of the data files for each cell line required to download data from ENCODE portal.

## get_data.sh bash script can be used to retrieve datasets
Given the ids, we can download the files using the get_data bash script provided. The bash script requires three arguments:
1. Cell Line: Either one of HSC, CMP, CFUE, Erythroblast 
2. Experiment Type: Either one of RNASeq, ATACSeq
3. Encode ID: Encode ID can be obtained from encodeIds.txt file

The script also takes an optional argument:
1. output_dir: The output directory pre-prefix where the dataset will be stored. Defaults to "../raw/"

The script creates a directory path <output_dir> + <cell_line>, if required, where the downloaded file will be stored. The file is named by its encode id followed by it tag, ".tsv" in case of RNASeq data and ".bigBED" in case of ATACSeq data.

## Other relevant project related information provided Dr. Li
Cell Line Mapping:
(1=HSC, 2=CMP, 3=CFUE, 4=Erythroblast)

Cell Line Comparison Assigned:
2-4; 1-2; 2-3; 1-4; 3-4; 1-3;

**Work only with ScriptSeq Data**
