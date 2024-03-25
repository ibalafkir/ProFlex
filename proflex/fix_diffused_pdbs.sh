#!/bin/bash

# Usage: runs a python script that corrects all diffused PDB models that stem from a PDB by using the last as reference

if [ $# -ne 1 ]; then
    echo "Usage: $0 solver python script and reference pdb are needed"
    exit 1
fi

solver=$1


for ref in *_proc.pdb; do

    start_time=$(date +%s)
    # diff="${ref%.*}"_diff_N*.pdb # To correct PDBs directly from RFdiffusion
    attdiffpack="${ref%.*}"_diff_N*attnpacked_cor.pdb # To correct PDB from RFd + AttnPacker (output of the last)

    for model in $attdiffpack; do # $diff for the first and $attdiffpack for the second
        python "$solver" "$ref" "$model"
    done
	
	end_time=$(date +%s)
	execution_time=$((end_time - start_time))
	echo "All models belonging to $ref were corrected in $execution_time seconds"
	echo ""
done





: '
ref=$2 

if [ ! -f "$ref" ]; then
    echo "'$ref' does not exist"
    exit 1
fi

diff="${ref%.*}"_diff_N*.pdb
start_time=$(date +%s)

for model in $diff; do
    python "$solver" "$ref" "$model"
done

end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "All models belonging to $ref were corrected in $execution_time seconds"

'
