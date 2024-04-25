#!/bin/bash

# Usage: runs a python script that corrects all diffused PDB models that stem from a PDB by using the last as reference

if [ $# -ne 1 ]; then
    echo "Usage: $0 solver python script and reference pdb are needed"
    exit 1
fi

solver=$1
files=$(ls *b.pdb && ls *ubnoT.pdb) # Reference PDBs, meaning the ones with which RFdiffusion was run and thus
                                    # the ones that serve as reference point for diffusion models correction
                                    # Edit this variable to change references e.g.: files=$(ls 6DBG_b.pdb)

for ref in $files; do

    start_time=$(date +%s)
    diff="${ref%.*}"_diff_N*.pdb # To correct PDBs directly from RFdiffusion
    for model in $diff; do # $diff for the first and $attdiffpack for the second
        python "$solver" "$ref" "$model"
    done
	end_time=$(date +%s)
	execution_time=$((end_time - start_time))
	echo "All models belonging to $ref were corrected in $execution_time seconds"
	echo ""
done