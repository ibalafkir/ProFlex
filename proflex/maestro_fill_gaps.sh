#!/bin/bash

# Run MAESTRO Schrodinger protein preparation wizard to fill residue loops
# Ideally provide an input PDB with gaps that have previous residue information:
# if existing, this information is always kept when PDBs are downloaded from
# https://files.rcsb.org/download/XXXX.pdb' where XXXX stands for the 4-letter PDB code

# The pipeline needs a PDB w/o gaps (blame RFdiffusion)
# When generating an ub structure, use this script in ab and ag per sepparate and then superpose with 
# the bound structure and minimize in the ppw GUI
# After this, use PDBProcessor



if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <SCH_PPW_PATH> <INPUT_PDB_PATH>"
    echo "Please provide the following arguments in the correct order:"
    echo "1. Path to Schrodinger Protein Preparation Wizard (SCH_PPW_PATH)"
    echo "2. Path to input PDB file (INPUT_PDB_PATH)"
    exit 1
fi

SCH_PPW_PATH="$1" # Path to Schrodinger Protein Preparation Wizard
INPUT_PDB_PATH="$2" # Path to input PDB file
OUTPUT_PDB_PATH="${INPUT_PDB_PATH%.*}_fill.pdb" # Output path for processed PDB file

"$SCH_PPW_PATH" -fillloops "$INPUT_PDB_PATH" "$OUTPUT_PDB_PATH"

echo "Job sent to $SCH_PPW_PATH, .log and output files expected in $PWD"
