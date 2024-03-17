"""
Processes a PDB to adjust it to the pipeline: 
- Keeps desired chains
- Fixes insertion codes
- Removes HETATOMS (ligands, metals, waters...)
- Removes other unnecessary lines

A PDB without any gaps is expected (if a Maestro Schrodinger license is available, maestro_fill_gaps.sh is advised to be used) 
"""

import argparse
from proflex.utils import PDBProcessor
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Creates a new PDB file fixed for insertion codes'
        )
    parser.add_argument('-pdb', type=str, help='Path to the PDB file')
    parser.add_argument('-chains', type=str, help='String of chains to keep')
    
    args = parser.parse_args()
    
    pdb_file = args.pdb
    chid = args.chains
    chidd = chid.split(",")
   
    print("Saving wanted chains...")
    PDBProcessor.pdb_keepchain(pdb_file,chidd)
    
    print("Fixing insertions...\n")
    PDBProcessor.fix_insertions(pdb_file[:-4]+'_ch.pdb')
    PDBProcessor.fix_ter_mistakes(pdb_file[:-4]+'_ch_insfixed.pdb')
    
    print("Removing HETATOMS...\n") 
    PDBProcessor.remove_hetatm(pdb_file[:-4]+'_ch_insfixed'+'_terfixed.pdb')

    print("Removing unnecessary lines...\n")
    PDBProcessor.remove_lines(pdb_file[:-4]+'_ch_insfixed'+'_terfixed_nohetatm.pdb', ['CRYST1', 'ANISOU', 'CONECT', 'MASTER' 'HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS', 'EXPDATA', 'AUTHOR', 'REVDAT', 'JRNL', 'REMARK'])
    # TODO more lines to be removed might be added soon

    print("Deleting temporal files...\n")
    pdb_name = pdb_file[:-4]
    pdb_name_ch = pdb_file[:-4]+'_ch.pdb'
    pdb_name_insfixed = pdb_name + '_ch_insfixed.pdb'
    pdb_name_terfixed = pdb_name_insfixed[:-4]+'_terfixed.pdb'
    pdb_name_nohetatm = pdb_name_terfixed[:-4]+'_nohetatm.pdb'
    pdb_name_removed = pdb_name_nohetatm[:-4]+'_removed.pdb'

    os.rename(pdb_name_removed, pdb_name+'_processed.pdb')
    
    os.remove(pdb_name_ch)
    os.remove(pdb_name_insfixed)
    os.remove(pdb_name_nohetatm)
    os.remove(pdb_name_terfixed)
   
    wd = os.getcwd()
    
    print(f"Output file in {wd} with the name {pdb_name+'_processed.pdb'}\n")
