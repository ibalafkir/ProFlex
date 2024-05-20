"""
Processes a PDB to adjust it to the pipeline: 
- Keeps desired chains
- Fixes insertion codes
- Removes HETATOMS (ligands, metals, waters...)
- Removes other unnecessary lines

A PDB without any structural gaps is expected
"""

import argparse
from proflex.pdb import PdbHandler
import os


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Creates a new PDB file fixed for insertion codes')
    parser.add_argument('-pdb', type=str, help='Path to the PDB file')
    parser.add_argument('-ab1', type=str, help='Antibody chain 1 or nanobody single chain')
    parser.add_argument('-ab2', type=str, help='Antibody chain 2', required=False)
    parser.add_argument('-ag', type=str, help='String of chains to keep')
    args = parser.parse_args()
    
    pdb_file = args.pdb
    ab1 = args.ab1
    ab2 = args.ab2
    ag = args.ag

    chains_to_keep = []
    if ab2:
        chains_to_keep = [ab1, ab2, ag]
        print("Keeping", chains_to_keep)
    else:
        chains_to_keep = [ab1, ag]
        print("Keeping", chains_to_keep)


    print("Saving atom lines from wanted chains and tidying...")
    PdbHandler.chain(pdb_file,chains_to_keep)
    PdbHandler.tidy(pdb_file[:-4]+'_ch.pdb')

    print("Fixing insertions...")
    PdbHandler.fix_insertions(pdb_file[:-4]+'_ch_tidied.pdb')
    PdbHandler.fix_ter_mistakes(pdb_file[:-4]+'_ch_tidied_insfixed.pdb')
    PdbHandler.tidy(pdb_file[:-4]+'_ch_tidied_insfixed_terfixed.pdb')
    
    """
    # NOT NECESSARY FOR NOW, DUE TO LINE 30 USE OF BIOPANDAS
    print("Removing HETATOMS...") 
    PdbHandler.remove_hetatm(pdb_file[:-4]+'_ch_insfixed'+'_terfixed.pdb')
    print("Removing unnecessary lines...")
    PdbHandler.custom_remove(pdb_file[:-4]+'_ch_insfixed'+'_terfixed_nohetatm.pdb', ['CRYST1', 'ANISOU', 'CONECT', 'MASTER' 'HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS', 'EXPDATA', 'AUTHOR', 'REVDAT', 'JRNL', 'REMARK'])
    # TODO more lines to be removed might be added soon

    """
    pdb_name = pdb_file[:-4]
    pdb_name_ch = pdb_file[:-4]+'_ch.pdb'
    pdb_name_ch_tidied = pdb_name_ch[:-4]+'_tidied.pdb'
    pdb_name_ch_tidied_insfixed = pdb_name_ch_tidied[:-4] + '_insfixed.pdb'
    pdb_name_ch_tidied_insfixed_terfixed = pdb_name_ch_tidied_insfixed[:-4]+'_terfixed.pdb'
    pdb_name_ch_tidied_insfixed_terfixed_tidied = pdb_name_ch_tidied_insfixed_terfixed[:-4]+'_tidied.pdb'
    
    os.rename(pdb_name_ch_tidied_insfixed_terfixed_tidied, pdb_name+'_proc.pdb')
    
    os.remove(pdb_name_ch)
    os.remove(pdb_name_ch_tidied)
    os.remove(pdb_name_ch_tidied_insfixed)
    os.remove(pdb_name_ch_tidied_insfixed_terfixed)
    
    print(f"Output file: {pdb_name+'_proc.pdb'}")
