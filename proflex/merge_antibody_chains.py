"""
Includes both antibody chains under the same chain_name

When DockQ is run over an antibody (H,L) -antigen (Ag) system in which 
the first has 2 chains, the stats are computed for the existing 
monomer-monomer interfaces (H-Ag, L-Ag, H-L)

This script merges the antibody chains so we get the stats for HL-A
"""
import argparse
from proflex.pdb import PdbDf, PdbHandler
from biopandas.pdb import PandasPdb
import pandas as pd
import os

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Includes both antibody chains under the same chain_name')
    parser.add_argument('-pdb', type=str, help='Path to the target PDB file')
    parser.add_argument('-ab_chs', type=str, help='Antibody chains to merge')
    parser.add_argument('-o_ch', type=str, help='Antibody final chain name')  
    args = parser.parse_args()
    pdb = args.pdb
    abchs = args.ab_chs
    och = args.o_ch
    ab_chs = abchs.split(",")
    
    # Merging chains
    atom_df = PdbDf.atoms(pdb)
    atom_df['chain_id'] = atom_df['chain_id'].replace({ab_chs[0]:och, ab_chs[1]: och})

    # Writing new PDB
    new_pdb = PandasPdb().read_pdb(pdb)
    new_pdb.df['ATOM'] = atom_df
    new_pdb.to_pdb(path=pdb[:-4]+'_ch.pdb', records=['ATOM'], gz = False)    
    
    # Renumbering residues
    PdbHandler.renres(pdb[:-4]+'_ch.pdb', 1)
    
    # Renumbering atoms
    # Warning it renumbers atoms since res_num 1
    PdbHandler.renatom(pdb[:-4]+'_ch_renum.pdb', 1)
    
    # Tidying
    PdbHandler.tidy(pdb[:-4]+'_ch_renum_atomsorted.pdb')
    
    # Removing the temp file
    os.remove(pdb[:-4]+'_ch.pdb')
    os.remove(pdb[:-4]+'_ch_renum.pdb')
    os.remove(pdb[:-4]+'_ch_renum_atomsorted.pdb')
    os.rename(pdb[:-4]+'_ch_renum_atomsorted_tidied.pdb', pdb[:-4]+'_toDockQ.pdb')