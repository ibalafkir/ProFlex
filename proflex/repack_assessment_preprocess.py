"""
Includes all atom lines of a PDB in the same chain identifier (A) and renumbered so that the output file serves as a valid input to rmsd + chi angles analysis
through AttnPacker scripts
"""

from proflex.pdb import PdbDf, PdbHandler
import argparse
from biopandas.pdb import PandasPdb
import os

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Remove side chains in a PDB')
    parser.add_argument('pdb', type=str, help='Path to the target PDB file (original)')  
    args = parser.parse_args()
    pdb = args.pdb
    
    # Opening as df
    atom_df = PdbDf.atoms(pdb)
    atom_df['chain_id'] = 'A'
    # Writing new PDB
    new_pdb = PandasPdb().read_pdb(pdb)
    new_pdb.df['ATOM'] = atom_df
    new_pdb.to_pdb(path=pdb[:-4]+'_ch.pdb', records=['ATOM'], gz = False)    
    
    # Renumbering residues
    PdbHandler.renres(pdb[:-4]+'_ch.pdb', 1)
    
    # Renumbering atoms
    PdbHandler.renatom(pdb[:-4]+'_ch_renum.pdb', 1)
    
    # Tidying
    PdbHandler.tidy(pdb[:-4]+'_ch_renum_atomsorted.pdb')
    
    # Removing the temp file
    os.remove(pdb[:-4]+'_ch.pdb')
    os.remove(pdb[:-4]+'_ch_renum.pdb')
    os.remove(pdb[:-4]+'_ch_renum_atomsorted.pdb')
    os.rename(pdb[:-4]+'_ch_renum_atomsorted_tidied.pdb', pdb[:-4]+'_toassess.pdb')
    
    
    
    

    
    