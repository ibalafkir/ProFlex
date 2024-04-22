from proflex.utils import *
import argparse
from biopandas.pdb import PandasPdb
import os

"""
Includes all atom lines of a PDB in the same chain identifier (A) and renumbered so that the output file serves as a valid input to rmsd + chi angles analysis
through AttnPacker scripts
"""

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Remove side chains in a PDB')
    parser.add_argument('pdb', type=str, help='Path to the target PDB file (original)')  
    args = parser.parse_args()
    pdb = args.pdb
    
    # Opening as df
    atom_df = PDBUtils.get_pdb_atoms_df(pdb)
    atom_df['chain_id'] = 'A'
    # Writing new PDB
    new_pdb = PandasPdb().read_pdb(pdb)
    new_pdb.df['ATOM'] = atom_df
    new_pdb.to_pdb(path=pdb[:-4]+'_ch.pdb', records=['ATOM'], gz = False)    
    
    # Renumbering residues
    PDBProcessor.pdb_resnum(pdb[:-4]+'_ch.pdb', 1)
    
    # Renumbering atoms
    PDBProcessor.pdb_atomrenumber(pdb[:-4]+'_ch_renum.pdb', 1)
    
    # Tidying
    RFDFixer.pdb_tidying(pdb[:-4]+'_ch_renum_atomsorted.pdb')
    
    # Removing the temp file
    os.remove(pdb[:-4]+'_ch.pdb')
    os.remove(pdb[:-4]+'_ch_renum.pdb')
    os.remove(pdb[:-4]+'_ch_renum_atomsorted.pdb')
    os.rename(pdb[:-4]+'_ch_renum_atomsorted_tidied.pdb', pdb[:-4]+'_toassess.pdb')
    
    
    
    

    
    