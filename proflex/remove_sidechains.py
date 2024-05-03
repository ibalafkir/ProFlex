"""
Side chain remover for Repacking method testing, will not be used in production of ProFlex    
"""

from proflex.pdb import PdbDf, PdbHandler
import argparse
from biopandas.pdb import PandasPdb
from rfdiffusion import RFdFix 
import os

def run(atom_df):
    atom_name_keep = ['N','CA', 'C', 'O']
    mask_delete_ch = (~atom_df['atom_name'].isin(atom_name_keep))
    atom_df_backbone = atom_df[~mask_delete_ch]
    return atom_df_backbone

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Remove side chains in a PDB')
    parser.add_argument('pdbfile', type=str, help='Path to the PDB file (input of RFD)')  
    args = parser.parse_args()
    pdb = args.pdbfile
    
    atom_df = PdbDf.atoms(pdb)
    atom_df_backbone = run(atom_df)
    
    output_pdb = PandasPdb().read_pdb(pdb)
    output_pdb.df['ATOM'] = atom_df_backbone
    output_pdb.to_pdb(path = pdb[:-4]+'_temp.pdb', records = ['ATOM'], gz=False) 
    
    PdbHandler.tidy(pdb[:-4]+'_temp.pdb')
    PdbHandler.renatom(pdb[:-4]+'_temp_tidied.pdb', 1)
    os.remove(pdb[:-4]+'_temp.pdb')
    os.remove(pdb[:-4]+'_temp_tidied.pdb')
    os.rename(pdb[:-4]+'_temp_tidied_atomsorted.pdb', pdb[:-4]+'_nosch.pdb')
    
    