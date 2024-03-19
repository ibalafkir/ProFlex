"""
Solves the PDB format of an RFdiffusion PDB output (which is all back atoms of all residues in a single chain)
"""

from proflex.utils import RFDFixer, PDBUtils, PDBProcessor
import argparse
import pandas as pd
from biopandas.pdb import PandasPdb 
from rfdiffusion import RFDSchains 
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Fixes PDB features missed in the output of RFDiffusion'
        )
    parser.add_argument('prerfd', type=str, help='Path to the PDB file (input of RFD)')
    parser.add_argument('postrfd', type=str, help='Path to the output PDB RFDiffusion file')
    # parser.add_argument('-contigs', type=str, help='Contigs string')
    args = parser.parse_args()
    
    pre = args.prerfd
    post = args.postrfd
    # contigs = args.contigs

    # Fixing the .pdb format in the RFD PDB output thanks to the PDB input in RFD

    print(f"Correcting {post}...")
    RFDFixer.pdb_backbone(pre) 
    RFDFixer.pdb_atom(pre[:len(pre)-4]+'_backbone.pdb')
    RFDFixer.correct_rfd_pdbs(post, pre[:len(pre)-4]+'_backbone.pdb')
    
    print("Tidying...")
    RFDFixer.pdb_tidying(post[:len(post)-4]+'_chainsfixed.pdb')
    
    print("Deleting temporal files...")
    pdb_name = pre[:len(pre)-4]
    pdb_file_backbone = pdb_name+'_backbone.pdb'
    pdb_file_backbone_atom = pdb_file_backbone[:-4]+'_atom.pdb'
    pdb_rfd_name = post[:len(post)-4]
    pdb_rfd_file_chainsfixed = pdb_rfd_name+'_chainsfixed.pdb'
    pdb_rfd_file_chainsfixed_tidied = pdb_rfd_file_chainsfixed[:-4]+'_tidied.pdb'
    os.rename(pdb_rfd_file_chainsfixed_tidied, pdb_rfd_name+'_rfdfixed.pdb')
    os.remove(pdb_file_backbone)
    os.remove(pdb_file_backbone_atom)
    os.remove(pdb_rfd_file_chainsfixed)
    
    wd = os.getcwd()
    
    print(f"Output: {pdb_rfd_name+'_rfdfixed.pdb'}")
    
    
    
    """
    
  
    
    ## HERE ONWARDS SUPERPOSITION AND OTHER PROCEDURES ARE DONE TO KEEP ONLY PREDICTIONS OF FLEXIBLE RESIDUES
    print("K")
    
    
    
    
    # Superposing the fixed PDB output with the PDB input in RFD
    RFDFixer.superpose(pre, post[:len(post)-4]+'_rfdfixed.pdb', post[:len(post)-4]+'_superposed.pdb')
    
    # Reads df of the preRFD PDB and postRFD PDB fixed
    pre_atom = PDBUtils.get_pdb_atoms_df(pre)
    postfixed_atom = PDBUtils.get_pdb_atoms_df(post[:len(post)-4]+'_superposed.pdb')
    
    # Keeps side chains of the preRFD PDB df
    pre_atom_sidechains = pre_atom[~pre_atom['atom_name'].isin(['N', 'CA', 'C', 'O'])] 

    # Mixes RFD backbone df with preRFD PDB df side chains
    mixed_df = pd.concat([postfixed_atom, pre_atom_sidechains]) 

    # Generates a temporal PDB file out of the previous df
    output_pdb = PandasPdb().read_pdb(post)
    output_pdb.df['ATOM'] = mixed_df
    output_pdb.to_pdb(path = post[:len(post)-4]+'_disordered.pdb', records = ['ATOM'], gz=False)  
    
    # Sorts the PDB according to chain id and residue number
    RFDFixer.pdb_sorting(post[:len(post)-4]+'_disordered.pdb')

    # Tidies the previous PDB (add TERs, END...). This PDB has all side chains from the prePDB
    RFDFixer.pdb_tidying(post[:len(post)-4]+'_disordered_sorted.pdb')

    os.rename(post[:len(post)-4]+'_disordered_sorted_tidied.pdb', post[:len(post)-4]+'_withschains.pdb')

    # Deletion of side chains in flexible zones
    target = post[:len(post)-4]+'_withschains.pdb'
    target_pdb = PandasPdb().read_pdb(target)
    target_pdb_df = PDBUtils.get_pdb_atoms_df(target)
  
    half_1, half_2 = RFDSchains.get_rigid_residues(contigs)



    # BUG de aqui abajo hace lo contrario de lo que yo quiero, recuperar solo las side chains flexibles del prepdb
    # y deberia hacer lo contrario. La cosa es que esta parte se pretende recuperar side chains rigidas aunque ahora
    # este el bug, y si esto va en la pipeline, el repacker que se tenga que utilizar deberia poder hacerlo solo
    # de las zonas que le digas
     
    df_noschains1 = RFDSchains.delete_sidechains(target_pdb_df, half_1)
    df_noschains1_noschains2 = RFDSchains.delete_sidechains(df_noschains1, half_2)
    
    
    # Side chain - removed df to PDB
    output_pdb = PandasPdb().read_pdb(post)
    output_pdb.df['ATOM'] = df_noschains1_noschains2
    output_pdb.to_pdb(path= post[:len(post)-4]+'_withrigidschains.pdb', records = ['ATOM'], gz=False)  

    # Sorting PDB for atom number
    PDBProcessor.pdb_atomrenumber(post[:len(post)-4]+'_withrigidschains.pdb', 1)
    
    # Removing unnecessary mid files
    
    os.remove(post[:len(post)-4]+'_disordered_sorted_tidied.pdb')
    os.remove(post[:len(post)-4]+'_disordered.pdb')
    os.remove(post[:len(post)-4]+'_withrigidschains.pdb')
    os.remove(post[:len(post)-4]+'_withschains.pdb') # it could be kept in case we want to readd all side chains
                                                    # an atom renumber would be necessary
    os.rename(post[:len(post)-4]+'_withrigidschains_atomsorted.pdb', post[:len(post)-4]+'_withrigidschains.pdb')
    """