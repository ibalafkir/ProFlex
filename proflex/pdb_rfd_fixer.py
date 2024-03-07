from proflex.utils import RFDFixer, PDBUtils, PDBProcessor
import argparse
import pandas as pd
from biopandas.pdb import PandasPdb 
from rfdiffusion import RFDSchains 
import os
"""
rfd = 'Ab_diff_N05_0_rfdfixed.pdb'
pre = '4pou.pdb'
post = 'hola.pdb'
contig_ex = '[A1-110/5-5/A116-124/0 B1-98/3-3/B102-105/4-4/B110-121]'

RFDFixer.superpose(pre, rfd, 'hola.pdb')

# Reads df of the preRFD PDB and postRFD fixed for format PDB
pre_atom = PDBUtils.get_pdb_atoms_df(pre)
post_atom = PDBUtils.get_pdb_atoms_df(post)

# Keeps side chains of the preRFD PDB df
pre_atom_sidechains = pre_atom[~pre_atom['atom_name'].isin(['N', 'CA', 'C', 'O'])] 

# Mixes RFD backbone df with preRFD PDB df side chains
mixed_df = pd.concat([post_atom, pre_atom_sidechains]) 

# Generates a temporal PDB file out of the previous df
output_pdb = PandasPdb().read_pdb('hola.pdb')
output_pdb.df['ATOM'] = mixed_df
output_pdb.to_pdb(path= 'hola_disordered.pdb', records = ['ATOM'], gz=False)  

# Sorts the PDB according to chain id and residue number
RFDFixer.pdb_sorting('hola_disordered.pdb')

# Tidies the previous PDB (add TERs, END...). This PDB has all side chains.
RFDFixer.pdb_tidying('hola_disordered_sorted.pdb')

# Deletion of side chains in flexible zones
prueba = 'hola_disordered_sorted_tidied.pdb'
prueba_pdb = PandasPdb().read_pdb(prueba)
prueba_df = PDBUtils.get_pdb_atoms_df(prueba)
half_1, half_2 = RFDSchains.get_rigid_residues(contig_ex)

pruebaa_df_noschains1 = RFDSchains.delete_sidechains(prueba_df, half_1)
pruebaa_df_noschains1_noschains2 = RFDSchains.delete_sidechains(pruebaa_df_noschains1, half_2)

# Side chain - removed df to PDB
output_pdb = PandasPdb().read_pdb('hola.pdb')
output_pdb.df['ATOM'] = pruebaa_df_noschains1_noschains2
output_pdb.to_pdb(path= 'hola_disordered_sorted_tidied_schree.pdb', records = ['ATOM'], gz=False)  

# Sorting PDB for atom number
PDBProcessor.pdb_atomrenumber('hola_disordered_sorted_tidied_schree.pdb', 1)
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Fixes PDB features missed in the output of RFDiffusion'
        )
    parser.add_argument('-prerfd', type=str, help='Path to the PDB file (input of RFD)')
    parser.add_argument('-postrfd', type=str, help='Path to the output PDB RFDiffusion file')
    parser.add_argument('-contigs', type=str, help='Contigs string')
    args = parser.parse_args()
    
    pre = args.prerfd
    post = args.postrfd
    contigs = args.contigs

    # Fixing the .pdb format in the RFD PDB output thanks to the PDB input in RFD
    n_atoms = RFDFixer.get_n_atoms(post) # RFD output PDBs have an ordered atom number so
                                            # the function can be well applied
    RFDFixer.pdb_backbone(pre) 
    RFDFixer.pdb_atom(pre[:len(pre)-4]+'_backbone.pdb')
    RFDFixer.correct_rfd_pdbs(post, pre[:len(pre)-4]+'_backbone.pdb')
    RFDFixer.pdb_tidying(post[:len(post)-4]+'_chainsfixed.pdb')
    RFDFixer.del_mid_files(pre, post)
    
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

    # Tidies the previous PDB (add TERs, END...). This PDB has all side chains.
    RFDFixer.pdb_tidying(post[:len(post)-4]+'_disordered_sorted.pdb')

    os.rename(post[:len(post)-4]+'_disordered_sorted.pdb', post[:len(post)-4]+'_withschains.pdb')

    # Deletion of side chains in flexible zones
    target = post[:len(post)-4]+'_withschains.pdb'
    target_pdb = PandasPdb().read_pdb(target)
    prueba_df = PDBUtils.get_pdb_atoms_df(target)
    half_1, half_2 = RFDSchains.get_rigid_residues(contigs)

    df_noschains1 = RFDSchains.delete_sidechains(prueba_df, half_1)
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