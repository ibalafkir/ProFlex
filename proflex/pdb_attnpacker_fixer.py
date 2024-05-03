"""
Version of pdb_rfd_fixer which solves the PDB format of an output PDB of RFd-AttnPacker (the first is
inserted in the second without processing) ot only AttnPacker or in general a PDB in which all residues
are mixed in the same chain

With respect to version 1, this program tries to keep the repacking prediction of only flexible regions,
that's why the contigs codes are needed, and recovering the rigid regions coordinates before diffunding
the system. As RFd does a brief translation, a superposition to the non-diffunded system (input in RFd,
aka prePDB in other logging statements)

PDB INPUT MUST BE IN A SINGLE CHAIN ID, OTHERWISE ATTNPACKER CAN'T REPACK ALL CHAINS AND ONLY REPACKS THE FIRST AND THUS THIS DOES NOT WORK!

$ pdb_gap -A <pdb>
# TODO probably delete due to no use of attnpacker
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
    parser.add_argument('postrfd', type=str, help='Path to the output PDB RFDiffusion+ATTNPacker file')
    #parser.add_argument('-contigs', type=str, help='Contigs string')
    args = parser.parse_args()

    pre = args.prerfd
    post = args.postrfd
    #contigs = args.contigs

                                ####### REFORMATTING RFD-ATTNPACKER OUTPUT ########

    # Strategy: reduceprePDB (the input in RFD-AttnPacker) to the same atom lines
    # that the postPDB (output of RFD-AttnPacker) has, i.e. backbone + side chains
    # without hydrogens and OXT atom names, and then make the chain, residue number
    # and atom number correction

    print(f"Correcting {post}")

    print("Editing reference")
    RFDFixer.pdb_delel(pre, ['H'], pre[:-4]+'_noh.pdb')
    RFDFixer.pdb_atom(pre[:len(pre)-4]+'_noh.pdb')
    # AttnPacker does not add 'OXT' atom types, they need to be removed in our pdb format correction strategy
    RFDFixer.pdb_delotherel(pre[:-4]+'_noh_atom.pdb', ['OXT'], pre[:-4]+'_noh_atom_noOXT.pdb')

    RFDFixer.pdb_atom(post)
    RFDFixer.correct_rfd_pdbs(post[:-4]+'_atom.pdb', pre[:-4]+'_noh_atom_noOXT.pdb')
    print("Tidying...")
    RFDFixer.pdb_tidying(post[:-4]+'_atom_chainsfixed.pdb')
    PDBProcessor.pdb_atomrenumber(post[:-4]+'_atom_chainsfixed_tidied.pdb', 1) # TODO TER line is counted as an atom

    # Removing temp files and rename the last to "name_of_pdb"_fixed.pdb
    os.remove(pre[:-4]+'_noh.pdb')
    os.remove(pre[:-4]+'_noh_atom.pdb')
    # os.remove(pre[:-4]+'_noh_atom_noOXT.pdb')
    os.remove(post[:-4]+'_atom.pdb')
    os.remove(post[:-4]+'_atom_chainsfixed.pdb')
    os.remove(post[:-4]+'_atom_chainsfixed_tidied.pdb')
    os.rename(post[:-4]+'_atom_chainsfixed_tidied_atomsorted.pdb', post[:-4]+'_attnfixed.pdb')



    os.remove(pre[:-4]+'_noh_atom_noOXT.pdb') # TODO IF THE SUPERPOSING PART IS DONE, DELETE THIS LINE


        ## HERE ONWARDS SUPERPOSITION AND OTHER PROCEDURES ARE DONE TO KEEP ONLY PREDICTIONS OF FLEXIBLE RESIDUES
        ## IT WORKS PRETTY WELL :)))







    """
                                    ####### SUPERPOSING AND RECOVERING RIGID SIDE CHAINS COORDINATES ########

    # Strategy: superpose the fixed prePDB with the postPDB fixed (#TODO with CA/backbone of rigid residues only if possible)
    # and recover the side chains coordinates of rigid residues from the prePDB

    # Processing the contigs
    half_1, half_2 = RFDSchains.get_rigid_residues(contigs)

    # Superposing #TODO superpose only with resnum+chain_id listed in half1 and half2

    # Como esta funcion requiere mismo n de atomos, superponemos no con prePDB original sino con este tras quitarle
    # hidrogenos y oxt, que son atom types que attnpacker no añade
    RFDFixer.superpose_v2(pre[:-4]+'_noh_atom_noOXT.pdb', post[:-4]+'_rfdattnfixed.pdb', post[:-4]+'_rfdattnfixed_superposed.pdb')
    os.remove(pre[:-4]+'_noh_atom_noOXT.pdb') # remove it as it is not needed anymore

    # Recovering side chains
    # Read both pdbs in df, remove rigid side chains in post_fixed, get rigid side chains of prePDB, merge prePDB
    # rigid side chains with post_fixed with
    # flexible side chains, sort, tidy, renumber atoms

    post_fixed_superposed_atomdf = PDBUtils.get_pdb_atoms_df(post[:-4]+'_rfdattnfixed_superposed.pdb')
    post_fixed_superposed_atomdf_norigidschains1 = RFDSchains.delete_sidechains(post_fixed_superposed_atomdf, half_1)
    post_fixed_superposed_atomdf_norigidschains1_norigidschains2 = RFDSchains.delete_sidechains(post_fixed_superposed_atomdf_norigidschains1, half_2)
    post_fixed_superposed_atomdf_norigidschains = post_fixed_superposed_atomdf_norigidschains1_norigidschains2
    post_fixed_superposed_atomdf_norigidschains1_norigidschains2_topdb = PandasPdb().read_pdb(post[:-4]+'_rfdattnfixed_superposed.pdb')
    post_fixed_superposed_atomdf_norigidschains1_norigidschains2_topdb.df['ATOM'] = post_fixed_superposed_atomdf_norigidschains1_norigidschains2
    post_fixed_superposed_atomdf_norigidschains1_norigidschains2_topdb.to_pdb(path= post[:-4]+'_rfdattnfixed_superposed_norigidschains.pdb', records = ['ATOM'], gz=False)

    pre_atomdf = PDBUtils.get_pdb_atoms_df(pre)
    pre_atomdf_rigidschains1 = RFDSchains.get_sidechains(pre_atomdf, half_1)
    pre_atomdf_rigidschains2 = RFDSchains.get_sidechains(pre_atomdf, half_2)
    pre_atomdf_rigidschains1_rigidschains2 = pd.concat([pre_atomdf_rigidschains1, pre_atomdf_rigidschains2])
    pre_atomdf_rigidschains1_rigidschains2_topdb = PandasPdb().read_pdb('4pou.pdb')
    pre_atomdf_rigidschains1_rigidschains2_topdb.df['ATOM'] = pre_atomdf_rigidschains1_rigidschains2
    pre_atomdf_rigidschains1_rigidschains2_topdb.to_pdb(pre[:-4]+'_rigidschains.pdb', records = ['ATOM'], gz=False)


    mixed_df = pd.concat([post_fixed_superposed_atomdf_norigidschains, pre_atomdf_rigidschains1_rigidschains2])
    output_pdb = PandasPdb().read_pdb(post)
    output_pdb.df['ATOM'] = mixed_df
    output_pdb.to_pdb(path = post[:len(post)-4]+'_disordered.pdb', records = ['ATOM'], gz=False)

    RFDFixer.pdb_sorting(post[:len(post)-4]+'_disordered.pdb')

    RFDFixer.pdb_tidying(post[:len(post)-4]+'_disordered_sorted.pdb')

    PDBProcessor.pdb_atomrenumber(post[:len(post)-4]+'_disordered_sorted_tidied.pdb', 1)

    os.remove(post[:-4]+'_rfdattnfixed_superposed.pdb')
    os.remove(post[:-4]+'_disordered_sorted_tidied.pdb')
    os.remove(post[:-4]+'_disordered.pdb')
    os.remove(post[:-4]+'_disordered_sorted.pdb')
    os.remove(pre[:-4]+'_rigidschains.pdb')
    os.remove(post[:-4]+'_rfdattnfixed_superposed_norigidschains.pdb')
    os.rename(post[:-4]+'_disordered_sorted_tidied_atomsorted.pdb', post[:-4]+ '_rfdattnprevrigidschainsfixed.pdb')

    # Being _rfdattnfixed.pdb the output of RFD+AttnPacker the format fixed version
    # and _rfdattnprevrigidschainsfixed.pdb the final version of post RFD+AttnPacker PDB well formatted
    # with prePDB's rigid side chains and predicted flexible side chains


    """
