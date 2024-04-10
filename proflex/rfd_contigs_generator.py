"""
Generates the contigs that indicate the regions to diffuse or fix in RFDiffusion
Accepts nanobodies (1 chain) or antibodies (2 chains) and antigens (1 chain for now) as inputs
Example:
python /path/to/rfd_contigs_generator.py --pdb XXXX.pdb --ab1 H --ab2 L --ag A --distance_threshold 7
python /path/to/rfd_contigs_generator.py --pdb XXXX.pdb --ab1 B --ag A --distance_threshold 7

"""

from proflex.utils import PDBUtils
from proflex.utils import InterfaceAnalyzer
from proflex.rfdiffusion import RFDContigs
import argparse
import pandas as pd

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Returns RFdiffusion contigs from a PDB file by taking chains ID and distance_threshold as inputs')
    parser.add_argument('--pdb', type=str, help='Path to the PDB file', required=True)
    parser.add_argument('--ab1', type=str, help='Antibody first chain or nanobody unique chain', required=True)
    parser.add_argument('--ab2', type=str, help='Antibody second chain', required=False)
    parser.add_argument('--ag', type=str, help='Antigen unique chain. If in different chains, merge in a single chain', required=False)
    parser.add_argument('--distance_threshold', type=float, help='Threshold to consider interaction: 7A default', default=7)
    args = parser.parse_args()
    pdb_file = args.pdb
    ab1 = args.ab1
    ab2 = args.ab2
    ag = args.ag
    d_thold = args.distance_threshold
        
    # Nanobody - antigen interaction
    if ab2 is None:

        # Loading data
        atom_df = PDBUtils.get_pdb_atoms_df(pdb_file)
        atom_df_ca = PDBUtils.get_ca(atom_df)
        chains_id_ordered = PDBUtils.get_chains_id(atom_df_ca)        

        # Setting the contigs order so that it coincides with the order of the chains in the PDB file
        # (compulsory for a success run in RFdiffusion)
        input_ch = [ab1, ag]
        if input_ch == chains_id_ordered:
            id1 = ab1
            id2 = ag
        else:
            id1 = ag
            id2 = ab1
        
        # Analyzing interface and getting contigs
        chain_1 = PDBUtils.get_chain(atom_df_ca, id1)
        chain_2 = PDBUtils.get_chain(atom_df_ca, id2)
        chain1_relevant = PDBUtils.get_relevant_columns(chain_1)
        chain2_relevant = PDBUtils.get_relevant_columns(chain_2)
        int1, int2, detected_interactions = InterfaceAnalyzer.get_interface_residues_by_chain(
            chain1_relevant, chain2_relevant, d_thold) # detected_interactions is unused
        intchain1additional = InterfaceAnalyzer.amplify_selection_residues(int1, chain1_relevant)
        intchain2additional = InterfaceAnalyzer.amplify_selection_residues(int2, chain2_relevant)
        contigs_1 = RFDContigs.get_contigs(chain_1, intchain1additional)
        contigs_2 = RFDContigs.get_contigs(chain_2, intchain2additional)
        contigs = f"[{contigs_1}0 {contigs_2[:-1]}]"
        print(contigs)
    
    # Antibody (2 chains) - antigen interaction
    if ab2 is not None:
        
        # Defining chains, loading data and analyzing interface of antibodies with respect to the antigen
        id1 = ab1
        id2 = ab2
        id3 = ag
        atom_df = PDBUtils.get_pdb_atoms_df(pdb_file)
        atom_df_ca = PDBUtils.get_ca(atom_df)
        chains_id_ordered = PDBUtils.get_chains_id(atom_df_ca)
        chain_1 = PDBUtils.get_chain(atom_df_ca, id1)
        chain_2 = PDBUtils.get_chain(atom_df_ca, id2)
        chain_3 = PDBUtils.get_chain(atom_df_ca, id3)
        chain1_relevant = PDBUtils.get_relevant_columns(chain_1)
        chain2_relevant = PDBUtils.get_relevant_columns(chain_2)
        chain3_relevant = PDBUtils.get_relevant_columns(chain_3)
        int1, int3_1, detected_interactions = InterfaceAnalyzer.get_interface_residues_by_chain(
            chain1_relevant, chain3_relevant, d_thold) # detected_interactions is unused
        intchain1additional = InterfaceAnalyzer.amplify_selection_residues(int1, chain1_relevant)
        int2, int3_2, detected_interactions = InterfaceAnalyzer.get_interface_residues_by_chain(
            chain2_relevant, chain3_relevant, d_thold) # detected_interactions is unused
        intchain2additional = InterfaceAnalyzer.amplify_selection_residues(int2, chain1_relevant)
        
        # Getting contigs for the antibody
        contigs_1 = RFDContigs.get_contigs(chain_1, intchain1additional)
        contigs_2 = RFDContigs.get_contigs(chain_2, intchain2additional)
        
        # Analyzing the interface of the antigen with respect to the antibody and getting contigs for the antigen
        merged_ag_int = pd.concat([int3_1, int3_2]) # Coges los seleccionados de cada interfaz y mezclas
        merged_ag_int_sorted = merged_ag_int.sort_values(by='residue_number') # Los ordenas por res num
        merged_ag_int_uniq = merged_ag_int_sorted.drop_duplicates(subset=['residue_number']) # Elimina duplicados en res num
        int3_mixed = merged_ag_int_uniq # Esto es todo lo seleccionado en la interfaz
        intchain3additional = InterfaceAnalyzer.amplify_selection_residues(int3_mixed, chain3_relevant) # Amplifica por vecinos
        contigs_3 = RFDContigs.get_contigs(chain_3, intchain3additional)
        
        # Setting the contigs order so that it coincides with the order of the chains in the PDB file
        # (compulsory for a success run in RFdiffusion) and building the contigs
        if chains_id_ordered[0] in contigs_1:
            first = contigs_1
        elif chains_id_ordered[0] in contigs_2:
            first = contigs_2
        elif chains_id_ordered[0] in contigs_3:
            first = contigs_3           
        if chains_id_ordered[1] in contigs_1:
            second = contigs_1
        elif chains_id_ordered[1] in contigs_2:
            second = contigs_2
        elif chains_id_ordered[1] in contigs_3:
            second = contigs_3    
        if chains_id_ordered[2] in contigs_1:
            third = contigs_1
        elif chains_id_ordered[2] in contigs_2:
            third = contigs_2
        elif chains_id_ordered[2] in contigs_3:
            third = contigs_3    
        contigs = f"[{first}0 {second}0 {third[:-1]}]"
        print(contigs)