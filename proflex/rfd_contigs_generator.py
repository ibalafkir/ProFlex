"""
Generates the contig codes that indicate the regions to diffuse or fix in RFDiffusion
Input from commandline: PDB with both proteins interacting and the chains to analyze interface

Example:

# TODO Check whether completely fixed chains are accepted by RFdiffusion

"""

from proflex.utils import PDBUtils
from proflex.utils import InterfaceAnalyzer
from proflex.rfdiffusion import RFDContigs
import argparse
import pandas as pd

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Returns RFdiffusion contig code from a PDB file')
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
        
        print(f"Nanobody - antigen system detected")
        #print(f"Antibody chain: {ab1}")
        #print(f"Antigen chain: {ag}")

        
        atom_df = PDBUtils.get_pdb_atoms_df(pdb_file)
        atom_df_ca = PDBUtils.get_ca(atom_df)
        chains_id_ordered = PDBUtils.get_chains_id(atom_df_ca)
        
        print(f"Order of chain appearence in {pdb_file}: first is {chains_id_ordered[0]}, second is {chains_id_ordered[1]}")
        

        input_ch = [ab1, ag]
        
        if input_ch == chains_id_ordered:
            id1 = ab1
            id2 = ag
        else:
            id1 = ag
            id2 = ab1
        
        print(f"In the contig output, the first id will be {id1} and the second {id2}, this MUST coincide with PDB order")

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
        
        print(f"Antibody - antigen system detected")
        #print(f"Antibody first input chain: {ab1}")
        #print(f"Antibody second input chain: {ab2}")
        #print(f"Antigen chain: {ag}")

        
        id1 = ab1
        id2 = ab2
        id3 = ag
        
        atom_df = PDBUtils.get_pdb_atoms_df(pdb_file)
        atom_df_ca = PDBUtils.get_ca(atom_df)
        chains_id_ordered = PDBUtils.get_chains_id(atom_df_ca)
        

        print(f"Order of chain appearence in {pdb_file}:"
            f" first is {chains_id_ordered[0]}"
            f", second is {chains_id_ordered[1]}"
            f", third is {chains_id_ordered[2]}")
        
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
        
        contigs_1 = RFDContigs.get_contigs(chain_1, intchain1additional)
        contigs_2 = RFDContigs.get_contigs(chain_2, intchain2additional)
        
        merged_ag_int = pd.concat([int3_1, int3_2]) # Coges los seleccionados de cada interfaz y mezclas
        merged_ag_int_sorted = merged_ag_int.sort_values(by='residue_number') # Los ordenas por res num
        merged_ag_int_uniq = merged_ag_int_sorted.drop_duplicates(subset=['residue_number']) # Elimina duplicados en res num
        int3_mixed = merged_ag_int_uniq # Esto es todo lo seleccionado en la interfaz
        intchain3additional = InterfaceAnalyzer.amplify_selection_residues(int3_mixed, chain3_relevant) # Amplifica por vecinos
        
        contigs_3 = RFDContigs.get_contigs(chain_3, intchain3additional)
        
        
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
        print(f"In principle, the contigs of {ag} should be correct, as a test run the same program twice for each ab chain")
        
        print("Int 3 con una interfaz")
        print(int3_1['residue_number'].tolist())
        print("Int 3 con una interfaz")
        print(int3_2['residue_number'].tolist())
        print("Int 3 mezclado interfaz y ordenado")
        print(int3_mixed['residue_number'].tolist())
        print("Int 3 con vecinos")
        print(intchain3additional['residue_number'].tolist())