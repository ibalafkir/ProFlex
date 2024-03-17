"""
Generates the contig codes that indicate the regions to diffuse in RFDiffusion
Input from commandline: PDB with both proteins interacting and the chains to analyze interface

Example:

python rfd_contigs_generator.py --pdb B_PD1_diff_8_prep_relax_fixed.pdb --id1 A --id2 B --distance_threshold 200

# TODO This program has to generate a contig of all chains, not just 2 of them
# TODO id1->ab1 should be an ab/nb chain, id2->ab2 should be another ab chain if existing, and ag should be the antigen
# TODO can the antigen have more than 1 chain?

"""

from proflex.utils import PDBUtils
from proflex.utils import InterfaceAnalyzer
from proflex.rfdiffusion import RFDContigs
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Returns RFdiffusion contig code from a PDB file'
        )
    parser.add_argument('--pdb', type=str, help='Path to the PDB file')
    parser.add_argument('--id1', type=str, help='First chain')
    parser.add_argument('--id2', type=str, help='Second chain')
    parser.add_argument('--distance_threshold', type=float,
                        help='Distance threshold for interface '
                             'residues (6 Angstrom is default)', default=6)
    args = parser.parse_args()
    pdb_file = args.pdb
    id1 = args.id1
    id2 = args.id2
    distance_threshold = args.distance_threshold

    atom_df = PDBUtils.get_pdb_atoms_df(pdb_file)
    atom_df_ca = PDBUtils.get_ca(atom_df)
    chain_1 = PDBUtils.get_chain(atom_df_ca, id1)
    chain_2 = PDBUtils.get_chain(atom_df_ca, id2)
    chain1_relevant = PDBUtils.get_relevant_columns(chain_1)
    chain2_relevant = PDBUtils.get_relevant_columns(chain_2)
    int1, int2, detected_interactions = InterfaceAnalyzer.get_interface_residues_by_chain(
        chain1_relevant, chain2_relevant, distance_threshold) # detected_interactions is unused
    intchain1additional = InterfaceAnalyzer.amplify_selection_residues(int1, chain1_relevant)
    intchain2additional = InterfaceAnalyzer.amplify_selection_residues(int2, chain2_relevant)
    print(RFDContigs.get_contigs(chain_1, chain_2, intchain1additional, intchain2additional))