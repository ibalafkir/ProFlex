"""
ProFlex interface analyzer of the interface of the interaction between 2 proteins (chain 1 and chain 2)
The execution will is done by inputting the PDB file, both chain IDs to study interactions and a maximum distance
threshold up to which 2 residues are considered to be interacting. The user needs to select the chains to study.
"""
import argparse
import pandas as pd
from proflex.pdb import PdbDf
from proflex.interface import InterfaceAnalyzer

def run(pdb_file, id1, id2, distance_threshold):

    atom_df = PdbDf.atoms(pdb_file)
    atom_df_ca = PdbDf.ca(atom_df)
    chain_1 = PdbDf.chain(atom_df_ca, id1)
    chain_2 = PdbDf.chain(atom_df_ca, id2)
    chain1_relevant = PdbDf.rel_col(chain_1)
    chain2_relevant = PdbDf.rel_col(chain_2)
    int1, int2, detected_interactions = InterfaceAnalyzer.calculate(chain1_relevant, chain2_relevant, distance_threshold)
    intchain1additional = InterfaceAnalyzer.extend_neighbourhood(int1, chain1_relevant)
    intchain2additional = InterfaceAnalyzer.extend_neighbourhood(int2, chain2_relevant)
    return int1, int2, detected_interactions, intchain1additional, intchain2additional

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculates interface residues between two chains of PDB file')
    parser.add_argument('--pdb_file', type=str, help='Path to the PDB file')
    parser.add_argument('--id1', type=str, help='First chain')
    parser.add_argument('--id2', type=str, help='Second chain')
    parser.add_argument('--distance_threshold', type=float,
                        help='Distance threshold for interface residues (6 Angstrom is default)', default=6)
    args = parser.parse_args()
    pdb_file = args.pdb_file
    id1 = args.id1
    id2 = args.id2
    distance_threshold = args.distance_threshold
    int1, int2, detected_interactions, intchain1additional, intchain2additional = run(pdb_file, args.id1, args.id2, distance_threshold)
    
    if not int1.empty and not int2.empty:
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        print("Interface residues from chain", id1)
        print(int1, "\n")
        print("Interface residues from chain", id2)
        print(int2, "\n")
        print("Interacting pairs")
        print(detected_interactions, "\n")
        print("Additional interacting residues (neighborhood = 2 residues) in chain", id1)
        print(intchain1additional, "\n")
        print("Additional interacting residues (neighborhood = 2 residues) in chain", id2)
        print(intchain2additional, "\n")
        print("Code for PELE com_distance metric \n")
        InterfaceAnalyzer.pele_com_distance(int1, int2, id1, id2)
    else:
        print("Interactions could not be detected at the distance threshold of", distance_threshold, "Angstroms")