"""
ProFlex interface analyser of the interface of the interaction between 2 proteins (chain 1 and chain 2)
"""
import argparse
from biopandas.pdb import PandasPdb
import numpy as np
import pandas as pd

def get_pdb_atoms_df(pdb_file):
    ppdb = PandasPdb().read_pdb(pdb_file)
    print("Reading PDB file")
    print(" ")
    return ppdb.df['ATOM'] # only keeps information related to atoms

def get_ca(atom_df):
    print("Filtering CA")
    print(" ")
    return atom_df[atom_df['atom_name'] == 'CA']

def get_chains_id(atom_df):
    result=[]
    array = pd.unique(atom_df["chain_id"])
    for i in array:
        result.append(i)
    print("Obtaining chains' identification")
    print(" ")
    return result

def get_chain(atom_df, chain_name):
    print("Obtaining different protein chains:", chain_name)
    print(" ")
    return atom_df[atom_df['chain_id'] == chain_name]

def get_relevant_columns(chain):
    desired_cols = ['chain_id', 'residue_number', 'residue_name', 'x_coord', 'y_coord', 'z_coord']
    x = chain[desired_cols]
    print("Filtering relevant columns of chain", chain['chain_id'].values[0])
    print(" ")
    return x

def get_interface_residues_by_chain(chain1: str, chain2: str, distance_threshold=6) -> (np.array, np.array, np.array):
    """
    Extraction of interface residues by a vectorised calculus of the distance between CA atoms
    :param chain1: str
    :param chain2: str
    :param distance_threshold: float
    :return: 2 np.arrays containing interface residues from each chain
             1 np.array containing interface residues and distance
    """
    coords_1 = chain1[['x_coord', 'y_coord', 'z_coord']].values
    coords_2 = chain2[['x_coord', 'y_coord', 'z_coord']].values
    distances = np.linalg.norm(coords_1[:, None] - coords_2, axis=2)
    close_residues = np.where(distances <= distance_threshold)
    residues_interface_1 = chain1.iloc[close_residues[0]]
    residues_interface_2 = chain2.iloc[close_residues[1]]
    detected_distances = pd.DataFrame({
        'First chain residue': residues_interface_1['residue_name'].values,
        'Residue number 1': residues_interface_1['residue_number'].values,
        'X_1': residues_interface_1['x_coord'].values,     # REMOVE?
        'Y_1': residues_interface_1['y_coord'].values,     # REMOVE?
        'Z_1': residues_interface_1['z_coord'].values,     # REMOVE?
        'Second chain residue': residues_interface_2['residue_name'].values,
        'Residue number 2': residues_interface_2['residue_number'].values,
        'X_2': residues_interface_2['x_coord'].values,     # REMOVE?
        'Y_2': residues_interface_2['y_coord'].values,     # REMOVE?
        'Z_2': residues_interface_2['z_coord'].values,     # REMOVE?
        'Distance': distances[close_residues]
    })
    print("Calculating distances")
    print(" ")
    return residues_interface_1, residues_interface_2, detected_distances

def amplify_selection_residues(int, chainrelevant):
    i = 0
    total_dif = []
    res_num = sorted(int['residue_number'].values)
    while i < len(res_num)-1:
        subs = res_num[i+1]-res_num[i]
        total_dif.append(subs)
        if subs == 2: # 2 means nb = 1 (|.|), 3 means nb = 2 (|..|), 4 means nb = 3 (|...|)
            res_num_to_select = res_num[i]+1
            residue_to_add = chainrelevant.loc[chainrelevant['residue_number'] == res_num_to_select]
            int = int._append(residue_to_add)

        if subs == 3:
            res_num_to_select = [res_num[i+1]-1, res_num[i]+1]
            residue_to_add1 = chainrelevant.loc[chainrelevant['residue_number'] == res_num_to_select[0]]
            residue_to_add2 = chainrelevant.loc[chainrelevant['residue_number'] == res_num_to_select[1]]
            int = int._append(residue_to_add1)
            int = int._append(residue_to_add2)

        if subs == 4:
            res_num_to_select = [res_num[i + 1] - 1, res_num[i] + 1, res_num[i] + 2]
            residue_to_add1 = chainrelevant.loc[chainrelevant['residue_number'] == res_num_to_select[0]]
            residue_to_add2 = chainrelevant.loc[chainrelevant['residue_number'] == res_num_to_select[1]]
            residue_to_add3 = chainrelevant.loc[chainrelevant['residue_number'] == res_num_to_select[2]]
            int = int._append(residue_to_add1)
            int = int._append(residue_to_add2)
            int = int._append(residue_to_add3)
        i += 1
    int_sorted = int.sort_values(by='residue_number')
    int_sorted_unique = int_sorted.drop_duplicates()
    print("Detecting additional interactions (neighbourhood = 2) residues for diffusion processes of chain", chainrelevant['chain_id'].values[0])
    print(" ")
    return int_sorted_unique

# Testing amplify_selection_residues
# pruebachainrelevant = chain2_relevant.head(5)
# pruebaint = pruebachainrelevant.drop(pruebachainrelevant.index[1])
# test_1_3 = amplify_selection_residues(pruebaint, pruebachainrelevant)
# pruebaint = pruebaint.drop(pruebaint.index[1])
# test_1_4 = amplify_selection_residues(pruebaint, pruebachainrelevant)
# pruebaint = pruebaint.drop(pruebaint.index[1])
# test_1_5 = amplify_selection_residues(pruebaint, pruebachainrelevant)

def get_interface_residues_by_each_chain(pdb_file, distance_threshold=6) -> (np.array, np.array, np.array, np.array):
    atom_df = get_pdb_atoms_df(pdb_file)
    atom_df_ca = get_ca(atom_df)
    chain_ids = get_chains_id(atom_df_ca)
    chain_1 = get_chain(atom_df_ca, chain_ids[0])
    chain_2 = get_chain(atom_df_ca, chain_ids[1])
    chain1_relevant = get_relevant_columns(chain_1)
    chain2_relevant = get_relevant_columns(chain_2)
    int1, int2, interactions = get_interface_residues_by_chain(chain1_relevant, chain2_relevant, distance_threshold)
    intchain1additional = amplify_selection_residues(int1, chain1_relevant)
    intchain2additional = amplify_selection_residues(int2, chain2_relevant)
    return int1, int2, interactions, intchain1additional, intchain2additional
###

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculates interface residues between protein chains.')
    parser.add_argument('pdb_file', type=str, help='Path to the PDB file')
    parser.add_argument('distance_threshold', type=float, help='Distance threshold for interface residues (6 A is default')
    args = parser.parse_args()
    pdb_file = args.pdb_file
    distance_threshold = args.distance_threshold
    int1, int2, detected_distances, intchain1additional, intchain2additional = get_interface_residues_by_each_chain(pdb_file, distance_threshold)

    if not int1.empty and not int2.empty and not detected_distances.empty:
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        print("INTERFACE RESIDUES FROM CHAIN 1")
        print(int1)
        print(" ")
        print("INTERFACE RESIDUES FROM CHAIN 2")
        print(int2)
        print(" ")
        print("DISTANCE VALUES IN INTERACTIONS")
        print(detected_distances)
        print(" ")
        print("INTERACTIONS WITH N = 2")
        print(intchain1additional)
        print(intchain2additional)
    else:
        print("Interactions could not be detected at the distance threshold of", distance_threshold, "Angstroms")