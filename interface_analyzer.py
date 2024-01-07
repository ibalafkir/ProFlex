"""
ProFlex interface analyzer of the interface of the interaction between 2 proteins (chain 1 and chain 2)
Code organization: libraries, functions* and executions**.
* Not all functions are used
** The execution will be done by inputting the PDB file, both chain IDs to study interactions and
   maximum distance threshold up to which 2 residues are considered to be interacting. The user
   needs to select which chains ID will study.
"""
import argparse
from biopandas.pdb import PandasPdb
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO
import warnings
from Bio import BiopythonWarning

def get_pdb_atoms_df(pdb_file):
    """
    Extracts ATOM information from a PDB file and returns a pandas DataFrame
    :param Path to PDB file: str
    :return Pandas DataFrame: pandas.DataFrame
    """
    ppdb = PandasPdb().read_pdb(pdb_file)
    print("Extracting atom information from PDB... \n")
    return ppdb.df['ATOM'] # only keeps information related to atoms

def get_ca(atom_df):
    """
    Extracts CA rows from a pandas DataFrame with ATOM information
    :param Pandas DataFrame with ATOM information: pandas.DataFrame
    :return Pandas DataFrame with CA information: pandas.DataFrame
    """
    print("Getting alpha carbons... \n")
    return atom_df[atom_df['atom_name'] == 'CA']

def get_chains_id(atom_df):
    """
    Extracts ID chains from a pandas DataFrame with ATOM information
    :param Pandas DataFrame with ATOM information: pandas.DataFrame
    :return List with chains ID: list
    """
    result=[]
    array = pd.unique(atom_df["chain_id"])
    for i in array:
        result.append(i)
    print("Obtaining chains' identification... \n")
    return result

def get_chain(df, chain_name):
    """
    Returns pandas DataFrame with the ATOM information of an input chain
    :param Pandas DataFrame with ATOM information (CA atom_names, all atom_names...):
           pandas.DataFrame
    :return Pandas DataFrame with ATOM information of an indicated chain:
            pandas.DataFrame
    """
    print("Obtaining different protein chains...:", chain_name, "\n")
    return df[df['chain_id'] == chain_name]

def extract_pdb_chains(pdbid):
    """
    Extracts PDB files with the different chains from a PDB ID
    :param PDB ID: str
    :return PDB files in working directory: .pdb
    """
    warnings.simplefilter('ignore', BiopythonWarning) # ignore discontinuity lines, revise
    io = PDBIO()
    pdb = PDBParser().get_structure(pdbid, pdbid+".pdb")
    for chain in pdb.get_chains():
        io.set_structure(chain)
        io.save(pdb.get_id() + "_" + chain.get_id() + ".pdb")

def get_relevant_columns(chain):
    """
    Extracts "relevant columns" from an ATOM PDB pandas DatFrame. These are:
    chain_id, residue_number, residue_name, x_coord, y_coord, z_coord
    :param Atom DataFrame: pandas.DataFrame
    :return Atom DataFrame with relevant columns
    """
    desired_cols = ['chain_id', 'residue_number', 'residue_name', 'x_coord', 'y_coord', 'z_coord']
    x = chain[desired_cols]
    print("Filtering relevant columns of chain", chain['chain_id'].values[0], "... \n")
    return x

def get_interface_residues_by_chain(chain1: str, chain2: str, distance_threshold=6):
    """
    Extraction of interface residues, meaning residues which CA are at a maximum given distance
    (default = 6 Angstrom), which is the maximum distance up to which interaction is considered
    to occur
    :param First chain ID: str
    :param Second chain ID: str
    :param Distance threshold: float
    :return 2 np.arrays containing interface residues from each chain
    """
    coords_1 = chain1[['x_coord', 'y_coord', 'z_coord']].values
    coords_2 = chain2[['x_coord', 'y_coord', 'z_coord']].values
    distances = np.linalg.norm(coords_1[:, None] - coords_2, axis=2)
    close_residues = np.where(distances <= distance_threshold)
    residues_interface_1 = chain1.iloc[close_residues[0]]
    residues_interface_2 = chain2.iloc[close_residues[1]]
    #detected_distances = pd.DataFrame({
    #    'First chain residue': residues_interface_1['residue_name'].values,
    #    'Residue number 1': residues_interface_1['residue_number'].values,
    #    'X_1': residues_interface_1['x_coord'].values,
    #    'Y_1': residues_interface_1['y_coord'].values,
    #    'Z_1': residues_interface_1['z_coord'].values,
    #    'Second chain residue': residues_interface_2['residue_name'].values,
    #    'Residue number 2': residues_interface_2['residue_number'].values,
    #    'X_2': residues_interface_2['x_coord'].values,
    #    'Y_2': residues_interface_2['y_coord'].values,
    #    'Z_2': residues_interface_2['z_coord'].values,
    #    'Distance': distances[close_residues]
    #})
    print("Calculating distances... \n")
    return residues_interface_1, residues_interface_2

def get_interaction_distances(chain1: str, chain2: str, distance_threshold=6):
    """
    Extraction of distances between interface residues, meaning residues which CA are at a maximum given
    distance (default = 6 Angstrom), which is the maximum distance up to which interaction is considered
    to occur
    :param First chain ID: str
    :param Second chain ID: str
    :param Distance threshold: float
    :return single np.array containing interface residues from each chain and the distance
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
        'X_1': residues_interface_1['x_coord'].values,
        'Y_1': residues_interface_1['y_coord'].values,
        'Z_1': residues_interface_1['z_coord'].values,
        'Second chain residue': residues_interface_2['residue_name'].values,
        'Residue number 2': residues_interface_2['residue_number'].values,
        'X_2': residues_interface_2['x_coord'].values,
        'Y_2': residues_interface_2['y_coord'].values,
        'Z_2': residues_interface_2['z_coord'].values,
        'Distance': distances[close_residues]
    })
    print("Calculating distances... \n")
    return detected_distances

def amplify_selection_residues(int, chainrelevant):
    """
    Amplifies the selected residues in a DataFrame by neighborhood = 2
    :param pandas.DataFrame of a chain with interacting residues: pandas.DataFrame
    :param pandas.DataFrame of the chain with which the interaction DataFrame was generated: pandas.DataFrame
    :return: pandas.DataFrame of a chain with interacting residues amplified by neighborhood = 2: pandas.DataFrame
    """
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
    print("Detecting additional interactions (neighborhood = 2) residues for diffusion processes of chain", chainrelevant['chain_id'].values[0], "... \n")
    return int_sorted_unique

def analyze_interface(pdb_file, id1, id2, distance_threshold):
    atom_df = get_pdb_atoms_df(pdb_file)
    atom_df_ca = get_ca(atom_df)
    chain_1 = get_chain(atom_df_ca, id1)
    chain_2 = get_chain(atom_df_ca, id2)
    chain1_relevant = get_relevant_columns(chain_1)
    chain2_relevant = get_relevant_columns(chain_2)
    int1, int2 = get_interface_residues_by_chain(chain1_relevant, chain2_relevant, distance_threshold)
    intchain1additional = amplify_selection_residues(int1, chain1_relevant)
    intchain2additional = amplify_selection_residues(int2, chain2_relevant)
    return int1, int2, intchain1additional, intchain2additional


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculates interface residues between two chains of PDB file')
    parser.add_argument('--pdb_file', type=str, help='Path to the PDB file')
    parser.add_argument('--id1', type=str, help='First chain')
    parser.add_argument('--id2', type=str, help='Second chain')
    parser.add_argument('--distance_threshold', type=float,
                        help='Distance threshold for interface '
                             'residues (6 Angstrom is default)', default=6)
    args = parser.parse_args()
    pdb_file = args.pdb_file
    id1 = args.id1
    id2 = args.id2
    distance_threshold = args.distance_threshold
    int1, int2, intchain1additional, intchain2additional = analyze_interface(pdb_file, args.id1, args.id2, distance_threshold)

    if not int1.empty and not int2.empty:
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        print("INTERFACE RESIDUES FROM CHAIN 1")
        print(int1, "\n")
        print("INTERFACE RESIDUES FROM CHAIN 2")
        print(int2, "\n")
        print("INTERACTIONS WITH NEIGHBORHOOD = 2 FROM CHAIN 1")
        print(intchain1additional, "\n")
        print("INTERACTIONS WITH NEIGHBORHOOD = 2 FROM CHAIN 2")
        print(intchain2additional)
    else:
        print("Interactions could not be detected at the distance threshold of", distance_threshold, "Angstroms")