"""
Generates the codes that indicate the regions to diffuse in RFDiffusion
Input from commandline: PDB with both proteins interacting and the chains to analyze interface
"""
# TODO Solve issue that happens when the first or last residue of a chain needs to be diffunded

## Functions from interface_analyzer 8de0854
import argparse
from biopandas.pdb import PandasPdb
import numpy as np
import pandas as pd
import copy

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

def get_resnum_list(df):
    result = []
    for i in df['residue_number']:
        result.append(int(i))
    return result

def get_chunks(lst):
    result = []
    subgroup = []
    for num in lst:
        if not subgroup or num == subgroup[-1] + 1:
            subgroup.append(num)
        else:
            result.append(subgroup)
            subgroup = [num]
    if subgroup:
        result.append(subgroup)
    for i in result: # quitar
        if len(i)==1: # si
            result.remove(i) # no quieres residuos solos en plan [1,2,3], [5], [8, 9, 10], no cuenta el 5
    return result

def chunk_filter(lst):
    """
    :param lst: list
    :return list without elements of length < 3 so as not to diffuse regions of less than 3 aas: list
    """
    result = []
    for i in lst:
        if len(i)>2:
            result.append(i)
    return result

def get_contigs_rfd(chain_1, chain_2, intchain1additional, intchain2additional):
    """
    id1 and id2 are indicated from terminal.
    :param chain_1:
    :param chain_2:
    :param intchain1additional:
    :param intchain2additional:
    :return:
    """

    chain1_start = get_resnum_list(chain_1)[0]
    chain1_end = get_resnum_list(chain_1)[len(chain_1)-1]
    chain2_start = get_resnum_list(chain_2)[0]
    chain2_end = get_resnum_list(chain_2)[len(chain_2)-1]
    intera1_resnum = get_resnum_list(intchain1additional)
    intera2_resnum = get_resnum_list(intchain2additional)
    intera1_lst = get_chunks(intera1_resnum)
    intera2_lst = get_chunks(intera2_resnum)
    intera1_lst = chunk_filter(intera1_lst)
    intera2_lst = chunk_filter(intera2_lst)

    contig = f"{id1}{chain1_start}-"
    for curr_list in intera1_lst:
        list_start = curr_list[0]
        list_end = curr_list[len(curr_list)-1]
        contig += f"{list_start-1}/{list_end-list_start+1}-{list_end-list_start+1}/{id1}{list_end+1}-"
    contig = contig[:len(contig)-1] + f"-{chain1_end}" + "/0 "

    contig += f"{id2}{chain2_start}-"
    for curr_list in intera2_lst:
        list_start = curr_list[0]
        list_end = curr_list[len(curr_list)-1]
        contig += f"{list_start-1}/{list_end-list_start+1}-{list_end-list_start+1}/{id2}{list_end+1}-"
    contig += f"{chain2_end}"
    return contig

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Returns diffusion code from a PDB file')
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

    ####### CODE
    atom_df = get_pdb_atoms_df(pdb_file)
    atom_df_ca = get_ca(atom_df)
    chain_ids = get_chains_id(atom_df_ca)
    chain_1 = get_chain(atom_df_ca, id1)
    chain_2 = get_chain(atom_df_ca, id2)
    chain1_relevant = get_relevant_columns(chain_1)
    chain2_relevant = get_relevant_columns(chain_2)
    int1, int2 = get_interface_residues_by_chain(chain1_relevant, chain2_relevant, distance_threshold=distance_threshold)
    intchain1additional = amplify_selection_residues(int1, chain1_relevant)
    intchain2additional = amplify_selection_residues(int2, chain2_relevant)
    print(get_contigs_rfd(chain_1, chain_2, intchain1additional, intchain2additional))
    #######




