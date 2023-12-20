"""
ProFlex PyDock poses filtration from ab-ag input proteins in different PDB files
"""
######################## INTERFACE ANALYSER CODE EXTRACTED FROM interface_analyzer.py 8de0854 ###########
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

###################################################################################################

def get_only_coords(chain):
    """
    Gets coordinates (x, y, z) from a PDB ATOM DataFrame
    :param Protein chain: pandas.DataFrame
    :return: Coordinates: pandas.DataFrame
    """
    desired_cols = ['x_coord', 'y_coord', 'z_coord']
    x = chain[desired_cols]
    return x

def from_new_coord_to_ca(newcoord, lig_coord):
    """
    Reconstructs an ATOM PDB DataFrame with the new coordinates
    :param newcoord: pandas.DataFrame
    :param lig_coord: pandas.DataFrame
    :return:
    """
    chain = copy.deepcopy(lig_coord) # so that the original chain is not affected
    chain.iloc[:, 3] = newcoord.iloc[:, 0]
    chain.iloc[:, 4] = newcoord.iloc[:, 1]
    chain.iloc[:, 5] = newcoord.iloc[:, 2]
    return chain

def rotate_and_translate(coords, rotation_translation):
    """
    Applies rotation and translation to a given set of coordinates using the parameters
    of a .rot PyDock file
    :param Coordinates dataframe: pandas.DataFrame
    :param ONE ROW of Rotation-translation dataframe. If more than 1 row are input, df will contain
    the new coordinates according to the last row of the .rot dataframe: pandas.DataFrame
    :return Dataframe with the new coordinates
    """
    new_coords = []
    i = 0
    while i <= len(rotation_translation)-1:
        for coord in coords.values:
            new_x = (rotation_translation.values[i][0] * coord[0] +
                     rotation_translation.values[i][1] * coord[1] +
                     rotation_translation.values[i][2] * coord[2] +
                     rotation_translation.values[i][9])
            new_y = (rotation_translation.values[i][3] * coord[0] +
                     rotation_translation.values[i][4] * coord[1] +
                     rotation_translation.values[i][5] * coord[2] +
                     rotation_translation.values[i][10])
            new_z = (rotation_translation.values[i][6] * coord[0] +
                     rotation_translation.values[i][7] * coord[1] +
                     rotation_translation.values[i][8] * coord[2] +
                     rotation_translation.values[i][11])
            new_coords.append([new_x, new_y, new_z])

        df = pd.DataFrame(new_coords, columns=['x_coord', 'y_coord', 'z_coord'])
        i += 1
        print("Advancing until the rotation+translation matrix row number:", i, "...")
        new_coords = []
    return df

def rotate_and_translate_and_filter(moving_coords, rotation_translation, rec_coord):


    rot_transl_result = []
    new_coords = []
    u = 0
    while u <= len(rotation_translation)-1:
        for coord in moving_coords.values:
            new_x = (rotation_translation.values[u][0] * coord[0] +
                     rotation_translation.values[u][1] * coord[1] +
                     rotation_translation.values[u][2] * coord[2] +
                     rotation_translation.values[u][9])
            new_y = (rotation_translation.values[u][3] * coord[0] +
                     rotation_translation.values[u][4] * coord[1] +
                     rotation_translation.values[u][5] * coord[2] +
                     rotation_translation.values[u][10])
            new_z = (rotation_translation.values[u][6] * coord[0] +
                     rotation_translation.values[u][7] * coord[1] +
                     rotation_translation.values[u][8] * coord[2] +
                     rotation_translation.values[u][11])
            new_coords.append([new_x, new_y, new_z]) # new_coords is a list of lists with new coordinates,
                                                     # next steps are executed when
                                                     # len(coords) == len(new_coords)

        df = pd.DataFrame(new_coords, columns=['x_coord', 'y_coord', 'z_coord'])
        new_lig = from_new_coord_to_ca(df, lig_coord)
        print("Generating new chain from new coordinates...")

        # Analyzing the interface
        int1, int2 = get_interface_residues_by_chain(rec_coord, new_lig, 6)
        print("Calculating interfaces...")

        # Starting filtration criteria: presence of any antibody CDR residue in the interface
        cdr_res = np.array(list(range(27,39)) + list(range(56,66)) + list(range(105,118)))
        for w in cdr_res:
            if w in int2.values[:, 1]:
                rot_transl_result.append(rotation_translation.values[u])
                break
        print("verified filtration criterion")
        # finished filtration criteria

        u += 1
        print("Advancing till the rot+transl matrix row number:", u) # logging
        new_coords = []

    column_names = [
        'Column1', 'Column2', 'Column3', 'Column4', 'Column5', 'Column6',
        'Column7', 'Column8', 'Column9', 'Column10', 'Column11', 'Column12',
        'Label'
    ]
    rot_transl_result_df = pd.DataFrame(rot_transl_result, columns=column_names)
    return rot_transl_result_df

def from_new_coord_to_wholedf(newcoord, chain):
    """
    :param newcoord: Pandas.DataFrame
    :param chain: Pandas.DataFrame
    :return dataframe containing all atom information with new coordinates: Pandas.DataFrame
    """
    chain.iloc[:, 11] = newcoord.iloc[:, 0]
    chain.iloc[:, 12] = newcoord.iloc[:, 1]
    chain.iloc[:, 13] = newcoord.iloc[:, 2]
    return chain

def generate_new_poses(pdb_file_lig, filtered_poses):
    """
    Generates new poses for a ligand by applying a (filtered or not) rot+transl matrix
    :param PDB path: str
    :param DataFrame with rotation+translation parameters: Pandas.DataFrame
    :return PDB files in path
    """
    f_poses = copy.deepcopy(filtered_poses)
    lig_atoms = get_pdb_atoms_df(pdb_file_lig)
    lig_atoms_coord = get_only_coords(lig_atoms)
    c = 0
    w = 0
    while len(f_poses)>0:
        row = f_poses.head(1)
        new_coords = rotate_and_translate(lig_atoms_coord, row)
        # Always "Log: advancing till the rot+transl matrix row number: 1" bc row = filtered_poses.head(1)
        # This can be solved by deleting the logging action in functions or the index "i"
        new_lig = from_new_coord_to_wholedf(new_coords, lig_atoms)
        ppdb = PandasPdb().read_pdb(pdb_file_lig)
        ppdb.df['ATOM'] = new_lig
        ppdb.to_pdb(path= "output"+str(w)+".pdb", records = ['ATOM', 'HETATM',
                                           'OTHERS'], gz=False) # w+1 to start output numbers from 1 and not from 0
        f_poses.drop(c, axis=0, inplace=True)
        c += 1
        w += 1
    print("Files generated in path")
###################################################################################################

# Loading files
pdb_file = 'B_PD1_diff_8_prep_relax_fixed.pdb'
pdb_file_rec = "apo_pd1_6umv_noMt_optH_minH-DOCK-Q15116_20_1.51_igfold_imgt_rec.pdb"
pdb_file_lig = "apo_pd1_6umv_noMt_optH_minH-DOCK-Q15116_20_1.51_igfold_imgt_lig.pdb"
rot_file_pydock = "apo_pd1_6umv_noMt_optH_minH-DOCK-Q15116_20_1.51_igfold_imgt.rot"
rotation_translation = pd.read_table(rot_file_pydock, delim_whitespace=True, header=None, names=[
    'Column1', 'Column2', 'Column3', 'Column4', 'Column5', 'Column6',
    'Column7', 'Column8', 'Column9', 'Column10', 'Column11', 'Column12',
    'Label'
])

## STARTING POINT: two PDB files with the ligand and receptor each
# Receptor
rec_atoms = get_pdb_atoms_df(pdb_file_rec)
rec_ca = get_ca(rec_atoms)
rec_coord = get_relevant_columns(rec_ca)

# Ligand
lig_atoms = get_pdb_atoms_df(pdb_file_lig)
lig_ca = get_ca(lig_atoms)
lig_coord = get_relevant_columns(lig_ca)
lig_coord_only = get_only_coords(lig_coord)

# Analyzing interface, rotating and translating the ligand, filtering poses and generating the corresponding PDB files
rttest = rotation_translation.head(20)
filtered_poses = rotate_and_translate_and_filter(lig_coord_only, rttest, rec_coord)
generate_new_poses(pdb_file_lig, filtered_poses)