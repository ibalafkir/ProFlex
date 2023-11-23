"""
ProFlex PyDock poses filtration
"""
######################## INTERFACE ANALYSER CODE EXTRACTED FROM interface_analyser.py eacbc3c ###########
import argparse
from biopandas.pdb import PandasPdb
import numpy as np
import pandas as pd
import copy
def get_pdb_atoms_df(pdb_file):
    ppdb = PandasPdb().read_pdb(pdb_file)
    return ppdb.df['ATOM'] # only keeps information related to atoms
def get_ca(atom_df):
    return atom_df[atom_df['atom_name'] == 'CA']
def get_chain(atom_df, chain_name):
    print("Obtaining different protein chains:", chain_name)
    print(" ")
    return atom_df[atom_df['chain_id'] == chain_name]
def get_relevant_columns(chain):
    desired_cols = ['chain_id', 'residue_number', 'residue_name', 'x_coord', 'y_coord', 'z_coord']
    x = chain[desired_cols]
    return x

def get_interface_residues_by_chain(chain1: str, chain2: str, distance_threshold=6) -> (np.array, np.array, np.array):
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
    return residues_interface_1, residues_interface_2

def amplify_selection_residues(int, chainrelevant):
    i = 0
    total_dif = []
    res_num = sorted(int['residue_number'].values)
    while i < len(res_num)-1:
        subs = res_num[i+1]-res_num[i]
        total_dif.append(subs)
        if subs == 2:
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
    return int_sorted_unique

def get_interface_residues_by_each_chain(pdb_file, distance_threshold=6) -> (np.array, np.array, np.array, np.array):
    atom_df = get_pdb_atoms_df(pdb_file)
    atom_df_ca = get_ca(atom_df)
    chain_ids = get_chains_id(atom_df_ca)
    chain_1 = get_chain(atom_df_ca, chain_ids[0])
    chain_2 = get_chain(atom_df_ca, chain_ids[1])
    chain1_relevant = get_relevant_columns(chain_1)
    chain2_relevant = get_relevant_columns(chain_2)
    int1, int2 = get_interface_residues_by_chain(chain1_relevant, chain2_relevant, distance_threshold)
    intchain1additional = amplify_selection_residues(int1, chain1_relevant)
    intchain2additional = amplify_selection_residues(int2, chain2_relevant)
    return intchain1additional, intchain2additional

###################################################################################################

def rotate_and_translate(coords, rotation_translation): ### This function works: gives new coordinates from a rot+tran
                                                        # matrix ###
    """Applies rotation and translation matrix to a set of coordinates."""
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

            if len(new_coords) == len(coords):
                df = pd.DataFrame(new_coords, columns=['x_coord', 'y_coord', 'z_coord'])
                i += 1
                print("Log: advancing till the rot+transl matrix row number:", i)
                new_coords = []
    return df

def rotate_and_translate_and_filter(coords, rotation_translation):
    rot_transl_result = []
    new_coords = []
    u = 0
    while u <= len(rotation_translation)-1:
        for coord in coords.values:
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
            new_coords.append([new_x, new_y, new_z])

        # new_coords is a list of lists with new coordinates, next steps are executed when
        # len(coords) == len(new_coords)

        df = pd.DataFrame(new_coords, columns=['x_coord', 'y_coord', 'z_coord'])

        new_lig = from_new_coord_to_ca(df, lig_coord)
        print("generating new chain from new_coord")

        int1, int2 = get_interface_residues_by_chain(rec_coord, new_lig, 6)

        print("calculated interfaces")
        # starting filtration criteria
        cdr3_res = [105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117]
                # res 105 --- res 118 belong to cd3
        for w in cdr3_res:
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

def get_only_coords(chain):
    desired_cols = ['x_coord', 'y_coord', 'z_coord']
    x = chain[desired_cols]
    print("Getting only coordinates")
    print(" ")
    return x

def from_new_coord_to_ca(newcoord, lig_coord):
    chain = copy.deepcopy(lig_coord) # so that the original chain is not affected (doesnt matter...)
    chain.iloc[:, 3] = newcoord.iloc[:, 0]
    chain.iloc[:, 4] = newcoord.iloc[:, 1]
    chain.iloc[:, 5] = newcoord.iloc[:, 2]
    return chain

def from_new_coord_to_wholedf(newcoord, lig):
    lig.iloc[:, 11] = newcoord.iloc[:, 0]
    lig.iloc[:, 12] = newcoord.iloc[:, 1]
    lig.iloc[:, 13] = newcoord.iloc[:, 2]
    return lig

def generate_new_poses(pdb_file_lig, filtered_poses):
    """
    Generates new poses for a ligand by applying a (filtered or not) rot+transl matrix
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

###### Testing rotate_and_translate and rotate_and_translate_and_filter #####

pdb_file_rec = "apo_pd1_6umv_noMt_optH_minH-DOCK-Q15116_20_1.51_igfold_imgt_rec.pdb"
pdb_file_lig = "apo_pd1_6umv_noMt_optH_minH-DOCK-Q15116_20_1.51_igfold_imgt_lig.pdb"
rot_file_pydock = "apo_pd1_6umv_noMt_optH_minH-DOCK-Q15116_20_1.51_igfold_imgt.rot"

rec_atoms = get_pdb_atoms_df(pdb_file_rec)
rec_ca = get_ca(rec_atoms)
rec_coord = get_relevant_columns(rec_ca)

lig_atoms = get_pdb_atoms_df(pdb_file_lig)
lig_ca = get_ca(lig_atoms)
lig_coord = get_relevant_columns(lig_ca)
lig_coord_only = get_only_coords(lig_coord)

rotation_translation = pd.read_table(rot_file_pydock, delim_whitespace=True, header=None, names=[
    'Column1', 'Column2', 'Column3', 'Column4', 'Column5', 'Column6',
    'Column7', 'Column8', 'Column9', 'Column10', 'Column11', 'Column12',
    'Label'
])
rttest = rotation_translation.head(50)
filtered_poses = rotate_and_translate_and_filter(lig_coord_only, rttest)
generate_new_poses(pdb_file_lig, filtered_poses)