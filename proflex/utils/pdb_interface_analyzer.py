import numpy as np
import pandas as pd

class InterfaceAnalyzer:
    """
    Functions to extract interface residues from PDB files
    """

    def get_interface_residues_by_chain(chain1: str, chain2: str, distance_threshold=6):
        """
        Extraction of interface residues, meaning residues which CA are at a maximum given distance
        (default = 6 Angstrom), which is the maximum distance up to which interaction is considered
        to occur
        :param First chain ID: str
        :param Second chain ID: str
        :param Distance threshold: float
        :return 2 np.arrays containing interface residues from each chain, 1 np.array containing interacting pairs
        """
        coords_1 = chain1[['x_coord', 'y_coord', 'z_coord']].values
        coords_2 = chain2[['x_coord', 'y_coord', 'z_coord']].values
        distances = np.linalg.norm(coords_1[:, None] - coords_2, axis=2)
        close_residues = np.where(distances <= distance_threshold)
        residues_interface_1 = chain1.iloc[close_residues[0]]
        residues_interface_2 = chain2.iloc[close_residues[1]]
        detected_interactions = pd.DataFrame({
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
        # print("Calculating distances... \n")
        return residues_interface_1, residues_interface_2, detected_interactions

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
        # print("Calculating distances... \n")
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
        # print("Detecting additional interactions (neighborhood = 2) residues for diffusion processes of chain", chainrelevant['chain_id'].values[0], "... \n")
        return int_sorted_unique
