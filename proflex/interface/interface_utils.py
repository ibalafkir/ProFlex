"""
ProFlex interface analyzing functionalities    
"""
import numpy as np
import pandas as pd
from proflex.pdb import PdbDf

class InterfaceAnalyzer:
    """
    PDB files interface residues extractor
    """
    def calculate(chain1: str, chain2: str, distance_threshold=8):
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

    def distances(chain1: str, chain2: str, distance_threshold=6):
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

    def extend_neighbourhood(int, chainrelevant):
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

    def run(pdb_file, id1, id2, distance_threshold):

        atom_df = PdbDf.atoms(pdb_file)
        atom_df_ca = PdbDf.ca(atom_df)
        chain_1 = PdbDf.chain(atom_df_ca, id1)
        chain_2 = PdbDf.chain(atom_df_ca, id2)
        chain1_relevant = PdbDf.rel_col(chain_1)
        chain2_relevant = PdbDf.rel_col(chain_2)
        int1, int2, detected_interactions = InterfaceAnalyzer.calculate(chain1_relevant, chain2_relevant, distance_threshold)
        intchain1_extended = InterfaceAnalyzer.extend_neighbourhood(int1, chain1_relevant)
        intchain2_extended = InterfaceAnalyzer.extend_neighbourhood(int2, chain2_relevant)
        return int1, int2, detected_interactions, intchain1_extended, intchain2_extended


    def pele_com_distance(int1, int2, id1, id2):
        """
        # TODO Create a PELETask builder in a class and include this
        From 2 pandas df inputs (each containing residues involved in interactions) and their chain IDs,
        it obtains the way this info must be introduced in the com_distance PELE metric
        """
        
        def del_repeated(lst):
            """
            Deletes repeated elements in a list
            """
            result = []
            for i in lst:
                if i not in result:
                    result.append(i)
                else:
                    continue
            return result
        int1_resnum = int1['residue_number'].to_list()
        int1_resnum_sorted_unique = sorted(del_repeated(int1_resnum))
        int1_resnum_sorted_unique_chname = [f"{id1}:" + str(element) for element in int1_resnum_sorted_unique]
        int1_resnum_sorted_unique_chname_mod = ['"' + element + '"' for element in int1_resnum_sorted_unique_chname]
        pele_code_1 = ', '.join(int1_resnum_sorted_unique_chname_mod)
        
        int2_resnum = int2['residue_number'].to_list()
        int2_resnum_sorted_unique = sorted(del_repeated(int2_resnum))
        int2_resnum_sorted_unique_chname = [f"{id2}:" + str(element) for element in int2_resnum_sorted_unique]
        int2_resnum_sorted_unique_chname_mod = ['"' + element + '"' for element in int2_resnum_sorted_unique_chname]
        pele_code_2 = ', '.join(int2_resnum_sorted_unique_chname_mod)
        
        print("{")
        print('"type":"com_distance",')
        print(f'"tag":"{id1}{id2}_distance",')
        print('\t"selection_group_1":{')
        print('\t\t"links": { "ids":[' + str(pele_code_1) + ']}')
        print('\t\t},')
        print('\t"selection_group_2":{')
        print('\t\t"links": { "ids":[' + str(pele_code_2) + ']}')
        print('\t\t}')
        print("},")  


    def pele_interacting_res(int1, id1, int2, id2, int3=None, id3=None, int4=None, id4=None):
        """
        Int1 should be the selected residues of the first ab or nb chain
        Int2 should be the selected residues of the antigen with respect to the previous chain
        Int3 should be the selected residues of the second ab 
        Int4 should be the selected residues of the antigen with respect to the previous chain
        
        Int2 and Int4 are merged to avoid repetitions
        
        """
        def del_repeated(lst):
            """
            Deletes repeated elements in a list
            """
            result = []
            for i in lst:
                if i not in result:
                    result.append(i)
                else:
                    continue
            return result
                
        def extend_neighbourhood_ag(lst):
            result = []
            for i in range(len(lst) - 1):
                result.append(lst[i])
                dif = lst[i+1] - lst[i]
                if 2 <= dif <= 4:
                    for j in range(1, dif):
                        result.append(lst[i] + j)
            result.append(lst[-1])  # Añadimos el último elemento de la lista original
            return result
        
        # Nb-Ag case (1 interface)
        if int3 is None and id3 is None and int4 is None and id4 is None:
            
            int1_resnum = int1['residue_number'].to_list()
            int1_resnum_sorted_unique = sorted(del_repeated(int1_resnum))
            int1_resnum_sorted_unique_chname = [f"{id1}:" + str(element) for element in int1_resnum_sorted_unique]
            int1_resnum_sorted_unique_chname_mod = ['"' + element + '"' for element in int1_resnum_sorted_unique_chname]
            pele_code_1 = ', '.join(int1_resnum_sorted_unique_chname_mod)
            
            int2_resnum = int2['residue_number'].to_list()
            int2_resnum_sorted_unique = sorted(del_repeated(int2_resnum))
            int2_resnum_sorted_unique_chname = [f"{id2}:" + str(element) for element in int2_resnum_sorted_unique]
            int2_resnum_sorted_unique_chname_mod = ['"' + element + '"' for element in int2_resnum_sorted_unique_chname]
            pele_code_2 = ', '.join(int2_resnum_sorted_unique_chname_mod)
            
            print('"links": { "ids":[' + str(pele_code_1) + ', ' + str(pele_code_2) + ']},')
            
        # Ab-Ag case (2 interfaces)
        else:
            int1_resnum = int1['residue_number'].to_list()
            int1_resnum_sorted_unique = sorted(del_repeated(int1_resnum))
            int1_resnum_sorted_unique_chname = [f"{id1}:" + str(element) for element in int1_resnum_sorted_unique]
            int1_resnum_sorted_unique_chname_mod = ['"' + element + '"' for element in int1_resnum_sorted_unique_chname]
            pele_code_1 = ', '.join(int1_resnum_sorted_unique_chname_mod)
            
            int2_resnum = int2['residue_number'].to_list()
            int2_resnum_sorted_unique = sorted(del_repeated(int2_resnum))
            int4_resnum = int4['residue_number'].to_list()
            int4_resnum_sorted_unique = sorted(del_repeated(int4_resnum))
            
            # Merging antigen selected residues
            newIntAg_resnum = []
            for i in int2_resnum_sorted_unique:
                if i not in newIntAg_resnum:
                    newIntAg_resnum.append(i)
                else:
                    continue
            for i in int4_resnum_sorted_unique:
                if i not in newIntAg_resnum:
                    newIntAg_resnum.append(i)
                else:
                    continue
            newIntAg_resnum_sorted_unique = sorted(newIntAg_resnum)
            
            ## Extended the antigen interface by neighbourhood = 2
            newIntAg_resnum_sorted_unique = extend_neighbourhood_ag(newIntAg_resnum_sorted_unique)
            ##
            
            newIntAg_resnum_sorted_unique_chname = [f"{id2}:" + str(element) for element in newIntAg_resnum_sorted_unique]
            newIntAg_resnum_sorted_unique_chname_mod = ['"' + element + '"' for element in newIntAg_resnum_sorted_unique_chname]
            pele_code_2 = ', '.join(newIntAg_resnum_sorted_unique_chname_mod)
            
            
            int3_resnum = int3['residue_number'].to_list()
            int3_resnum_sorted_unique = sorted(del_repeated(int3_resnum))
            int3_resnum_sorted_unique_chname = [f"{id3}:" + str(element) for element in int3_resnum_sorted_unique]
            int3_resnum_sorted_unique_chname_mod = ['"' + element + '"' for element in int3_resnum_sorted_unique_chname]
            pele_code_3 = ', '.join(int3_resnum_sorted_unique_chname_mod)
            
            
            
            
            print('"links": { "ids":[' + str(pele_code_1) + ', ' + str(pele_code_2) + ', ' + str(pele_code_3) + ']}')