"""
ProFlex PyDock poses filtration from protein-protein input proteins in different PDB files
"""

from proflex.utils import PDBUtils, InterfaceAnalyzer
import pandas as pd
import numpy as np
from proflex.pydock import PDockDf
import argparse

"""
1) "_coords_only" is a pandas df of 3 columns with the coordinates of the system that PyDock moves
    
2) "_coords" is a pandas df that contains the relevant columns considered in PDBUtils: 
                        chain_id, residue_number, residue_name, x_coord, y_coord, z_coord

3) "rec" corresponds to the static protein

4) "lig" corresponds to the moving protein

Example:

python pydock_transformation.py --rec_pdb apo_pd1_6umv_noMt_optH_minH-DOCK-Q15116_20_1.51_igfold_imgt_rec.pdb 
--lig_pdb apo_pd1_6umv_noMt_optH_minH-DOCK-Q15116_20_1.51_igfold_imgt_lig.pdb 
--transf_matrix apo_pd1_6umv_noMt_optH_minH-DOCK-Q15116_20_1.51_igfold_imgt.rot

"""

def rotate_and_translate_and_filter(moving_coords, moving_coords_only, rotation_translation, static_coords):
    
    rot_transl_result = []
    new_moving_coords_only = []
    u = 0
    while u <= len(rotation_translation)-1:
        for coord in moving_coords_only.values:
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
            new_moving_coords_only.append([new_x, new_y, new_z]) # new_coords is a list of lists with new coordinates,
                                                        # next steps are executed when
                                                        # len(coords) == len(new_coords)

        df = pd.DataFrame(new_moving_coords_only, columns=['x_coord', 'y_coord', 'z_coord'])
        new_moving_coords = PDBUtils.from_new_coord_to_ca(df, moving_coords)
        print("Generating new chain from new coordinates...")

        # Analyzing the interface (int1 and detected_interactions dataframes are unused)
        # The "int2" pandas dataframe contains the interacting residues from the moving protein
        int1, int2, detected_interactions= InterfaceAnalyzer.get_interface_residues_by_chain(static_coords, new_moving_coords, d)
        
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
        print("Advancing till the rot+transl matrix row number:", u)
        new_moving_coords_only = []

    column_names = [
        'Column1', 'Column2', 'Column3', 'Column4', 'Column5', 'Column6',
        'Column7', 'Column8', 'Column9', 'Column10', 'Column11', 'Column12',
        'Label'
    ]
    rot_transl_result_df = pd.DataFrame(rot_transl_result, columns=column_names)
    return rot_transl_result_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Returns new poses of a protein in a protein-protein interaction'
        )
    parser.add_argument('--rec_pdb', type=str, help='Path to the static protein PDB file')
    parser.add_argument('--lig_pdb', type=str, help='Path to the moving protein PDB file')
    parser.add_argument('--transf_matrix', type=str, help='Path to the PyDock transformation matrix')
    parser.add_argument('--distance_threshold', type=float,
                        help='Distance threshold for interface '
                             'residues (6 Angstrom is default)', default=6)
    args = parser.parse_args()
    pdb_file_rec = args.rec_pdb
    pdb_file_lig = args.lig_pdb
    rot_file_pydock = args.transf_matrix
    d = args.distance_threshold 

    rotation_translation = PDockDf.read_transf_matr(rot_file_pydock)

    rec_atoms = PDBUtils.get_pdb_atoms_df(pdb_file_rec)
    rec_ca = PDBUtils.get_ca(rec_atoms)
    rec_coord = PDBUtils.get_relevant_columns(rec_ca)

    lig_atoms = PDBUtils.get_pdb_atoms_df(pdb_file_lig)
    lig_ca = PDBUtils.get_ca(lig_atoms)
    lig_coord = PDBUtils.get_relevant_columns(lig_ca)
    lig_coord_only = PDBUtils.get_only_coords(lig_coord)

    # Analyzing interface, rotating and translating the ligand, filtering poses and generating the corresponding PDB files
    
    rttest = rotation_translation.head(20) # TODO must be removed when integrated into the iteration process, now it's
                                           # here for testing purposes

    filtered_poses = rotate_and_translate_and_filter(lig_coord, lig_coord_only, rttest, rec_coord)
    PDockDf.generate_new_poses(pdb_file_lig, filtered_poses)