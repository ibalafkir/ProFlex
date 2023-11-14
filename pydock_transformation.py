"""
ProFlex PyDock poses filtration
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
def get_relevant_columns(chain):
    desired_cols = ['chain_id', 'residue_number', 'residue_name', 'x_coord', 'y_coord', 'z_coord']
    x = chain[desired_cols]
    print("Filtering relevant columns of chain", chain['chain_id'].values[0])
    print(" ")
    return x
def get_only_coords(chain):
    desired_cols = ['x_coord', 'y_coord', 'z_coord']
    x = chain[desired_cols]
    print("Getting only coordinates")
    print(" ")
    return x

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

############ good till here
def rotate_and_translate(coords, rotation_translation): # FALLA ESTA FUNCION
                                                        # Me quede por ver como sacar columna entera
    """Applies rotation and translation matrix to a set of coordinates."""
    new_coords = []
    i = 0
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
        new_coords.append((new_x, new_y, new_z))
        i += 1
        df = pd.DataFrame(new_coords, columns=['x_coord', 'y_coord', 'z_coord'])
    return df
lig_coord_after_rot_trans = rotate_and_translate(lig_coord_only, rotation_translation)