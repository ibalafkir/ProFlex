from proflex.utils import PDBUtils
import pandas as pd
import copy
from biopandas.pdb import PandasPdb

class PDockDf:
    
    """
    1) "_coords_only" is a pandas df of 3 columns with the coordinates of the system that PyDock moves
    
    2) "_coords" is a pandas df that contains the relevant columns considered in PDBUtils: 
                        chain_id, residue_number, residue_name, x_coord, y_coord, z_coord

    3) "rec" corresponds to the static protein

    4) "lig" corresponds to the moving protein
    """
    
    def read_transf_matr(path_to_matrix):
        """
        Reads a pydock .rot file
        """
        return pd.read_table(path_to_matrix, delim_whitespace=True, header=None, names=[
        'Column1', 'Column2', 'Column3', 'Column4', 'Column5', 'Column6',
        'Column7', 'Column8', 'Column9', 'Column10', 'Column11', 'Column12',
        'Label'])
    
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
            # print("Advancing until the rotation+translation matrix row number:", i, "...")
            new_coords = []
        return df
    
    def generate_new_poses(pdb, rotation_translation):
        """
        Generates new poses for a ligand (moving protein) by applying a (filtered or not) rot+transl matrix
        :param PDB path: str
        :param DataFrame with rotation+translation parameters: Pandas.DataFrame
        :return PDB files in path
        """
        poses = copy.deepcopy(rotation_translation)
        atom_df = PDBUtils.get_pdb_atoms_df(pdb)
        atoms_coord_only = PDBUtils.get_only_coords(atom_df)
        c = 0
        w = 0
        while len(poses)>0:
            row = poses.head(1)
            new_coords = PDockDf.rotate_and_translate(atoms_coord_only, row)
            # Logging statement "advancing till the rot+transl matrix row number: 1" because row = filtered_poses.head(1)

            new_atoms_df = PDBUtils.from_new_coord_to_wholedf(new_coords, atom_df)
            ppdb = PandasPdb().read_pdb(pdb)
            ppdb.df['ATOM'] = new_atoms_df
            ppdb.to_pdb(path= "output"+str(w)+".pdb", records = ['ATOM', 'HETATM',
                                            'OTHERS'], gz=False) # w+1 to start output numbers from 1 and not from 0
            poses.drop(c, axis=0, inplace=True)
            c += 1
            w += 1
        # print("Files generated in path")

class PDockFilter:
    
    # TODO proflex.pydock: class to filter interactions (PDockFilter)

    pass