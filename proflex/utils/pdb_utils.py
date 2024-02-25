from biopandas.pdb import PandasPdb
import pandas as pd

class PDBUtils:
    @staticmethod
    def get_pdb_atoms_df(pdb_file):
        """
        Extracts ATOM information from a PDB file and returns a pandas DataFrame
        :param Path to PDB file: str
        :return Pandas DataFrame: pandas.DataFrame
        """
        ppdb = PandasPdb().read_pdb(pdb_file)
        # print("Extracting atom information from PDB... \n")
        return ppdb.df['ATOM'] # only keeps information related to atoms
    
    @staticmethod
    def get_chains_id(atom_df):
        result=[]
        array = pd.unique(atom_df["chain_id"])
        for i in array:
            result.append(i)
        # print("Obtaining chains' identification")
        # print(" ")
        return result

    @staticmethod
    def get_ca(atom_df):
        """
        Extracts CA rows from a pandas DataFrame with ATOM information
        :param Pandas DataFrame with ATOM information: pandas.DataFrame
        :return Pandas DataFrame with CA information: pandas.DataFrame
        """
        # print("Getting alpha carbons... \n")
        return atom_df[atom_df['atom_name'] == 'CA']

    @staticmethod
    def get_chain(df, chain_name):
        """
        Returns pandas DataFrame with the ATOM information of an input chain
        :param Pandas DataFrame with ATOM information (CA atom_names, all atom_names...):
               pandas.DataFrame
        :return Pandas DataFrame with ATOM information of an indicated chain:
                pandas.DataFrame
        """
        # print("Obtaining different protein chains...:", chain_name, "\n")
        result = df[df['chain_id'] == chain_name]
        if len(result) == 0: # if the input chain_name is absent, the df will have size 0
            raise ValueError('Error: the input chain name is absent')
        return result

    @staticmethod
    def get_relevant_columns(chain):
        """
        Extracts "relevant columns" from an ATOM PDB pandas DatFrame. These are:
        chain_id, residue_number, residue_name, x_coord, y_coord, z_coord
        :param Atom DataFrame: pandas.DataFrame
        :return Atom DataFrame with relevant columns
        """
        desired_cols = ['chain_id', 'residue_number', 'residue_name', 'x_coord', 'y_coord', 'z_coord']
        x = chain[desired_cols]
        # print("Filtering relevant columns of chain", chain['chain_id'].values[0], "... \n")
        return x
