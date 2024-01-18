from biopandas.pdb import PandasPdb


class PDBUtils:
    def get_pdb_atoms_df(pdb_file):
        """
        Extracts ATOM information from a PDB file and returns a pandas DataFrame
        :param Path to PDB file: str
        :return Pandas DataFrame: pandas.DataFrame
        """
        ppdb = PandasPdb().read_pdb(pdb_file)
        print("Extracting atom information from PDB... \n")
        return ppdb.df['ATOM'] # only keeps information related to atoms

    @staticmethod
    def get_ca(atom_df):
        """
        Extracts CA rows from a pandas DataFrame with ATOM information
        :param Pandas DataFrame with ATOM information: pandas.DataFrame
        :return Pandas DataFrame with CA information: pandas.DataFrame
        """
        print("Getting alpha carbons... \n")
        return atom_df[atom_df['atom_name'] == 'CA']

    @staticmethod
    def delete_res_tag(df):
        def delete_numbers_in_string(string):
            if isinstance(string, str):
                return ''.join(c for c in string if c.isdigit())
            else:
                return str(string)

        # Uses .loc to avoid SettingWithCopyWarning
        df.loc[:, 'residue_number'] = df['residue_number'].apply(delete_numbers_in_string)
        df['residue_number'] = df['residue_number'].astype(int)
        return df

    @staticmethod
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
        print("Filtering relevant columns of chain", chain['chain_id'].values[0], "... \n")
        return x