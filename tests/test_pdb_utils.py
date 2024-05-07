import os
import pytest
import numpy as np
from biopandas.pdb import PandasPdb
from proflex.pdb import PdbDf


@pytest.fixture
def sample_pdb_file():
    this_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_file = os.path.join(this_dir, 'data', 'sample_pdb.pdb')
    return pdb_file


@pytest.fixture
def sample_atom_df(sample_pdb_file):
    return PandasPdb().read_pdb(sample_pdb_file).df['ATOM']


class TestPdbDf:
    def TestPdbDf_atoms(self, sample_pdb_file, sample_atom_df):
        result_df = PdbDf.atoms(sample_pdb_file)
        expected_df = sample_atom_df
        assert result_df.equals(expected_df)

    def TestPdbDf_ca(self, sample_atom_df):
        result_df = PdbDf.ca(sample_atom_df)
        expected_df = sample_atom_df[sample_atom_df['atom_name'] == 'CA']
        unique_atoms = np.unique(result_df['atom_name'].values)
        assert len(unique_atoms) == 1, 'More than one atom in result_df'
        assert 'CA' == unique_atoms[0]
        assert result_df.equals(expected_df)

    @pytest.mark.parametrize("sel_chain, expected_chain",
                             [('A', 'A'),
                              ('B', 'B')])
    def TestPdbDf_chains_id_IfInPDB(self, sample_atom_df, sel_chain, expected_chain):
        result_df = PdbDf.chains_id(sample_atom_df, sel_chain)
        unique_chains = np.unique(result_df['chain_id'].values)
        assert len(unique_chains) == 1, 'More than one chain in result_df'
        assert expected_chain in unique_chains[0]

    @pytest.mark.parametrize("missing_chain", ['C', 'D'])
    def TestPdbDf_chains_id_IfNotInPDB(self, sample_atom_df, missing_chain):
        with pytest.raises(ValueError):
            PdbDf.chains(sample_atom_df, missing_chain)

    def test_get_relevant_columns(self, sample_atom_df):
        pdb_utils = PdbDf()
        result_df = pdb_utils.rel_col(sample_atom_df)
        relevant_cols = ['chain_id', 'residue_number', 'residue_name', 'x_coord', 'y_coord', 'z_coord']
        assert all(col in result_df.columns for col in relevant_cols)
        assert len(result_df.columns) == len(relevant_cols)
