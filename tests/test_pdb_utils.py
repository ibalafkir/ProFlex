import os
import pytest
import numpy as np
from biopandas.pdb import PandasPdb
from proflex.utils import PDBUtils


@pytest.fixture
def sample_pdb_file():
    this_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_file = os.path.join(this_dir, 'data', 'sample_pdb.pdb')
    return pdb_file


@pytest.fixture
def sample_atom_df(sample_pdb_file):
    return PandasPdb().read_pdb(sample_pdb_file).df['ATOM']


class TestPDBUtils:
    def test_get_pdb_atoms_df(self, sample_pdb_file, sample_atom_df):
        result_df = PDBUtils.get_pdb_atoms_df(sample_pdb_file)
        expected_df = sample_atom_df
        assert result_df.equals(expected_df)

    def test_get_ca(self, sample_atom_df):
        result_df = PDBUtils.get_ca(sample_atom_df)
        expected_df = sample_atom_df[sample_atom_df['atom_name'] == 'CA']
        unique_atoms = np.unique(result_df['atom_name'].values)
        assert len(unique_atoms) == 1, 'More than one atom in result_df'
        assert 'CA' == unique_atoms[0]
        assert result_df.equals(expected_df)

    def test_delete_res_tag(self, sample_atom_df):
        # TODO implement test_delete_res_tag
        pass

    @pytest.mark.parametrize("sel_chain, expected_chain",
                             [('A', 'A'),
                              ('B', 'B')])
    def test_get_chain_when_its_in_pdb(self, sample_atom_df, sel_chain, expected_chain):
        result_df = PDBUtils.get_chain(sample_atom_df, sel_chain)
        unique_chains = np.unique(result_df['chain_id'].values)
        assert len(unique_chains) == 1, 'More than one chain in result_df'
        assert expected_chain in unique_chains[0]

    @pytest.mark.parametrize("missing_chain", ['C', 'D'])
    def test_get_chain_when_its_not_in_pdb(self, sample_atom_df, missing_chain):
        # TODO Update PDBUtils.get_chain to raise a ValueError when the chain is not in the PDB
        with pytest.raises(ValueError):
            PDBUtils.get_chain(sample_atom_df, missing_chain)

    def test_get_relevant_columns(self, sample_atom_df):
        pdb_utils = PDBUtils()
        result_df = pdb_utils.get_relevant_columns(sample_atom_df)
        relevant_cols = ['chain_id', 'residue_number', 'residue_name', 'x_coord', 'y_coord', 'z_coord']
        assert all(col in result_df.columns for col in relevant_cols)
        assert len(result_df.columns) == len(relevant_cols)
