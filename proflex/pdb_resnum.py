from proflex.utils import PDBProcessor, PDBUtils
from rfdiffusion import RFDContigs
from biopandas.pdb import PandasPdb
import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Returns total number of residues')
    parser.add_argument('pdb', type=str, help='Path to the PDB file')
    args = parser.parse_args()
    pd = args.pdb
    PDBProcessor.pdb_resnum(pd, 1)
    atom_df = PDBUtils.get_pdb_atoms_df(pd[:-4]+'_renum.pdb')
    atom_df_relevant = PDBUtils.get_relevant_columns(atom_df)
    res_num_list = RFDContigs.get_resnum_list(atom_df_relevant)
    print(res_num_list[-1])
    os.remove(pd[:-4]+'_renum.pdb')