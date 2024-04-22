from proflex.utils import PDBProcessor, PDBUtils
from rfdiffusion import RFDContigs
from biopandas.pdb import PandasPdb
import argparse
import os

"""
Gets the number of residues in a PDB
Assumes the input PDB has no alt_locations 
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Returns total number of residues')
    parser.add_argument('pdb', type=str, help='Path to the PDB file')
    args = parser.parse_args()
    pd = args.pdb
    atom_df = PDBUtils.get_pdb_atoms_df(pd)
    atom_df_ca = atom_df[atom_df['atom_name']=='CA']
    print(f'Number of residues assuming the inexistence of alternative locations and that there is a CA per residue: {len(atom_df_ca)}')
