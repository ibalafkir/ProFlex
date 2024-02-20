import argparse
from proflex.utils import PDBAdaptor

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Creates a new PDB file fixed for insertion codes'
        )
    parser.add_argument('-pdb', type=str, help='Path to the PDB file')

    args = parser.parse_args()
    
    pdb_file = args.pdb 
    PDBAdaptor.fix_insertions(pdb_file)
    PDBAdaptor.fix_ter_mistakes(PDBAdaptor.extract_without_extension(pdb_file)+'_fixed.pdb')