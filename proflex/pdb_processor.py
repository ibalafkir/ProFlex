import argparse
from proflex.utils import PDBProcessor

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Creates a new PDB file fixed for insertion codes'
        )
    parser.add_argument('-pdb', type=str, help='Path to the PDB file')

    args = parser.parse_args()
    
    pdb_file = args.pdb
    PDBProcessor.fix_insertions(pdb_file)
    PDBProcessor.fix_ter_mistakes(PDBProcessor.extract_without_extension(pdb_file)+'_insfixed.pdb')
    PDBProcessor.remove_other_lines(PDBProcessor.extract_without_extension(pdb_file)+'_insfixed'+'_terfixed.pdb')
    PDBProcessor.del_mid_files(pdb_file)