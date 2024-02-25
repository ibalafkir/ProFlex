from proflex.utils import RFDFixer
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Fixes PDB features missed in the output of RFDiffusion'
        )
    parser.add_argument('-prerfd', type=str, help='Path to the PDB file (input of RFD)')
    parser.add_argument('-postrfd', type=str, help='Path to the output PDB RFDiffusion file')
    
    args = parser.parse_args()
    
    pre = args.prerfd
    post = args.postrfd

    n_atoms = RFDFixer.get_n_atoms(post) # RFD output PDBs have an ordered atom number so
                                            # the function can be well applied
    RFDFixer.pdb_backbone(pre) 
    RFDFixer.pdb_atom(pre[:len(pre)-4]+'_backbone.pdb')
    RFDFixer.correct_rfd_pdbs(post, pre[:len(pre)-4]+'_backbone.pdb')
    RFDFixer.pdb_tidying(post[:len(post)-4]+'_chainsfixed.pdb')
    RFDFixer.del_mid_files(pre, post)