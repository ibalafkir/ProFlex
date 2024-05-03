"""
Solves the PDB format of an RFdiffusion model output (which is all backbone atoms of all residues in a single chain)
using reference input PDB to RFdiffusion and the diffusion models starting by {pdb_code}_diff

"""

from proflex.pdb import PdbHandler
from proflex.rfdiffusion import RFdFix
from biopandas.pdb import PandasPdb 
import os
import sys

if __name__ == "__main__":
    
    w_dir = sys.argv[1]
    items = os.listdir(w_dir)
    reference = sys.argv[2]
    pdbs_to_correct = [os.path.join(w_dir,pdb) for pdb in items if pdb.startswith(os.path.basename(reference)[:-4]+'_diff')]
    # Fixing the .pdb format in the RFD PDB output thanks to the PDB input in RFD
    print(f"Correcting diffusion models from {reference}...")
    PdbHandler.backbone(reference) 
    PdbHandler.atom(reference[:len(reference)-4]+'_backbone.pdb')
    
    for mdl in pdbs_to_correct:
        print(f"Correcting {mdl}...")
        RFdFix.correct(mdl, reference[:-4]+'_backbone.pdb')
        
        print("Tidying...")
        PdbHandler.tidy(mdl[:len(mdl)-4]+'_chainsfixed.pdb')
        
        pdb_rfd_name = mdl[:len(mdl)-4]
        pdb_rfd_file_chainsfixed = pdb_rfd_name+'_chainsfixed.pdb'
        pdb_rfd_file_chainsfixed_tidied = pdb_rfd_file_chainsfixed[:-4]+'_tidied.pdb'
        
        print("Deleting and renaming temporal files...")
        os.rename(pdb_rfd_file_chainsfixed_tidied, pdb_rfd_name+'_rfdfixed.pdb')
        os.remove(pdb_rfd_file_chainsfixed)
        
    print("Deleting reference files... \n")
    pdb_name = reference[:-4]
    pdb_file_backbone = pdb_name+'_backbone.pdb'
    pdb_file_backbone_atom = pdb_file_backbone[:-4]+'_atom.pdb'
    os.remove(pdb_file_backbone)
    os.remove(pdb_file_backbone_atom)