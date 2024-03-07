""""
Preprocesses a PDB file before entering the pipeline, mainly by removing
insertion codes that antibodies usually have
"""

from pdbtools import pdb_fixinsert
import os
import pandas as pd
from biopandas.pdb import PandasPdb
from pdbtools import pdb_tidy, pdb_selatom, pdb_sort, pdb_reatom
from proflex.utils import PDBUtils
import mdtraj as md


class PDBProcessor:
    
    def extract_without_extension(n):
        """
        Extracts path of PDB without the extension
        :param Path to PDB file
            str
        :return Path to PDB file without extension
            str
        """
        point_positions = [pos for pos, char in enumerate(n) if char == "."]
        return n[0:point_positions[-1]]
     
        
    def fix_insertions(n):
        """
        Does a renumbering deleting antibody insertion codes using pdb-tools
        TER lines tend to encounter errors (easily solved afterwards)
        :param Path to PDB
        :return Solved PDB file for insertion codes
        """
        f = open(n, "rt")
        f_fixed = open(PDBProcessor.extract_without_extension(n)+'_insfixed.pdb', "wt")
        lines = f.readlines()
        #f.seek(0)
        for modified_line in pdb_fixinsert.run(lines, []):
            f_fixed.write(modified_line)
        f.close()
        f_fixed.close()


    def fix_ter_mistakes(pdb_insfixed):
        """
        Solves TER lines mistakes from the pdb-tools pdb_fixinsert tool
        :param Path to processed PDB
        :return Solved PDB file for TER lines
        """
        with open(pdb_insfixed, "r") as f_input, open(PDBProcessor.extract_without_extension(pdb_insfixed)+'_terfixed.pdb', "w") as f_output:
            for line in f_input:
                if line.startswith("TER"):
                    f_output.write("TER\n")
                    atom_index = line.find("ATOM")
                    if atom_index != -1:
                        f_output.write(line[atom_index:])  
                    
                    # Conect lines are not included in the fixed PDB because the bonds between atoms
                    # are infered by the protein residue numbers
                    # CONECT_index = line.find("CONECT")
                    # if CONECT_index != -1:
                    #     f_output.write(line[CONECT_index:])
                    
                else:
                    f_output.write(line)
      
        
    def remove_other_lines(pdb_insfixed_terfixed):
        """
        Removes HETATM lines and their ANISOU lines (for now, other removals might be added here)
        :param Path to PDB
        :return PDB with previous lines removed
        """
        with open(pdb_insfixed_terfixed, 'r') as f_input, open(PDBProcessor.extract_without_extension(pdb_insfixed_terfixed)+'_removed.pdb', "w") as f_output:
            c = 0 # Counter
            for line in f_input: # Start reading
                if line.startswith('HETATM'): 
                    c += 1 # If bumping into a HETATM, counter stops being 0 -> this serves to know we've entered
                           # the HETATM area
                    continue # We don't write HETATM lines
                if c != 0 and line.startswith('ANISOU'): # If we are in the HETATM area (c = 0) and we bump into
                                                         # ANISOU lines
                    continue # We don't write them
                f_output.write(line)
    

    def pdb_atomrenumber(pdb, atom_number):
        """
        Renumbers a PDB from the desired atom_number
        """
        pdb_atomsorted = pdb[:len(pdb)-4]+'_atomsorted.pdb'
        f = open(pdb, 'rt')
        f_atomsorted = open(pdb_atomsorted, 'wt')
        lines = f.readlines()        
        for modified_line in pdb_reatom.run(lines, atom_number):
            f_atomsorted.write(modified_line)
        f.close()
        f_atomsorted.close()


        
    def del_mid_files(pdb):
        """
        Deletes all the files that previous functions but the last generate and renames the last
        to <pdb>_processed.pdb
        :param Path to PDB

        """
        pdb_name = PDBProcessor.extract_without_extension(pdb)
        pdb_name_insfixed = pdb_name + '_insfixed.pdb'
        pdb_name_terfixed = PDBProcessor.extract_without_extension(pdb_name_insfixed)+'_terfixed.pdb'
        pdb_name_removed = PDBProcessor.extract_without_extension(pdb_name_terfixed)+'_removed.pdb'
        os.rename(pdb_name_removed, pdb_name+'_processed.pdb')
        os.remove(pdb_name_insfixed)
        os.remove(pdb_name_terfixed)

class RFDFixer:
    
    def get_n_atoms(pdb):
        """
        Gets number of atoms assuming the last ATOM line has the highest atom number and that
        all atoms are followed: e.g. 1-2-3... and not 1-3-5... with jumps
        Input: PDB file, a RFD output pdb fulfills previous conditions
        
        All atom renumber problems can be easily solved with pdbtools
        
        """
        pdb_df = PDBUtils.get_pdb_atoms_df(pdb)
        n_atoms = pdb_df['atom_number'].iloc[-1]
        return n_atoms
    
    def pdb_backbone(pdb):
        """
        Deletes all atoms but the ones belonging to the backbone
        Does not touch other lines (headers...)
        In the fixing of RFD output PDBs approach, serves with pdb_atom to modify
        the preRFD-PDB to make possible the comparison with postRFD-PDB
        """
        pdb_backbone = pdb[:len(pdb)-4]+'_backbone.pdb'
        f = open(pdb, 'rt')
        f_backbone = open(pdb_backbone, 'wt')
        lines = f.readlines()        
        for modified_line in pdb_selatom.run(lines, ['CA','C','N','O']):
            f_backbone.write(modified_line)
        f.close()
        f_backbone.close()


    def pdb_atom(pdb):
        """
        Keeps only atom lines
        """
        pdb_atom = pdb[:len(pdb)-4]+'_atom.pdb'
        ppdb = PandasPdb().read_pdb(pdb)
        ppdb.to_pdb(path= pdb_atom, records = ['ATOM'], gz=False)


    def correct_rfd_pdbs(pdb_rfd, pdb_before_rfd_backbone_atom):
        """
        Assigns chain ID and residues ID from a the input PDB to RFD in
        the output PDB of RFD
        """
        pdb_rfd_df = PDBUtils.get_pdb_atoms_df(pdb_rfd)
        pdb_before_rfd_backbone_atom_df = PDBUtils.get_pdb_atoms_df(pdb_before_rfd_backbone_atom)
        pdb_rfd_df['chain_id'] = pdb_before_rfd_backbone_atom_df['chain_id']
        pdb_rfd_df['residue_number'] = pdb_before_rfd_backbone_atom_df['residue_number']
        
        pdb_rfd_fixed_name = pdb_rfd[:len(pdb_rfd)-4]+ "_chainsfixed.pdb"
        pdb_rfd_chainsfixed = PandasPdb().read_pdb(pdb_rfd)
        pdb_rfd_chainsfixed.df['ATOM'] = pdb_rfd_df
        pdb_rfd_chainsfixed.to_pdb(path= pdb_rfd_fixed_name, records = ['ATOM'], gz=False)   


    def pdb_tidying(pdb):
        """
        Detects chain ID changes to assign TER and END lines
        """
        pdb_tidied = pdb[:len(pdb)-4]+'_tidied.pdb'
        f = open(pdb, 'rt')
        f_tidied = open(pdb_tidied, 'wt')
        lines = f.readlines()        
        for modified_line in pdb_tidy.run(lines, False):
            f_tidied.write(modified_line)
        f.close()
        f_tidied.close()

 
    def del_mid_files(pdb_before_rfd, pdb_rfd):
        """
        Deletes temporal files that previous functions create to
        correct the PDB
        """
        pdb_name = pdb_before_rfd[:len(pdb_before_rfd)-4]
        pdb_file_backbone = pdb_name+'_backbone.pdb'
        pdb_file_backbone_atom = pdb_file_backbone[:len(pdb_file_backbone)-4]+'_atom.pdb'

        pdb_rfd_name = pdb_rfd[:len(pdb_rfd)-4]
        pdb_rfd_file_chainsfixed = pdb_rfd_name+'_chainsfixed.pdb'
        pdb_rfd_file_chainsfixed_tidied = pdb_rfd_file_chainsfixed[:len(pdb_rfd_file_chainsfixed)-4]+'_tidied.pdb'

        os.rename(pdb_rfd_file_chainsfixed_tidied, pdb_rfd_name+'_rfdfixed.pdb')
        os.remove(pdb_file_backbone)
        os.remove(pdb_file_backbone_atom)
        os.remove(pdb_rfd_file_chainsfixed)
    
    def superpose(pdb1, pdb2, output):
        """
        pdb1 contains reference coordinates
        pdb2 contains the moving coordinates in the superimposition
        # TODO Superpose only backbone/CA of rigid residues (not the designed ones by RFdiffusion)
        """
        traj1 = md.load_pdb(pdb1)
        traj2 = md.load_pdb(pdb2)
        backbone_atoms_1 = traj1.topology.select("backbone") # 
        backbone_atoms_2 = traj2.topology.select("backbone") # 

        traj1_backbone = traj1.atom_slice(backbone_atoms_1)
        traj2_backbone = traj2.atom_slice(backbone_atoms_2)

        superposed_traj2 = traj2_backbone.superpose(traj1_backbone)
        superposed_traj2.save(output)
    
    def pdb_sorting(pdb):
        """ Sorts a PDB according to chain ID and residue number
        """
        pdb_sorted = pdb[:len(pdb)-4]+'_sorted.pdb'
        f = open(pdb, 'rt')
        f_sorted = open(pdb_sorted, 'wt')
        lines = f.readlines()        
        for modified_line in pdb_sort.run(lines, 'CR'): # 'CR' orders first by chain and then by residues
            f_sorted.write(modified_line)
        f.close()
        f_sorted.close()
    
    
