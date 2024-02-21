""""
Preprocesses a PDB file before entering the pipeline, mainly by removing
insertion codes that antibodies usually have
"""

from pdbtools import pdb_fixinsert
import os

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
                    CONECT_index = line.find("CONECT")
                    if CONECT_index != -1:
                        f_output.write(line[CONECT_index:])
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