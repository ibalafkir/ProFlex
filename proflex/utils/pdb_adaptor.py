""""
Preprocesses a PDB file before entering the pipeline, mainly by removing
insertion codes that antibodies usually have
"""

from pdbtools import pdb_fixinsert

class PDBAdaptor:
    
    def extract_without_extension(n):
        """
        Extracts name of PDB without the extension
        :param Name of PDB file
            str
        :return Name of PDB file without extension
            str
        """
        point_positions = [pos for pos, char in enumerate(n) if char == "."]
        return n[0:point_positions[-1]]
        
    def fix_insertions(n):
        """
        Does a renumbering deleting antibody insertion codes using pdb-tools
        TER lines tend to encounter errors
        :param Path to PDB
        :return Solved PDB file
        """
        f = open(n, "rt")
        f_fixed = open(PDBAdaptor.extract_without_extension(n)+'_fixed.pdb', "wt")
        lines = f.readlines()
        #f.seek(0)
        for modified_line in pdb_fixinsert.run(lines, []):
            f_fixed.write(modified_line)
        f.close()
        f_fixed.close()

    def fix_ter_mistakes(input_pdb):
        """
        Solves TER lines mistakes
        :param Path to fixed PDB
        :return Solved PDB file for TER lines
        """
        temp_file = input_pdb + '.tmp'
        with open(input_pdb, "r") as f_input, open(temp_file, "w") as f_output:
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
        
        # Replace original file with the modified one
        import os
        os.replace(temp_file, input_pdb)