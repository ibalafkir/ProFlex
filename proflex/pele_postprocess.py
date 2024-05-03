from proflex import *
import argparse
"""
Processes PELE outputs so that they are well read by RFdiffusion 
(changes alternative res_names to common ones for now)
#TODO In development for adding other res_name alternatives
"""

def run(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                residue_name = line[17:20]
                if residue_name.strip() == "CYT" or residue_name.strip() == "CYX":
                    line = line[:17] + "CYS" + line[20:]
                if residue_name.strip() == "HIE" or residue_name.strip() == "HIP" or residue_name.strip() == "HID":
                    line = line[:17] + "HIS" + line[20:]      
            outfile.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PELE poses processor')
    parser.add_argument('pdb', type=str, help='Path to the PDB file')
    args = parser.parse_args()
    pdb = args.pdb
    run(pdb, pdb[:-4]+'_processed.pdb')