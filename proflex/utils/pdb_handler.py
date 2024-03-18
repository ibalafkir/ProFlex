""""
Preprocesses a PDB file before entering the pipeline, mainly by removing
insertion codes that antibodies usually have
"""

from pdbtools import pdb_fixinsert
import os
import pandas as pd
from biopandas.pdb import PandasPdb
from pdbtools import pdb_tidy, pdb_selatom, pdb_sort, pdb_reatom, pdb_delchain, pdb_delelem
from proflex.utils import PDBUtils
import mdtraj as md


class PDBProcessor:
        
    def fix_insertions(pdb):
        """
        Does a renumbering deleting antibody insertion codes using pdb-tools
        TER lines tend to encounter errors (easily solved afterwards)
        :param Path to PDB
        :return Solved PDB file for insertion codes
        """
        f = open(pdb, "rt")
        f_fixed = open(pdb[:-4]+'_insfixed.pdb', "wt")
        lines = f.readlines()
        #f.seek(0)
        for modified_line in pdb_fixinsert.run(lines, []):
            f_fixed.write(modified_line)
        f.close()
        f_fixed.close()


    def fix_ter_mistakes(pdb):
        """
        Solves TER lines mistakes from the pdb-tools pdb_fixinsert tool
        :param Path to processed PDB
        :return Solved PDB file for TER lines
        """
        with open(pdb, "r") as f_input, open(pdb[:-4]+'_terfixed.pdb', "w") as f_output:
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
      
    def remove_hetatm(pdb):
        """
        Removes HETATM lines and their ANISOU lines (for now, other removals might be added here)
        :param Path to PDB
        :return PDB with previous lines removed
        """
        with open(pdb, 'r') as f_input, open(pdb[:-4]+'_nohetatm.pdb', "w") as f_output:
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
    
    def remove_lines(pdb, lst):
        with open(pdb, 'r') as f_input, open(pdb[:-4]+'_removed.pdb', "w") as f_output:
            for line in f_input:
                if not any(line.startswith(element) for element in lst):
                    f_output.write(line)
    
    def pdb_keepchain(pdb, keep):
        name = pdb[:-4]+'_ch.pdb'
        atom_df = PDBUtils.get_pdb_atoms_df(pdb)
        atom_df_ch = atom_df[atom_df['chain_id'].isin(keep)]
        new_pdb = PandasPdb().read_pdb(pdb)
        new_pdb.df['ATOM'] = atom_df_ch
        new_pdb.to_pdb(path=name, records=['ATOM'], gz = False)
    
    def pdb_atomrenumber(pdb, atom_number=1):
        """
        Renumbers a PDB starting from the desired atom_number
        """
        pdb_atomsorted = pdb[:len(pdb)-4]+'_atomsorted.pdb'
        f = open(pdb, 'rt')
        f_atomsorted = open(pdb_atomsorted, 'wt')
        lines = f.readlines()        
        for modified_line in pdb_reatom.run(lines, atom_number):
            f_atomsorted.write(modified_line)
        f.close()
        f_atomsorted.close()


class RFDFixer:
    
    def pdb_delel(pdb, lst_el, output_name):
        """
        Deletes all atomic elements indicated in the list
        """
        f = open(pdb, 'rt')
        f_deleted = open(output_name, 'wt')
        lines = f.readlines()        
        for modified_line in pdb_delelem.run(lines, lst_el):
            f_deleted.write(modified_line)
        f.close()
        f_deleted.close()
        
    
    def pdb_delotherel(pdb, lst_el, output_name):
        """
        Deletes atom types indicated in the list apart from
        the conventional ones (N, C, O, H; manageable by pdb_delel)
        Used when the atom type is not manageable by pdb_delel
        """
        
        pdb_atom_df = PDBUtils.get_pdb_atoms_df(pdb)
        
        for i in lst_el:
            pdb_atom_df = pdb_atom_df[pdb_atom_df['atom_name'] != i]
        
        pdb_file = PandasPdb().read_pdb(pdb)
        pdb_file.df['ATOM'] = pdb_atom_df
        pdb_file.to_pdb(path= output_name, records = ['ATOM'], gz=False) 
    
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
        Assigns chain ID and residues ID (using the input PDB in RFD) 
        to the output of RFD
        Both must have the same atom lines
        
        """
        pdb_rfd_df = PDBUtils.get_pdb_atoms_df(pdb_rfd)
        pdb_before_rfd_backbone_atom_df = PDBUtils.get_pdb_atoms_df(pdb_before_rfd_backbone_atom)
        
        if len(pdb_rfd_df) == len(pdb_before_rfd_backbone_atom_df):
            print("Number of backbone atoms coincide thus the correction can be done\n")
        else:
            print("Number of backbone atoms do not coincide, check if any backbone atom is repeated and delete lines manually") 
            # As certain choices like b-factors needs to be managed by the user, we cannot guess
            exit(1)

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
    
    
    def superpose_v2(pdb1, pdb2, output):
            """ 
            pdb1 contains reference coordinates
            pdb2 contains the moving coordinates in the superimposition
            PDB1 AND PDB2, THROUGH THIS WAY, MUST HAVE THE SAME NUMBER OF ATOMS!!! OTHERWISE IT DOESN'T WORK
            # TODO Superpose only backbone/CA of rigid residues (not the designed ones by RFdiffusion)
            """
            pdb1_traj = md.load(pdb1)
            pdb2_traj = md.load(pdb2)

            # Superponer las estructuras minimizando el RMSD del pdb2 sobre el pdb1
            aligned_traj = pdb2_traj.superpose(pdb1_traj)

            # Guardar la estructura superpuesta
            aligned_traj.save(output)
    
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
    
    


"""
 def pdb_keepchainbad(pdb,keep): 
    
        pdb_del = pdb[:-4]+'_ch.pdb'
        atom_df = PDBUtils.get_pdb_atoms_df(pdb)
        chains_id = PDBUtils.get_chains_id(atom_df)
        delet = [] 
        for i in chains_id:
            if i not in keep:
                delet.append(i)  
        f = open(pdb, 'rt')
        f_del = open(pdb_del, 'wt')
        lines = f.readlines()
        for modified_line in pdb_delchain.run(lines, delet):
            f_del.write(modified_line)   
        f.close()
        f_del.close()
        return

"""