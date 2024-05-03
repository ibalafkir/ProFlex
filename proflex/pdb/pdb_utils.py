"""
General PDB manipulating functionalities
"""
from biopandas.pdb import PandasPdb
import copy
import pandas as pd
from pdbtools import pdb_fixinsert, pdb_tidy, pdb_selatom, pdb_sort, pdb_reatom, pdb_delchain, pdb_delelem, pdb_reres

class PdbDf:
    """
    PDB manipulation through dataframes
    """
    
    def atoms(pdb):
        """
        Extracts ATOM information from a PDB file and returns a pandas DataFrame
        """
        ppdb = PandasPdb().read_pdb(pdb)
        return ppdb.df['ATOM']
    
    def chains_id(atom_df):
        result=[]
        array = pd.unique(atom_df["chain_id"])
        for i in array:
            result.append(i)
        return result

    def ca(atom_df):
        """
        Extracts CA rows from a pandas DataFrame with ATOM information
        """
        return atom_df[atom_df['atom_name'] == 'CA']

    def chain(atom_df, chain_id):
        """
        Returns pandas DataFrame with the ATOM information of an input chain
        """
        result = atom_df[atom_df['chain_id'] == chain_id]
        if len(result) == 0:
            raise ValueError('Error: the input chain name is absent')
        return result

    def rel_col(chain_atom_df):
        """
        Extracts "relevant columns" from an ATOM PDB pandas DataFrame. 
        ProFlex considers relevant: chain_id, residue_number, residue_name, x_coord, y_coord, z_coord
        """
        desired_cols = ['chain_id', 'residue_number', 'residue_name', 'x_coord', 'y_coord', 'z_coord']
        x = chain_atom_df[desired_cols]
        return x

    def coords(chain):
        """
        Gets ONLY coordinates (x, y, z) from a PDB ATOM DataFrame
        """
        desired_cols = ['x_coord', 'y_coord', 'z_coord']
        x = chain[desired_cols]
        return x

    def to_ca(new_coords, atom_df):
        """
        Reconstructs an ATOM PDB DataFrame with the new coordinates
        """
        chain = copy.deepcopy(new_coords) # so that the original chain is not affected
        chain.iloc[:, 3] = atom_df.iloc[:, 0]
        chain.iloc[:, 4] = atom_df.iloc[:, 1]
        chain.iloc[:, 5] = atom_df.iloc[:, 2]
        return chain

    def update_coords(new_coords, atom_df):
        """
        Reconstructs an ATOM PDB DataFrame with the new coordinates
        """
        atom_df.iloc[:, 11] = new_coords.iloc[:, 0]
        atom_df.iloc[:, 12] = new_coords.iloc[:, 1]
        atom_df.iloc[:, 13] = new_coords.iloc[:, 2]
        return atom_df


class PdbHandler:
    """
    General PDB processing functions that generate new PDBs by applying these functions
    """
        
    def fix_insertions(pdb):
        """
        Does a renumbering deleting antibody insertion codes using pdb-tools
        TER lines tend to encounter errors (easily solved afterwards)
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
        """
        with open(pdb, "r") as f_input, open(pdb[:-4]+'_terfixed.pdb', "w") as f_output:
            for line in f_input:
                if line.startswith("TER"):
                    f_output.write("TER\n")
                    atom_index = line.find("ATOM")
                    if atom_index != -1:
                        f_output.write(line[atom_index:])  
                else:
                    f_output.write(line)

    def remove_hetatm(pdb):
        """
        Removes HETATM lines and their ANISOU lines (for now, other removals might be added here)
        """
        with open(pdb, 'r') as f_input, open(pdb[:-4]+'_nohetatm.pdb', "w") as f_output:
            c = 0 # Counter
            for line in f_input: # Start reading
                if line.startswith('HETATM'): 
                    c += 1 # If bumping into a HETATM, counter stops being 0 -> this serves to know we've entered the HETATM area
                    continue # We don't write HETATM lines
                if c != 0 and line.startswith('ANISOU'): # If we are in the HETATM area (c = 0) and we bump into ANISOU lines
                    continue # We don't write them
                f_output.write(line)
    
    def custom_remove(pdb, lst):
        """
        Removes all lines starting by all the elements indicated in the list (e.g. ['ATOM', 'HETATM']) in the output PDB
        """
        with open(pdb, 'r') as f_input, open(pdb[:-4]+'_removed.pdb', "w") as f_output:
            for line in f_input:
                if not any(line.startswith(element) for element in lst):
                    f_output.write(line)
    
    def chain(pdb, chains_lst):
        """
        Keeps the atom names from the desired chain names
        """
        name = pdb[:-4]+'_ch.pdb'
        atom_df = PdbDf.atoms(pdb)
        atom_df_ch = atom_df[atom_df['chain_id'].isin(chains_lst)]
        atom_df_ch['segment_id']=''
        new_pdb = PandasPdb().read_pdb(pdb)
        new_pdb.df['ATOM'] = atom_df_ch
        new_pdb.to_pdb(path=name, records=['ATOM'], gz = False)
    
    def renatom(pdb, atom_num=1):
        """
        Renumbers a PDB starting from the desired atom number
        """
        pdb_atomsorted = pdb[:len(pdb)-4]+'_atomsorted.pdb'
        f = open(pdb, 'rt')
        f_atomsorted = open(pdb_atomsorted, 'wt')
        lines = f.readlines()        
        for modified_line in pdb_reatom.run(lines, atom_num):
            f_atomsorted.write(modified_line)
        f.close()
        f_atomsorted.close()
    
    def renres(pdb, resnum):
        """
        Renumbers a PDB starting from the desired residue number
        """
        f = open(pdb, 'rt')
        f_resnum = open(pdb[:-4]+'_renum.pdb', 'wt')
        lines = f.readlines()        
        for modified_line in pdb_reres.run(lines, resnum):
            f_resnum.write(modified_line)
        f.close()
        f_resnum.close()

    def delel(pdb, lst_el, output_name):
        """
        Deletes all atomic elements indicated in the list
        #TODO quitar lo de output
        """
        f = open(pdb, 'rt')
        f_deleted = open(output_name, 'wt')
        lines = f.readlines()        
        for modified_line in pdb_delelem.run(lines, lst_el):
            f_deleted.write(modified_line)
        f.close()
        f_deleted.close()

    def delotherel(pdb, lst_el, output_name):
        # TODO QUITAR LO DE OUTPUT
        """
        Deletes atom types indicated in the list apart from
        the conventional ones (N, C, O, H; manageable by pdb_delel)
        Also some H are not manageable by pdb_delel
        Used when the atom type is not manageable by pdb_delel
        """
        pdb_atom_df = PdbDf.atom(pdb)
        lst_elh = [term for term in pdb_atom_df['atom_name'] if term.startswith('H')]
        lst_elh.append('OXT')
        pdb_atom_df = pdb_atom_df[~pdb_atom_df['atom_name'].isin(lst_elh)]
        pdb_file = PandasPdb().read_pdb(pdb)
        pdb_file.df['ATOM'] = pdb_atom_df
        pdb_file.to_pdb(path= output_name, records = ['ATOM'], gz=False) 

    def backbone(pdb):
        """
        Deletes all atoms but the ones belonging to the backbone
        Does not touch other lines (headers...)
        """
        pdb_backbone = pdb[:len(pdb)-4]+'_backbone.pdb'
        f = open(pdb, 'rt')
        f_backbone = open(pdb_backbone, 'wt')
        lines = f.readlines()        
        for modified_line in pdb_selatom.run(lines, ['CA','C','N','O']):
            f_backbone.write(modified_line)
        f.close()
        f_backbone.close()
        
    def atom(pdb):
        """
        Keeps only atom lines
        """
        pdb_atom = pdb[:len(pdb)-4]+'_atom.pdb'
        ppdb = PandasPdb().read_pdb(pdb)
        ppdb.to_pdb(path= pdb_atom, records = ['ATOM'], gz=False)

    def tidy(pdb):
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

    def sort(pdb):
        """ 
        Sorts a PDB according to chain ID and residue number
        """
        pdb_sorted = pdb[:len(pdb)-4]+'_sorted.pdb'
        f = open(pdb, 'rt')
        f_sorted = open(pdb_sorted, 'wt')
        lines = f.readlines()        
        for modified_line in pdb_sort.run(lines, 'CR'): # 'CR' orders first by chain and then by residues
            f_sorted.write(modified_line)
        f.close()
        f_sorted.close()