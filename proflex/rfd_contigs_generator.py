"""
Generates the contig codes that indicate the regions to diffuse in RFDiffusion
Input from commandline: PDB with both proteins interacting and the chains to analyze interface

Example:

python rfd_contigs_generator.py --pdb B_PD1_diff_8_prep_relax_fixed.pdb --id1 A --id2 B --distance_threshold 200
"""

from proflex.utils import PDBUtils
from proflex.utils.pdb_interface_analyzer import InterfaceAnalyzer
import argparse

def get_resnum_list(df):
    result = []
    for i in df['residue_number']:
        result.append(int(i))
    return result

def get_chunks(lst):
    result = []
    subgroup = []
    for num in lst:
        if not subgroup or num == subgroup[-1] + 1:
            subgroup.append(num)
        else:
            result.append(subgroup)
            subgroup = [num]
    if subgroup:
        result.append(subgroup)
    
    # Deletes isolate residue numbers like 5 in [1,2,3], [5], [8, 9, 10] 
    
    for i in result: 
        if len(i)==1: 
            result.remove(i) 
    return result

def chunk_filter(lst):
    """
    :param lst: list
    :return list without elements of length < 3 so as not to diffuse regions of less than 3 aas: list
    """
    result = []
    for i in lst:
        if len(i)>2:
            result.append(i)
    return result

def get_contigs(chain_1, chain_2, intchain1additional, intchain2additional):
    """
    id1 and id2 are indicated from terminal.
    :param chain_1:
    :param chain_2:
    :param intchain1additional:
    :param intchain2additional:
    :return:
    """
    chain1_start = get_resnum_list(chain_1)[0]
    chain1_end = get_resnum_list(chain_1)[len(chain_1)-1]
    chain2_start = get_resnum_list(chain_2)[0]
    chain2_end = get_resnum_list(chain_2)[len(chain_2)-1]
    intera1_resnum = get_resnum_list(intchain1additional)
    intera2_resnum = get_resnum_list(intchain2additional)
    intera1_lst = get_chunks(intera1_resnum)
    intera2_lst = get_chunks(intera2_resnum)
    intera1_lst = chunk_filter(intera1_lst)
    intera2_lst = chunk_filter(intera2_lst)

    # First chain

    if chain1_start in intera1_resnum and chain1_end in intera1_resnum:
        contig = f"{chain1_end-chain1_start+1}-{chain1_end-chain1_start+1}"
      
    elif chain1_start in intera1_resnum and chain1_end not in intera1_resnum:
        
        # Checks if only the last residue is considered rigid
        c = 0
        if len(intera1_resnum) == chain1_end-chain1_start: # only chain1_end residue is excluded
            c += 1
            contig = f"{chain1_end-chain1_start+1}-{chain1_end-chain1_start+1}/{id1}{chain1_end}"
            
        if c == 0:
            contig = f"{intera1_lst[0][-1]-intera1_lst[0][0]+1}-{intera1_lst[0][-1]-intera1_lst[0][0]+1}/{id1}"
            contig = contig + f"{intera1_lst[0][-1]+1}-"
            intera1_lst = intera1_lst[1:]
            for curr_list in intera1_lst:
                list_start = curr_list[0]
                list_end = curr_list[len(curr_list)-1]
                contig += f"{list_start-1}/{list_end-list_start+1}-{list_end-list_start+1}/{id1}{list_end+1}-"
            contig = contig[:len(contig)-1] + f"-{chain1_end}"
    
    elif chain1_start not in intera1_lst and chain1_start+1 in intera1_resnum:
        
        contig = f"{id1}"
        for curr_list in intera1_lst:
            list_start = curr_list[0]
            list_end = curr_list[len(curr_list)-1]
            contig += f"{list_start-1}/{list_end-list_start+1}-{list_end-list_start+1}/{id1}{list_end+1}-"
        contig = contig[:len(contig)-1] + f"-{chain1_end}"
      
    else: 
        contig = f"{id1}{chain1_start}-"
        for curr_list in intera1_lst:
            list_start = curr_list[0]
            list_end = curr_list[len(curr_list)-1]
            contig += f"{list_start-1}/{list_end-list_start+1}-{list_end-list_start+1}/{id1}{list_end+1}-"
        contig = contig[:len(contig)-1] + f"-{chain1_end}"
    
    ## Change of chain ######
    
    contig += "/0 "
    
    #########################
    
    # Second chain
    
    if chain2_start in intera2_resnum and chain2_end in intera2_resnum:
        contig += f"{chain2_end-chain2_start+1}-{chain2_end-chain2_start+1}"

    elif chain2_start in intera2_resnum and chain2_end not in intera2_resnum:
        c = 0
        if len(intera2_resnum) == chain2_end-chain2_start: # only chain1_end residue is excluded
            c += 1
            contig += f"{chain2_end-chain2_start+1}-{chain2_end-chain2_start+1}"
                
        if c == 0:
            contig += f"{intera2_lst[0][-1]-intera2_lst[0][0]+1}-{intera2_lst[0][-1]-intera2_lst[0][0]+1}/{id2}"
            contig = contig + f"{intera2_lst[0][-1]+1}-"
            intera2_lst = intera2_lst[1:]
            for curr_list in intera2_lst:
                list_start = curr_list[0]
                list_end = curr_list[len(curr_list)-1]
                contig += f"{list_start-1}/{list_end-list_start+1}-{list_end-list_start+1}/{id2}{list_end+1}-"
            contig = contig[:len(contig)-1] + f"-{chain1_end}"
                
    else:
        contig += f"{id2}{chain2_start}-"
        for curr_list in intera2_lst:
            list_start = curr_list[0]
            list_end = curr_list[len(curr_list)-1]
            contig += f"{list_start-1}/{list_end-list_start+1}-{list_end-list_start+1}/{id2}{list_end+1}-"
        contig += f"{chain2_end}"
    
    contig = '['+contig+']'    
    return contig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Returns RFdiffusion contig code from a PDB file'
        )
    parser.add_argument('--pdb', type=str, help='Path to the PDB file')
    parser.add_argument('--id1', type=str, help='First chain')
    parser.add_argument('--id2', type=str, help='Second chain')
    parser.add_argument('--distance_threshold', type=float,
                        help='Distance threshold for interface '
                             'residues (6 Angstrom is default)', default=6)
    args = parser.parse_args()
    pdb_file = args.pdb
    id1 = args.id1
    id2 = args.id2
    distance_threshold = args.distance_threshold

    atom_df = PDBUtils.get_pdb_atoms_df(pdb_file)
    atom_df_ca = PDBUtils.get_ca(atom_df)
    chain_1 = PDBUtils.get_chain(atom_df_ca, id1)
    chain_2 = PDBUtils.get_chain(atom_df_ca, id2)
    chain1_relevant = PDBUtils.get_relevant_columns(chain_1)
    chain2_relevant = PDBUtils.get_relevant_columns(chain_2)
    int1, int2, detected_interactions = InterfaceAnalyzer.get_interface_residues_by_chain(
        chain1_relevant, chain2_relevant, distance_threshold) # detected_interactions is unused
    intchain1additional = InterfaceAnalyzer.amplify_selection_residues(int1, chain1_relevant)
    intchain2additional = InterfaceAnalyzer.amplify_selection_residues(int2, chain2_relevant)
    print(get_contigs(chain_1, chain_2, intchain1additional, intchain2additional))