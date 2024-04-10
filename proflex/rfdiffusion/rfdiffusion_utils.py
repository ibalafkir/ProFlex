from proflex.utils import PDBUtils

class RFDContigs:
    """
    Generation of contigs in RFdiffusion  
    """
    def get_resnum_list(df):
        """
        Extracts list of residue numbers from the 
        residue_number column in an ATOM pandas dataframe        
        """
        result = []
        for i in df['residue_number']:
            result.append(int(i))
        return result

    def get_chunks(lst):
        """
        Makes a list of lists with chunks of followed residues
        BUG when it gets two same followed numbers eg [..., 220, 220,...] bc it starts a new chunk, example in 5HGG_b,
            this happens when the first 220 corresponds to a conformatin (eg A) and the second to another (eg B)
            No problem if the list it gets has no repeated res numbers
        TODO the removal of lists with only one element is not ideal here and should be in chunk_filter
        """
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
        
        # Deletes isolate residue numbers like 5 in [1,2,3], [5], [8, 9, 10] / THIS FITS BETTER IN chunk_filter (it works here though)
        for i in result: 
            if len(i)==1: 
                result.remove(i) 
        return result

    def chunk_filter(lst):
        """
        Returns a list without elements of length < 3 
        so as not to diffuse regions of less than 3 aas
        """
        result = []
        for i in lst:
            if len(i)>2:
                result.append(i)
        return result

    def get_contigs(chain, intchainadditional, include_extremes=True):
        """
        Generates contigs
        Inputs. PDB df of a SINGLE chain and PDB df with selected as interactive residues with the other chain
        Tip. The last df should be amplified so that the diffused regions altogether are not too segmented (that's why the parameter is 'additional')
        """
        
        # Loading data and defining variables
        id = PDBUtils.get_chains_id(chain)[0]
        chain_start = RFDContigs.get_resnum_list(chain)[0]
        chain_end = RFDContigs.get_resnum_list(chain)[-1]
        chain_allresnum = RFDContigs.get_resnum_list(chain)
        intera_resnum =  RFDContigs.get_resnum_list(intchainadditional)
        intera_lst = RFDContigs.get_chunks(intera_resnum)
        intera_lst = RFDContigs.chunk_filter(intera_lst)

        # RFdiffusion cannot deal with rigid extremes in a chain if they are just next to a flexible zone
        # this nuance is written here because it has to do with how the program works, but it might be
        # more logical to transfer and adapt it to amplify in interface_analyzer
        if include_extremes == True:
            if intera_lst:
                if chain_start not in intera_lst[0] and chain_start+1 in intera_lst[0]:
                    intera_lst[0].insert(0, chain_start)
                if chain_end not in intera_lst[-1] and chain_end-1 in intera_lst[-1]:
                    intera_lst[-1].append(chain_end)
        
        # Flattens the list
        intera_lst_extended = []
        for sublst in intera_lst:
            intera_lst_extended.extend(sublst)
        
        # Generates contigs
        if intera_lst:
            print(intera_lst)       
            if len(intera_lst_extended) == len(chain_allresnum):
                contig = f"{len(chain_allresnum)}-{len(chain_allresnum)}/"
                
            elif chain_start in intera_lst[0] and chain_end not in intera_lst[-1]:
                contig = f"{len(intera_lst[0])}-{len(intera_lst[0])}/"
                contig += f"{id}{intera_lst[0][-1]+1}-"
                for curr_list in intera_lst[1:]:
                    curr_list_start = curr_list[0]
                    curr_list_end = curr_list[-1]
                    contig += f"{curr_list_start-1}/{len(curr_list)}-{len(curr_list)}/"
                    contig += f"{id}{curr_list_end+1}-"
                contig += f"{chain_end}"+"/"
                
            elif chain_start in intera_lst[0] and chain_end in intera_lst[-1]:
                contig = f"{len(intera_lst[0])}-{len(intera_lst[0])}/"
                contig += f"{id}{intera_lst[0][-1]+1}-"
                for curr_list in intera_lst[1:]:
                    curr_list_start = curr_list[0]
                    curr_list_end = curr_list[-1]
                    contig += f"{curr_list_start-1}/{len(curr_list)}-{len(curr_list)}/"
                    if curr_list == intera_lst[-1]:
                        pass
                    else:
                        contig += f"{id}{curr_list_end+1}-"
            
            elif chain_start not in intera_lst[0] and chain_end in intera_lst[-1]:
                contig = f"{id}{chain_start}-"
                for curr_list in intera_lst:
                    curr_list_start = curr_list[0]
                    curr_list_end = curr_list[-1]
                    contig += f"{curr_list_start-1}/{len(curr_list)}-{len(curr_list)}/"
                    if curr_list == intera_lst[-1]:
                        pass
                    else:
                        contig += f"{id}{curr_list_end+1}-"
                
            
            else:
                contig = f"{id}{chain_start}-"
                for curr_list in intera_lst:
                    curr_list_start = curr_list[0]
                    curr_list_end = curr_list[-1]
                    contig += f"{curr_list_start-1}/{len(curr_list)}-{len(curr_list)}/"
                    contig += f"{id}{curr_list_end+1}-"
                contig += f"{chain_end}"+"/"


        else:
            contig = f"{id}{chain_start}-{chain_end}"+"/"

        return contig
class RFDSchains:
    
    def get_contig_chid(contiglist):
        """
        Gets a list of chains ID contained in a contig code.
        E.g. [A1-A30/3-3/A34-37/0 B1-B3/2-2/B6-B10] results in the lis [A, B]
        """
        result = []
        for i in contiglist:
            if i.isalpha() and not i in result:
                result.append(i)
        return result

    def get_rigid_ranges(lst):
        """
        Gets rigid ranges (contigs) from a list of chunks of the contig
        In ['A1-110', '5-5', 'A116-124'], where the 1st and 3rd elements are contigs and
        the second is a flexible zone, the function would keep the 1st and 3rd elements
        """
        result = []
        for item in lst:
            if any(char.isalpha() for char in item):
                result.append(item)
        return result

    def expand_ranges(lst):
        """
        Creates a list of (residue numbers) from rigid ranges in contigs code
        [A1-3, A5-9] results in [1, 2, 3, 5, 6, 7, 8, 9]
        """
        result = []
        for item in lst:
            parts = item.split('-')  
            start = int(parts[-2][1:])
            end = int(parts[-1]) 
            result.extend(range(start, end + 1))
        return result

    def del_spaces(lst):
        """
        Deletes spaces in elements of a list if any, often happening in previous results of functions    
        """
        result = []
        for i in lst:
            result.append(i.replace(" ", ""))
        return result
            
    def get_rigid_residues(contig):
        """
        Combines functions to get 2 lists, each with the chain ID (first element) and
        residue numbers of rigid residues (rest of elements)
        """

        chain_id_1 = RFDSchains.get_contig_chid(contig)[0]
        first_chain = contig.split('/0')[0]
        first_chain = first_chain[1:]
        first_chain_chunks = first_chain.split('/')
        first_chain_rigid_ranges = RFDSchains.get_rigid_ranges(first_chain_chunks)
        first_chain_rigid_ranges = RFDSchains.del_spaces(first_chain_rigid_ranges)
        first_chain_rigid_resnumbers = RFDSchains.expand_ranges(first_chain_rigid_ranges)

        chain_id_2 = RFDSchains.get_contig_chid(contig)[1]
        second_chain = contig.split('/0')[1]
        second_chain = second_chain[:len(second_chain)-1]
        second_chain_chunks = second_chain.split('/')
        second_chain_rigid_ranges = RFDSchains.get_rigid_ranges(second_chain_chunks)
        second_chain_rigid_ranges = RFDSchains.del_spaces(second_chain_rigid_ranges)
        second_chain_rigid_resnumbers = RFDSchains.expand_ranges(second_chain_rigid_ranges)

        return list(chain_id_1) + first_chain_rigid_resnumbers, list(chain_id_2) + second_chain_rigid_resnumbers

    def delete_sidechains(df, half):
        """ 
        Deletes side chains of the given chain ID and residue numbers (contained in a list)
        in an ATOM pandas dataframe of a PDB
        """
        atom_name_keep = ['N','CA', 'C', 'O']
            
        # Boolean mask strategy to know which rows will be deleted
        # i.e. the ones that are not included in the mask
            
        mask_delete_ch = (df['chain_id'].isin(half)) & \
                        (df['residue_number'].isin(half)) & \
                        (~df['atom_name'].isin(atom_name_keep))
            
        df_filt = df[~mask_delete_ch]
            
        return df_filt

    def get_sidechains(df, half):
        """ 
        Gets side chains of the given chain ID and residue numbers (contained in a list)
        in an ATOM pandas dataframe of a PDB
        """
        atom_name_nokeep = ['N','CA', 'C', 'O']
            
        # Boolean mask strategy to know which rows will be deleted
        # i.e. the ones that are not included in the mask
        mask_get_ch = (df['chain_id'].isin(half)) & \
                        (df['residue_number'].isin(half)) & \
                        (~df['atom_name'].isin(atom_name_nokeep))
            
        df_filt = df[mask_get_ch]
            
        return df_filt