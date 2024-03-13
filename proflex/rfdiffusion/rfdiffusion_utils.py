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
        
        # Deletes isolate residue numbers like 5 in [1,2,3], [5], [8, 9, 10] 
        
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

    def get_contigs(chain_1, chain_2, intchain1additional, intchain2additional):
        """
        #TODO Continuously revise possible mistakes in the code generation when diffusing first
        and last residues
        """
        id1 = PDBUtils.get_chains_id(chain_1)[0]
        id2 = PDBUtils.get_chains_id(chain_2)[0]
        chain1_start = RFDContigs.get_resnum_list(chain_1)[0]
        chain1_end = RFDContigs.get_resnum_list(chain_1)[len(chain_1)-1]
        chain2_start = RFDContigs.get_resnum_list(chain_2)[0]
        chain2_end = RFDContigs.get_resnum_list(chain_2)[len(chain_2)-1]
        intera1_resnum = RFDContigs.get_resnum_list(intchain1additional)
        intera2_resnum = RFDContigs.get_resnum_list(intchain2additional)
        intera1_lst = RFDContigs.get_chunks(intera1_resnum)
        intera2_lst = RFDContigs.get_chunks(intera2_resnum)
        intera1_lst = RFDContigs.chunk_filter(intera1_lst)
        intera2_lst = RFDContigs.chunk_filter(intera2_lst)

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
        """ Deletes side chains of the given chain ID and residue numbers (contained in a list)
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
        """ Gets side chains of the given chain ID and residue numbers (contained in a list)
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