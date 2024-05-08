"""
ProFlex interface analyzer of the interface of the interaction between 2 proteins in which:
- The antigen needs to be in a single chain
- The antibody can have 2 chains (ab1 and ab2) or 1 chain (ab1) in case of nanobodies
The distance threshold is the distance up to which 2 residues are considered to be interacting.
"""
import argparse
import pandas as pd
from proflex.pdb import PdbDf
from proflex.interface import InterfaceAnalyzer

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculates interface residues between two chains of PDB file')
    parser.add_argument('--pdb', type=str, help='Path to the PDB file')
    parser.add_argument('--ab1', type=str, help='First chain')
    parser.add_argument('--ab2', type=str, help='Second chain if existing')
    parser.add_argument('--ag', type=str, help='An chain')
    parser.add_argument('--distance_threshold', type=float, help='Distance threshold for interface', default=8)
    args = parser.parse_args()
    pdb = args.pdb
    ab1 = args.ab1
    ab2 = args.ab2
    ag = args.ag
    
    distance_threshold = args.distance_threshold
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    

    # For NANOBODIES
    
    if ab1 and not ab2:
    
        int1, int2, detected_interactions, int1_extended, int2_extended = InterfaceAnalyzer.run(pdb, args.ab1, args.ag, distance_threshold)
        print("----------------------------------------------------------------------\n")
        print("[INFO] Running ProFlex InterfaceAnalyzer for a nanobody-antigen system\n")
        print("----------------------------------------------------------------------\n")
        print(f"[INFO] 1) CHAINS {ab1} AND {ag}\n")    
        
        #print("[INFO]Interface residues from chain", ab1)
        #print(int1, "\n")
        #print("[INFO] Interface residues from chain", ag)
        #print(int2, "\n")
        # print("[INFO] Interacting pairs")
        # print(detected_interactions, "\n")
        #print("[INFO] Additional interacting residues (neighborhood = 2 residues) in chain", ab1)
        #print(intchain1_extended, "\n")
        #print("[INFO] Additional interacting residues (neighborhood = 2 residues) in chain", ag)
        # print(intchain2_extended, "\n")
        print("[INFO] Code for PELE com_distance metric \n")
        InterfaceAnalyzer.pele_com_distance(int1, int2, ab1, ag)

    # For ANTIBODIES

    if ab1 and ab2:
        print("----------------------------------------------------------------------\n")
        print("[INFO] ProFlex InterfaceAnalyzer for an antibody-antigen system\n")
        print("----------------------------------------------------------------------\n") 

        int1_1, int2_1, detected_interactions, int1_1_extended, int2_1_extended = InterfaceAnalyzer.run(pdb, args.ab1, args.ag, distance_threshold)      
        print(f"[INFO] 1) CHAINS {ab1} AND {ag}")
        #print("[INFO]Interface residues from chain", ab1)
        #print(int1_1, "\n")
        #print("[INFO] Interface residues from chain", ag)
        #print(int2_1, "\n")
        #print("[INFO] Interacting pairs")
        #print(detected_interactions, "\n")
        #print("[INFO] Additional interacting residues (neighborhood = 2 residues) in chain", ab1)
        #print(intchain1_extended, "\n")
        #print("[INFO] Additional interacting residues (neighborhood = 2 residues) in chain", ag)
        #print(intchain2_extended, "\n")
        print("[INFO] Code for PELE com_distance metric")
        InterfaceAnalyzer.pele_com_distance(int1_1, int2_1, ab1, ag)

        print(f"\n[INFO] 2) CHAINS {ab2} AND {ag}")  
        int1_2, int2_2, detected_interactions, int1_2_extended, int2_2_extended = InterfaceAnalyzer.run(pdb, args.ab2, args.ag, distance_threshold)
        #print("\n[INFO] Interface residues from chain", ab2)
        #print(int1_2, "\n")
        #print("[INFO] Interface residues from chain", ag)
        #print(int2_2, "\n")
        #print("[INFO] Interacting pairs")
        #print(detected_interactions, "\n")
        #print("[INFO] Additional interacting residues (neighborhood = 2 residues) in chain", ab2)
        #print(intchain1_extended, "\n")
        #print("[INFO] Additional interacting residues (neighborhood = 2 residues) in chain", ag)
        #print(intchain2_extended, "\n")
        print("[INFO] Code for PELE com_distance metric")
        InterfaceAnalyzer.pele_com_distance(int1_2, int2_2, ab2, ag)
        
        
        
        
                            ####### Residues to omit ######## TODO Create a PELETask builder in a class and include this
        
        print(f"\n[INFO] SELECTIONS TO OMIT (residues with CA distance < 8 Angstrom)")
        def del_repeated(lst):
            """
            Deletes repeated elements in a list
            """
            result = []
            for i in lst:
                if i not in result:
                    result.append(i)
                else:
                    continue
            return result
        
        # Analyze antibodies interface
        ab1, ab2, d, int1_ex, int2_ex = InterfaceAnalyzer.run(pdb, args.ab1, args.ab2, distance_threshold=8)
        
        # Selected residues with the antigen
        ab1Antigen_lst = int1_1['residue_number'].to_list()
        ab2Antigen_lst = int1_2['residue_number'].to_list()     
        ab1Antigen_lst = sorted(del_repeated(ab1Antigen_lst))
        ab2Antigen_lst = sorted(del_repeated(ab2Antigen_lst))
        
        # Making the selection to omit tag
        
        ab1_lst = ab1['residue_number'].to_list()
        ab2_lst = ab2['residue_number'].to_list()
        ab1_lst = sorted(del_repeated(ab1_lst))
        ab2_lst = sorted(del_repeated(ab2_lst))
        #print("Esto es ab1_lst", ab1_lst)
        #print("Esto es ab2_lst", ab2_lst)
        
        ab1_lst_updated = []
        for element in ab1_lst:
            if element in ab1Antigen_lst:
                continue
            else:
                ab1_lst_updated.append(element)
        ab2_lst_updated = []
        for element in ab2_lst:
            if element in ab2Antigen_lst:
                continue
            else:
                ab2_lst_updated.append(element)
        
        ab1_lst_updated = sorted(del_repeated(ab1_lst_updated))
        ab2_lst_updated = sorted(del_repeated(ab2_lst_updated))
        
        ab1_lst_updated_chname = [f"{args.ab1}:" + str(element) for element in ab1_lst_updated]
        ab1_lst_updated_chname_mod = ['"' + element + '"' for element in ab1_lst_updated_chname]
        pele_code_1 = ', '.join(ab1_lst_updated_chname_mod)


        ab2_lst_updated_chname = [f"{args.ab2}:" + str(element) for element in ab2_lst_updated]
        ab2_lst_updated_chname_mod = ['"' + element + '"' for element in ab2_lst_updated_chname]
        pele_code_2 = ', '.join(ab2_lst_updated_chname_mod)
        
        print('"selectionToOmit": { "links": { "ids":[' + str(pele_code_1) + ', ' + str(pele_code_2) + ']}},')