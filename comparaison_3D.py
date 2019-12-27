'''
Created on 2 mai 2019

@author: coline

Calcul de la RMSD sur notre jeu de donnees
'''

from networkx.algorithms import isomorphism
import Bio.PDB
from Bio.PDB import MMCIFParser
import os
import math
import copy
import csv
from scipy.integrate import quad     
import numpy as np
from recup_data.constantes import PATH_MMCIF, EXTENSION_PATH_TAILLE,\
    GROUPES_TOUTES_ARETES_MAX_4_10_07, EXTENSION_PATH,\
    CLUSTERING_PEREZ_VERSION_NON_CAN_2, GROUPE_GNRA_ETENDU, GROUPE_GNRA,\
    GROUPE_ARICH, GROUPE_ARICH_ETENDU, NEW_EXTENSION_PATH_TAILLE
import pickle
import urllib
import networkx as nx
from networkx.classes.function import create_empty_copy
import matplotlib.pyplot as plt
import re
import seaborn as sns

''' Recupere les positions des extensions dans le graphe de la molecule '''
def recup_motif_et_chaines(graphe, taille_ext, taille_mol):
    
    chaines = []
    motif = []
    nb_nts_par_chaine = [taille_ext]*5
    for i in range(1,6) :
        chaines.append(graphe.nodes[i]["position"][0])
        motif.append(graphe.nodes[i]["position"][0])
        
#         print(graphe.nodes[i]["position"][0])
#         print(chaines)
        if i < 5 :
            for j in range(1,taille_ext) :
                if i == 1 or i == 4 :
                    if graphe.nodes[i]["position"][0]+j <= taille_mol and graphe.nodes[i]["position"][0]+j not in chaines :
                        chaines.append(graphe.nodes[i]["position"][0]+j) 
                else :
                    if graphe.nodes[i]["position"][0]-j > 0 and graphe.nodes[i]["position"][0]-j not in chaines :
                        chaines.append(graphe.nodes[i]["position"][0]-j)
        

    if min(graphe.nodes[1]["position"][0],graphe.nodes[2]["position"][0]) + 2*(taille_ext-1) >=  max(graphe.nodes[1]["position"][0],graphe.nodes[2]["position"][0]) :
        print("ramous")
    
    if min(graphe.nodes[3]["position"][0],graphe.nodes[4]["position"][0]) + 2*(taille_ext-1) >=  max(graphe.nodes[3]["position"][0],graphe.nodes[4]["position"][0]) :
        print("ramous")
    
    return motif, chaines, nb_nts_par_chaine

''' Recupere les positions des extensions dans le graphe de la molecule 

modif pour le nouveau format de graphe'''
def recup_motif_et_chaines_new_data(graphe, taille_ext, taille_mol, graphes):
    
    chaines = []
    motif = []
    nb_nts_par_chaine = [taille_ext]*5
    for i in range(1,6) :
        chaines.append((graphe.nodes[i]["num_ch"], graphe.nodes[i]["position"][0]))
        motif.append((graphe.nodes[i]["num_ch"], graphe.nodes[i]["position"][0]))
        
#         print(graphe.nodes[i]["position"][0])
#         print(chaines)
        if i < 5 :
            for j in range(1,taille_ext) :
                if i == 1 or i == 4 :
                    if graphe.nodes[i]["position"][0]+j <= taille_mol and graphe.nodes[i]["position"][0]+j not in chaines :
                        chaines.append((graphe.nodes[i]["num_ch"], graphe.nodes[i]["position"][0]+j)) 
                else :
                    if graphe.nodes[i]["position"][0]-j > 0 and graphe.nodes[i]["position"][0]-j not in chaines :
                        chaines.append((graphe.nodes[i]["num_ch"], graphe.nodes[i]["position"][0]-j))
#         if i == 5 :
#             for j in range(1, taille_ext+1) :
#                 print((graphe.nodes[i]["num_ch"], graphe.nodes[i]["position"][0]+j))
#                 if graphe.nodes[i]["position"][0]-j <= taille_mol and (graphe.nodes[i]["num_ch"], graphe.nodes[i]["position"][0]-j) not in chaines :
#                         chaines.append((graphe.nodes[i]["num_ch"],graphe.nodes[i]["position"][0]-j))
#                         print("petit rat")

    if min(graphe.nodes[1]["position"][0],graphe.nodes[2]["position"][0]) + 2*(taille_ext-1) >=  max(graphe.nodes[1]["position"][0],graphe.nodes[2]["position"][0]) :
        print("ramous")
    
    if min(graphe.nodes[3]["position"][0],graphe.nodes[4]["position"][0]) + 2*(taille_ext-1) >=  max(graphe.nodes[3]["position"][0],graphe.nodes[4]["position"][0]) :
        print("ramous")
    
    return motif, chaines, nb_nts_par_chaine

'''issu du code de Vladimir
deux aretes sont considerees comme identiques si elles ont le meme label et le meme type de distance
'''
def edge_match(d_g1, d_g2):
    if (d_g1['label'] == d_g2['label'] and
        d_g1['long_range'] == d_g2['long_range']):
        return True
    return False


'''issu du code de Vladimir
renvoie le mapping des deux graphes en entree
cad la correspondance entre les sommets des deux graphes quand on les superpose selon les aretes de meme label et de meme type de distance
les deux graphes doivent etre isomorphes
'''
def graphs_mapping(g1, g2):#
    GM = isomorphism.DiGraphMatcher(g1, g2, edge_match=edge_match)
    if not GM.is_isomorphic():
        print("Try to map non-isomorphic graphs!!!")
        return -1
    mapping = GM.mapping
    return mapping


''' modifie a partir du code de Vladimir
creation d'un fichier PDB ne comprenant que les positions des atomes des nucleotides
de tab_positions a partir du fichier PDB complet 
'''
def extract_pdb(tab_positions, name, graphes): #name = numPDB_chaine, tab_positions = ensemble des positions a considerer, graphes = graphes.2.92 
    rna = name.split("_")[0] ## num pdb
    chain = name.split("_")[1] ## num chaine
    
    p = MMCIFParser()
    file = PATH_MMCIF + rna + ".cif"
    print(PATH_MMCIF)
    if rna+".cif" not in os.listdir(PATH_MMCIF) : ## si le fichier n est pas encore stocke en interne, on va le chercher sur la PDB et on le stocke ou on veut
        print("petit rat")
        url = 'http://www.rcsb.org/pdb/files/%s.cif' % rna
        urllib.request.urlretrieve(url, file)
    print(rna)
    try:
        mmcifd = p.get_structure(rna, file)
    except:
        print("Error in RNA %s" % file)
        exit()
 
    mmcifd = mmcifd[0][chain]
    nodes = [graphes[(rna, chain)].node[x]['fr3d'] for x in graphes[(rna,chain)].nodes() if x in tab_positions]
    real_nodes = []
    print(nodes)
    p = re.compile('[-]*[0-9]*')
    for elt in nodes :
        m = p.match(elt)
        if m :
            num = m.group()
        else :
            print("probleme5")
        real_nodes.append(int(num))
        
    nb = 0
    l_atoms = []
    print(mmcifd.get_residues())
    for res in mmcifd.get_residues():
        print(res.id)
        if res.id[1] not in real_nodes:
            continue
        for atom in res.get_atom():
            nb += 1
            #number name resname chainid(only first char) seqnb coord occupancy bfactor element
            coords = list(atom.get_coord())
            string = 'ATOM  {:>5d}  {:<4}{:>3} {}{:>4}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.1f}{:>6.2f}         {:>2}  '.format(
                        nb, atom.name, res.resname, chain[0], res.id[1], coords[0], coords[1], coords[2],
                        atom.occupancy, atom.bfactor, atom.element)
            l_atoms.append(string)
    return '\n'.join(l_atoms)

''' modifie a partir du code de Vladimir
creation d'un fichier PDB ne comprenant que les positions des atomes des nucleotides
de tab_positions a partir du fichier PDB complet 

modif a cause du nouveau format de fichiers de graphes
'''
def extract_pdb_new_data(tab_positions, rna, chains, graphe): #name = numPDB_chaine, tab_positions = ensemble des positions a considerer, graphes = graphes.2.92 
    ## rna = num_pdb
    #chain = name.split("_")[1] ## num chaine
    
    print(chains)
    
    p = MMCIFParser()
    file = PATH_MMCIF + rna + ".cif"
    print(PATH_MMCIF)
    if rna+".cif" not in os.listdir(PATH_MMCIF) : ## si le fichier n est pas encore stocke en interne, on va le chercher sur la PDB et on le stocke ou on veut
        print("petit rat")
        url = 'http://www.rcsb.org/pdb/files/%s.cif' % rna
        urllib.request.urlretrieve(url, file)
    print(rna)
    try:
        mmcifd = p.get_structure(rna, file)
    except:
        print("Error in RNA %s" % file)
        exit()
    
    nb = 0
    l_atoms = []
    compteur = 1
    for chain in chains :
        #mmcifds.append(mmcifd[0][chain])
        mmcif_ch = mmcifd[0][chain]
        nodes = [graphe.node[x]['fr3d'] for x in graphe.nodes() if x in tab_positions and x[0] == chain]
        #test_nodes = [x for x in graphe.nodes() if x in tab_positions and x[0] == chain]
        print("ramousnif")
        print(nodes)
        real_nodes = []
        print(nodes)
        p = re.compile('[-]*[0-9]*')
        for elt in nodes :
            m = p.match(elt)
            if m :
                num = m.group()
                if num == elt :
                    real_nodes.append((' ',int(num),' '))
                else :
                    
                    c = elt[0]
                    i = 0
                    while i < len(num) and c == num[i] :
                        c = elt[i+1]
                        i = i+1
                    if i != len(num) :
                        print("bizarre")
                        exit(0)
                    num_en_plus = elt[i:]   
                    
                    real_nodes.append((' ', int(num), num_en_plus))
            else :
                print("probleme5")
            #real_nodes.append(int(num))
            
        
        print(mmcifd.get_residues())
        print(real_nodes)
        
        
        #for mmcifd in mmcifds :
        for res in mmcif_ch.get_residues():
            print(res.id)
            if res.id not in real_nodes:
                continue
            for atom in res.get_atom():
                nb += 1
                #number name resname chainid(only first char) seqnb coord occupancy bfactor element
                coords = list(atom.get_coord())
                string = 'ATOM  {:>5d}  {:<4}{:>3} {}{:>4}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.1f}{:>6.2f}         {:>2}  '.format(
                            nb, atom.name, res.resname, compteur, res.id[1], coords[0], coords[1], coords[2],
                            atom.occupancy, atom.bfactor, atom.element)
                l_atoms.append(string)
        compteur += 1
    return '\n'.join(l_atoms)

'''issu de github/charnley/rmsd'''
def calcul_rmsd(V, W): ## 
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    Parameters
    ----------
    V : array
        (N,D) matrix, where N is points and D is dimension.
    W : array
        (N,D) matrix, where N is points and D is dimension.
    Returns
    -------
    rmsd : float
        Root-mean-square deviation between the two vectors
    """
    if len(V) == 0 or len(W) == 0 :
        return -1;
    
    D = len(V[0])
    N = len(V)
    result = 0.0
    for v, w in zip(V, W):
        result += sum([(v[i] - w[i])**2.0 for i in range(D)])
    res = np.sqrt(result/N)   
    if res == None :
        print("chabou")
        print(result)
        print(N)
        exit()
    return res

''' modifie a partir du code de Vladimir
A partir des graphes des molecules et de la position des extensions,
recupere les graphes composes du motif et des 4 chaines (sans les liens non covalents, sauf pour le motif)
superpose 2à2 les motifs en minimisant la RMSD (version nt represente par son C1')
calcule la RMSD associee a la comparaison des motifs + 4 chaines 2à2 selon la rotation/translation qu'on a trouvee a l'etape precedente (version get_rms_2)
 '''
def rmsd(tab_fichiers, tab_fichier_positions, tab_motif_positions, tab_nb_nts_par_chaine, graphes):
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
  
    def get_rms(s1, s2, graph1, graph2, subgraph1, subgraph2):
        print(s1)
        print(s2)
        s1_rep = PATH_MMCIF + s1
        s2_rep = PATH_MMCIF + s2
        print(s1)
        print(s2)
        
        ref_structure = pdb_parser.get_structure(s1_rep, s1_rep)
        sample_structure = pdb_parser.get_structure(s2_rep, s2_rep)
          
        ref_model = ref_structure[0]
        sample_model = sample_structure[0]
  
        #ref_atoms = [at for chain in ref_model for res in chain for at in res]
        #sample_atoms = [at for chain in sample_model for res in chain for at in res]
  
        ref_res = [res for chain in ref_model for res in chain]
        sample_res = [res for chain in sample_model for res in chain]
        m = graphs_mapping(graph1, graph2)
        
        for chain in ref_model :
            num_chaine_ref = chain.id
            
        for chain in sample_model : 
            num_chaine_sample = chain.id

#         print(m)
        #print s1, s2
        #print i, j
        #print cluster[CLASS-1]['names'][i], cluster[CLASS-1]['names'][j]
        
        print(len(ref_res))
        print(ref_res)
        print(len(sample_res))
        print(sample_res)
        if len(sample_res) != len(ref_res) :
            print("probleme3")
            return
        
        for x in sample_res :
            print(x.id)
        m = {graphes[(s1.split("_")[1], s1.split("_")[2])].node[x]['fr3d']:graphes[(s2.split("_")[1], s2.split("_")[2])].node[m[x]]['fr3d']
             for x in m}
        print(m)
        ref_atoms = [] 
        sample_atoms = []
        
        new_sample_res = Bio.PDB.Chain.Chain(num_chaine_sample)
        new_ref_res = Bio.PDB.Chain.Chain(num_chaine_ref)
        print(new_sample_res)
        print(new_ref_res)
        
        for res1 in ref_res:
            pos1 =  str(res1.id[1])
            pos2 = m[pos1]
#             print(pos2)
            
            res2 = [x for x in sample_res if str(x.id[1]) == pos2][0]
            
            print("petit rat")
            print(res1)
            
            names_1 = [x.name for x in res1]
            names_2 = [x.name for x in res2]
            print(names_1)
            
            ref_atoms.extend([atom for atom in res1 if atom.name in names_2 ]) #and atom.name=="C1'"])
            sample_atoms.extend([atom for atom in res2 if atom.name in names_1]) #and atom.name=="C1'"])
            
            new_ref_res.add(res1)
            print(res2.id)
            res2.id = (' ', int(pos1), ' ')
            print(res2.id)
            new_sample_res.add(res2) 
            
            
        print(m)
        print(new_ref_res)
        print(new_sample_res)
        print("gros rat")
        print(ref_atoms)
        print(sample_atoms)
        for res in new_ref_res :
            print("ramou")
            print(res)
#         print(len(ref_atoms))
#         print(len(sample_atoms))
#         print(ref_atoms)
#         print(sample_atoms)
#         if (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_1U9S_A_58_11_2.pdb" and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_5J7L_DA_48_15_2.pdb") or (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_5J7L_DA_48_15_2.pdb" and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_1U9S_A_58_11_2.pdb") :
#             exit()
#         if (s2 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_1U9S_A_58_11_2.pdb' and s1 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_5J7L_DA_48_15_2.pdb') or (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_1U9S_A_58_11_2.pdb' and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_5J7L_DA_48_15_2.pdb'):
#             print("gros rat /n /n /n /n")
#             exit()
        super_imposer = Bio.PDB.Superimposer()
        super_imposer.set_atoms(sample_atoms, ref_atoms)
        super_imposer.apply(ref_model.get_atoms())
        print(super_imposer.rms)
        print(len(ref_atoms))
        print(len(sample_atoms))
        
        new_sample_model = Bio.PDB.Model.Model(0)
        new_ref_model = Bio.PDB.Model.Model(0)
        new_sample_model.add(new_sample_res)
        new_ref_model.add(new_ref_res)
        
        new_structure_sample = Bio.PDB.Structure.Structure(0)
        new_structure_ref = Bio.PDB.Structure.Structure(0)
        new_structure_ref.add(new_ref_model)
        new_structure_sample.add(new_sample_model)
        
        io_sample = Bio.PDB.PDBIO()
        io_sample.set_structure(new_structure_sample) 
        io_sample.save(s2_rep[:len(s2_rep)-4]+"_aligned_with_"+s1[:len(s1)-4]+".pdb")
        
        io_ref = Bio.PDB.PDBIO()
        io_ref.set_structure(new_structure_ref) 
        io_ref.save(s1_rep[:len(s1_rep)-4]+"_aligned_with_"+s2[:len(s2)-4]+".pdb")
        
        return super_imposer.rms 
  
    def get_rms_2(s1, s2, graph1, graph2, subgraph1, subgraph2):
        with open(PATH_MMCIF+"fichiers_problemes.txt", 'a') as fichier :
            print(s1)
            print(s2)
            s1_rep = PATH_MMCIF + s1
            s2_rep = PATH_MMCIF + s2
            print(s1)
            print(s2)
            
            ref_structure = pdb_parser.get_structure(s1_rep, s1_rep)
            sample_structure = pdb_parser.get_structure(s2_rep, s2_rep)
              
            ref_model = ref_structure[0]
            sample_model = sample_structure[0]
      
            #ref_atoms = [at for chain in ref_model for res in chain for at in res]
            #sample_atoms = [at for chain in sample_model for res in chain for at in res]
      
            ref_res = [res for chain in ref_model for res in chain]
            sample_res = [res for chain in sample_model for res in chain]
            
            
            ##On commence par chercher a superposer les motifs seulement en minimisant la RMSD
            m = graphs_mapping(subgraph1, subgraph2)
            
            for chain in ref_model :
                num_chaine_ref = chain.id
                
            for chain in sample_model : 
                num_chaine_sample = chain.id
    
    #         print(m)
            #print s1, s2
            #print i, j
            #print cluster[CLASS-1]['names'][i], cluster[CLASS-1]['names'][j]
            
            print(len(ref_res))
            print(ref_res)
            print(len(sample_res))
            print(sample_res)
            if len(sample_res) != len(ref_res) :
                print("probleme")
                return
            
            for x in sample_res :
                print(x.id)
            m = {graphes[(s1.split("_")[1], s1.split("_")[2])].node[x]['fr3d']:graphes[(s2.split("_")[1], s2.split("_")[2])].node[m[x]]['fr3d']
                 for x in m}
            print(m)
            ref_atoms = [] 
            sample_atoms = []
        
            
            for res1 in ref_res:
                pos1 =  str(res1.id[1])
                if pos1 in m.keys() :
                    pos2 = m[pos1]
        #             print(pos2)
                    
                    res2 = [x for x in sample_res if str(x.id[1]) == pos2][0]
                    
                    print("petit rat")
                    print(res1)
                    
                    names_1 = [x.name for x in res1]
                    names_2 = [x.name for x in res2]
                    print(names_1)
                    
                    ref_atoms.extend([atom for atom in res1 if atom.name in names_2 and atom.name=="C1'"])
                    sample_atoms.extend([atom for atom in res2 if atom.name in names_1 and atom.name=="C1'"])
                
                
            print(m)
    #         print(new_ref_res)
    #         print(new_sample_res)
            print("gros rat")
            print(ref_atoms)
            print(sample_atoms)
    #         print(len(ref_atoms))
    #         print(len(sample_atoms))
    #         print(ref_atoms)
    #         print(sample_atoms)
    #         if (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_1U9S_A_58_11_2.pdb" and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_5J7L_DA_48_15_2.pdb") or (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_5J7L_DA_48_15_2.pdb" and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_1U9S_A_58_11_2.pdb") :
    #             exit()
    #         if (s2 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_1U9S_A_58_11_2.pdb' and s1 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_5J7L_DA_48_15_2.pdb') or (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_1U9S_A_58_11_2.pdb' and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_5J7L_DA_48_15_2.pdb'):
    #             print("gros rat /n /n /n /n")
    #             exit()
            print("avant")
            for atom in ref_atoms :
                print(atom.get_vector())
    
            super_imposer = Bio.PDB.Superimposer()
            super_imposer.set_atoms(sample_atoms, ref_atoms)
            #super_imposer.apply(ref_model.get_atoms())
            print(super_imposer.rms)
            print(len(ref_atoms))
            print(len(sample_atoms))
            print(super_imposer.rotran)
            
            #On recupere la rotation qui correspond a cette meilleure RMSD
            rotation,translation = super_imposer.rotran
            
            
            ## On fait la superposition sur toute l'extension avec la rotation des atomes qu'on vient de recuperer et on calcule la RMSD de cette superposition-là
    
    #         ref_structure_new = pdb_parser.get_structure(s1_rep, s1_rep)
    #         sample_structure_new = pdb_parser.get_structure(s2_rep, s2_rep)
              
    #         ref_model_new = ref_structure_new[0]
    #         sample_model_new = sample_structure_new[0]
      
            #ref_atoms = [at for chain in ref_model for res in chain for at in res]
            #sample_atoms = [at for chain in sample_model for res in chain for at in res]
      
    #         ref_res_new = [res for chain in ref_model_new for res in chain]
    #         sample_res_new = [res for chain in sample_model_new for res in chain]
    #         
    
             
            m_pymol = graphs_mapping(graph1, graph2) ## superposition avec le motif pour la visualisation
            if m_pymol == -1 :
                print("probleme2 : non isomorphic graphs (peut etre branches qui se rejoignent)")
                return 
            a_enlever = []
            for elt in m_pymol.keys() :
                if elt in subgraph1.nodes() :
                    a_enlever.append(elt)
              
            m_calcul = copy.deepcopy(m_pymol)   ##superposition sans le motif pour le calcul de RMSD
            for elt in a_enlever :
                del(m_calcul[elt])
            print("tout petit rat")
            print(m_pymol)
            print(m_calcul)
             
              
    #         for chain in ref_model_new :
    #             num_chaine_ref = chain.id
    #               
    #         for chain in sample_model_new : 
    #             num_chaine_sample = chain.id
    #  
    # #         print(m)
    #         #print s1, s2
    #         #print i, j
    #         #print cluster[CLASS-1]['names'][i], cluster[CLASS-1]['names'][j]
    #          
    #         print(len(ref_res_new))
    #         print(ref_res_new)
    #         print(len(sample_res_new))
    #         print(sample_res_new)
    #         if len(sample_res_new) != len(ref_res_new) :calcule la RMSD  
    #             print("probleme")
    #             return
    #           
    #         for x in sample_res_new :
    #             print(x.id)
            m_pymol = {graphes[(s1.split("_")[1], s1.split("_")[2])].node[x]['fr3d']:graphes[(s2.split("_")[1], s2.split("_")[2])].node[m_pymol[x]]['fr3d']
                 for x in m_pymol}
            
            m_calcul = {graphes[(s1.split("_")[1], s1.split("_")[2])].node[x]['fr3d']:graphes[(s2.split("_")[1], s2.split("_")[2])].node[m_calcul[x]]['fr3d']
                 for x in m_calcul}
            
            ref_atoms_calcul = [] 
            sample_atoms_calcul = []
              
            sample_res_pymol = Bio.PDB.Chain.Chain(num_chaine_sample)
            ref_res_pymol = Bio.PDB.Chain.Chain(num_chaine_ref)
            print(sample_res_pymol)
            print(ref_res_pymol)
            
            print(m_pymol)
            print(m_calcul)
            
            compteur = 1
            for res1 in ref_res:
                pos1 =  str(res1.id[1])
                print(pos1)
                if pos1 in m_pymol.keys() :
                    pos2 = m_pymol[pos1]
                    print(pos2)
                    
                    for x in sample_res :
                        print(x.id)
                       
                    tab_res2 = [x for x in sample_res if str(x.id[1]) == pos2]
                    if len(tab_res2) > 0 :
                        res2 = tab_res2[0]
                      
                        print("petit rat")
                        print(res1)
                          
                        names_1 = [x.name for x in res1]
                        names_2 = [x.name for x in res2]
                        print(names_1)
                        
                        
                        
                        res1_pymol = Bio.PDB.Residue.Residue((' ',compteur , ' '), res1.get_resname(), res1.get_segid())
                        for atom in res1 :
                            if atom.name in names_2 :
                                print(atom.get_vector())
                                atom.transform(rotation, translation)
                                
                                res1_pymol.add(atom)
                        #res1.id = (' ',compteur , ' ')  
                        ref_res_pymol.add(res1_pymol)
                        print(res2.id)
                        #res2.id = (' ', compteur, ' ')
                        res2_pymol = Bio.PDB.Residue.Residue((' ',compteur , ' '), res2.get_resname(), res2.get_segid())
                        for atom in res2 :
                            if atom.name in names_1 :
                                print(atom.get_vector())
                                res2_pymol.add(atom)
                        print(res2.id)
                        sample_res_pymol.add(res2_pymol) 
                        
                        compteur += 1
                        
                        if pos1 in m_calcul.keys() :
                            if m_calcul[pos1] != pos2 : ## pas meme superposition entre visualisation et calcul => probleme
                                print("probleme")
                                print(m_pymol)
                                print(m_calcul)
                                print(m_pymol[pos1])
                                print(pos2)
                                exit()
                            else :
                                ref_atoms_calcul.extend([atom for atom in res1 if atom.name in names_2 and atom.name=="C1'"])
                                sample_atoms_calcul.extend([atom for atom in res2 if atom.name in names_1 and atom.name=="C1'"])
                    else :
                        fichier.write(str(s2_rep)+ " " + str(pos2) +"\n")
            
            liste_atomes_ref = []
            for atom in ref_atoms_calcul :
                #atom.transform(rotation, translation)   
                liste_atomes_ref.append(list(atom.get_vector()))
                print(list(atom.get_vector()))
            
            liste_atomes_sample = []    
            for atom in sample_atoms_calcul :
    #             atom.transform(rotation, translation)
                liste_atomes_sample.append(list(atom.get_vector()))
            
            print(liste_atomes_ref)
            
            print("avant")
            for res in ref_res_pymol :
                for atom in res :
                    print(atom.get_vector())
                    break
            
    #         new_new_ref_res = Bio.PDB.Chain.Chain(num_chaine_ref)
    #         for res in new_ref_res :
    #             #res2 = [x for x in new_sample_res if str(x.id[1]) == str(res.id[1])][0]
    # #             new_res = Bio.PDB.Residue.Residue(res.id[0], res.id[1], res.id[2])
    #             for atom in res :
    #                 #if atom.name in [x.name for x in res2] :
    #                     atom.transform(rotation, translation)  
    #                     new_res.add(atom)
    #             new_new_ref_res.add(new_res)       
                        
                        
            print("avpres")
            for res in ref_res_pymol :
                for atom in res :
                    print(atom.get_vector())
                    break
                    
            for res in sample_res_pymol :
                for atom in res :
                    print(atom.get_vector())
                    break
                        
    #         for res in new_sample_res :
    #             res2 = [x for x in new_ref_res if str(x.id[1]) == str(res.id[1])][0]
    #             for atom in res :
    #                 if atom.name in [x.name for x in res2] :
    #                     atom.transform(rotation, translation)  
                  
            print(m)
            print(ref_res_pymol)
            print(sample_res_pymol)
            print("gros rat")
            print(ref_atoms)
            print(sample_atoms)
            for res in ref_res_pymol :
                print("ramou")
                print(res)
                
                
            sample_model_pymol = Bio.PDB.Model.Model(0)
            ref_model_pymol = Bio.PDB.Model.Model(0)
            sample_model_pymol.add(sample_res_pymol)
            ref_model_pymol.add(ref_res_pymol)
               
            structure_sample_pymol = Bio.PDB.Structure.Structure(0)
            structure_ref_pymol = Bio.PDB.Structure.Structure(0)
            structure_ref_pymol.add(ref_model_pymol)
            structure_sample_pymol.add(sample_model_pymol)
    
            io_sample = Bio.PDB.PDBIO()
            io_sample.set_structure(structure_sample_pymol) 
            io_sample.save(s2_rep[:len(s2_rep)-4]+"_aligned_with_"+s1[:len(s1)-4]+".pdb")
             
            io_ref = Bio.PDB.PDBIO()
            io_ref.set_structure(structure_ref_pymol) 
            io_ref.save(s1_rep[:len(s1_rep)-4]+"_aligned_with_"+s2[:len(s2)-4]+".pdb")#   
    
    
    #         print(len(ref_atoms))
    #         print(len(sample_atoms))
    #         print(ref_atoms)
    #         print(sample_atoms)
    #         if (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_1U9S_A_58_11_2.pdb" and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_5J7L_DA_48_15_2.pdb") or (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_5J7L_DA_48_15_2.pdb" and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_1U9S_A_58_11_2.pdb") :
    #             exit()
    #         if (s2 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_1U9S_A_58_11_2.pdb' and s1 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_5J7L_DA_48_15_2.pdb') or (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_1U9S_A_58_11_2.pdb' and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_5J7L_DA_48_15_2.pdb'):
    #             print("gros rat /n /n /n /n")
    #             exit()
            
            
            return calcul_rmsd(liste_atomes_ref, liste_atomes_sample)
    
    rms = {}
    for i in range(len(tab_fichiers)) :
        print(tab_motif_positions[i])
        subgraph1_motif = nx.DiGraph(graphes[(tab_fichiers[i].split("_")[1], tab_fichiers[i].split("_")[2])].subgraph(tab_motif_positions[i]))
        graphe1_motif = create_empty_copy(subgraph1_motif)
        
        for noeud in graphe1_motif.nodes() :
                if noeud+1 in graphe1_motif.nodes() and (noeud,noeud+1) not in graphe1_motif.edges() and noeud != tab_motif_positions[i][4] and noeud+1 != tab_motif_positions[i][4]  :
                    graphe1_motif.add_edge(noeud, noeud+1, label='B53', long_range=False)
        
        graphe1_motif.add_edge(tab_motif_positions[i][0], tab_motif_positions[i][1], label='CSS', long_range=True)
        graphe1_motif.add_edge(tab_motif_positions[i][1], tab_motif_positions[i][0], label='CSS', long_range=True)
        
        graphe1_motif.add_edge(tab_motif_positions[i][2], tab_motif_positions[i][3], label='CSS', long_range=True)
        graphe1_motif.add_edge(tab_motif_positions[i][3], tab_motif_positions[i][2], label='CSS', long_range=True)
        
        graphe1_motif.add_edge(tab_motif_positions[i][0], tab_motif_positions[i][4], label='TSS', long_range=True)
        graphe1_motif.add_edge(tab_motif_positions[i][4], tab_motif_positions[i][0], label='TSS', long_range=True)
    
        graphe1_motif.add_edge(tab_motif_positions[i][1], tab_motif_positions[i][4], label='CWW', long_range=False)
        graphe1_motif.add_edge(tab_motif_positions[i][4], tab_motif_positions[i][1], label='CWW', long_range=False)
        
        subgraph = nx.DiGraph(graphes[(tab_fichiers[i].split("_")[1], tab_fichiers[i].split("_")[2])].subgraph(tab_fichier_positions[i]))
        graph1 = create_empty_copy(subgraph)
            
        for noeud in graph1.nodes() :
                if noeud+1 in graph1.nodes() and (noeud,noeud+1) not in graph1.edges() and noeud != tab_motif_positions[i][4] and noeud+1 != tab_motif_positions[i][4]  :
                    graph1.add_edge(noeud, noeud+1, label='B53', long_range=False)
        
        graph1.add_edge(tab_motif_positions[i][0], tab_motif_positions[i][1], label='CSS', long_range=True)
        graph1.add_edge(tab_motif_positions[i][1], tab_motif_positions[i][0], label='CSS', long_range=True)
        
        graph1.add_edge(tab_motif_positions[i][2], tab_motif_positions[i][3], label='CSS', long_range=True)
        graph1.add_edge(tab_motif_positions[i][3], tab_motif_positions[i][2], label='CSS', long_range=True)
        
        graph1.add_edge(tab_motif_positions[i][0], tab_motif_positions[i][4], label='TSS', long_range=True)
        graph1.add_edge(tab_motif_positions[i][4], tab_motif_positions[i][0], label='TSS', long_range=True)
    
        graph1.add_edge(tab_motif_positions[i][1], tab_motif_positions[i][4], label='CWW', long_range=False)
        graph1.add_edge(tab_motif_positions[i][4], tab_motif_positions[i][1], label='CWW', long_range=False)

                
        for j in range(i+1, len(tab_fichiers)) :
            print(tab_motif_positions[j])
            subgraph2_motif = nx.DiGraph(graphes[(tab_fichiers[j].split("_")[1], tab_fichiers[j].split("_")[2])].subgraph(tab_motif_positions[j]))
           
            graphe2_motif = create_empty_copy(subgraph2_motif)
            
            for noeud in graphe2_motif.nodes() :
                if noeud+1 in graphe2_motif.nodes() and (noeud,noeud+1) not in graphe2_motif.edges() and noeud != tab_motif_positions[j][4] and noeud+1 != tab_motif_positions[j][4] :
                    graphe2_motif.add_edge(noeud, noeud+1, label='B53', long_range=False)
        
            graphe2_motif.add_edge(tab_motif_positions[j][0], tab_motif_positions[j][1], label='CSS', long_range=True)
            graphe2_motif.add_edge(tab_motif_positions[j][1], tab_motif_positions[j][0], label='CSS', long_range=True)
            
            graphe2_motif.add_edge(tab_motif_positions[j][2], tab_motif_positions[j][3], label='CSS', long_range=True)
            graphe2_motif.add_edge(tab_motif_positions[j][3], tab_motif_positions[j][2], label='CSS', long_range=True)
            
            graphe2_motif.add_edge(tab_motif_positions[j][0], tab_motif_positions[j][4], label='TSS', long_range=True)
            graphe2_motif.add_edge(tab_motif_positions[j][4], tab_motif_positions[j][0], label='TSS', long_range=True)
        
            graphe2_motif.add_edge(tab_motif_positions[j][1], tab_motif_positions[j][4], label='CWW', long_range=False)
            graphe2_motif.add_edge(tab_motif_positions[j][4], tab_motif_positions[j][1], label='CWW', long_range=False)
            
            subgraph = nx.DiGraph(graphes[(tab_fichiers[j].split("_")[1], tab_fichiers[j].split("_")[2])].subgraph(tab_fichier_positions[j]))
            graph2 = create_empty_copy(subgraph)
            
            for noeud in graph2.nodes() :
                if noeud+1 in graph2.nodes() and (noeud,noeud+1) not in graph2.edges() and noeud != tab_motif_positions[j][4] and noeud+1 != tab_motif_positions[j][4]  :
                    graph2.add_edge(noeud, noeud+1, label='B53', long_range=False)
                    
            graph2.add_edge(tab_motif_positions[j][0], tab_motif_positions[j][1], label='CSS', long_range=True)
            graph2.add_edge(tab_motif_positions[j][1], tab_motif_positions[j][0], label='CSS', long_range=True)
            
            graph2.add_edge(tab_motif_positions[j][2], tab_motif_positions[j][3], label='CSS', long_range=True)
            graph2.add_edge(tab_motif_positions[j][3], tab_motif_positions[j][2], label='CSS', long_range=True)
            
            graph2.add_edge(tab_motif_positions[j][0], tab_motif_positions[j][4], label='TSS', long_range=True)
            graph2.add_edge(tab_motif_positions[j][4], tab_motif_positions[j][0], label='TSS', long_range=True)
        
            graph2.add_edge(tab_motif_positions[j][1], tab_motif_positions[j][4], label='CWW', long_range=False)
            graph2.add_edge(tab_motif_positions[j][4], tab_motif_positions[j][1], label='CWW', long_range=False)
                
#             print(tab_nb_nts_par_chaine[i])
#             print(graph1.number_of_edges())
#             print(graph1.edges.data())
#             print(tab_motif_positions[i])
#             print(tab_nb_nts_par_chaine[j])
#             print(graph2.number_of_edges())
#             print(graph2.edges.data())
#             print(tab_motif_positions[j])
            #print(graphes[(tab_fichiers[j].split("_")[1], tab_fichiers[j].split("_")[2])].nodes.data())
            
            if graph1.number_of_edges() != graph2.number_of_edges() :
                print(tab_fichiers[i])
                print(tab_fichiers[j])
                print(graph1.number_of_edges())
                print(graph1.edges.data())
                print(graph2.number_of_edges())
                print("probleme1")
                
            
            
            print(graphe1_motif.nodes.data())
            print(graphe2_motif.nodes.data())
            
            rms.update({(tab_fichiers[i], tab_fichiers[j]) : get_rms_2(tab_fichiers[i], tab_fichiers[j], graph1, graph2, graphe1_motif, graphe2_motif) })
            
#             if (tab_fichiers[i], tab_fichiers[j]) == ('fichier_3JCS_1_282_1_taille_4.pdb', 'fichier_4V9F_0_25_56_taille_4.pdb') :
#                 print(rms[(tab_fichiers[i], tab_fichiers[j])])
#                 exit()
    
    return rms


''' modifie a partir du code de Vladimir
A partir des graphes des molecules et de la position des extensions,
recupere les graphes composes du motif et des 4 chaines (sans les liens non covalents, sauf pour le motif)
superpose 2à2 les motifs en minimisant la RMSD (version nt represente par son C1')
calcule la RMSD associee a la comparaison des motifs + 4 chaines 2à2 selon la rotation/translation qu'on a trouvee a l'etape precedente (version get_rms_2)
 modif pour nouveau format de graphe
 '''
def rmsd_new_data(tab_fichiers, tab_fichier_positions, tab_motif_positions, tab_nb_nts_par_chaine, graphes, taille):
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
  
    def get_rms(s1, s2, graph1, graph2, subgraph1, subgraph2):
        print(s1)
        print(s2)
        s1_rep = PATH_MMCIF + s1
        s2_rep = PATH_MMCIF  + s2
        print(s1)
        print(s2)
        
        ref_structure = pdb_parser.get_structure(s1_rep, s1_rep)
        sample_structure = pdb_parser.get_structure(s2_rep, s2_rep)
          
        ref_model = ref_structure[0]
        sample_model = sample_structure[0]
  
        #ref_atoms = [at for chain in ref_model for res in chain for at in res]
        #sample_atoms = [at for chain in sample_model for res in chain for at in res]
  
        ref_res = [res for chain in ref_model for res in chain]
        sample_res = [res for chain in sample_model for res in chain]
        m = graphs_mapping(graph1, graph2)
        
        for chain in ref_model :
            num_chaine_ref = chain.id
            
        for chain in sample_model : 
            num_chaine_sample = chain.id

#         print(m)
        #print s1, s2
        #print i, j
        #print cluster[CLASS-1]['names'][i], cluster[CLASS-1]['names'][j]
        
        print(len(ref_res))
        print(ref_res)
        print(len(sample_res))
        print(sample_res)
        if len(sample_res) != len(ref_res) :
            print("probleme3")
            return
        
        for x in sample_res :
            print(x.id)
        m = {graphes[(s1.split("_")[1], s1.split("_")[2])].node[x]['fr3d']:graphes[(s2.split("_")[1], s2.split("_")[2])].node[m[x]]['fr3d']
             for x in m}
        print(m)
        ref_atoms = [] 
        sample_atoms = []
        
        new_sample_res = Bio.PDB.Chain.Chain(num_chaine_sample)
        new_ref_res = Bio.PDB.Chain.Chain(num_chaine_ref)
        print(new_sample_res)
        print(new_ref_res)
        
        for res1 in ref_res:
            pos1 =  str(res1.id[1])
            pos2 = m[pos1]
#             print(pos2)
            
            res2 = [x for x in sample_res if str(x.id[1]) == pos2][0]
            
            print("petit rat")
            print(res1)
            
            names_1 = [x.name for x in res1]
            names_2 = [x.name for x in res2]
            print(names_1)
            
            ref_atoms.extend([atom for atom in res1 if atom.name in names_2 ]) #and atom.name=="C1'"])
            sample_atoms.extend([atom for atom in res2 if atom.name in names_1]) #and atom.name=="C1'"])
            
            new_ref_res.add(res1)
            print(res2.id)
            res2.id = (' ', int(pos1), ' ')
            print(res2.id)
            new_sample_res.add(res2) 
            
            
        print(m)
        print(new_ref_res)
        print(new_sample_res)
        print("gros rat")
        print(ref_atoms)
        print(sample_atoms)
        for res in new_ref_res :
            print("ramou")
            print(res)
#         print(len(ref_atoms))
#         print(len(sample_atoms))
#         print(ref_atoms)
#         print(sample_atoms)
#         if (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_1U9S_A_58_11_2.pdb" and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_5J7L_DA_48_15_2.pdb") or (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_5J7L_DA_48_15_2.pdb" and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_1U9S_A_58_11_2.pdb") :
#             exit()
#         if (s2 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_1U9S_A_58_11_2.pdb' and s1 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_5J7L_DA_48_15_2.pdb') or (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_1U9S_A_58_11_2.pdb' and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_5J7L_DA_48_15_2.pdb'):
#             print("gros rat /n /n /n /n")
#             exit()
        super_imposer = Bio.PDB.Superimposer()
        super_imposer.set_atoms(sample_atoms, ref_atoms)
        super_imposer.apply(ref_model.get_atoms())
        print(super_imposer.rms)
        print(len(ref_atoms))
        print(len(sample_atoms))
        
        new_sample_model = Bio.PDB.Model.Model(0)
        new_ref_model = Bio.PDB.Model.Model(0)
        new_sample_model.add(new_sample_res)
        new_ref_model.add(new_ref_res)
        
        new_structure_sample = Bio.PDB.Structure.Structure(0)
        new_structure_ref = Bio.PDB.Structure.Structure(0)
        new_structure_ref.add(new_ref_model)
        new_structure_sample.add(new_sample_model)
        
        io_sample = Bio.PDB.PDBIO()
        io_sample.set_structure(new_structure_sample) 
        io_sample.save(s2_rep[:len(s2_rep)-4]+"_aligned_with_"+s1[:len(s1)-4]+".pdb")
        
        io_ref = Bio.PDB.PDBIO()
        io_ref.set_structure(new_structure_ref) 
        io_ref.save(s1_rep[:len(s1_rep)-4]+"_aligned_with_"+s2[:len(s2)-4]+".pdb")
        
        return super_imposer.rms 
  
    def get_rms_2(s1, s2, graph1, graph2, subgraph1, subgraph2, graphs1, graphs2, taille):
        with open(PATH_MMCIF+"fichiers_problemes.txt", 'a') as fichier :
            print(s1)
            print(s2)
            s1_rep = "/media/coline/Maxtor/Fichiers_mmcif/taille_%d/"%taille + s1
            s2_rep = "/media/coline/Maxtor/Fichiers_mmcif/taille_%d/"%taille  + s2
            print(s1_rep)
            print(s2)
            
            ref_structure = pdb_parser.get_structure(s1_rep, s1_rep)
            sample_structure = pdb_parser.get_structure(s2_rep, s2_rep)
              
            print(ref_structure)  
              
            ref_model = ref_structure[0]
            sample_model = sample_structure[0]
      
            #ref_atoms = [at for chain in ref_model for res in chain for at in res]
            #sample_atoms = [at for chain in sample_model for res in chain for at in res]
      
            ref_res = [res for chain in ref_model for res in chain]
            sample_res = [res for chain in sample_model for res in chain]
            
            print(ref_res)
            print(sample_res)
            print(len(ref_res))
            
            ##On commence par chercher a superposer les motifs seulement en minimisant la RMSD
            m = graphs_mapping(subgraph1, subgraph2)
            
            print(m)
            
            for chain in ref_model :
                num_chaine_ref = chain.id
                
            for chain in sample_model : 
                num_chaine_sample = chain.id
    
    #         print(m)
            #print s1, s2
            #print i, j
            #print cluster[CLASS-1]['names'][i], cluster[CLASS-1]['names'][j]
            
            print(len(ref_res))
            print(ref_res)
            print(len(sample_res))
            print(sample_res)
            
            
            if len(sample_res) != len(ref_res) :
                print("probleme3")
                return
            
            for x in sample_res :
                print(x.id)
            m = {graphs1.node[x]['fr3d']:graphs2.node[m[x]]['fr3d']
                 for x in m}
            print(m)
            ref_atoms = [] 
            sample_atoms = []
        
            print(ref_res)
            
            for res1 in ref_res:
                pos1 =  str(res1.id[1])
                print(pos1)
                if pos1 in m.keys() :
                    pos2 = m[pos1]
        #             print(pos2)
                    
                    res2 = [x for x in sample_res if str(x.id[1]) == pos2][0]
                    
                    print("petit rat")
                    print(res1)
                    
                    names_1 = [x.name for x in res1]
                    names_2 = [x.name for x in res2]
                    print(names_1)
                    
                    ref_atoms.extend([atom for atom in res1 if atom.name in names_2 and atom.name=="C3'"])
                    sample_atoms.extend([atom for atom in res2 if atom.name in names_1 and atom.name=="C3'"])
                
                
            print(m)
            
    #         print(new_ref_res)
    #         print(new_sample_res)
            print("gros rat")
            print(ref_atoms)
            print(sample_atoms)
            
    #         print(len(ref_atoms))
    #         print(len(sample_atoms))
    #         print(ref_atoms)
    #         print(sample_atoms)
    #         if (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_1U9S_A_58_11_2.pdb" and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_5J7L_DA_48_15_2.pdb") or (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_5J7L_DA_48_15_2.pdb" and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_1U9S_A_58_11_2.pdb") :
    #             exit()
    #         if (s2 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_1U9S_A_58_11_2.pdb' and s1 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_5J7L_DA_48_15_2.pdb') or (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_1U9S_A_58_11_2.pdb' and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_5J7L_DA_48_15_2.pdb'):
    #             print("gros rat /n /n /n /n")
    #             exit()
            print("avant")
            for atom in ref_atoms :
                print(atom.get_vector())
            for atom in sample_atoms :
                print(atom.get_vector())   
            
    
            super_imposer = Bio.PDB.Superimposer()
            super_imposer.set_atoms(sample_atoms, ref_atoms)
            
            print("apres")
            for atom in ref_atoms :
                print(atom.get_vector())
            for atom in sample_atoms :
                print(atom.get_vector())
            #super_imposer.apply(ref_model.get_atoms())
            print(super_imposer.rms)
            print(len(ref_atoms))
            print(len(sample_atoms))
            print(super_imposer.rotran)
            
            if taille == 1 :
#                 sample_model_pymol = Bio.PDB.Model.Model(0)
#                 ref_model_pymol = Bio.PDB.Model.Model(0)
#                 sample_model_pymol.add(sample_res)
#                 ref_model_pymol.add(ref_res)
#                    
#                 structure_sample_pymol = Bio.PDB.Structure.Structure(0)
#                 structure_ref_pymol = Bio.PDB.Structure.Structure(0)
#                 structure_ref_pymol.add(ref_model_pymol)
#                 structure_sample_pymol.add(sample_model_pymol)
#         
#                 io_sample = Bio.PDB.PDBIO()
#                 io_sample.set_structure(structure_sample_pymol) 
#                 io_sample.save(s2_rep[:len(s2_rep)-4]+"_aligned_with_"+s1[:len(s1)-4]+".pdb")
#                  
#                 io_ref = Bio.PDB.PDBIO()
#                 io_ref.set_structure(structure_ref_pymol) 
#                 io_ref.save(s1_rep[:len(s1_rep)-4]+"_aligned_with_"+s2[:len(s2)-4]+".pdb")# 
                
                return super_imposer.rms
            
            else :
            
                #On recupere la rotation qui correspond a cette meilleure RMSD
                rotation,translation = super_imposer.rotran
                
                
                ## On fait la superposition sur toute l'extension avec la rotation des atomes qu'on vient de recuperer et on calcule la RMSD de cette superposition-là
        
        #         ref_structure_new = pdb_parser.get_structure(s1_rep, s1_rep)
        #         sample_structure_new = pdb_parser.get_structure(s2_rep, s2_rep)
                  
        #         ref_model_new = ref_structure_new[0]
        #         sample_model_new = sample_structure_new[0]
          
                #ref_atoms = [at for chain in ref_model for res in chain for at in res]
                #sample_atoms = [at for chain in sample_model for res in chain for at in res]
          
        #         ref_res_new = [res for chain in ref_model_new for res in chain]
        #         sample_res_new = [res for chain in sample_model_new for res in chain]
        #         
        
                 
                m_pymol = graphs_mapping(graph1, graph2) ## superposition avec le motif pour la visualisation
                if m_pymol == -1 :
                    print("probleme2 : non isomorphic graphs (peut etre branches qui se rejoignent)")
                    return 
                a_enlever = []
                for elt in m_pymol.keys() :
                    if elt in subgraph1.nodes() :
                        a_enlever.append(elt)
                  
                m_calcul = copy.deepcopy(m_pymol)   ##superposition sans le motif pour le calcul de RMSD
                for elt in a_enlever :
                    del(m_calcul[elt])
                print("tout petit rat")
                print(m_pymol)
                print(m_calcul)
                  
        #         for chain in ref_model_new :
        #             num_chaine_ref = chain.id
        #               
        #         for chain in sample_model_new : 
        #             num_chaine_sample = chain.id
        #  
        # #         print(m)
        #         #print s1, s2
        #         #print i, j
        #         #print cluster[CLASS-1]['names'][i], cluster[CLASS-1]['names'][j]
        #          
        #         print(len(ref_res_new))
        #         print(ref_res_new)
        #         print(len(sample_res_new))
        #         print(sample_res_new)
        #         if len(sample_res_new) != len(ref_res_new) :calcule la RMSD  
        #             print("probleme")
        #             return
        #           
        #         for x in sample_res_new :
        #             print(x.id)
                m_pymol = {graphs1.node[x]['fr3d']:graphs2.node[m_pymol[x]]['fr3d']
                     for x in m_pymol}
                
                m_calcul = {graphs1.node[x]['fr3d']:graphs2.node[m_calcul[x]]['fr3d']
                     for x in m_calcul}
                
                ref_atoms_calcul = [] 
                sample_atoms_calcul = []
                  
                sample_res_pymol = Bio.PDB.Chain.Chain(num_chaine_sample)
                ref_res_pymol = Bio.PDB.Chain.Chain(num_chaine_ref)
                #print(sample_res_pymol)
                #print(ref_res_pymol)
                
                #print(m_pymol)
                #print(m_calcul)
                
                #print(ref_res)
                
                compteur = 1
                for res1 in ref_res:
                    pos1 =  str(res1.id[1])
                    #print(pos1)
                    if pos1 in m_pymol.keys() :
                        pos2 = m_pymol[pos1]
                        #print(pos2)
                        
                        #for x in sample_res :
                        #    print(x.id)
                           
                        tab_res2 = [x for x in sample_res if str(x.id[1]) == pos2]
                        if len(tab_res2) > 0 :
                            res2 = tab_res2[0]
                          
                            #print("petit rat")
                            #print(res1)
                              
                            names_1 = [x.name for x in res1]
                            names_2 = [x.name for x in res2]
                            print(names_1)
                            
                            
                            res1_pymol = Bio.PDB.Residue.Residue((' ',compteur , ' '), res1.get_resname(), res1.get_segid())
                            for atom in res1 :
                                #if atom.name in names_2 :
                                #print(atom.get_vector())
                                atom.transform(rotation, translation)
                                    
                                res1_pymol.add(atom)
                                
                    
                            #res1.id = (' ',compteur , ' ')  
                            ref_res_pymol.add(res1_pymol)
                            #print(res2.id)
                            #res2.id = (' ', compteur, ' ')
                            res2_pymol = Bio.PDB.Residue.Residue((' ',compteur , ' '), res2.get_resname(), res2.get_segid())
                            for atom in res2 :
                                #if atom.name in names_1 :
                                print(atom.get_vector())
                                res2_pymol.add(atom)
                            #print(res2.id)
                            sample_res_pymol.add(res2_pymol) 
                            
                            compteur += 1
                            
                            if pos1 in m_calcul.keys() :
                                if m_calcul[pos1] != pos2 : ## pas meme superposition entre visualisation et calcul => probleme
                                    print("probleme4")
                                    print(m_pymol)
                                    print(m_calcul)
                                    print(m_pymol[pos1])
                                    print(pos2)
                                    exit()
                                else :
                                    ref_atoms_calcul.extend([atom for atom in res1 if atom.name in names_2 and atom.name=="C3'"])
                                    sample_atoms_calcul.extend([atom for atom in res2 if atom.name in names_1 and atom.name=="C3'"])
                        else :
                            fichier.write(str(s2_rep)+ " " + str(pos2) +"\n")
                
                liste_atomes_ref = []
                for atom in ref_atoms_calcul :
                    #atom.transform(rotation, translation)   
                    liste_atomes_ref.append(list(atom.get_vector()))
                    #print(list(atom.get_vector()))
                
                liste_atomes_sample = []    
                for atom in sample_atoms_calcul :
        #             atom.transform(rotation, translation)
                    liste_atomes_sample.append(list(atom.get_vector()))
                
                #print(liste_atomes_ref)
                
                print("avant")
                for res in ref_res_pymol :
                    for atom in res :
                        print(atom.get_vector())
                        #break
                
        #         new_new_ref_res = Bio.PDB.Chain.Chain(num_chaine_ref)
        #         for res in new_ref_res :
        #             #res2 = [x for x in new_sample_res if str(x.id[1]) == str(res.id[1])][0]
        # #             new_res = Bio.PDB.Residue.Residue(res.id[0], res.id[1], res.id[2])
        #             for atom in res :
        #                 #if atom.name in [x.name for x in res2] :
        #                     atom.transform(rotation, translation)  
        #                     new_res.add(atom)
        #             new_new_ref_res.add(new_res)       
                            
                            
                print("avpres")
                for res in ref_res_pymol :
                    for atom in res :
                        print(atom.get_vector())
                        #break
                        
                for res in sample_res_pymol :
                    for atom in res :
                        print(atom.get_vector())
                        #break
                            
        #         for res in new_sample_res :
        #             res2 = [x for x in new_ref_res if str(x.id[1]) == str(res.id[1])][0]
        #             for atom in res :
        #                 if atom.name in [x.name for x in res2] :
        #                     atom.transform(rotation, translation)  
                      
                print(m)
                print(ref_res_pymol)
                print(sample_res_pymol)
                print("gros rat")
                print(ref_atoms)
                print(sample_atoms)
                for res in ref_res_pymol :
                    print("ramou")
                    print(res)
                    
                    
                sample_model_pymol = Bio.PDB.Model.Model(0)
                ref_model_pymol = Bio.PDB.Model.Model(0)
                sample_model_pymol.add(sample_res_pymol)
                ref_model_pymol.add(ref_res_pymol)
                   
                structure_sample_pymol = Bio.PDB.Structure.Structure(0)
                structure_ref_pymol = Bio.PDB.Structure.Structure(0)
                structure_ref_pymol.add(ref_model_pymol)
                structure_sample_pymol.add(sample_model_pymol)
        
                io_sample = Bio.PDB.PDBIO()
                io_sample.set_structure(structure_sample_pymol) 
                io_sample.save(s2_rep[:len(s2_rep)-4]+"_aligned_with_"+s1[:len(s1)-4]+".pdb")
                 
                io_ref = Bio.PDB.PDBIO()
                io_ref.set_structure(structure_ref_pymol) 
                io_ref.save(s1_rep[:len(s1_rep)-4]+"_aligned_with_"+s2[:len(s2)-4]+".pdb")#   
                
                print(liste_atomes_ref)
                print(liste_atomes_sample)
        
        #         print(len(ref_atoms))
        #         print(len(sample_atoms))
        #         print(ref_atoms)
        #         print(sample_atoms)
        #         if (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_1U9S_A_58_11_2.pdb" and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_5J7L_DA_48_15_2.pdb") or (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_5J7L_DA_48_15_2.pdb" and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6] + "fichier_1U9S_A_58_11_2.pdb") :
        #             exit()
        #         if (s2 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_1U9S_A_58_11_2.pdb' and s1 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_5J7L_DA_48_15_2.pdb') or (s1 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_1U9S_A_58_11_2.pdb' and s2 == PATH_MMCIF[:len(PATH_MMCIF)-6]+'fichier_5J7L_DA_48_15_2.pdb'):
        #             print("gros rat /n /n /n /n")
        #             exit()
                
                print(calcul_rmsd(liste_atomes_ref, liste_atomes_sample))
                return calcul_rmsd(liste_atomes_ref, liste_atomes_sample)
       
    compteur = 0
    rms = {}
    for i in range(len(tab_fichiers)) :
        with open("Graphs/"+tab_fichiers[i].split("_")[1]+".pickle", 'rb') as fichier_graphs1 :
            mon_depickler_gg1 = pickle.Unpickler(fichier_graphs1)
            graphs1 = mon_depickler_gg1.load()
        
            print(tab_motif_positions[i])
            subgraph1_motif = nx.DiGraph(graphs1.subgraph(tab_motif_positions[i]))
            graphe1_motif = create_empty_copy(subgraph1_motif)
            
            for noeud in graphe1_motif.nodes() :
                    if (noeud[0], noeud[1]+1) in graphe1_motif.nodes() and (noeud,(noeud[0], noeud[1]+1)) not in graphe1_motif.edges() and noeud != tab_motif_positions[i][4] and (noeud[0], noeud[1]+1) != tab_motif_positions[i][4]  :
                        graphe1_motif.add_edge(noeud, (noeud[0], noeud[1]+1), label='B53', long_range=False)
            
            graphe1_motif.add_edge(tab_motif_positions[i][0], tab_motif_positions[i][1], label='CSS', long_range=True)
            graphe1_motif.add_edge(tab_motif_positions[i][1], tab_motif_positions[i][0], label='CSS', long_range=True)
            
            graphe1_motif.add_edge(tab_motif_positions[i][2], tab_motif_positions[i][3], label='CSS', long_range=True)
            graphe1_motif.add_edge(tab_motif_positions[i][3], tab_motif_positions[i][2], label='CSS', long_range=True)
            
            graphe1_motif.add_edge(tab_motif_positions[i][0], tab_motif_positions[i][4], label='TSS', long_range=True)
            graphe1_motif.add_edge(tab_motif_positions[i][4], tab_motif_positions[i][0], label='TSS', long_range=True)
        
            graphe1_motif.add_edge(tab_motif_positions[i][1], tab_motif_positions[i][4], label='CWW', long_range=False)
            graphe1_motif.add_edge(tab_motif_positions[i][4], tab_motif_positions[i][1], label='CWW', long_range=False)
            
            subgraph = nx.DiGraph(graphs1.subgraph(tab_fichier_positions[i]))
            graph1 = create_empty_copy(subgraph)
                
            for noeud in graph1.nodes() :
                    if (noeud[0], noeud[1]+1) in graph1.nodes() and (noeud,(noeud[0], noeud[1]+1)) not in graph1.edges() and noeud != tab_motif_positions[i][4] and (noeud[0], noeud[1]+1) != tab_motif_positions[i][4]  :
                        graph1.add_edge(noeud, (noeud[0], noeud[1]+1), label='B53', long_range=False)
            
            graph1.add_edge(tab_motif_positions[i][0], tab_motif_positions[i][1], label='CSS', long_range=True)
            graph1.add_edge(tab_motif_positions[i][1], tab_motif_positions[i][0], label='CSS', long_range=True)
            
            graph1.add_edge(tab_motif_positions[i][2], tab_motif_positions[i][3], label='CSS', long_range=True)
            graph1.add_edge(tab_motif_positions[i][3], tab_motif_positions[i][2], label='CSS', long_range=True)
            
            graph1.add_edge(tab_motif_positions[i][0], tab_motif_positions[i][4], label='TSS', long_range=True)
            graph1.add_edge(tab_motif_positions[i][4], tab_motif_positions[i][0], label='TSS', long_range=True)
        
            graph1.add_edge(tab_motif_positions[i][1], tab_motif_positions[i][4], label='CWW', long_range=False)
            graph1.add_edge(tab_motif_positions[i][4], tab_motif_positions[i][1], label='CWW', long_range=False)
    
                    
            for j in range(i+1, len(tab_fichiers)) :
                print(compteur)
                with open("Graphs/"+tab_fichiers[j].split("_")[1]+".pickle", 'rb') as fichier_graphs2 :
                    mon_depickler_gg2 = pickle.Unpickler(fichier_graphs2)
                    graphs2 = mon_depickler_gg2.load()
                    
                    print(tab_motif_positions[j])
                    subgraph2_motif = nx.DiGraph(graphs2.subgraph(tab_motif_positions[j]))
                   
                    graphe2_motif = create_empty_copy(subgraph2_motif)
                    
                    for noeud in graphe2_motif.nodes() :
                        if (noeud[0], noeud[1]+1) in graphe2_motif.nodes() and (noeud,(noeud[0], noeud[1]+1)) not in graphe2_motif.edges() and noeud != tab_motif_positions[j][4] and (noeud[0], noeud[1]+1) != tab_motif_positions[j][4] :
                            graphe2_motif.add_edge(noeud, (noeud[0], noeud[1]+1), label='B53', long_range=False)
                
                    graphe2_motif.add_edge(tab_motif_positions[j][0], tab_motif_positions[j][1], label='CSS', long_range=True)
                    graphe2_motif.add_edge(tab_motif_positions[j][1], tab_motif_positions[j][0], label='CSS', long_range=True)
                    
                    graphe2_motif.add_edge(tab_motif_positions[j][2], tab_motif_positions[j][3], label='CSS', long_range=True)
                    graphe2_motif.add_edge(tab_motif_positions[j][3], tab_motif_positions[j][2], label='CSS', long_range=True)
                    
                    graphe2_motif.add_edge(tab_motif_positions[j][0], tab_motif_positions[j][4], label='TSS', long_range=True)
                    graphe2_motif.add_edge(tab_motif_positions[j][4], tab_motif_positions[j][0], label='TSS', long_range=True)
                
                    graphe2_motif.add_edge(tab_motif_positions[j][1], tab_motif_positions[j][4], label='CWW', long_range=False)
                    graphe2_motif.add_edge(tab_motif_positions[j][4], tab_motif_positions[j][1], label='CWW', long_range=False)
                    
                    subgraph = nx.DiGraph(graphs2.subgraph(tab_fichier_positions[j]))
                    graph2 = create_empty_copy(subgraph)
                    
                    for noeud in graph2.nodes() :
                        if (noeud[0], noeud[1]+1) in graph2.nodes() and (noeud,(noeud[0], noeud[1]+1)) not in graph2.edges() and noeud != tab_motif_positions[j][4] and (noeud[0], noeud[1]+1) != tab_motif_positions[j][4]  :
                            graph2.add_edge(noeud, (noeud[0], noeud[1]+1), label='B53', long_range=False)
                            
                    graph2.add_edge(tab_motif_positions[j][0], tab_motif_positions[j][1], label='CSS', long_range=True)
                    graph2.add_edge(tab_motif_positions[j][1], tab_motif_positions[j][0], label='CSS', long_range=True)
                    
                    graph2.add_edge(tab_motif_positions[j][2], tab_motif_positions[j][3], label='CSS', long_range=True)
                    graph2.add_edge(tab_motif_positions[j][3], tab_motif_positions[j][2], label='CSS', long_range=True)
                    
                    graph2.add_edge(tab_motif_positions[j][0], tab_motif_positions[j][4], label='TSS', long_range=True)
                    graph2.add_edge(tab_motif_positions[j][4], tab_motif_positions[j][0], label='TSS', long_range=True)
                
                    graph2.add_edge(tab_motif_positions[j][1], tab_motif_positions[j][4], label='CWW', long_range=False)
                    graph2.add_edge(tab_motif_positions[j][4], tab_motif_positions[j][1], label='CWW', long_range=False)
                        
        #             print(tab_nb_nts_par_chaine[i])
        #             print(graph1.number_of_edges())
        #             print(graph1.edges.data())
        #             print(tab_motif_positions[i])
        #             print(tab_nb_nts_par_chaine[j])
        #             print(graph2.number_of_edges())
        #             print(graph2.edges.data())
        #             print(tab_motif_positions[j])
                    #print(graphes[(tab_fichiers[j].split("_")[1], tab_fichiers[j].split("_")[2])].nodes.data())
                    
                    if graph1.number_of_edges() != graph2.number_of_edges() :
                        print(tab_fichiers[i])
                        print(tab_fichiers[j])
                        print(graph1.number_of_edges())
                        print(graph1.edges.data())
                        print(graph2.number_of_edges())
                        print("probleme1")
                        
                    
                    
                    print(graphe1_motif.nodes.data())
                    print(graphe2_motif.nodes.data())
                    
                    compteur +=1
                    
                    rms.update({(tab_fichiers[i], tab_fichiers[j]) : get_rms_2(tab_fichiers[i], tab_fichiers[j], graph1, graph2, graphe1_motif, graphe2_motif, graphs1, graphs2,taille) })
            
#             if (tab_fichiers[i], tab_fichiers[j]) == ('fichier_3JCS_1_282_1_taille_4.pdb', 'fichier_4V9F_0_25_56_taille_4.pdb') :
#                 print(rms[(tab_fichiers[i], tab_fichiers[j])])
#                 exit()
    
    return rms

'''calcule l'exp de -x^2 '''
def function(x): 
    return math.exp((-x)**2) 

''' affiche la distribution des valeurs de RMSD selon une certaine taille d'extension
et calcule la p-value pour chaque valeur (mais en fait, c'est pas terrible) '''
def distrib_rmsd(taille_ext, i):
    
    homologues = [['4V9F_0_62_12', '5J7L_DA_62_5', '5DM6_X_328_2', '5FDU_1A_62_14'],
                  ['4V9F_0_30_4', '5J7L_DA_30_15', '5FDU_1A_30_17'],
                  ['4V9F_0_48_13', '5J7L_DA_48_20', '5DM6_X_48_28', '5FDU_1A_48_25'],
                  ['4V9F_0_207_3', '5J7L_DA_272_2', '5FDU_1A_272_1'],
                  ['4V9F_0_25_56', '5J7L_DA_25_12', '5FDU_1A_25_68', '4V88_A5_25_30'],
                  ['4V9F_0_137_5', '5J7L_DA_48_1', '5FDU_1A_137_6', '5DM6_X_48_10', '4V88_A5_48_3'],
                  ['4V9F_0_134_5', '5FDU_1A_74_7', '5DM6_X_127_7'],
                  ['4V9F_0_48_21', '5J7L_DA_197_4', '5FDU_1A_197_3', '5DM6_X_48_9'],
                  ['4V9F_0_127_6', '5J7L_DA_134_1', '5FDU_1A_134_3', '5DM6_X_134_2'],
                  ['4V9F_0_287_2', '5DM6_X_25_15'],
                  ['5J7L_DA_25_10', '5DM6_X_25_34', '4V88_A5_25_47', '5FDU_1A_25_78'],
                  ['1FJG_A_48_17', '5J5B_BA_48_14'],
                  ['1FJG_A_48_8', '5J5B_BA_48_23'],
                  ['1FJG_A_294_1', '5J5B_BA_294_2'],
                  ['1FJG_A_58_23', '5J5B_BA_58_3'],
                  ['1FJG_A_138_3', '5J5B_BA_138_2'],
                  ['5J5B_BA_48_7', '4V88_A6_48_12'],
                  ['4YAZ_R_36_25', '3UCZ_R_62_15']]
    
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_1.pickle"%taille_ext, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        dico_rmsd = mon_depickler.load()
        
        with open(EXTENSION_PATH%i+"sim_extensions_toutes_aretes_coeff_all1_taille_%i.pickle"%i, 'rb') as fichier_sim :
            mon_depickler_2 = pickle.Unpickler(fichier_sim)
            dico_sim = mon_depickler_2.load()
            
            tab_sim_1 = []
            tab_sim_2 = []
            tab_rmsd_1 = []
            tab_rmsd_2 = []
            
            tab_rmsd_int_g11 = []
            tab_sim_int_g11 = []
            tab_rmsd_int_g12 = []
            tab_sim_int_g12 = []
            tab_sim_ext = []
            tab_rmsd_ext = []
            
            array_dico_rmsd = []
            for elt in dico_rmsd.values() :
                if elt != None :
                    array_dico_rmsd.append(elt)
            
            print(array_dico_rmsd)
            sns.distplot(array_dico_rmsd, kde=False)
            
            axes = plt.gca()
# #             plt.scatter(tab_rmsd_1, tab_sim_1, color='blue')
# #             plt.scatter(tab_rmsd_2, tab_sim_2, color='red')
            plt.title("Distribution des RMSD pour une taille d'extension de 4")
#             #plt.plot([0,max(array_rmsd)], [0.7,0.7], color='red')
            axes.set_xlabel("Valeurs de RMSD")
            axes.set_ylabel("Nombre")
            
            plt.show()
            
            
            #groupe = ['5FDU_1A_48_19', '5J7L_DA_30_15', '4V9F_0_30_4', '1FJG_A_48_8', '5J7L_DA_48_15', '5J5B_BA_48_23', '5FDU_1A_30_17', '1U9S_A_58_11']
            groupe_1 = GROUPE_GNRA
            groupe_1.extend(GROUPE_GNRA_ETENDU)
            
            groupe_1 = [x for x in groupe_1 if x not in ['5J7L_DA_25_10', '5DM6_X_25_34', '3JCS_1_25_46', '5J7L_DA_25_12', '4V9F_0_25_56', '4V88_A5_25_30']]

            
            groupe_2 = CLUSTERING_PEREZ_VERSION_NON_CAN_2[12]
            for elt in CLUSTERING_PEREZ_VERSION_NON_CAN_2[11] :
                if elt not in groupe_2 :
                    groupe_2.append(elt)
            for elt in groupe_2 :
                print(elt)
                groupe_1.remove(elt)
            #groupe_2.extend(GROUPE_ARICH_ETENDU)
            #groupe = ['5FDU_1A_30_17', '1U9S_A_58_11']
            for cle in dico_sim.keys() :
                nom_1 = "fichier_"+cle[0]+"_taille_%d.pdb"%i
                nom_2 = "fichier_"+cle[1]+"_taille_%d.pdb"%i
                p_value = 0
                if nom_1 != 'fichier_2XD0_V_36_21_taille_%d.pdb'%i and nom_2 != "fichier_2XD0_V_36_21_taille_%d.pdb"%i:
                    #if dico_sim[cle] < 0.4 : #and dico_sim[cle] < 0.6 : 
                        if cle[0] in groupe_1 and cle[1] in groupe_1 :
                            print(cle)
                            print(dico_sim[cle])
                            
                            rmsd_standard = 5.1*(29**0.41)-19.8
                            
                            if dico_sim[cle] < 0.7 :
                                print('ramousnif')
                            #print(dico_rmsd.keys())
                            
                            if (nom_1,nom_2) in dico_rmsd.keys() :
                                if dico_rmsd[(nom_1, nom_2)] != None :
                                    print(dico_rmsd[(nom_1, nom_2)])
                                    tab_rmsd_int_g11.append(dico_rmsd[(nom_1, nom_2)])
                                    tab_sim_int_g11.append(dico_sim[cle])
                                    print(dico_rmsd[(nom_1, nom_2)])
                                    z = (dico_rmsd[(nom_1, nom_2)] - rmsd_standard)/1.8
#                                     print(rmsd_standard)
#                                     print(z)
                                    integrale, err = quad(function, 0, z/math.sqrt(2))
#                                     print(integrale)
#                                     print(err)
                                    fct_erreur = (2/math.sqrt(math.pi))*integrale
                                    #print(fct_erreur)
                                    p_value = (1+fct_erreur)/2
                                else :
                                    print("ramous")    
                                    
#                                     if (cle[0] in ['5J7L_DA_30_15', '5FDU_1A_30_17', '4V9F_0_30_4'] and cle[1] in ['5FDU_1A_48_19', '5J7L_DA_48_15', '1U9S_A_58_11']) or (cle[1] in ['5J7L_DA_30_15', '5FDU_1A_30_17', '4V9F_0_30_4'] and cle[0] in ['5FDU_1A_48_19', '5J7L_DA_48_15', '1U9S_A_58_11']) :
#                                         tab_rmsd_2.append(dico_rmsd[(nom_1, nom_2)])
#                                         tab_sim_2.append(dico_sim[cle])
#                                     else :
#                                         tab_rmsd_1.append(dico_rmsd[(nom_1, nom_2)])
#                                         tab_sim_1.append(dico_sim[cle])
                                    
                            else :
                                if dico_rmsd[(nom_2, nom_1)] != None :
                                    print(dico_rmsd[(nom_2, nom_1)])
                                    tab_rmsd_int_g11.append(dico_rmsd[(nom_2, nom_1)])
                                    tab_sim_int_g11.append(dico_sim[cle])
                                    z = (dico_rmsd[(nom_2, nom_1)] - rmsd_standard)/1.8
                                    integrale, err = quad(function, 0, z/math.sqrt(2))
                                    fct_erreur = (2/math.sqrt(math.pi))*integrale
                                    p_value = (1+fct_erreur)/2
                                else :
                                    print("ramous")
                                    
                                    
                        if cle[0] in groupe_2 and cle[1] in groupe_2 :
                            
                            print(cle)
                            print(dico_sim[cle])
                            
                            rmsd_standard = 5.1*(12**0.41)-19.8
                            
                            if dico_sim[cle] < 0.7 :
                                print('ramous')
                            #print(dico_rmsd.keys())
                            
                            if (nom_1,nom_2) in dico_rmsd.keys() :
                                if dico_rmsd[(nom_1, nom_2)] != None :
                                    print(dico_rmsd[(nom_1, nom_2)])
                                    tab_rmsd_int_g12.append(dico_rmsd[(nom_1, nom_2)])
                                    tab_sim_int_g12.append(dico_sim[cle])
                                    print(dico_rmsd[(nom_1, nom_2)])
                                    z = (dico_rmsd[(nom_1, nom_2)] - rmsd_standard)/1.8
#                                     print(rmsd_standard)
#                                     print(z)
                                    integrale, err = quad(function, 0, z/math.sqrt(2))
#                                     print(integrale)
#                                     print(err)
                                    fct_erreur = (2/math.sqrt(math.pi))*integrale
                                    #print(fct_erreur)
                                    p_value = (1+fct_erreur)/2
                                    
                                    
#                                     if (cle[0] in ['5J7L_DA_30_15', '5FDU_1A_30_17', '4V9F_0_30_4'] and cle[1] in ['5FDU_1A_48_19', '5J7L_DA_48_15', '1U9S_A_58_11']) or (cle[1] in ['5J7L_DA_30_15', '5FDU_1A_30_17', '4V9F_0_30_4'] and cle[0] in ['5FDU_1A_48_19', '5J7L_DA_48_15', '1U9S_A_58_11']) :
#                                         tab_rmsd_2.append(dico_rmsd[(nom_1, nom_2)])
#                                         tab_sim_2.append(dico_sim[cle])
#                                     else :
#                                         tab_rmsd_1.append(dico_rmsd[(nom_1, nom_2)])
#                                         tab_sim_1.append(dico_sim[cle])
                                    
                            else :
                                if dico_rmsd[(nom_2, nom_1)] != None :
                                    print(dico_rmsd[(nom_2, nom_1)])
                                    tab_rmsd_int_g12.append(dico_rmsd[(nom_2, nom_1)])
                                    tab_sim_int_g12.append(dico_sim[cle])
                                    z = (dico_rmsd[(nom_2, nom_1)] - rmsd_standard)/1.8
                                    integrale, err = quad(function, 0, z/math.sqrt(2))
                                    fct_erreur = (2/math.sqrt(math.pi))*integrale
                                    p_value = (1+fct_erreur)/2
#                                     if (cle[0] in ['5J7L_DA_30_15', '5FDU_1A_30_17', '4V9F_0_30_4'] and cle[1] in ['5FDU_1A_48_19', '5J7L_DA_48_15', '1U9S_A_58_11']) or (cle[1] in ['5J7L_DA_30_15', '5FDU_1A_30_17', '4V9F_0_30_4'] and cle[0] in ['5FDU_1A_48_19', '5J7L_DA_48_15', '1U9S_A_58_11']) :
#                                         tab_rmsd_2.append(dico_rmsd[(nom_2, nom_1)])
#                                         tab_sim_2.append(dico_sim[cle])
#                                     else :
#                                         tab_rmsd_1.append(dico_rmsd[(nom_2, nom_1)])
#                                         tab_sim_1.append(dico_sim[cle])
                                    
                
                            print(p_value)
                        if (cle[0] in groupe_1 and not cle[0] in groupe_2 and cle[1] in groupe_2 and not cle[1] in groupe_1) or (cle[0] in groupe_2 and not cle[0] in groupe_1 and cle[1] in groupe_1 and not cle[1] in groupe_2) :
                            print(cle)
                            print(dico_sim[cle])
                            
                            rmsd_standard = 5.1*(12**0.41)-19.8
                            
                            if dico_sim[cle] < 0.7 :
                                print('ramousnif')
                            #print(dico_rmsd.keys())
                            
                            if (nom_1,nom_2) in dico_rmsd.keys() :
                                if dico_rmsd[(nom_1, nom_2)] != None :
                                    print(dico_rmsd[(nom_1, nom_2)])
                                    tab_rmsd_ext.append(dico_rmsd[(nom_1, nom_2)])
                                    tab_sim_ext.append(dico_sim[cle])
                                    print(dico_rmsd[(nom_1, nom_2)])
                                    z = (dico_rmsd[(nom_1, nom_2)] - rmsd_standard)/1.8
#                                     print(rmsd_standard)
#                                     print(z)
                                    integrale, err = quad(function, 0, z/math.sqrt(2))
#                                     print(integrale)
#                                     print(err)
                                    fct_erreur = (2/math.sqrt(math.pi))*integrale
                                    #print(fct_erreur)
                                    p_value = (1+fct_erreur)/2
                                    
                                    
#                                     if (cle[0] in ['5J7L_DA_30_15', '5FDU_1A_30_17', '4V9F_0_30_4'] and cle[1] in ['5FDU_1A_48_19', '5J7L_DA_48_15', '1U9S_A_58_11']) or (cle[1] in ['5J7L_DA_30_15', '5FDU_1A_30_17', '4V9F_0_30_4'] and cle[0] in ['5FDU_1A_48_19', '5J7L_DA_48_15', '1U9S_A_58_11']) :
#                                         tab_rmsd_2.append(dico_rmsd[(nom_1, nom_2)])
#                                         tab_sim_2.append(dico_sim[cle])
#                                     else :
#                                         tab_rmsd_1.append(dico_rmsd[(nom_1, nom_2)])
#                                         tab_sim_1.append(dico_sim[cle])
                                    
                            else :
                                if dico_rmsd[(nom_2, nom_1)] != None :
                                    print(dico_rmsd[(nom_2, nom_1)])
                                    tab_rmsd_ext.append(dico_rmsd[(nom_2, nom_1)])
                                    tab_sim_ext.append(dico_sim[cle])
                                    z = (dico_rmsd[(nom_2, nom_1)] - rmsd_standard)/1.8
                                    integrale, err = quad(function, 0, z/math.sqrt(2))
                                    fct_erreur = (2/math.sqrt(math.pi))*integrale
                                    p_value = (1+fct_erreur)/2
#                                     if (cle[0] in ['5J7L_DA_30_15', '5FDU_1A_30_17', '4V9F_0_30_4'] and cle[1] in ['5FDU_1A_48_19', '5J7L_DA_48_15', '1U9S_A_58_11']) or (cle[1] in ['5J7L_DA_30_15', '5FDU_1A_30_17', '4V9F_0_30_4'] and cle[0] in ['5FDU_1A_48_19', '5J7L_DA_48_15', '1U9S_A_58_11']) :
#                                         tab_rmsd_2.append(dico_rmsd[(nom_2, nom_1)])
#                                         tab_sim_2.append(dico_sim[cle])
#                                     else :
#                                         tab_rmsd_1.append(dico_rmsd[(nom_2, nom_1)])
#                                         tab_sim_1.append(dico_sim[cle])
                                    
                
                            print(p_value)
                            #exit()
            print(type(dico_rmsd))
            #print(dico_rmsd.values())
            print(type(dico_rmsd.values()))
            
            
#             print(dico_rmsd[('fichier_5FDU_1A_30_17_2.pdb', 'fichier_5J5B_BA_48_23_2.pdb')])
#             #print(dico_rmsd.items())
#             dico_trie = sorted((elt for elt in dico_rmsd.items() if elt[1] != None), key=lambda t: t[1])
#             
#             compteur = 0
#             for cle in dico_trie :
#                 nom_1 = cle[0][0][8:len(cle[0][0])-6]
#                 nom_2 = cle[0][1][8:len(cle[0][1])-6]
#     #             print(nom_1)
#     #             print(nom_2)
#             
#                         
#                 if nom_1 in GROUPES_TOUTES_ARETES_MAX_4_10_07[2]  and nom_2  in GROUPES_TOUTES_ARETES_MAX_4_10_07[2] :
#                     print(cle)
#                     print(compteur)
#                 compteur += 1
            print(min([elt for elt in dico_rmsd.values() if elt != None]))
            print(max([elt for elt in dico_rmsd.values() if elt != None]))
 
            print(tab_sim_int_g11)
            print(tab_rmsd_int_g11)
            print(tab_sim_int_g12)
            print(tab_rmsd_int_g12)
            print(tab_sim_ext)
            print(tab_rmsd_ext)
    #         for elt in dico_rmsd.keys() :
    #             if dico_rmsd[elt] != None and dico_rmsd[elt] > 4.0 :
    #                 print(elt)
    #                 print(dico_rmsd[elt])
            
            #plt.hist([elt for elt in dico_rmsd.values() if elt != None], range = (1,20), bins=20,  edgecolor = 'black')
#             array_sim = np.array(tab_sim)
#             array_rmsd = np.array(tab_rmsd)
#             print(tab_rmsd)
#             print(len(array_sim))
#             print(array_sim)
#             print(len(array_rmsd))
#             print(array_rmsd)
#             moy_sim = array_sim.mean()
#             moy_rmsd = array_rmsd.mean()
#             
#             covariance = 0
#             for i in range(len(tab_sim)) :
#                 covariance += float(tab_sim[i])*float(tab_rmsd[i]) - moy_sim*moy_rmsd
#             covariance = covariance/len(tab_sim)
#             print(moy_sim)
#             print(moy_rmsd)
#             print(covariance)
#             print(array_rmsd.std())
#             print(array_sim.std())
#             print(covariance/(array_rmsd.std()*array_sim.std()))
            
            axes = plt.gca()
# #             plt.scatter(tab_rmsd_1, tab_sim_1, color='blue')
# #             plt.scatter(tab_rmsd_2, tab_sim_2, color='red')
            plt.scatter(tab_rmsd_ext, tab_sim_ext)
            plt.title("Comparaison RMSD et similarite pour les groupes clustering perez 11 et 12 \n par rapport aux autres GNRA (taille 4)")
#             #plt.plot([0,max(array_rmsd)], [0.7,0.7], color='red')
            axes.set_xlabel("RMSD")
            axes.set_ylabel("Similarite")
#             
#             #plt.plot(tab_sim)
#             
#             
            plt.savefig(PATH_MMCIF+"distrib_rmsd_sim_groupe_clustering_perez_11_12_autres_taille_4.png")
            #plt.show()
            
            print(len(tab_sim_2))
            print(len(tab_sim_1))
            print(len(groupe_1))
            
            with open(PATH_MMCIF+"fichier_csv_groupes_gnra_12_autres.csv", 'w', newline="") as fichier_csv :
                csvwriter = csv.writer(fichier_csv)
                csvwriter.writerow(["groupe gnra ", "groupe autres", "entre les deux"])
                compteur = 0
                
                
                
                for elt in tab_rmsd_int_g11 :
                    if compteur < len(tab_rmsd_ext) and compteur < len(tab_rmsd_int_g12) :
                        csvwriter.writerow([elt, tab_rmsd_int_g12[compteur], tab_rmsd_ext[compteur]])
                    elif compteur < len(tab_rmsd_int_g12) :
                        csvwriter.writerow([elt, tab_rmsd_int_g12[compteur], ""])
                    elif compteur < len(tab_rmsd_ext) :
                        csvwriter.writerow([elt,"", tab_rmsd_ext[compteur]])
                    else :
                        csvwriter.writerow([elt, "", ""])
                    compteur += 1
            
#             print(dico_rmsd.keys())
#             print(dico_rmsd[('fichier_4YAZ_R_36_25_taille_4.pdb', 'fichier_5DM6_X_25_34_taille_4.pdb')])
#             print(dico_sim[( '4YAZ_R_36_25', '5DM6_X_25_34')])
            #plt.plot()
        
def distrib_rmsd_new_data_par_cluster(num_ARN):
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load()
            
        with open("/media/coline/Maxtor/clustering_perez_%s_new_data_res_inf_3.pickle"%num_ARN, 'rb') as fichier_sortie :
            mon_depickler = pickle.Unpickler(fichier_sortie)
            clusters = mon_depickler.load()
            
            liste_arn = []
            for cluster in clusters :
                for elt in cluster :
                    if elt not in liste_arn :
                        liste_arn.append(elt)
            
            #print(rmsd)
            print(num_ARN)
            print(len(clusters))
            
            liste_paires = []
            liste_val_rmsd_par_cluster = []
            liste_val_rmsd_inter_cluster = []
            
#             for cle in rmsd.keys() :
#                 meme_cluster = False
#                 noeud1 = (cle[0].split("_")[1], int(cle[0].split("_")[2]))
#                 noeud2 = (cle[1].split("_")[1], int(cle[1].split("_")[2]))
#                 for cluster in clusters :
#                     if noeud1 in cluster and noeud2 in cluster :
#                         meme_cluster = True
#                         if rmsd[cle] != None and (noeud1, noeud2) not in liste_paires :
#                                     #print((noeud1, noeud2)) 
#                             liste_val_rmsd_par_cluster.append(rmsd[cle])
#                             liste_paires.append((noeud1, noeud2))
#                             
#                 if not meme_cluster and noeud1 in liste_arn and noeud2 in liste_arn: 
#                     if rmsd[cle] != None and (noeud1, noeud2) not in liste_paires :
#                                     #print((noeud1, noeud2))
#                             liste_val_rmsd_inter_cluster.append(rmsd[cle])
#                             liste_paires.append((noeud1, noeud2)) 

                
            for cluster in clusters :
                #print(cluster)
                #del(liste_val_rmsd_par_cluster[:])
                if len(cluster) > 1 :
                    for i in range(len(cluster)) :
                        for j in range(i+1, len(cluster)) :
                            noeud1 = "fichier_%s_%d_taille_4.pdb"%(cluster[i][0], cluster[i][1])
                            noeud2 = "fichier_%s_%d_taille_4.pdb"%(cluster[j][0], cluster[j][1])
                             
                            if (noeud1, noeud2) in rmsd.keys() :
                                #print("petit rat")
                                if rmsd[(noeud1, noeud2)] != None and (noeud1, noeud2) not in liste_paires :
                                    #print((noeud1, noeud2))
                                    #if rmsd[(noeud1, noeud2)] < 40.0 :
                                        #print(noeud1, noeud2) 
                                        liste_val_rmsd_par_cluster.append(rmsd[(noeud1, noeud2)])
                                        liste_paires.append((noeud1, noeud2))
                                else :
                                    print((noeud2, noeud1))
                            elif (noeud2, noeud1) in rmsd.keys() :
                                #print("petit rat")
                                if rmsd[(noeud2, noeud1)] != None and (noeud1, noeud2) not in liste_paires :
                                    #print((noeud1, noeud2))
                                    #if rmsd[(noeud2, noeud1)] < 40.0 :
                                        liste_val_rmsd_par_cluster.append(rmsd[(noeud2, noeud1)])
                                        liste_paires.append((noeud2, noeud1))
                                else :
                                    print((noeud2, noeud1))
                                     
            print(len(liste_val_rmsd_par_cluster))
            #print(liste_val_rmsd_par_cluster)
                 
              
            liste_val_rmsd_inter_cluster = []
            for i in range(len(clusters)) :
                for j in range(i+1, len(clusters)) :
                    for elt1 in clusters[i] :
                        for elt2 in clusters[j] :
                            noeud1 = "fichier_%s_%d_taille_4.pdb"%(elt1[0], elt1[1])
                            noeud2 = "fichier_%s_%d_taille_4.pdb"%(elt2[0], elt2[1])
                         
                            if (noeud1, noeud2) in rmsd.keys() and (noeud1, noeud2) not in liste_paires :
                                #print("petit rat")
                                if rmsd[(noeud1, noeud2)] != None :
                                    #if rmsd[(noeud1, noeud2)] < 40.0 :
                                    #print((noeud1, noeud2))
                                        liste_val_rmsd_inter_cluster.append(rmsd[(noeud1, noeud2)])
                                    #else :
                                    #    print((noeud1, noeud2))
                                else :
                                    print((noeud2, noeud1))
                            elif (noeud2, noeud1) in rmsd.keys() and (noeud2, noeud1) not in liste_paires  :
                                #print("petit rat")
                                if rmsd[(noeud2, noeud1)] != None :
                                    #print((noeud1, noeud2))
                                    #if rmsd[(noeud2, noeud1)] < 40.0 :
                                        liste_val_rmsd_inter_cluster.append(rmsd[(noeud2, noeud1)])
                                    #else :
                                    #    print((noeud2, noeud1))
                                else :
                                    print((noeud2, noeud1))
            print(len(liste_val_rmsd_inter_cluster))
            #print(liste_val_rmsd_inter_cluster)    
            
            figure = plt.figure(figsize = (6,7))
            axs1 = sns.distplot(liste_val_rmsd_par_cluster, kde=False, norm_hist = True, label="intra-cluster")
            axs1.set_xlabel("RMSD C3'")
            axs1.set_ylabel("Pair proportion")
            
            
            axs2 = sns.distplot(liste_val_rmsd_inter_cluster, kde=False, norm_hist = True, label="inter-cluster")
            plt.title("Distribution of RMSD values within clusters and between clusters \n total number of intra-cluster pairs : %d \n total number of inter-cluster pairs : %d"%(len(liste_val_rmsd_par_cluster), len(liste_val_rmsd_inter_cluster)))
            #plt.title("Distribution des RMSD pour tous les types \n nombre total de paires intracluster : %d"%( len(liste_val_rmsd_par_cluster)))

            plt.legend()
            plt.savefig("/media/coline/Maxtor/Distribution_rmsd_taille_4/distrib_rmsd_tout_taile_4_carbone_3'_new.svg", format='svg', transparent=True)
            
            plt.show()
            
            compteur_none = 0
            for cle in rmsd.keys() :
                if rmsd[cle] == None :
                    compteur_none += 1
                    
            print(compteur_none)
            print(len(rmsd))
            #print(min(liste_val_rmsd_par_cluster))
                
#                 if len(liste_val_rmsd_par_cluster) > 0 :
#                     sns.distplot(liste_val_rmsd_par_cluster, kde = False)
#                     plt.show()


if __name__ == '__main__':
    types_arn = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16s_arnm", "arnt_16s"]
#     groupe = CLUSTERING_PEREZ_VERSION_NON_CAN_2[12]
#     #groupe = ['5FDU_1A_30_17', '1U9S_A_58_11']
#       
#     #['5FDU_1A_48_19', '5J7L_DA_30_15', '4V9F_0_30_4', '1FJG_A_48_8', '5J7L_DA_48_15', '5J5B_BA_48_23', '5FDU_1A_30_17', '1U9S_A_58_11']

    ''' Version 90 molécules '''
#     with open("graphs_2.92.pickle", 'rb') as fichier_tout :
#         i = 5
#         mon_depickler_graphes = pickle.Unpickler(fichier_tout)
#         graphes = mon_depickler_graphes.load()    
#         compteur = 0
#         tab_fichiers = []
#         tab_fichier_positions = []
#         tab_motif_positions = []
#         tab_nb_nts_par_chaine = []
#         for fichier in os.listdir(EXTENSION_PATH_TAILLE%i) :
#             if  ".pickle" in fichier and len(fichier.split("_")) == 5 and "couples_possibles" not in fichier :
#                           
#                 with open(EXTENSION_PATH_TAILLE%i +fichier,'rb') as fichier_graphe :
#                     mon_depickler = pickle.Unpickler(fichier_graphe)
#                     graphe = mon_depickler.load()
#                     print(fichier)
#                     motif, tab_positions, nb_nts_par_chaines = recup_motif_et_chaines(graphe, i, graphes[(fichier.split("_")[1], fichier.split("_")[2])].number_of_nodes())
#                     print(tab_positions)
#                     cle  = fichier.split("_")[1] + "_" + fichier.split("_")[2] + "_" + fichier.split("_")[3] + "_" + fichier.split("_")[4][:len(fichier.split("_")[4])-7]
#                     print(cle)
#                     if fichier != "fichier_2XD0_V_36_21.pickle" : #and cle in groupe :
#                         tab_fichier_positions.append(tab_positions)
#                         tab_motif_positions.append(motif)
#                         tab_nb_nts_par_chaine.append(nb_nts_par_chaines)
#                         print(tab_fichier_positions)
#                         #print(graphes[(fichier.split("_")[1], fichier.split("_")[2])].nodes.data())
#                         if fichier[:len(fichier)-7]+"_taille_%d.pdb"%i not in os.listdir(PATH_MMCIF) :
#                             res = extract_pdb(tab_positions, fichier.split("_")[1] + "_"+fichier.split("_")[2], graphes)
#                             print(res)
#                                            
#                             with open(PATH_MMCIF+fichier[:len(fichier)-7]+"_taille_%d.pdb"%i, 'w') as fichier_pdb:
#                                 fichier_pdb.write(res)
#                         tab_fichiers.append(fichier[:len(fichier)-7]+"_taille_%d.pdb"%i)
#                               
#                     compteur += 1
#     
#         print(compteur)
#         print(tab_fichiers)
#         print(tab_fichier_positions)
#         rmsd = rmsd(tab_fichiers, tab_fichier_positions, tab_motif_positions, tab_nb_nts_par_chaine, graphes)
#         print(rmsd)
#                     
#         with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_1.pickle"%i, 'wb') as fichier_rmsd :
#             mon_pickler = pickle.Pickler(fichier_rmsd)
#             mon_pickler.dump(rmsd)  
            
            
    ''' Version new data '''
          
    tab_fichiers = []
    tab_fichier_positions = []
    tab_motif_positions = []
    tab_nb_nts_par_chaine = []
              
    with open("resolutions.pickle", 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        resolutions = mon_depickler.load()
                  
        proportion_inf_3 = 0
        total = 0
                  
        liste_a_garder = []
                  
        for typ in types_arn :
            with open("Nouvelles_donnees/liste_representant_%s.pickle"%typ, 'rb') as fichier_repr :
                mon_depickler = pickle.Unpickler(fichier_repr)
                liste_representant = mon_depickler.load() 
                          
                #print(resolutions)
                for elt in liste_representant :
                    if resolutions[elt[0]] <= 3.0 :
                        liste_a_garder.append(elt)
          
    with open("occ_multi_chaine.pickle", 'rb') as fichier_multi_chaines :
        mon_depickler = pickle.Unpickler(fichier_multi_chaines)
        liste_plusieurs_chaines = mon_depickler.load() 
#      
#      
    for i in range(4,5) :
        for elt in liste_a_garder :
            if elt == ('1vq8', 11) :# or elt == ('4w2f', 40) :
                with open(NEW_EXTENSION_PATH_TAILLE + "fichier_%s_%s.pickle"%(elt[0],elt[1]),'rb') as fichier_graphe :
                        mon_depickler = pickle.Unpickler(fichier_graphe)
                        graphe = mon_depickler.load() 
                                  
                        with open("Graphs/"+elt[0]+".pickle", 'rb') as fichier_graphes :
                            mon_depickler_g  = pickle.Unpickler(fichier_graphes)
                            graphes = mon_depickler_g.load()
                                  
                            motif, tab_positions, nb_nts_par_chaines = recup_motif_et_chaines_new_data(graphe, i, graphes.number_of_nodes(), graphes)
                                  
                            print(len(tab_positions))
                                  
                            tab_fichier_positions.append(tab_positions)
                            tab_motif_positions.append(motif)
                            tab_nb_nts_par_chaine.append(nb_nts_par_chaines)
                            print(tab_positions)
                            print(elt)
                                #print(graphes[(fichier.split("_")[1], fichier.split("_")[2])].nodes.data())
                            if "fichier_"+elt[0]+"_"+str(elt[1])+"_taille_%d_test_.pdb"%i not in os.listdir("/media/coline/Maxtor/Fichiers_mmcif/"):
                                           
                                if elt in liste_plusieurs_chaines :
                                    chains = []
                                    for noeud, data in graphe.nodes(data=True) :
                                        if data["num_ch"] not in chains :
                                            chains.append(data["num_ch"])
                                else :            
                                    chains = [graphe.nodes[1]["num_ch"]]
                                res = extract_pdb_new_data(tab_positions, elt[0], chains , graphes)
                                print(res)
                                                              
                                with open("/media/coline/Maxtor/Fichiers_mmcif/fichier_%s_%s_taille_%d_test.pdb"%(elt[0], str(elt[1]),i), 'w') as fichier_pdb:
                                    fichier_pdb.write(res)
                                    print("gros rat")
                            tab_fichiers.append("fichier_"+elt[0]+"_"+str(elt[1])+"_taille_%d.pdb"%i)
# # # #              
        #rmsd = rmsd_new_data(tab_fichiers, tab_fichier_positions, tab_motif_positions, tab_nb_nts_par_chaine, graphes, i)
#         #print(rmsd)
# #                          
#         with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%i, 'wb') as fichier_rmsd :
#             mon_pickler = pickle.Pickler(fichier_rmsd)
#             mon_pickler.dump(rmsd) 
#          
#     with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_1_new_data.pickle"%4, 'rb') as fichier_rmsd :
#             mon_depickler = pickle.Unpickler(fichier_rmsd)
#             rmsd = mon_depickler.load()
#               
#             dico_none = {}
#             compter = 0
#             for elt in rmsd.keys() :
#                 if rmsd[elt] == None :
#                     print(elt)
#                     dico_none.update({elt : rmsd[elt]})
#                     compter += 1
#                 #print(rmsd[elt])
#               
#             print(compter)
                
            #sns.distplot([x for x in rmsd.values() if x != None], kde=False, norm_hist=True)
            #plt.show()
    #for typ in types_arn :          
    #distrib_rmsd_new_data_par_cluster(types_arn)
            
            
    
            
            
#         print(graphes[('1U9S', 'A')].nodes.data())
        #print(graphes[('5J7L', 'DA')].nodes.data())
            
            
          