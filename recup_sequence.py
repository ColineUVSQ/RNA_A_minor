'''
Created on 7 mai 2019

@author: coline
'''

import pickle
import os
from recup_data.constantes import EXTENSION_PATH_TAILLE,\
    GROUPES_TOUTES_ARETES_MAX_4_10_07, EXTENSION_PATH, GROUPE_ARICH_DE_07

def recup_sequence(graphe_extension, graphe, taille_extension, num_chaine):
    sequence = ""
    position_depart = graphe_extension.nodes[num_chaine]["position"][0]
    
    compteur = position_depart
    for i in range(taille_extension) :
        if compteur > 0 and compteur < graphe.number_of_nodes() :
            sequence += graphe.nodes[compteur]["nt"]
            if num_chaine == 1 or num_chaine == 4 :
                compteur += 1
            else :
                compteur -= 1
    
    print(sequence)    
    return sequence


    
    

def recup_sequence_tige_groupe3(graphe_extension, graphe):   
    sequence1 = ""
    sequence2 = ""
    
    position_depart = graphe_extension.nodes[4]["position"][0]
    
    compteur = position_depart + 1
    for _ in range(3) :
        if compteur > 0 and compteur < graphe.number_of_nodes() :
            sequence1 += graphe.nodes[compteur]["nt"]
            for voisin in graphe[compteur] :
                if graphe.edges[compteur, voisin]["label"] == 'CWW' :
                    sequence2 += graphe.nodes[voisin]["nt"]
            compteur += 1
    
#     sequence2 =  ""   
#     position_depart = graphe_extension.nodes[3]["position"][0] - 3
#     compteur = position_depart
#     for _ in range(3) :
#         if compteur > 0 and compteur < graphe.number_of_nodes() :
#             sequence2 += graphe.nodes[compteur]["nt"]
#             compteur -= 1    
    
    print(sequence1)
    print(sequence2)
    return sequence1, sequence2
            
if __name__ == '__main__':
    with open("graphs_2.92.pickle", 'rb') as fichier_tout :
        mon_depickler_graphes = pickle.Unpickler(fichier_tout)
        graphes = mon_depickler_graphes.load()
        
        sequences = []
        for fichier in os.listdir(EXTENSION_PATH_TAILLE%8) :
            if "_2.pickle" in fichier and "couples_possibles" not in fichier and len(fichier.split("_")) == 6  :
                #print(fichier[8:len(fichier)-9])
                if  fichier[8:len(fichier)-9] in GROUPES_TOUTES_ARETES_MAX_4_10_07[2] :
                    with open(EXTENSION_PATH_TAILLE%8 +fichier,'rb') as fichier_graphe :
                        mon_depickler = pickle.Unpickler(fichier_graphe)
                        graphe = mon_depickler.load()
                        
                        print(">"+fichier)
                        sequences.append(recup_sequence_tige_groupe3(graphe, graphes[(fichier.split("_")[1], fichier.split("_")[2])]))
        
#         with open(EXTENSION_PATH%"taille_max"+"result_k_max_4_10_toutes_aretes/fichier_sequence_groupe_0.7_5.2.txt", 'a') as fichier_ecriture :
#             fichier_ecriture.write("chaine 4 : \n")
#             fichier_ecriture.write("sequence = c(")            
#             for seq in sequences :
#                 fichier_ecriture.write(', "'+seq+'" ')
#             fichier_ecriture.write(")\n")