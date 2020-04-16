'''
Created on 19 nov. 2018

@author: coline
'''
import os
import pickle
from recup_data.calcul_sim import calcul_sim_avec_poids

# with open("graphs_2.92.pickle", 'rb') as fichier :
#     mon_depickler_graphes = pickle.Unpickler(fichier)
#     graphes = mon_depickler_graphes.load()
#     compteur = 0
#     with open("fichier_pas_liaison_cov_2.txt", 'w') as fichier_sortie :
#         for element in os.listdir('graphes_extension/'):
#             if "pickle" in element : 
#                 with open("graphes_extension/"+element, 'rb') as fichier_entree :
#                     mon_depickler = pickle.Unpickler(fichier_entree)
#                     G_ext = mon_depickler.load()
#                     
#                     occ = (element.split("_")[1], element.split("_")[2])
#                     
#                     for noeud_ext in G_ext.nodes() :
#                         liaison_cov = False
#                         elt = G_ext.nodes[noeud_ext]["position"][0]
#                         while elt <= G_ext.nodes[noeud_ext]["position"][1] :
#                             #for voisin in graphes[occ][elt] :
#                             #if graphes[occ][elt][elt+1]["label"] == 'B53' :
#                             #    liaison_cov = True
#                             if elt+1 < graphes[occ].number_of_nodes() :
#                                 if graphes[occ][elt][elt+1]["label"] != 'B53' and G_ext.nodes[noeud_ext]["type"] not in [None,11,12,13,14] :
#                                     print(element)
#                                     print(elt)
#                                     fichier_sortie.write(str(occ) + " num RIN : " + str(int(element.split('_')[3]) - 1)+ " num occ : " + str(int(element.split('_')[4][:len(element.split('_')[4]) - 7])-1) + " num nucleotide " + str(elt) + " " + str(graphes[occ][elt][elt+1]["label"]) + "\n")
#                                     compteur = compteur+1
#                             elt = elt + 1
#         
#     print(compteur)
   
   
# with open("dico_comp_complet_nouvelle_metrique.pickle", 'rb') as fichier :
#     mon_depickler_graphes = pickle.Unpickler(fichier)
#     graphes = mon_depickler_graphes.load()
#     print(graphes.keys())
#     
#     for cle in graphes.keys() :
#         if (cle[0] == 'fichier_3UCZ_R_62_15' and cle[1] == 'fichier_1U9S_A_58_11') or (cle[1] == 'fichier_3UCZ_R_62_15' and cle[0] == 'fichier_1U9S_A_58_11'):   
#             print(graphes[cle].edges.data())    
#             
#     with open("result_graphes_comp_serveur_metrique_sommets_aretes/graphe_comp_couples_possibles_fichier_1U9S_A_58_11_fichier_3UCZ_R_62_15.pickle", 'rb') as fichier_2 :
#         mon_depickler_graphe = pickle.Unpickler(fichier_2)
#         graphe = mon_depickler_graphe.load()    
#         
#         for cle in graphe.keys() :
#             print(graphe[cle].edges.data()) 
#             
#             
#             with open("graphes_extension/fichier_1U9S_A_58_11.pickle", 'rb') as fichier_graphe1 :
#                 mon_depickler_graphe1 = pickle.Unpickler(fichier_graphe1)
#                 graphe1 = mon_depickler_graphe1.load()
#                 
#                 with open("graphes_extension/fichier_3UCZ_R_62_15.pickle", 'rb') as fichier_graphe2 :
#                     mon_depickler_graphe2 = pickle.Unpickler(fichier_graphe2)
#                     graphe2 = mon_depickler_graphe2.load()
#             
#                     sim = calcul_sim_avec_poids(graphe1, graphe2, graphe[cle], cle)
#                     print(sim)

with open("grands_graphes.pickle", 'rb') as fichier_ecriture :
    mon_pickler = pickle.Unpickler(fichier_ecriture)
    dico_graphes = mon_pickler.load()
    
    pas_bon = []
    for cle in dico_graphes.keys() :
#         print(cle)
#         print(dico_graphes[cle].nodes.data())
#         print(dico_graphes[cle].edges.data())
        cww_non_can = False
        compteur = 0
        for u,v,data in dico_graphes[cle].edges(data=True) :
            if data["label"] == 'CWW' :
                if (dico_graphes[cle].nodes[u]["nt"] == 'A' and dico_graphes[cle].nodes[v]["nt"] == 'U') or (dico_graphes[cle].nodes[u]["nt"] == 'U' and dico_graphes[cle].nodes[v]["nt"] == 'A') or (dico_graphes[cle].nodes[u]["nt"] == 'C' and dico_graphes[cle].nodes[v]["nt"] == 'G') or (dico_graphes[cle].nodes[u]["nt"] == 'G' and dico_graphes[cle].nodes[v]["nt"] == 'C') or (dico_graphes[cle].nodes[u]["nt"] == 'G' and dico_graphes[cle].nodes[v]["nt"] == 'U') or (dico_graphes[cle].nodes[u]["nt"] == 'U' and dico_graphes[cle].nodes[v]["nt"] == 'G') :
                    compteur += 1
                else :
                    print(cle)
                    print(dico_graphes[cle].nodes[u]["nt"])
                    print(dico_graphes[cle].nodes[v]["nt"])
                    print("non_ok")
                    cww_non_can = True
        if cww_non_can :
            nom = "fichier_"+ cle[0] + "_" + cle[1] + "_" + str(cle[2]) + "_" + str(cle[3]) + ".pickle"
            pas_bon.append(nom)
            if nom == "fichier_1GID_B_25_20.pickle" :
                print("ramou")
    print(len(pas_bon))
    print(pas_bon)
    print(len(dico_graphes))
            
                     
    