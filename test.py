'''
Created on 11 oct. 2018

@author: Coline Gi
'''
import pickle
import os
import networkx as nx
import numpy as np
import math
import matplotlib.pyplot as plt
from networkx.algorithms import isomorphism
from recup_data.constantes import EXTENSION_PATH, EXTENSION_PATH_TAILLE,\
    CLUSTERING_PEREZ_VERSION_NON_CAN_2, PATH_MMCIF, GROUPE_GNRA,\
    GROUPE_GNRA_ETENDU, HOMOLOGUES, NEW_EXTENSION_PATH_TAILLE
from recup_data import calcul_sim
from recup_data.calcul_sim import calcul_aretes_communes_avec_coeff,\
    calcul_sim_aretes_avec_coeff
from recup_data.draw_like_carnaval import rang_sim
import Bio.PDB
from Bio.PDB import MMCIF2Dict
from Bio.PDB import MMCIFParser
from Bio import SeqIO
from gemmi import cif
import tabnanny
import time
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import csv
from recup_data.etude_composantes_connexes import construction_matrice_distance
from sklearn.cluster import KMeans
from recup_data.new_algo_comparaison import comparaison, test_compatibilite

# with open("graphs_2.92.pickle", 'rb') as fichier_tout :
#         mon_depickler_graphes = pickle.Unpickler(fichier_tout)
#         graphes = mon_depickler_graphes.load()
#         with open("fichier_struct_second_2XD0_V.txt","w") as fichier :
#             fichier.write(str(graphes[('2XD0', 'V')].nodes.data()))
#             fichier.write(str(graphes[('2XD0', 'V')].edges.data()))
#             print("ramousnif")
#             print(graphes[('2XD0', 'V')].nodes.data())
#             print(graphes[('2XD0', 'V')].edges.data())
            
#         with open("graphes_extension/fichier_1FJG_A_48_8.pickle", 'rb') as fichier_extension :
#             mon_depickler_extension = pickle.Unpickler(fichier_extension)
#             extension = mon_depickler_extension.load()
#             print(extension.nodes.data())

# with open("fichier_max.pickle", 'rb') as fichier :
#     mon_depickler = pickle.Unpickler(fichier)
#     tab = mon_depickler.load()
#     print(len(tab))
#     print(tab)


# for fic in os.listdir("graphes_extension/fichiers_couples/") :
#         element1 = fic.split('_')[2] + '_' + fic.split('_')[3] + '_' + fic.split('_')[4] + '_' + fic.split('_')[5] + '_' + fic.split('_')[6]
#         element2 = fic.split('_')[7] + '_' + fic.split('_')[8] + '_' + fic.split('_')[9] + '_' + fic.split('_')[10] + '_' + fic.split('_')[11][:len(fic.split('_')[11])-7]
#         #element1 = "fichier_1FJG_A_48_11"
#         #element2 = "fichier_1FJG_A_48_8"
#         #fic = "couples_possibles_fichier_1FJG_A_48_11_fichier_1FJG_A_48_8.pickle"
#         with open("graphes_extension/"+element1+".pickle", 'rb') as fichier1 :
#                 mon_depickler1 = pickle.Unpickler(fichier1)
#                 graphe1 = mon_depickler1.load()     
#                 with open("graphes_extension/"+element2+".pickle", 'rb') as fichier2 :
#                     mon_depickler2 = pickle.Unpickler(fichier2)
#                     graphe2 = mon_depickler2.load()
#      
#                     with open("graphes_extension/fichiers_couples/" + fic, 'rb') as fichier_pickle :
#                                     memory_error = False
#                                     graphe_motif = nx.MultiGraph()
#                                     for i in range(1,6) :
#                                         graphe_motif.add_node((i,i))
#                                     graphe_motif.add_edge((1,1),(2,2), type="NON_CAN")
#                                     graphe_motif.add_edge((1,1),(3,3), type="COV")
#                                     graphe_motif.add_edge((1,1),(5,5), type="NON_CAN")
#                                     graphe_motif.add_edge((2,2),(4,4), type="COV")
#                                     graphe_motif.add_edge((2,2),(5,5), type="CAN")
#                                     graphe_motif.add_edge((3,3),(4,4), type="NON_CAN")
#                                     
#                                                             
#                                     mon_depickler = pickle.Unpickler(fichier_pickle)
#                                     couples_possibles = mon_depickler.load()
#                                     print(fic)
#                                     #print(couples_possibles)
#                                     
#                                     new_couples_possibles = []
#                                     for i in range(4) :
#                                         new_couples = []
#                                         if 'memory error' in couples_possibles[i] :
#                                             memory_error = True
#                                             break
#                                         for chaine in couples_possibles[i] :
#                                             new_chaine = [(i+1, i+1)]
#                                             for couple in chaine :
#                                                 new_chaine.append(couple)
#                                             new_couples.append(new_chaine)
#                                         new_couples_possibles.append(new_couples)
#                                     print(new_couples_possibles)

# with open("tab_nexiste_pas.pickle", 'rb') as fichier :
#     mon_depickler = pickle.Unpickler(fichier)
#     tab = mon_depickler.load()
#     print(tab)
#     
#     for elt in tab :
#         if '3JCS_1_282_1' in elt[0] or '3JCS_1_282_1' in elt[1] :
#             print(elt)
            
    
# with open("graphes_extension/fichiers_couples_qui_manquent/couples_possibles_fichier_4PRF_B_25_69_fichier_5J7L_DA_272_2.pickle", 'rb') as fichier :
#     mon_depickler = pickle.Unpickler(fichier)
#     couples_possibles = mon_depickler.load()
#     print(couples_possibles[3])       
            

# with open("fichier_comp_grands_graphes_test.pickle", 'rb') as fichier:
#     mon_depickler_comp = pickle.Unpickler(fichier)
#     dico_comp_4v9f = mon_depickler_comp.load()
#     
#     with open("fichier_comp_grands_graphes_V2.pickle", 'rb') as fichier_comp:
#         mon_depickler_tout = pickle.Unpickler(fichier_comp)
#         dico_comp = mon_depickler_tout.load()
#         
#         print(len(dico_comp))
#         print(len(dico_comp_4v9f))
#         print(dico_comp_4v9f.keys())
#         
# #         for elt in dico_comp.keys() :
# #             if elt[0] == ('4V9F', '0', 25, 4) or elt[1] == ('4V9F', '0', 25, 4) :
# #                 print("ramousnif")
#         
#         compteur = 0
#         for elt in dico_comp_4v9f.keys() :
#             if elt[0] == ('4V9F', '0', 25, 4) or elt[1] == ('4V9F', '0', 25, 4) :
#                 dico_comp.update({elt : dico_comp_4v9f[elt]})
#                 compteur += 1
#         print(compteur)
# 
#         print(len(dico_comp))          

# with open("grands_graphes.pickle", 'rb') as fichier:
#         mon_depickler_tout = pickle.Unpickler(fichier)
#         dico_graphe = mon_depickler_tout.load()
# #         for u,v, dict in dico_graphe[('1FJG', 'A', 294, 1)].edges(data=True) :
# #             print((u,v))
# #             print(dict)
#             
#         compteur_arc = 0
#         compteur = 0
#         for (u, v, keys, t) in dico_graphe[('1FJG', 'A', 294, 1)].edges(data="label", keys = True) :
#             if (v,u) in dico_graphe[('1FJG', 'A', 294, 1)].edges() :
#                 print("ramou")
#                 compteur += 1
#         print(compteur)
#                 
#         print(compteur_arc)
#         print(dico_graphe[('1FJG', 'A', 294, 1)].number_of_edges())
    
#     
# with open("graphs_2.92.pickle", 'rb') as fichier_tout :
#     mon_depickler_graphes = pickle.Unpickler(fichier_tout)
#     graphes = mon_depickler_graphes.load() 
#       
#     with open("fichier_struct_second_4L81_A.txt", 'w') as fichier :
#           
#         for noeud, data in graphes[('4L81', 'A')].nodes(data=True) :
#             fichier.write(str(noeud) + " "+ str(data)+ "\n")
#           
#         for u,v,data in graphes[('4L81', 'A')].edges(data=True) :
#             fichier.write(str((u,v)) + " "+ str(data) + "\n")
# 
# pas_bon = []         
# for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension") :
#     if "pickle" in fic :
#         for fic_2 in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_test") :
#             if fic == fic_2 :
#                 with open("Extensions/Metrique_toutes_aretes/graphes_extension/"+fic, 'rb') as fichier_1 :
#                     mon_depickler = pickle.Unpickler(fichier_1)
#                     graphe1 = mon_depickler.load()
#                        
#                     with open("Extensions/Metrique_toutes_aretes/graphes_extension_test/"+fic_2, 'rb') as fichier_2 :
#                         mon_depickler_2 = pickle.Unpickler(fichier_2)
#                         graphe2 = mon_depickler_2.load()
#                            
#                            
#                         print(fic)
#                         for noeud, data in graphe1.nodes(data=True) :
#                             if noeud not in graphe2.nodes() and data["type"] != -1 :
#                                 print("pb1")
#                                 print(noeud)
#                                 if fic not in pas_bon :
#                                     pas_bon.append(fic)
#                            
#                         for noeud, data in graphe2.nodes(data=True) :
#                             if noeud not in graphe1.nodes() :
#                                 print("pb2")
#                                 print(noeud)
#                                 if fic not in pas_bon :
#                                     pas_bon.append(fic)            
#                            
#                                
#                         for u,v, data in graphe1.edges(data=True) :
#                             if (u,v) not in graphe2.edges() and data["label"] != '0' :
#                                 print("pb3")
#                                 print((u,v))
#                                 if fic not in pas_bon :
#                                     pas_bon.append(fic) 
#                                     
#                             elif (u,v) in graphe2.edges() :
#                                 compteur_bon = 0
#                                 for edge_1 in graphe1[u][v] :
#                                     for edge_2 in graphe2[u][v] :
#                                         if graphe1[u][v][edge_1]["label"] == graphe2[u][v][edge_2]["label"] :
#                                             compteur_bon += 1
#                                 if compteur_bon != len(graphe1[u][v]) :
#                                     print("pb5")
#                                     print((u,v))
#                                     pas_bon.append(fic)
#                                        
#                         for u,v, data in graphe2.edges(data=True) :
#                             if (u,v) not in graphe1.edges() :
#                                 print("pb4")
#                                 print((u,v))
#                                 if fic not in pas_bon :
#                                     pas_bon.append(fic) 
#                             
#                             elif (u,v) in graphe1.edges() :
#                                 compteur_bon = 0
#                                 for edge_1 in graphe1[u][v] :
#                                     for edge_2 in graphe2[u][v] :
#                                         if graphe1[u][v][edge_1]["label"] == graphe2[u][v][edge_2]["label"] :
#                                             compteur_bon += 1
#                                 if compteur_bon != len(graphe2[u][v]) :
#                                     print("pb6")
#                                     print((u,v))
#                                     pas_bon.append(fic[:len(fic)-7])
#                        
# print(pas_bon) 
# print(len(pas_bon))


# with open("graphs_2.92.pickle", 'rb') as fichier_tout :
#     mon_depickler_graphes = pickle.Unpickler(fichier_tout)
#     graphes = mon_depickler_graphes.load() 
#     for fic in os.listdir("graphes_extension/graphes_extension_test") :
#         if "pickle" in fic :
#             with open("graphes_extension/graphes_extension_test/"+fic, 'rb') as fichier_1 :
#                 mon_depickler = pickle.Unpickler(fichier_1)
#                 graphe1 = mon_depickler.load()
#                 print(fic)
#                 nom_cle = (fic.split("_")[1], fic.split("_")[2])
#                 ## verification de pas de boucle dans le graphe
#                 for noeud in graphe1.nodes() : 
#                     if (noeud, noeud) in graphe1.edges() :
#                         print("pb")
#                  
#                 compte = [0,0,0,0,0]
#                 for noeud, data in graphe1.nodes(data=True) :
#                     if data["type"] != None and data["type"] != -1 :
#                         for elt in data["chaine"]:
#                             if 1 in data["position"] :
#                                 compte[elt-1] = 10
#                             elif graphes[nom_cle].number_of_nodes() in data["position"] :
#                                 compte[elt-1] = 10
#                             else :
#                                 compte[elt-1] += data["poids"]
#                          
#                 for i in range(4) :
#                     if compte[i] != 10 :
#                          
#                         print('pb')
#                         print(compte)
#                         print(graphe1.nodes.data())
                     


# with open("Extensions/Metrique_toutes_aretes/couples_possibles_bonne_version/couples_possibles_fichier_5FJC_A_138_1_fichier_1FJG_A_48_8.pickle", 'rb') as fichier :
#     mon_depickler = pickle.Unpickler(fichier)
#     couples = mon_depickler.load()
# 
#     print(couples[3])
    
# with open("Extensions/Metrique_raymond_non_cov/sim_extensions_longue_distance.pickle", 'rb') as fichier :
#     mon_depickler = pickle.Unpickler(fichier)
#     dico_sim = mon_depickler.load()
#     
#     
#     print(dico_sim.keys())
#     print(len(dico_sim))
#     print(dico_sim[('3JCS_1_25_16', '5DM6_X_328_2')])


# with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_4/fichier_1FJG_A_48_8.pickle", 'rb') as fichier :
#     mon_depickler = pickle.Unpickler(fichier)
#     graphe = mon_depickler.load()
#     print(graphe.nodes.data())
    
# with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_1/graphe_comp_test_couples_possibles_test_fichier_3JCS_1_48_18_fichier_3JCS_1_25_46.pickle", 'rb') as fichier :
#     mon_depickler = pickle.Unpickler(fichier)
#     graphe_comp = mon_depickler.load()
#     
#     for cle in graphe_comp.keys() :
#         print(graphe_comp[cle][(4,4)])
#     
# with open("graphs_2.92.pickle", 'rb') as fichier_tout :
#         mon_depickler_graphes = pickle.Unpickler(fichier_tout)
#         graphes = mon_depickler_graphes.load()
#         
#         
#         print(graphes[('3JCS', '1')][865])


## Verif difference avec version 3 (13 juillet 2019)
# for j in range(3, 11) :
#     pas_bon = []
#     for fic in os.listdir(EXTENSION_PATH_TAILLE%j) :
#             #if "graphe_comp_test" in fic :
#             if "couples_possibles" not in fic and ".pickle" in fic and len(fic.split("_")) == 5  :
#                 #fic_2 = fic[:11]+"_"+fic[17:34]+"_"+fic[40:]
#                 fic_2 = fic[:len(fic)-7] + "_3.pickle"
#                 #print(fic)
#                 #print(fic_2)
#                 if fic_2 in os.listdir(EXTENSION_PATH_TAILLE%j) :
#                     with open(EXTENSION_PATH_TAILLE%j+fic, 'rb') as fichier_1 :
#                             mon_depickler = pickle.Unpickler(fichier_1)
#                             graphe1 = mon_depickler.load()
#                                        
#                                            
#                             with open(EXTENSION_PATH_TAILLE%j+fic_2, 'rb') as fichier_2 :
#                                     mon_depickler_2 = pickle.Unpickler(fichier_2)
#                                     graphe2 = mon_depickler_2.load()
#                                                
#                                 #for cle in graphe1.keys() :     
#                                     #print(fic)
#                                     #print(cle)
#                                     for noeud, data in graphe1.nodes(data=True) :
#                                         if noeud not in graphe2.nodes() :
#                                             print("pb1")
#                                             print(noeud)
#                                             if fic[:len(fic)-7] not in pas_bon :
#                                                 pas_bon.append(fic[:len(fic)-7])
#                                         else :
#                                             if data["type"] != graphe2.nodes[noeud]["type"] :
#                                                 if fic[:len(fic)-7] not in pas_bon :
#                                                     pas_bon.append(fic[:len(fic)-7])
#                                                    
#                                     for noeud, data in graphe2.nodes(data=True) :
#                                         if noeud not in graphe1.nodes() :
#                                             print("pb2")
#                                             print(noeud)
#                                             if fic[:len(fic)-7] not in pas_bon :
#                                                 pas_bon.append(fic[:len(fic)-7]) 
#                                         else :
#                                             if data["type"] != graphe2.nodes[noeud]["type"] :
#                                                 if fic[:len(fic)-7] not in pas_bon :
#                                                     pas_bon.append(fic[:len(fic)-7])           
#                                                    
#                                                        
#                                     for u,v, data in graphe1.edges(data=True) :
#                                         if (u,v) not in graphe2.edges() :
#                                             print("pb3")
#                                             print((u,v))
#                                             if fic[:len(fic)-7] not in pas_bon :
#                                                          
#                                                 pas_bon.append(fic[:len(fic)-7]) 
#                                                             
#                                         elif (u,v) in graphe2.edges() :
#                                             compteur_bon = 0
#                                             for edge_1 in graphe1[u][v] :
#                                                 for edge_2 in graphe2[u][v] :
#                                                     if graphe1[u][v][edge_1]["label"] == graphe2[u][v][edge_2]["label"] or ( len(graphe1[u][v][edge_1]["label"]) == 3 and  len(graphe2[u][v][edge_2]["label"]) == 3 and graphe1[u][v][edge_1]["label"][1] == graphe2[u][v][edge_2]["label"][2] and graphe1[u][v][edge_1]["label"][2] == graphe2[u][v][edge_2]["label"][1])    :
#                                                         compteur_bon += 1
#                                             if compteur_bon != len(graphe1[u][v]) :
#                                                 print("pb5")
#                                                 print((u,v))
#                                                 print(graphe2[u][v][edge_2]["label"])
#                                                 if fic[:len(fic)-7] not in pas_bon :
#                                                     pas_bon.append(fic[:len(fic)-7])
#                                                                
#                                     for u,v, data in graphe2.edges(data=True) :
#                                         if (u,v) not in graphe1.edges() :
#                                             print("pb4")
#                                             print((u,v))
#                                             if fic[:len(fic)-7] not in pas_bon :
#                                                 pas_bon.append(fic[:len(fic)-7]) 
#                                                     
#                                         elif (u,v) in graphe1.edges() :
#                                             compteur_bon = 0
#                                             for edge_1 in graphe1[u][v] :
#                                                 for edge_2 in graphe2[u][v] :
#                 
#                                                     if graphe1[u][v][edge_1]["label"] == graphe2[u][v][edge_2]["label"] or ( len(graphe1[u][v][edge_1]["label"]) >= 3 and  len(graphe2[u][v][edge_2]["label"]) >= 3 and graphe1[u][v][edge_1]["label"][1] == graphe2[u][v][edge_2]["label"][2] and graphe1[u][v][edge_1]["label"][2] == graphe2[u][v][edge_2]["label"][1])    :
#                                                         compteur_bon += 1
#                                             if compteur_bon != len(graphe2[u][v]) :
#                                                 print("pb6")
#                                                 print((u,v))
#                                                 print(graphe2[u][v][edge_2]["label"])
#                                                 if fic[:len(fic)-7] not in pas_bon :
#                                                     pas_bon.append(fic[:len(fic)-7])          
#     print(j)
#     print(pas_bon) 
#     print(len(pas_bon))
    
#    
# compteur = 0
# with open("script_enlever_a_changer.sh", 'w') as fichier :
#     for fic in os.listdir("/home/coline/Bureau/graphes_extension_autres_tailles_new/graphes_extension_autres_tailles_new/graphes_extension_taille_%d"%j) :
#         if "pickle" in fic and "couples_possibles" in fic :
#             for elt in pas_bon :
#                 if elt in fic :
#                     fichier.write("rm /home/coline/Bureau/graphes_extension_autres_tailles_new/graphes_extension_autres_tailles_new/graphes_extension_taille_%d/%s\n"%(j, fic))
#                     compteur += 1
# print(compteur)

# compteur = 0
# with open("script_renommer.sh", 'w') as fichier :
#     for fic in os.listdir("/home/coline/Bureau/graphes_extension_autres_tailles_new/graphes_extension_autres_tailles_new/graphes_extension_taille_%d"%j) :
#         if "pickle" in fic and "couples_possibles" in fic :
#             if "graphe_comp" not in fic :
#                 #fic_avec_3 = fic.split("_")[0] + "_" + fic.split("_")[1] + "_" + fic.split("_")[2] + "_" + fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6] + "_3_" +fic.split("_")[8] + "_" +fic.split("_")[9] + "_" + fic.split("_")[10] + "_" + fic.split("_")[11] + "_" + fic.split("_")[12] + "_3.pickle"
#                 fic_avec_3 = fic.split("_")[0] + "_" + fic.split("_")[1] + "_" + fic.split("_")[2] + "_" + fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6] + "_3_" +fic.split("_")[7] + "_" +fic.split("_")[8] + "_" + fic.split("_")[9] + "_" + fic.split("_")[10] + "_" + fic.split("_")[11][:len(fic.split("_")[11])-7] +  "_3.pickle"
#                 fichier.write("mv /home/coline/Bureau/graphes_extension_autres_tailles_new/graphes_extension_autres_tailles_new/graphes_extension_taille_%d/%s /home/coline/Bureau/graphes_extension_autres_tailles_new/graphes_extension_autres_tailles_new/graphes_extension_taille_%d/%s\n"%(j, fic, j, fic_avec_3))
#                 print(fic_avec_3)
#                 compteur += 1
# print(compteur)


# pas_bon_avec_2 = []
# 
# for elt in pas_bon :
#     pas_bon_avec_2.append(elt[:len(elt)-7]+"_2.pickle")
# 
# print(pas_bon_avec_2)
# autre_pas_bon = ['fichier_4V88_A6_17_55.pickle', 'fichier_5J7L_DA_25_10.pickle', 'fichier_5DM6_X_25_15.pickle', 'fichier_1GID_B_25_20.pickle', 'fichier_5DM6_X_25_34.pickle', 'fichier_3JCS_1_25_46.pickle', 'fichier_4V88_A5_25_47.pickle', 'fichier_4PRF_B_25_69.pickle', 'fichier_1U9S_A_25_74.pickle', 'fichier_4L81_A_25_77.pickle', 'fichier_5FDU_1A_25_78.pickle', 'fichier_4V9F_0_30_4.pickle', 'fichier_5J7L_DA_30_15.pickle', 'fichier_5FDU_1A_30_17.pickle', 'fichier_5J7L_DA_30_36.pickle', 'fichier_5J7L_DA_36_3.pickle', 'fichier_2XD0_V_36_21.pickle', 'fichier_1FJG_A_48_8.pickle', 'fichier_5DM6_X_48_9.pickle', 'fichier_4V9F_0_48_13.pickle', 'fichier_5J5B_BA_48_14.pickle', 'fichier_4V9F_0_48_16.pickle', 'fichier_5J7L_DA_48_20.pickle', 'fichier_4V9F_0_48_21.pickle', 'fichier_5J5B_BA_48_23.pickle', 'fichier_5FDU_1A_48_25.pickle', 'fichier_4V9F_0_48_26.pickle', 'fichier_5DM6_X_48_28.pickle', 'fichier_5J5B_BA_58_3.pickle', 'fichier_1FJG_A_58_23.pickle', 'fichier_1FJG_A_109_6.pickle', 'fichier_4V9F_0_127_6.pickle', 'fichier_5J7L_DA_134_1.pickle', 'fichier_5DM6_X_134_2.pickle', 'fichier_5FDU_1A_134_3.pickle', 'fichier_5FDU_1A_197_3.pickle', 'fichier_5J7L_DA_197_4.pickle', 'fichier_4V9F_0_207_3.pickle', 'fichier_1FJG_A_271_1.pickle', 'fichier_5FDU_1A_272_1.pickle', 'fichier_5J7L_DA_272_2.pickle', 'fichier_4V88_A5_287_1.pickle', 'fichier_4V9F_0_287_2.pickle', 'fichier_5FDU_1A_290_2.pickle']
#  
# for elt in autre_pas_bon :
#     if elt not in pas_bon :
#         print(elt)

# compteur = 0
# with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_3/script_renommage.py", 'w') as fichier :
#     for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_3") :
#         if "couples_possibles" in fic and "test" not in fic and "graphe_comp" not in fic :
#             fic_2 = fic[:17]+"_test"+fic[17:]
#             #print(fic_2)
#             if fic_2 in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_3") :
#                 fichier.write("rm %s\n"%(fic))
#                 compteur += 1
#     print(compteur)
#     for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_3") :
#         if "couples_possibles" in fic and "test" in fic and "graphe_comp" not in fic :
#             fic_2 = fic[:17]+fic[22:]
#             print(fic_2)
#             fichier.write("mv %s %s\n"%(fic, fic_2))
#             compteur += 1
#     print(compteur)
#       
#     for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_3") :
#         if "graphe_comp" in fic and "test" not in fic :
#             fic_2 = fic[:11]+"_test"+fic[11:30]+ "test_" + fic[30:]
#             print(fic_2)
#             if fic_2 in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_3") :
#                 fichier.write("rm %s\n"%(fic))
#                 compteur += 1
#     print(compteur)
#     compteur = 0
#     for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_3") :
#         if "graphe_comp" in fic and "test" in fic :
#             fic_2 = fic[:11]+fic[16:35]+ fic[40:]
#             print(fic_2)
#             fichier.write("mv %s %s\n"%(fic, fic_2))
#             compteur += 1
#     print(compteur)





# 
# with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/motif_taille_4_0.7_4/digraphe_commun_0.7_4_taille_4.pickle", 'rb') as fichier :
#     mon_depickler = pickle.Unpickler(fichier)
#     crible = mon_depickler.load()
#      
#     print(crible.nodes.data())
#     print(crible.number_of_nodes())
#     print(crible.number_of_edges())
#  

# with open("fichier_ajout_dossiers.sh", 'w') as fichier :
#     for i in range(1,10) :
#         fichier.write("mkdir result_graphes_extension_%d/composantes_connexes_extensions_toutes_aretes_coeff_all1_taille_%s/\n"%(i,i))
#    


# with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/result_graphes_extension_6/sim_extensions_toutes_aretes_coeff_all1_par_k_taille_6.pickle", 'rb') as fichier :
#     mon_depickler = pickle.Unpickler(fichier)
#     dico_sim = mon_depickler.load()
#     print(dico_sim)



# G1 = nx.Graph()
# G2 = nx.Graph()
# G1.add_nodes_from([1,2,3,4])
# G2.add_nodes_from([1,2,3,4])
# G1.add_edges_from([(1,2), (3,4)])
# G2.add_edges_from([(2,3), (1,4)])
# 
# GM = isomorphism.GraphMatcher(G1,G2)
# print(GM.is_isomorphic())
# print(GM.mapping)

# with open(EXTENSION_PATH%"taille_max"+"result_k_max_4_10_toutes_aretes/digraphe_commun_0.7_55_0_taille_5.pickle", 'rb') as fichier :
#     mon_depickler = pickle.Unpickler(fichier)
#     graphe_commun = mon_depickler.load()
#     print(graphe_commun.nodes.data())
#     print(graphe_commun.edges.data())

# with open("graphs_2.92.pickle", 'rb') as fichier_tout :
#     mon_depickler_graphes = pickle.Unpickler(fichier_tout)
#     graphes = mon_depickler_graphes.load()
#     
#     print(graphes[('5J7L', 'DA')].nodes[685])
#     print(graphes[('5J7L', 'DA')][685])

# compte_gd_chgt = 0
# compte_changement = 0
# compte_plus_petit = 0
# compte_plus_grand = 0
# with open(EXTENSION_PATH%10+"sim_extensions_toutes_aretes_coeff_all1_taille_10.pickle", 'rb') as fichier_ecriture :
#     mon_depickler = pickle.Unpickler(fichier_ecriture)
#     dico_sim = mon_depickler.load()
#  
#     compteur = 0
#     liste_vus = {}
#     for fic in os.listdir("/home/coline/Bureau/graphes_extension_autres_tailles_new/graphes_extension_autres_tailles_new") :
#         if os.path.isfile("/home/coline/Bureau/graphes_extension_autres_tailles_new/graphes_extension_autres_tailles_new/"+fic) :
#              
#             with open("/home/coline/Bureau/graphes_extension_autres_tailles_new/graphes_extension_autres_tailles_new/"+fic, 'rb') as fichier_crible :
#                 mon_depickler = pickle.Unpickler(fichier_crible)
#                 crible = mon_depickler.load()
#              
#                 elt1 = fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6]
#                 elt2 = fic.split("_")[7] + "_" + fic.split("_")[8] + "_" + fic.split("_")[9] + "_" + fic.split("_")[10]
#                 if fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6] not in liste_vus.keys() :
#                     liste_vus.update({fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6] : [fic.split("_")[7] + "_" + fic.split("_")[8] + "_" + fic.split("_")[9] + "_" + fic.split("_")[10]]} )
#                 else :
#                     liste_vus[fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6]].append(fic.split("_")[7] + "_" + fic.split("_")[8] + "_" + fic.split("_")[9] + "_" + fic.split("_")[10])
#                   
#                 if fic.split("_")[7] + "_" + fic.split("_")[8] + "_" + fic.split("_")[9] + "_" + fic.split("_")[10] not in liste_vus.keys() :
#                     liste_vus.update({fic.split("_")[7] + "_" + fic.split("_")[8] + "_" + fic.split("_")[9] + "_" + fic.split("_")[10] : [fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6]]} )
#                 else :
#                     liste_vus[fic.split("_")[7] + "_" + fic.split("_")[8] + "_" + fic.split("_")[9] + "_" + fic.split("_")[10]].append(fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6])
#              
#                 compteur += 1
#                  
# #                 print(crible)
#                 for cle in crible.keys() :
#                     if (elt1, elt2) in dico_sim.keys() :
#                         sim1 = dico_sim[(elt1, elt2)] 
#                         sim2 = crible[cle][9] 
#                         diff = dico_sim[(elt1, elt2)] - crible[cle][9] 
#                     else :
#                         sim1 = dico_sim[(elt2, elt1)]
#                         sim2 = crible[cle][9] 
#                         diff = dico_sim[(elt2, elt1)] - crible[cle][9]
#                  
#                  
#                 if diff != 0 :
#                     #if max(sim1, sim2) > 0.7 :
# #                         print(elt1)
# #                         print(elt2)
# #                         print(sim1)
# #                         print(sim2)
# #                         print(diff)
#                     compte_changement += 1
#                 if diff < -0.1 or diff > 0.1 :
#                     #if max(sim1, sim2) > 0.7 :
# #                         print(elt1)
# #                         print(elt2)
# #                         print(sim1)
# #                         print(sim2)
# #                         print(diff)
#                     compte_gd_chgt += 1
#                         
#                 if diff < 0 :
#                     compte_plus_grand += 1
#                     
#                 elif diff > 0 :
#                     if max(sim1, sim2) > 0.6 :
#                         print(elt1)
#                         print(elt2)
#                         print(sim1)
#                         print(sim2)
#                         print(diff)
#                     compte_plus_petit += 1
#     print(compte_gd_chgt)
#     print(compte_changement)
#     print(compte_plus_petit)
#     print(compte_plus_grand)
#     print(liste_vus)        
#     print(compteur)


        
# compteur = 0
# liste_complete = []
# for fic in os.listdir(EXTENSION_PATH_TAILLE%10) :
#     if "pickle" in fic and "couples_possibles" not in fic and len(fic.split("_")) == 5 :
#         liste_complete.append(fic.split("_")[1]+ "_" + fic.split("_")[2]+ "_" + fic.split("_")[3]+ "_" + fic.split("_")[4][:len(fic.split("_")[4])-7])
#         compteur += 1
# 
# print(compteur)
# print(liste_complete)
# 
# with open(EXTENSION_PATH%"taille_max"+"result_k_max_4_10_toutes_aretes/sim_extensions_toutes_aretes_coeff_all1_max_taille_taille_max_avec_val_k.pickle", 'rb') as fichier_ecriture :
#     mon_depickler = pickle.Unpickler(fichier_ecriture)
#     dico_sim= mon_depickler.load()
# 
#     for elt in liste_vus.keys() :
#         if len(liste_vus[elt]) < 89 :
#             print("gros rat")
#             print(elt)
#             print(len(liste_vus[elt]))
#             for elt2 in liste_complete :
#                 if elt2 not in liste_vus[elt] and elt != elt2 :
#                     print(elt2)
#                     if (elt, elt2) in dico_sim.keys() :
#                         print(dico_sim[(elt, elt2)])
#                     else :
#                         print(dico_sim[(elt2, elt)])
        


# with open("script_suppr.sh", 'w') as fichier_sup :
#     for i in range(1,11) :
#         for fic in os.listdir(EXTENSION_PATH_TAILLE%i) :
#             if "graphe_comp" in fic and len(fic.split("_")) == 16 :
#                 new_fic = fic.split("_")[0] + "_" + fic.split("_")[1] + "_" + fic.split("_")[2] + "_" + fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6] + "_" + fic.split("_")[7] + "_" + fic.split("_")[8] + "_" + fic.split("_")[10] + "_" + fic.split("_")[11] + "_" + fic.split("_")[12] + "_" + fic.split("_")[13] + "_" + fic.split("_")[14] + ".pickle" 
#                 fichier_sup.write("mv graphes_extension_taille_%s/%s graphes_extension_taille_%s/%s\n"%(i, fic, i, new_fic))
#             elif "graphe_comp" not in fic and "couples_possibles" in fic and len(fic.split("_")) == 14 :
#                 new_fic = fic.split("_")[0] + "_" + fic.split("_")[1] + "_" + fic.split("_")[2] + "_" + fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6] + "_" + fic.split("_")[8] + "_" + fic.split("_")[9] + "_" + fic.split("_")[10] + "_" + fic.split("_")[11] + "_" + fic.split("_")[12] + ".pickle" 
#                 fichier_sup.write("mv graphes_extension_taille_%s/%s graphes_extension_taille_%s/%s\n"%(i, fic, i, new_fic))
#    
# with open("script_suppr_normal.sh", 'w') as fichier_sup :
#     for i in range(1,11) :
#         for fic in os.listdir(EXTENSION_PATH_TAILLE%i) :
#             if "graphe_comp" not in fic and "couples_possibles" not in fic :
#                 new_fic = fic.split("_")[0] + "_" + fic.split("_")[1] + "_" + fic.split("_")[2] + "_" + fic.split("_")[3] + "_" + fic.split("_")[4] + ".pickle" 
#                 fichier_sup.write("mv graphes_extension_taille_%s/%s graphes_extension_taille_%s/%s\n"%(i, fic, i, new_fic))


# with open("/home/coline/site_django/aminorsite/static/graphes_comparaison/toutes_aretes_non_can_int/script_suppr_nom_png.sh", 'w') as fichier_sup :
#     for i in range(1,6) :
#         for fic in os.listdir("/home/coline/site_django/aminorsite/static/graphes_comparaison/toutes_aretes_non_can_int/taille_%s"%i) :
# #             print(fic)
# #             print(len(fic.split("_")))
#             if "png" in fic :
#                 print(len(fic.split("_")))
#                 print(fic.split("_"))
#                 new_fic = fic.split("_")[2] + "_" + fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6] + "_" + fic.split("_")[7] + "_" + fic.split("_")[8] + "_" + fic.split("_")[9] + "_" + fic.split("_")[10] + "_" + fic.split("_")[11]
#                 fichier_sup.write("mv taille_%s/%s taille_%s/%s\n"%(i, fic, i, new_fic))
#     
#     for i in range(5,9) :        
#         for fic in os.listdir("/home/coline/site_django/aminorsite/static/graphes_comparaison/toutes_aretes_non_can_2/taille_%s"%i) :
# #             print(fic)
# #             print(len(fic.split("_")))
#             print(len(fic.split("_")))
#             print(fic.split("_"))
#             new_fic = fic.split("_")[2] + "_" + fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6] + "_" + fic.split("_")[7] + "_" + fic.split("_")[8] + "_" + fic.split("_")[9] + "_" + fic.split("_")[10] + "_" + fic.split("_")[11]
#             fichier_sup.write("mv taille_%s/%s taille_%s/%s\n"%(i, fic, i, new_fic))


# liste_a_faire = []
# for i in range(10,11) :
#     for fic in os.listdir(EXTENSION_PATH_TAILLE%i) :
#         if "graphe_comp" in fic and fic not in os.listdir('Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles_new/graphes_extension_taille_%s/'%i) :
#             cle1 = fic.split("_")[5] + "_" + fic.split("_")[6] + "_" + fic.split("_")[7] + "_" + fic.split("_")[8]
#             cle2 = fic.split("_")[10] + "_" + fic.split("_")[11] + "_" + fic.split("_")[12] + "_" + fic.split("_")[13][:len(fic.split("_")[13])-7]
#             
#             liste_a_faire.append((cle1, cle2))
# print(len(liste_a_faire))
# print(liste_a_faire)

# with open(EXTENSION_PATH%"taille_max"+"result_k_max_4_10_toutes_aretes/sim_extensions_toutes_aretes_coeff_all1_max_taille_taille_max_avec_val_k.pickle", 'rb') as fichier_ecriture :
#     mon_depickler = pickle.Unpickler(fichier_ecriture)
#     dico_sim = mon_depickler.load()
#     
#     if ("5J5B_BA_48_23", "5DM6_X_25_15") in dico_sim.keys() :
#         print(dico_sim[("5J5B_BA_48_23", "5DM6_X_25_15")])
#     else :
#         print(dico_sim[("5DM6_X_25_15", "5J5B_BA_48_23")])

# with open("graphs_2.92.pickle", 'rb') as fichier_tout :
#     mon_depickler_graphes = pickle.Unpickler(fichier_tout)
#     graphes = mon_depickler_graphes.load()
#     for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_10") :
#         #if "graphe_comp_test" in fic :
#         if "couples_possibles" not in fic and "_3.pickle" in fic and len(fic.split("_")) == 6  :
#             nom_mol = (fic.split("_")[1], fic.split("_")[2])
#             print(nom_mol)
#             with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_10/"+fic, 'rb') as fichier_1 :
#                 mon_depickler = pickle.Unpickler(fichier_1)
#                 graphe1 = mon_depickler.load()
#                 
#                 for u,v,key,data in graphe1.edges(data=True, keys=True) :
#                     if data["label"] not in ['B53', None, '0']:
#                         ok = False
#                         for edge in graphe1[v][u] :
#                             #print(data["label"])
#                             if data["label"][1] == graphe1[v][u][edge]["label"][2] and data["label"][2] == graphe1[v][u][edge]["label"][1] :
#                                 ok = True
#                         if not ok :
#                             print(fic)
#                             print(data["label"])
#                             print(graphe1[u][v])
#                             print(graphe1[v][u])
#                         
#                         if data["label"] not in ['B53', 'CWW', 'CWWn', '0', None] and (graphe1.nodes[u]["position"][0], graphe1.nodes[v]["position"][0]) in graphes[nom_mol].edges() and  data["label"] != graphes[nom_mol].edges[graphe1.nodes[u]["position"][0],graphe1.nodes[v]["position"][0]]["label"] :
#                             print(fic)
#                             #print(graphes[nom_mol].edges[graphe1.nodes[u]["position"][0],graphe1.nodes[v]["position"][0]]["label"])
#                             print(data["label"])
# #                         elif data["label"] not in ['B53', 'CWW', 'CWWn', '0', None] :
# #                             print(data["label"])
# #                             print("ok")
    
  
## Verif pas de doublon dans les aretes  (29 mai 2019)
# for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_10") :
#         #if "graphe_comp_test" in fic :
#         if "couples_possibles" not in fic and "_3.pickle" in fic and len(fic.split("_")) == 6  :
#             with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_10/"+fic, 'rb') as fichier_1 :
#                 mon_depickler = pickle.Unpickler(fichier_1)
#                 graphe1 = mon_depickler.load()
#                 
#                 compteur = 0
#                 for u, dicto in graphe1.adjacency():
#                     for cle1 in dicto.keys() :
#                         compteur = 0
#                         for cle2 in dicto.keys() :
#                             if cle1 == cle2 :
#                                 idem = 0
#                                 for edge1 in dicto[cle1].keys() :
#                                     for edge2 in dicto[cle2].keys() :
#                                         if dicto[cle1][edge1]["label"] == dicto[cle2][edge2]["label"] and dicto[cle1][edge1]["long_range"] == dicto[cle2][edge2]["long_range"] :
#                                             idem += 1
#                                 if idem == len(dicto[cle1]) :
#                                     compteur += 1
#                         print(compteur)
#                         if compteur > 1 :
#                             print("bizarre")

# with open("liste_extensions.pickle", 'wb') as fichier_pickle :
#     mon_pickler = pickle.Pickler(fichier_pickle)
#     liste = []
#     for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_10") :
#         if "couples_possibles" not in fic and "_3.pickle" in fic and len(fic.split("_")) == 6  :  
#             elt = fic.split("_")[1] + "_" + fic.split("_")[2] + "_" + fic.split("_")[3] + "_" + fic.split("_")[4]
#             liste.append(elt)
#             
#             
#     print(len(liste))
#     print(liste)
#     mon_pickler.dump(liste)
            
# with open("/home/coline/Bureau/graphes_extension_autres_tailles_new/graphes_extension_autres_tailles_new/graphes_extension_taille_10/couples_possibles_fichier_5J7L_DA_30_15_3_fichier_5DM6_X_127_7_3.pickle", 'rb') as fichier_1 :
#                 mon_depickler = pickle.Unpickler(fichier_1)
#                 graphe1 = mon_depickler.load()
#                 
#                 print(graphe1[2])

# 
# j = 5
# compteur_not_idem = 0
# compteur_tot = 0
# cles = []
# diff_moyenne = 0
# mini = 1.0
# maxi = 0.0
# for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles_new/graphes_extension_taille_%s"%j) :
#     if "graphe_comp" in fic :
#         with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles_new/graphes_extension_taille_%s/"%j+fic, 'rb') as fichier :
#             mon_depickler = pickle.Unpickler(fichier)
#             dico_graphe = mon_depickler.load()
#              
#             with open(EXTENSION_PATH%j+"sim_extensions_toutes_aretes_coeff_all1_taille_%d.pickle"%j, 'rb') as fichier_sim :
#                 mon_depickler_sim = pickle.Unpickler(fichier_sim)
#                 dico_sim_new = mon_depickler_sim.load()
#                 
#                 with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles_new_dernier/result_graphes_extension_%s/sim_extensions_toutes_aretes_coeff_all1_taille_%d.pickle"%(j,j), 'rb') as fichier_sim_dernier :
#                     mon_depickler_sim_dernier = pickle.Unpickler(fichier_sim_dernier)
#                     dico_sim_new_dernier = mon_depickler_sim_dernier.load()
# #                     print(dico_sim_new.keys())
#                     for cle in dico_graphe.keys() :
#                         
#                         cle1 = cle[0].split("_")[0] + "_" + cle[0].split("_")[1] + "_" + cle[0].split("_")[2] + "_" + cle[0].split("_")[3] + "_" + cle[0].split("_")[4]
#                         cle2 = cle[1].split("_")[0] + "_" + cle[1].split("_")[1] + "_" + cle[1].split("_")[2] + "_" + cle[1].split("_")[3] + "_" + cle[1].split("_")[4]
#                          
#         #                 if cle1 in ["fichier_5FDU_1A_272_1", "fichier_5J7L_DA_197_4"] and cle2 in ["fichier_5FDU_1A_272_1", "fichier_5J7L_DA_197_4"] :
#         #                     print(dico_graphe[cle].nodes.data())
#                         
#                         with open(EXTENSION_PATH_TAILLE%j+cle1+".pickle", 'rb') as fichier_graphe1 :
#                             mon_depickler_graphe1 = pickle.Unpickler(fichier_graphe1)
#                             graphe1 = mon_depickler_graphe1.load()
#                              
#                             with open(EXTENSION_PATH_TAILLE%j+cle2+".pickle", 'rb') as fichier_graphe2 :
#                                 mon_depickler_graphe2 = pickle.Unpickler(fichier_graphe2)
#                                 graphe2 = mon_depickler_graphe2.load()
#                                  
#                                 for u,v,data in dico_graphe[cle].edges(data=True) :
#                                     if data["type"] == 'NON_CAN' :
#                                         compteur_tot += 1
#                                         idem = False
#                                         #print(u)
#                                         #print(v)
#                                         for edge1 in graphe1[u[0]][v[0]] :
#                                             for edge2 in graphe2[u[1]][v[1]] :                             
#                                                 if graphe1[u[0]][v[0]][edge1]["label"] == graphe2[u[1]][v[1]][edge2]["label"] :
#                                                     idem = True
#                                         if not idem :
#                                             
#                                             
#                                             cle_sim_1 = cle[0].split("_")[1] + "_" + cle[0].split("_")[2] + "_" + cle[0].split("_")[3] + "_" + cle[0].split("_")[4]
#                                             cle_sim_2 = cle[1].split("_")[1] + "_" + cle[1].split("_")[2] + "_" + cle[1].split("_")[3] + "_" + cle[1].split("_")[4]
#                                             
#                                              
#                                             
#                                             if (cle_sim_1,cle_sim_2) in dico_sim_new.keys() :
#                                                 sim_new = dico_sim_new[(cle_sim_1,cle_sim_2)]
#                                             else :
#                                                 sim_new = dico_sim_new[(cle_sim_2,cle_sim_1)]
#                                             
#                                             if (cle_sim_1,cle_sim_2) in dico_sim_new_dernier.keys() :
#                                                 sim_new_dernier = dico_sim_new_dernier[(cle_sim_1,cle_sim_2)]
#                                             else :
#                                                 sim_new_dernier = dico_sim_new_dernier[(cle_sim_2,cle_sim_1)]
#                                                 
# #                                             if sim_new > 0.5 and sim_new_dernier > 0.5 :
# #                                                 print(cle)
# #                                                 print(graphe1[u[0]][v[0]][edge1]["label"])
# #                                                 print(graphe2[u[1]][v[1]][edge2]["label"])
# #                                                 print(sim_new)
# #                                                 print(sim_new_dernier)
#                                             if abs(sim_new - sim_new_dernier) != 0 :
#                                                 if cle not in cles :
#                                                     cles.append(cle)
#     
#                                                 diff_moyenne += abs(sim_new - sim_new_dernier)   
#                                                 if abs(sim_new - sim_new_dernier) > maxi :
#                                                     maxi = abs(sim_new - sim_new_dernier)  
#                                                     print(sim_new)
#                                                     print(sim_new_dernier) 
#                                                 if abs(sim_new - sim_new_dernier) < mini :
#                                                     mini = abs(sim_new - sim_new_dernier)   
#                                                 compteur_not_idem += 1
#                                 compteur_tot -= 3   
# print(compteur_tot)
# print(compteur_not_idem)    
# print(len(cles))       
# print(diff_moyenne/compteur_not_idem)  
# print(mini)
# print(maxi)                   
                                


# with open("graphs_2.92.pickle", 'rb') as fichier_tout :
#         mon_depickler_graphes = pickle.Unpickler(fichier_tout)
#         graphes = mon_depickler_graphes.load()    
#         
#         for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_10") :
#             if "couples_possibles" not in fic and "_3.pickle" in fic and len(fic.split("_")) == 6  :  
#                 
#                 with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_10/"+fic, 'rb') as fichier_graphe :
#                     mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
#                     graphe = mon_depickler_graphe.load()
#                     nom = (fic.split("_")[1], fic.split("_")[2])
#                     
#                     for i in range(1,5) :
#                         pos = graphe.nodes[i]["position"][0]
#                         if i in [1,4] :
#                             int_t = 1
#                         else :
#                             int_t = -1
#                         for j in range(10) :
#                             if pos+j*int_t > 0 and pos+j*int_t < graphes[nom].number_of_nodes() :
#                                 for voisin in graphes[nom][pos+j*int_t] :
#     #                                 print(graphes[nom][pos+j*int_t])
#     #                                 print(graphes[nom].nodes[pos+j*int_t])
#     #                                 print(graphes[nom].nodes[voisin])
#                                     if (pos+j*int_t, voisin) in graphes[nom].edges() and graphes[nom].edges[pos+j*int_t, voisin]["label"] == 'CWW' \
#                                     and graphes[nom].edges[pos+j*int_t, voisin]["long_range"] == False and \
#                                     not ((graphes[nom].nodes[pos+j*int_t]["nt"] == 'A' and graphes[nom].nodes[voisin]["nt"] == 'U') \
#                                     or (graphes[nom].nodes[pos+j*int_t]["nt"] == 'U' and graphes[nom].nodes[voisin]["nt"] == 'A') \
#                                     or (graphes[nom].nodes[pos+j*int_t]["nt"] == 'C' and graphes[nom].nodes[voisin]["nt"] == 'G') \
#                                     or (graphes[nom].nodes[pos+j*int_t]["nt"] == 'G' and graphes[nom].nodes[voisin]["nt"] == 'C') \
#                                     or (graphes[nom].nodes[pos+j*int_t]["nt"] == 'G' and graphes[nom].nodes[voisin]["nt"] == 'U') \
#                                     or (graphes[nom].nodes[pos+j*int_t]["nt"] == 'U' and graphes[nom].nodes[voisin]["nt"] == 'G')) :
#                                         print(fic)
#                                         print("%s %s"%(graphes[nom].nodes[pos+j*int_t]["nt"], graphes[nom].nodes[voisin]["nt"]))
                                    
                                    
# with open(EXTENSION_PATH%1+"dico_comp_complet_metrique_toutes_aretes_coeff_all1_taille_%s.pickle"%1, 'rb') as fichier_graphe :
#     mon_depickler = pickle.Unpickler(fichier_graphe)
#     dico_graphe = mon_depickler.load()
#     with open("liste_extensions.pickle", 'rb') as fichier_pickle :
#         mon_depickler = pickle.Unpickler(fichier_pickle)
#         liste_ext = mon_depickler.load()
#         compteur = 0
#         for i in range(len(liste_ext)) :
#             for j in range(i+1, len(liste_ext)) :
#                 ok = False
#                 for cle in dico_graphe.keys() :
#                     cle1 = cle[0].split("_")[1] + "_" + cle[0].split("_")[2] + "_" + cle[0].split("_")[3] + "_" + cle[0].split("_")[4]
#                     cle2 = cle[1].split("_")[1] + "_" + cle[1].split("_")[2] + "_" + cle[1].split("_")[3] + "_" + cle[1].split("_")[4]
#                      
#                     if liste_ext[i] in [cle1,cle2] and liste_ext[j] in [cle1,cle2] :
#                         ok = True
#                 if not ok :
#                     print("petit rat")
#                     print(liste_ext[i])
#                     print(liste_ext[j])
#                     compteur += 1
#                      
# print(compteur)
                    

# with open(EXTENSION_PATH%8+"sim_extensions_toutes_aretes_coeff_all1_taille_8.pickle", 'rb') as fichier_graphe :
#     mon_depickler = pickle.Unpickler(fichier_graphe)
#     dico_graphe = mon_depickler.load()
#      
#     print(dico_graphe[('5DM6_X_127_7', '4V9F_0_134_5')])
#     
#     with open(EXTENSION_PATH%8+"dico_comp_complet_metrique_toutes_aretes_coeff_all1_taille_8.pickle", 'rb') as fichier_graphe_complet :
#         mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
#         dico_graphe_complet = mon_depickler_complet.load()
#         
#         print(dico_graphe_complet.keys())
#         
#         with open(EXTENSION_PATH_TAILLE%8+"fichier_5DM6_X_127_7.pickle", 'rb') as fichier_graphe_1 :
#             mon_depickler_1 = pickle.Unpickler(fichier_graphe_1)
#             graphe1= mon_depickler_1.load()
#             
#             with open(EXTENSION_PATH_TAILLE%8+"fichier_4V9F_0_134_5.pickle", 'rb') as fichier_graphe_2 :
#                 mon_depickler_2 = pickle.Unpickler(fichier_graphe_2)
#                 graphe2 = mon_depickler_2.load()
#             
#                 calcul_sim_aretes_avec_coeff(dico_graphe_complet[('fichier_5DM6_X_127_7', 'fichier_4V9F_0_134_5')], graphe1, graphe2, ('5DM6_X_127_7', '4V9F_0_134_5'), 1,1,1)
#     
# with open('Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles_new/result_graphes_extension_8/sim_extensions_toutes_aretes_coeff_all1_taille_8.pickle', 'rb') as fichier_graphe :
#     mon_depickler = pickle.Unpickler(fichier_graphe)
#     dico_graphe = mon_depickler.load()
#      
#     print(dico_graphe[('5DM6_X_127_7', '4V9F_0_134_5')])
#         
#     with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles_new/result_graphes_extension_8/dico_comp_complet_metrique_toutes_aretes_coeff_all1_taille_8.pickle", 'rb') as fichier_graphe_complet :
#             mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
#             dico_graphe_complet = mon_depickler_complet.load()
#             
#             with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles_new/graphes_extension_taille_8/fichier_5DM6_X_127_7.pickle", 'rb') as fichier_graphe_1 :
#                 mon_depickler_1 = pickle.Unpickler(fichier_graphe_1)
#                 graphe1= mon_depickler_1.load()
#                 
#                 with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles_new/graphes_extension_taille_8/fichier_4V9F_0_134_5.pickle", 'rb') as fichier_graphe_2 :
#                     mon_depickler_2 = pickle.Unpickler(fichier_graphe_2)
#                     graphe2 = mon_depickler_2.load()
#                 
#                     calcul_sim_aretes_avec_coeff(dico_graphe_complet[('5DM6_X_127_7', '4V9F_0_134_5')], graphe1, graphe2, ('5DM6_X_127_7', '4V9F_0_134_5'), 1,1,1)

# with open("/home/coline/site_django/aminorsite/static/graphes_comparaison/toutes_aretes_non_can_int/script_copie.sh", 'w') as fichier_script :                              
#     for i in range(1,6) :
#         for fic in os.listdir(EXTENSION_PATH_TAILLE%i) :
#             if "pickle" not in fic and len(fic.split("_")) > 5 :
#                 fichier_script.write("cp %s %s\n"%("/home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/"+EXTENSION_PATH_TAILLE%i+fic, "/home/coline/site_django/aminorsite/static/graphes_comparaison/toutes_aretes_non_can_int/taille_%d"%i ))
# with open("/home/coline/site_django/aminorsite/static/graphes_comparaison/toutes_aretes_non_can_2/taille_9/script_enlever_test.sh", 'w') as fichier_script :                              
#         for fic in os.listdir("/home/coline/site_django/aminorsite/static/graphes_comparaison/toutes_aretes_non_can_2/taille_9") :
#             new_fic = fic[:len(fic)-9]+".png"
#             print(new_fic)
#             fichier_script.write("mv %s %s\n"%(fic, new_fic))
                           
# with open("fichier_script_copie.sh", 'w') as fichier_script :
#     for i in range(1,11) :
#         for fic in os.listdir("/home/coline/Bureau/graphes_extension_autres_tailles_new/graphes_extension_taille_%d"%i) :
#             if "couples_possibles" in fic :
#                 fichier_script.write("cp /home/coline/Bureau/graphes_extension_autres_tailles_new/graphes_extension_taille_%d/"%i+fic+" /home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles_new/graphes_extension_taille_%d\n"%i)

##Stats sur les sommets isoles dans les comparaisons 
# for i in range(1,11) : 
#     compte_sommet_isole = 0
#     compte_sommet_tot = 0
#     compte_sommet_isole_par_mol = 0
#     compte_liaisons_non_can = 0
#     compte_liaisons_can = 0
#     compte_liaisons_2 = 0
#     compte_bizarre = 0
#     compte_mol_bizarre = 0
#     for fic in os.listdir(EXTENSION_PATH_TAILLE%i) :
#         if "graphe_comp" in fic :
#             with open(EXTENSION_PATH_TAILLE%i+fic, 'rb') as fichier_pickle :
#                 mon_depickler = pickle.Unpickler(fichier_pickle)
#                 graphe_comp = mon_depickler.load()
#                 
#                 
#                 
#                 dans_mol = False
#                 bizarre = False
#                 for cle in graphe_comp.keys() :
#                     if len(cle[0].split("_")) != 5 :
#                         cle_new = cle[0].split("_")[0] + "_" + cle[0].split("_")[1] + "_" +cle[0].split("_")[2] + "_" +cle[0].split("_")[3] + "_" + cle[0].split("_")[4]
#                     else :
#                         cle_new = cle[0]
#                     with open(EXTENSION_PATH_TAILLE%i+cle_new+".pickle", 'rb') as fichier_graphe1 :
#                         mon_depickler_graphe1 = pickle.Unpickler(fichier_graphe1)
#                         graphe1 = mon_depickler_graphe1.load()
#                         for noeud in graphe_comp[cle].nodes() :
#                             if len(graphe_comp[cle][noeud]) == 0 :
#                                 compte_sommet_isole += 1
#                                 dans_mol = True
#                                 
#                                 if graphe1.nodes[noeud[0]]["type"] == 3 :
#                                     compte_liaisons_non_can += 1
#                                 elif graphe1.nodes[noeud[0]]["type"] == 1 : 
#                                     compte_liaisons_can += 1
#                                 else : 
#                                     if graphe1.nodes[noeud[0]]["type"] != 2 : 
# #                                         print("bizarre")
# #                                         print(fic)
# #                                         print(graphe1.nodes[noeud[0]]["type"])
#                                         compte_bizarre += 1
#                                         bizarre = True
#                                     else :
#                                         compte_liaisons_2 += 1 
#                 compte_sommet_tot += graphe_comp[cle].number_of_nodes() - 5
#                 if dans_mol :
#                     compte_sommet_isole_par_mol += 1
#                 if bizarre :
#                     compte_mol_bizarre += 1
#                 
#     print("i : %d"%i)
#     print("nb de sommets isoles : %d "%compte_sommet_isole)
#     print("nb de sommets isoles de type 3 : %d"%compte_liaisons_non_can)
#     print("nb de sommets isoles de type 1 : %d"%compte_liaisons_can)
#     print("nb de sommets isoles de type 2 : %d"%compte_liaisons_2)
#     
#     print("nb de sommets au total : %d"%compte_sommet_tot)
#     print("proportion des sommets isoles de type 3 sur le total des sommets isoles : %f"%(compte_liaisons_non_can/max(1,compte_sommet_isole)))
#     print("proportion des sommets isoles de type 1 sur le total des sommets isoles : %f"%(compte_liaisons_can/max(1,compte_sommet_isole)))
#     print("proportion des sommets isoles de type 2 sur le total des sommets isoles : %f"%(compte_liaisons_2/max(1,compte_sommet_isole)))
#     print("proportion des sommets isoles sur le total des sommets : %f"%(compte_sommet_isole/compte_sommet_tot))
#     print("nb moyen de sommets isoles par comparaison : %f"%(compte_sommet_isole/4005))
#     
#     print("nb de comparaison avec des sommets isoles (un ou plusieurs) : %d"%compte_sommet_isole_par_mol)
#             
#                             
#     print("nb de sommets isoles bizarres : %d"%compte_bizarre)                         
#     print("nb de comparaison avec des sommets isoles bizarres (un ou plusieurs) : %d"%compte_mol_bizarre)
                            
    
    
# with open(EXTENSION_PATH_TAILLE%4+"graphe_comp_couples_possibles_fichier_5FDU_1A_272_1_fichier_5J7L_DA_272_2.pickle", 'rb') as fichier_pickle :                        
#         mon_depickler = pickle.Unpickler(fichier_pickle)
#         comp_sim = mon_depickler.load()
#         
#         for cle in comp_sim.keys() :
#             print(cle)
#             print(comp_sim[cle].edges.data())
#             with open(EXTENSION_PATH_TAILLE%4+cle[0][:len(cle[0])-2]+".pickle", 'rb') as fichier_graphe1 :
#                 mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
#                 graphe1 = mon_depickler_1.load()
#                 
#                 with open(EXTENSION_PATH_TAILLE%4+cle[1][:len(cle[1])-2]+".pickle", 'rb') as fichier_graphe2 :
#                     mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
#                     graphe2 = mon_depickler_2.load()  
#                     
#                     sim = calcul_sim_aretes_avec_coeff(graphe1, graphe2, comp_sim[cle], cle, 1, 1, 1)
#                     print(sim)
                    
                    
# with open(EXTENSION_PATH_TAILLE%7+"graphe_comp_couples_possibles_fichier_5DM6_X_48_10_fichier_5J7L_DA_62_5.pickle", 'rb') as fichier_pickle_1 :                        
#         mon_depickler = pickle.Unpickler(fichier_pickle_1)
#         comp_sim = mon_depickler.load()     
#         for cle in comp_sim.keys() :          
#             print(comp_sim[cle].edges.data())
#         
#         with open(EXTENSION_PATH_TAILLE%7+"graphe_comp_test_couples_possibles_fichier_5DM6_X_48_10_fichier_5J7L_DA_62_5.pickle", 'rb') as fichier_pickle_2 :                        
#             mon_depickler = pickle.Unpickler(fichier_pickle_2)
#             comp_sim_new = mon_depickler.load() 
#             for cle in comp_sim_new.keys() :
#                 print(comp_sim_new[cle].edges.data())
#                 print(comp_sim_new[cle].nodes.data())
#                 
#                                     
#             with open(EXTENSION_PATH_TAILLE%7+"fichier_5DM6_X_48_10.pickle", 'rb') as fichier_pickle_3 :                        
#                 mon_depickler = pickle.Unpickler(fichier_pickle_3)
#                 graphe1 = mon_depickler.load() 
#                 print(graphe1.edges.data())               
#                 
#                 with open(EXTENSION_PATH_TAILLE%7+"fichier_5J7L_DA_62_5.pickle", 'rb') as fichier_pickle_4 :                        
#                     mon_depickler = pickle.Unpickler(fichier_pickle_4)
#                     graphe2 = mon_depickler.load() 
#                     print(graphe2.edges.data())            

# with open("script_copie_serveur.sh", 'w') as fichier :  
#     with open("script_copie_serveur_couples.sh", 'w') as fichier_couples :                        
#         for i in range(1,11) :
#             for fic in os.listdir(EXTENSION_PATH_TAILLE%i) :
#                 if "pickle" in fic and "couples_possibles" not in fic :
#                     fichier.write("cp /home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/"+EXTENSION_PATH_TAILLE%i + fic + " /home/coline/Bureau/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_taille_%d/"%i+fic + "\n")
#                 elif "graphe_comp" not in fic : 
#                     fichier_couples.write("cp /home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/"+EXTENSION_PATH_TAILLE%i + fic + " /home/coline/Bureau/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_taille_%d/"%i+fic + "\n")
                                                         

# with open(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes/"+"digraphe_commun_0.7_clustering_perez_groupe_11_taille_8.pickle", 'rb') as fichier_pickle_1 :                        
#         mon_depickler = pickle.Unpickler(fichier_pickle_1)
#         comp_sim = mon_depickler.load()    
#         
#         print(comp_sim.edges.data())    
#         
#         for i in range(len(CLUSTERING_PEREZ_VERSION_NON_CAN_2[11])) :
#             for j in range(i+1, len(CLUSTERING_PEREZ_VERSION_NON_CAN_2[11])) :
#                 with open(EXTENSION_PATH_TAILLE%8+"graphe_comp_couples_possibles_fichier_%s_fichier_%s.pickle"%(CLUSTERING_PEREZ_VERSION_NON_CAN_2[11][i], CLUSTERING_PEREZ_VERSION_NON_CAN_2[11][j]), 'rb') as fichier_pickle_2 :                        
#                     mon_depickler = pickle.Unpickler(fichier_pickle_2)
#                     graphe1 = mon_depickler.load() 
#                     
#                     for cle in graphe1.keys() :
#                         print(graphe1[cle].edges.data())   
                        
                        

# with open("graphs_2.92.pickle", 'rb') as fichier_tout :
#         mon_depickler_graphes = pickle.Unpickler(fichier_tout)
#         graphes = mon_depickler_graphes.load()
#         
#         print(graphes[('3JCS', '1')][867])                         

## Verif difference version non_can_2 1 et 2 (01 juillet 2019)
# j = 4
# pas_bon = []
# for fic in os.listdir(EXTENSION_PATH_TAILLE%j) :
#         #if "graphe_comp_test" in fic :
#         if "couples_possibles" not in fic and ".pickle" in fic and len(fic.split("_")) == 5  :
#             #fic_2 = fic[:11]+"_"+fic[17:34]+"_"+fic[40:]
#             if fic in os.listdir("/home/coline/Bureau/graphes_extension_autres_tailles_new_version_avec_non_can_result/graphes_extension_taille_%d"%j) :
#                 with open(EXTENSION_PATH_TAILLE%j+fic, 'rb') as fichier_1 :
#                         mon_depickler = pickle.Unpickler(fichier_1)
#                         graphe1 = mon_depickler.load()
#                                   
#                                       
#                         with open("/home/coline/Bureau/graphes_extension_autres_tailles_new_version_avec_non_can_result/graphes_extension_taille_%d/"%j+fic, 'rb') as fichier_2 :
#                                 mon_depickler_2 = pickle.Unpickler(fichier_2)
#                                 graphe2 = mon_depickler_2.load()
#                                 print(type(graphe1))
#                                 print(type(graphe2))
#                                 
#                                 print(graphe1.edges.data())
#                                 print(graphe2.edges.data())
#                             #for cle in graphe1.keys() :     
#                                 print(fic)
#                                 #print(cle)
#                                 for noeud, data in graphe1.nodes(data=True) :
#                                     if noeud not in graphe2.nodes() :
#                                         print("pb1")
#                                         print(noeud)
#                                         if fic[:len(fic)-7] not in pas_bon :
#                                             pas_bon.append(fic[:len(fic)-7])
#                                     else :
#                                         if data["type"] != graphe2.nodes[noeud]["type"] :
#                                             if fic[:len(fic)-7] not in pas_bon :
#                                                 pas_bon.append(fic[:len(fic)-7])
#                                               
#                                 for noeud, data in graphe2.nodes(data=True) :
#                                     if noeud not in graphe1.nodes() :
#                                         print("pb2")
#                                         print(noeud)
#                                         if fic[:len(fic)-7] not in pas_bon :
#                                             pas_bon.append(fic[:len(fic)-7]) 
#                                     else :
#                                         if data["type"] != graphe2.nodes[noeud]["type"] :
#                                             if fic[:len(fic)-7] not in pas_bon :
#                                                 pas_bon.append(fic[:len(fic)-7])           
#                                               
#                                                   
#                                 for u,v, data in graphe1.edges(data=True) :
#                                     if (u,v) not in graphe2.edges() :
#                                         print("pb3")
#                                         print((u,v))
#                                         if fic[:len(fic)-7] not in pas_bon :
#                                                     
#                                             pas_bon.append(fic[:len(fic)-7]) 
#                                                        
#                                     elif (u,v) in graphe2.edges() :
#                                         compteur_bon = 0
#                                         for edge_1 in graphe1[u][v] :
#                                             for edge_2 in graphe2[u][v] :
#                                                 if graphe1[u][v][edge_1]["label"] == graphe2[u][v][edge_2]["label"] or ( len(graphe1[u][v][edge_1]["label"]) == 3 and  len(graphe2[u][v][edge_2]["label"]) == 3 and graphe1[u][v][edge_1]["label"][1] == graphe2[u][v][edge_2]["label"][2] and graphe1[u][v][edge_1]["label"][2] == graphe2[u][v][edge_2]["label"][1])    :
#                                                     compteur_bon += 1
#                                         if compteur_bon != len(graphe1[u][v]) :
#                                             print("pb5")
#                                             print((u,v))
#                                             print(graphe2[u][v][edge_2]["label"])
#                                             if fic[:len(fic)-7] not in pas_bon :
#                                                 pas_bon.append(fic[:len(fic)-7])
#                                                           
#                                 for u,v, data in graphe2.edges(data=True) :
#                                     if (u,v) not in graphe1.edges() :
#                                         print("pb4")
#                                         print((u,v))
#                                         if fic[:len(fic)-7] not in pas_bon :
#                                             pas_bon.append(fic[:len(fic)-7]) 
#                                                
#                                     elif (u,v) in graphe1.edges() :
#                                         compteur_bon = 0
#                                         for edge_1 in graphe1[u][v] :
#                                             for edge_2 in graphe2[u][v] :
#            
#                                                 if graphe1[u][v][edge_1]["label"] == graphe2[u][v][edge_2]["label"] or ( len(graphe1[u][v][edge_1]["label"]) >= 3 and  len(graphe2[u][v][edge_2]["label"]) >= 3 and graphe1[u][v][edge_1]["label"][1] == graphe2[u][v][edge_2]["label"][2] and graphe1[u][v][edge_1]["label"][2] == graphe2[u][v][edge_2]["label"][1])    :
#                                                     compteur_bon += 1
#                                         if compteur_bon != len(graphe2[u][v]) :
#                                             print("pb6")
#                                             print((u,v))
#                                             print(graphe2[u][v][edge_2]["label"])
#                                             if fic[:len(fic)-7] not in pas_bon :
#                                                 pas_bon.append(fic[:len(fic)-7])          
# print(pas_bon) 
# print(len(pas_bon))


# with open(EXTENSION_PATH_TAILLE%8+"graphe_comp_couples_possibles_fichier_5DM6_X_48_10_fichier_5FDU_1A_197_3.pickle", 'rb') as fichier :
#     mon_depickler_3 = pickle.Unpickler(fichier)
#     graphe_comp = mon_depickler_3.load()
#     
#     for cle in graphe_comp.keys() :
#         print(graphe_comp[cle].edges.data())

# with open(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes"+"digraphe_commun_0.7_clustering_perez_groupe_12_taille_4.pickle", "rb") as fichier :
#     mon_depickler = pickle.Unpickler(fichier)
#     digraphe = mon_depickler.load()
#     
#     print(digraphe.nodes.data())

# with open("graphs_2.92.pickle", 'rb') as fichier_tout :
#     mon_depickler_graphes = pickle.Unpickler(fichier_tout)
#     graphes = mon_depickler_graphes.load()      
#      
#     print(graphes[('4L81', 'A')][83])                


## 04 juillet 2019
## Deplacer les fichiers representation 2D sur Django et les transformer en png
# with open("/home/coline/site_django/aminorsite/static/graphes_comparaison/graphes_sequence/script_copie.sh", 'w') as fichier_script :                              
#     for i in range(4,11) :
#         for fic in os.listdir("/home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/Graphes_globaux/graphes_sequence/taille_%d/"%i) :
#             if ".pdf" in fic and '_' in fic :
#                 fichier_script.write("cp %s %s\n"%("/home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/Graphes_globaux/graphes_sequence/taille_%d/"%i+fic, "/home/coline/site_django/aminorsite/static/graphes_comparaison/graphes_sequence/taille_%d"%i ))


# with open("/home/coline/site_django/aminorsite/static/graphes_comparaison/graphes_sequence/script_copie.sh", 'w') as fichier_script :                              
#     for i in range(4,11) :
#         for fic in os.listdir("/home/coline/site_django/aminorsite/static/graphes_comparaison/graphes_sequence/taille_%d"%i) :
#             if ".pdf" in fic and len(fic.split("-")) == 7:
#                 print(fic)
#                 fic_new = fic.split("-")[0] +"_" + fic.split("-")[1] + "_" + fic.split("-")[2] + "_" + fic.split("-")[3] + "_"+ fic.split("-")[4] + "_"+ fic.split("-")[5] + "_"+ fic.split("-")[6]
#                 fichier_script.write("mv %s %s\n"%("/home/coline/site_django/aminorsite/static/graphes_comparaison/graphes_sequence/taille_%d/"%i+fic, "/home/coline/site_django/aminorsite/static/graphes_comparaison/graphes_sequence/taille_%d/"%i+fic_new ))

# with open("/home/coline/site_django/aminorsite/static/graphes_comparaison/graphes_sequence/script_conversion.sh", 'w') as fichier_script :                              
#     for i in range(4,11) :
#         for fic in os.listdir("/home/coline/site_django/aminorsite/static/graphes_comparaison/graphes_sequence/taille_%d"%i) :
#             if ".pdf" in fic :
#                 fichier_script.write("pdftoppm %s %s -png\n"%("taille_%d/"%i+fic, "taille_%d/"%i+fic[:len(fic)-4] ))



def distrib_sim_range_par_taille_ext(groupe, nom_groupe):
    os.makedirs(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes/"+"distrib_groupe_%s"%nom_groupe, exist_ok=True)
    for i in range(len(groupe)) :
        for j in range(i+1,len(groupe)) : 
            print(groupe[i])
            print(groupe[j])
            liste_sim = []
            for taille in range(4,11) :
                
                with open(EXTENSION_PATH%taille+"sim_extensions_toutes_aretes_coeff_all1_taille_%d.pickle"%taille, 'rb') as fichier :
                    mon_depickler = pickle.Unpickler(fichier)
                    dico_sim = mon_depickler.load()
                    if (groupe[i], groupe[j]) in dico_sim.keys() :
                        liste_sim.append(dico_sim[(groupe[i], groupe[j])])
                    else :
                        liste_sim.append(dico_sim[(groupe[j], groupe[i])])    
                        
            print(liste_sim)
            
            ax1 = plt.gca()
            plt.plot(range(4,11),liste_sim) 
            ax1.xaxis.set_ticks(range(4,11))
            ax1.yaxis.set_ticks(np.arange(0,1,0.1))
            ax1.set_xlabel("Taille d'extension")
            ax1.set_ylabel("Similarite")
    plt.title("Distribution des similarites dans le groupe \n %s par taille d'extension"%nom_groupe)
    plt.savefig(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes/"+"distrib_groupe_%s/distrib_sim.png"%nom_groupe)
    #plt.show()
    
    for i in range(len(groupe)) :
        for j in range(i+1,len(groupe)) : 
            print(groupe[i])
            print(groupe[j])
            liste_rang = []
            for taille in range(4,11) :
                
                with open(EXTENSION_PATH%taille+"sim_extensions_toutes_aretes_coeff_all1_taille_%d.pickle"%taille, 'rb') as fichier :
                    mon_depickler = pickle.Unpickler(fichier)
                    dico_sim = mon_depickler.load()
                    if (groupe[i], groupe[j]) in dico_sim.keys() :
                        rang = rang_sim(dico_sim, dico_sim[(groupe[i], groupe[j])])
                        liste_rang.append(rang)
                    else :
                        rang = rang_sim(dico_sim, dico_sim[(groupe[j], groupe[i])])
                        liste_rang.append(rang)      
                        
            print(liste_rang)
            
            ax1 = plt.gca()
            plt.plot(range(4,11),liste_rang) 
            ax1.xaxis.set_ticks(range(4,11))
            ax1.yaxis.set_ticks(np.arange(1,4005,500))
            ax1.set_xlabel("Taille d'extension")
            ax1.set_ylabel("Rang de similarité")
    plt.title("Distribution des rangs de similarites dans le groupe \n %s par taille d'extension"%nom_groupe)
    plt.savefig(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes/"+"distrib_groupe_%s/distrib_rang.png"%nom_groupe)
    #plt.show()        

'''04/10/19'''
def test_combien_reste():
    dico_sim = {}
    #with open("script_copie_disque.sh", 'w') as fichier_script :
    with open("groupes_23S_homologues.pickle", "rb") as fichier_homologues :
            mon_depickler_1 = pickle.Unpickler(fichier_homologues)
            groupes_homologues = mon_depickler_1.load()
                  
            liste_manque = []
            liste_manque_complet = []
            for i in range(0, len(groupes_homologues)) :
                for j in range(i+1, len(groupes_homologues)) :
                    if "dico_sim_new_algo_23S_homologues_groupe_%s_groupe_%s_part_1.pickle"%(i,j) in os.listdir("/media/coline/Maxtor/Resultats/23S") : #or "dico_sim_new_algo_23S_homologues_groupe_%s_groupe_%s_part_1.pickle"%(i,j) in os.listdir("/media/coline/Maxtor/Resultats_dans_dico_sim") :
    #                     date = int(time.ctime(os.path.getmtime("/home/coline/Bureau/Resultats/Resultats/dico_sim_new_algo_23S_homologues_groupe_%s_groupe_%s_part_1.pickle"%(i,j))).split(" ")[2])
    #                     temps = time.ctime(os.path.getmtime("/home/coline/Bureau/Resultats/Resultats/dico_sim_new_algo_23S_homologues_groupe_%s_groupe_%s_part_1.pickle"%(i,j))).split(" ")[3]
    #                     temps = int(temps.split(":")[0])*60+int(temps.split(":")[1])
    #                     if not ((date < 18 ) or (date == 18 and temps <= 832)) :
                            print(i)
                            print(j)
                            taille = math.ceil((len(groupes_homologues[i])*len(groupes_homologues[j]))/8000)
                            ok = True
                            for k in range(1, taille+1) :
                                if "dico_sim_new_algo_23S_homologues_groupe_%s_groupe_%s_part_%s.pickle"%(i,j,k) not in os.listdir("/media/coline/Maxtor/Resultats/23S") :# and "dico_sim_new_algo_23S_homologues_groupe_%s_groupe_%s_part_%d.pickle"%(i,j, k) not in os.listdir("/media/coline/Maxtor/Resultats_dans_dico_sim") :
                                    if (i,j,k) not in liste_manque :
                                        liste_manque.append((i,j,k))
                                        ok = False
                    else :
                        liste_manque_complet.append((i,j))
                                         
                                 
#                             print(ok)
#                             if ok :
#                                 for k in range(1, taille+1) :
#                                     fichier_script.write("mv /home/coline/Bureau/Resultats/Resultats/dico_sim_new_algo_23S_homologues_groupe_%s_groupe_%s_part_%d.pickle /media/coline/Maxtor/Resultats/dico_sim_new_algo_23S_homologues_groupe_%s_groupe_%s_part_%d.pickle\n"%(i,j, k,i,j,k))
#                                     with open("/home/coline/Bureau/Resultats/Resultats/dico_sim_new_algo_23S_homologues_groupe_%s_groupe_%s_part_%d.pickle"%(i,j, k), 'rb') as fichier :
#                                         mon_depickler = pickle.Unpickler(fichier)
#                                         dico = mon_depickler.load()  
#                                         for cle in dico.keys() :
#                                             if cle not in dico_sim.keys() :
#                                                 dico_sim.update({cle : dico[cle]["sim"]})
#                                             else :
#                                                 print("gros rat")
#                                                 exit(0)
                          
    print(liste_manque)
    print(len(liste_manque))
    print(liste_manque_complet)
    print(len(liste_manque_complet))
    
    
    
def script_a_faire(numARN):

    with open(NEW_EXTENSION_PATH_TAILLE+"/liste_representant_%s.pickle"%numARN, 'rb') as fichier_vraiment_id :
        mon_depickler = pickle.Unpickler(fichier_vraiment_id)
        liste_representant = mon_depickler.load()
# 
#             print(liste_representant)
#             print(len(liste_representant))
#             #exit(0)
#     big_dico_sans_doublons = {}
#     for fic in os.listdir("/media/coline/Maxtor/Resultats/%s"%numARN) :
#         
#         with open("/media/coline/Maxtor/Resultats/%s/"%numARN+fic, 'rb') as fichier_dico_sim :
#             mon_depickler = pickle.Unpickler(fichier_dico_sim)
#             big_dico = mon_depickler.load() 
#             print(len(big_dico))  
#             
#             compteur = 0
#             for cle in big_dico.keys() :
#                 #print(compteur)
#                 #print((cle[0].split("_")[0], cle[0].split("_")[1]))
#                 if (cle[0].split("_")[0], int(cle[0].split("_")[1])) in liste_representant or (cle[1].split("_")[0], int(cle[1].split("_")[1])) in liste_representant :
#                     if isinstance(big_dico[cle], dict) :
#                         big_dico_sans_doublons.update({cle : big_dico[cle]["sim"]})
#                     else :
#                         big_dico_sans_doublons.update({cle : big_dico[cle]})
#                 compteur += 1
         
#     
#     with open("/media/coline/Maxtor/big_dico_sim_sans_doublons_%s.pickle"%numARN, 'wb') as fichier_dico_sim_sans_doublons :
#         mon_pickler = pickle.Pickler(fichier_dico_sim_sans_doublons)
#         mon_pickler.dump(big_dico_sans_doublons) 
    
#     with open("/media/coline/Maxtor/big_dico_sim_sans_doublons.pickle", 'rb') as fichier_dico_sim_sans_doublons :
#         mon_pickler = pickle.Unpickler(fichier_dico_sim_sans_doublons)
#         big_dico_sans_doublons =  mon_pickler.load() 
#     
#         #print(list(big_dico.keys())[0])
#         #print(len(big_dico))
#         print(len(big_dico_sans_doublons))   
    
#     with open("Nouvelles_donnees/groupe_25S.pickle", 'rb') as fichier :
#         mon_depickler = pickle.Unpickler(fichier)
#         groupes = mon_depickler.load()
     
        graphe_complet = nx.Graph()
           
        compteur = 1
        #for groupe in groupes_homologues :
        for elt in liste_representant :
            graphe_complet.add_node(compteur, nom=elt)
            compteur += 1
        print(graphe_complet.nodes.data())
        print(graphe_complet.number_of_nodes())
        #exit(0)
        compteur = 0
        
        #for cle in big_dico_sans_doublons.keys() :
#             print(compteur)
#             noeud1 = -1
#             noeud2 = -1
#             print(cle)
#            # print(len(big_dico_sans_doublons))
#             for noeud, data in graphe_complet.nodes(data=True) :
#                 
#                     #print(noeud)
#                     #print(data)
#                 if data["nom"] ==  (cle[0].split("_")[0], int(cle[0].split("_")[1])) :
#                         noeud1 = noeud
#                  
#                 elif data["nom"] ==  (cle[1].split("_")[0], int(cle[1].split("_")[1])) :
#                         noeud2 = noeud
#                 if noeud1 != -1 and noeud2 != -1 :
#                     break
#                 
#             if noeud1 != -1 and noeud2 != -1 :
# #                 print(cle)
# #                 print(liste_representant.index((cle[0].split("_")[0], int(cle[0].split("_")[1]))))
# #                 print(liste_representant.index((cle[1].split("_")[0], int(cle[1].split("_")[1]))))         
#                 graphe_complet.add_edge(noeud1, noeud2, sim=big_dico_sans_doublons[cle])
#             compteur +=1
           
        print(graphe_complet.nodes.data())
        print(graphe_complet.edges.data())
        print(graphe_complet.number_of_nodes())
        #print(graphe_complet.number_of_edges())
           
#         a_enlever = []
#         for noeud in graphe_complet.nodes() :
#             if len(graphe_complet[noeud]) == 0 :
#                 a_enlever.append(noeud) 
#                    
#         for elt in a_enlever :
#             graphe_complet.remove_node(elt)
               
        print(graphe_complet.number_of_nodes())
        print(graphe_complet.number_of_edges())
          
        with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_%s.pickle"%numARN, 'wb') as fichier_graphe_sim_sans_doublons :
            mon_pickler = pickle.Pickler(fichier_graphe_sim_sans_doublons)
            mon_pickler.dump(graphe_complet) 
        
#         for noeud in graphe_complet.nodes() :
#             print(len(graphe_complet[noeud]))

#         with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_25s.pickle", 'rb') as fichier_graphe_sim_sans_doublons :
#             mon_depickler = pickle.Unpickler(fichier_graphe_sim_sans_doublons)
#             graphe_complet = mon_depickler.load() 
            
            
#             print(graphe_complet.nodes[602]["nom"])
#             print(len(graphe_complet[602]))
#             print(graphe_complet[602])
            
            dico_manque = {}
#             for voisin in graphe_complet[602] :
#                 if voisin == 940 :
#                     print("trouve")
            #print(graphe_complet.nodes[940]["nom"])
            compter = 0
            for noeud1,data1 in graphe_complet.nodes(data=True) :
                for noeud2,data2 in graphe_complet.nodes(data=True) :
                    if (noeud1,noeud2) not in graphe_complet.edges() :
                        print(compter)
                        with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+data1["nom"][0]+"_"+str(data1["nom"][1])+".pickle", 'rb') as fichier1 :
                            mon_depickler_1 = pickle.Unpickler(fichier1)
                            graphe1 = mon_depickler_1.load()
                                    
                            with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+data2["nom"][0]+"_"+str(data2["nom"][1])+".pickle", 'rb') as fichier2 :
                                    mon_depickler = pickle.Unpickler(fichier2)
                                    graphe2 = mon_depickler.load()
                                    graphe_commun_max, sim_max = comparaison(graphe1, graphe2, "petit rat") 
                                    graphe_complet.add_edge(noeud1, noeud2, sim=sim_max)
                                    dico_manque.update({(data1["nom"][0]+"_"+str(data1["nom"][1]), data2["nom"][0]+"_"+str(data2["nom"][1])) : {"graphe" : graphe_commun_max, "sim" : sim_max} })
                        compter += 1
#                     print(len(graphe_complet[noeud]))
#                     print(graphe_complet[noeud])
#                          
#             print(compter)
#             print(graphe_complet.number_of_nodes())
#             print(graphe_complet.number_of_edges())
#             
# #             compteur =0
# #             for groupe in groupes_homologues :
# #                 if ('4v7s',20) in groupe and ('4tud', 16) in groupe :
# #                     print("meme_groupe")
# #                 if (('4v7s',20) in groupe and not ('4tud', 16) in groupe) or (not ('4v7s',20) in groupe and ('4tud', 16) in groupe) :
# #                     print(compteur)
# #                     print("pas meme groupe")
# #                 compteur += 1
# #                 
# #             if ('4v7s_20', '4tud_16') in big_dico.keys() or  ('4tud_16', '4v7s_20') in big_dico.keys() :
# #                 print("trouve")
        print(graphe_complet.number_of_nodes())
        print(graphe_complet.number_of_edges())         
        with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_avec_manque_%s.pickle"%numARN, 'wb') as fichier_graphe_sim_sans_doublons_avec_manque :
                mon_pickler = pickle.Pickler(fichier_graphe_sim_sans_doublons_avec_manque)
                mon_pickler.dump(graphe_complet) 
#                 
#                 print(graphe_complet.number_of_nodes())
#                 print(graphe_complet.number_of_edges())
#                 graphe_complet = nx.convert_node_labels_to_integers(graphe_complet, first_label=0, ordering='default', label_attribute=None)
#                 
#                 matrice = [[0] *graphe_complet.number_of_nodes() for _ in range(graphe_complet.number_of_nodes())]
#                 print(matrice)
#                 
#                 compteur = 0
#                 for i in range(graphe_complet.number_of_nodes()) :
#                     for j in range(graphe_complet.number_of_nodes()) :
#                         if i != j :
#                             if (i,j) in graphe_complet.edges() or (j,i) in graphe_complet.edges() :
#                                 matrice[i][j] = 1 - graphe_complet.edges[i,j]["sim"]
#                                 compteur += 1
#                             else :
#                                 matrice[i][j] = 50.0
#                         
#                 #print(matrice)
#                 #print(compteur)
#                 
#                 clustering = KMeans(n_clusters=10, n_init=50).fit(matrice)#init=init, n_init=1).fit(matrice)
#                     #print(kmeans.labels_)
#                 print(graphe_complet.nodes())
#                 tab_clustering = [['']]
#                 print(len(tab_clustering))
#                 print(len(clustering.labels_))
#                 print(clustering.labels_)
#                 compteur = 0
#                 for elt in clustering.labels_ :
#                     if elt != -1 :
#                         for _ in range(len(tab_clustering), elt+1) :
#                                 tab_clustering.append([''])
#                         print(elt)
#                         print(compteur)
#                         tab_clustering[elt].append(graphe_complet.nodes[compteur]["nom"])
#                     compteur += 1
#                 
#                 
#                 a_enlever = []
#                 for u,v,data in graphe_complet.edges(data=True) :
#                     if data["sim"] < 0.6 :
#                         a_enlever.append((u,v))
#         #     print(a_enlever)
#                 for elt in a_enlever :
#                     graphe_complet.remove_edge(elt[0], elt[1])
#                 print(graphe_complet.number_of_nodes())
#                 print(graphe_complet.number_of_edges())
#                 
#                 
#                 
#                 with open("/media/coline/Maxtor/fichier_csv_grandes_composantes_test_seuil_0.6.csv",'w') as fichier_csv:
#                     csvwriter = csv.writer(fichier_csv)
#                     csvwriter.writerow(["source", "target", "label"])
#                     
#                     seuil = 0.6
#                     for u,v,data in graphe_complet.edges(data=True) :
#                         if data["sim"] > seuil :
#                             csvwriter.writerow([u,v,round(data["sim"],2)])
#                         
#                         with open("/media/coline/Maxtor/fichier_csv_grandes_composantes_noeuds_test_seuil_0.6.csv",'w') as fichier_csv:
#                             csvwriter2 = csv.writer(fichier_csv)
#                             csvwriter2.writerow(["id", "label", "homologues", "clustering_kmeans"])
#                                 
#                             for noeud,data in graphe_complet.nodes(data=True) :
#                                 num_homologues = -1
#                                 compteur = 1
#                                 for groupe in groupes_homologues :
#                                     if data["nom"] in groupe :
#                                         num_homologues = compteur
#                                     compteur += 1
#                                 if num_homologues == -1 :
#                                     print("gros rat")
#                                  
#                                 num_clusters = -1   
#                                 for groupe in tab_clustering :
#                                     if data["nom"] in groupe :
#                                         num_clusters = compteur
#                                     compteur += 1
#                                 if num_clusters == -1 :
#                                     print("gros rat")
#                                 if data["nom"] in groupe :
#                                     numero_homologues = 0
#                                     for j in range(len(HOMOLOGUES)) :
#                                         
#                                         if data["nom"] in HOMOLOGUES[j] :
#                                             numero_homologues = j+1
#                                             
#                                             
#                                     if data["nom"] in groupe_base :
#                                         dans_groupe = 1
#                                     else :
#                                         dans_groupe = 0
                                #csvwriter2.writerow([noeud, data["nom"][0] + "_"+str(data["nom"][1]), num_homologues, num_clusters])
#                   #         with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_25s.pickle", 'rb') as fichier_graphe_sim_sans_doublons :
#             mon_depickler = pickle.Unpickler(fichier_graphe_sim_sans_doublons)
#             graphe_complet = mon_depickler.load() 
        with open("/media/coline/Maxtor/dico_sim_manque_%s.pickle"%numARN, 'wb') as fichier_sim_sans_doublons_avec_manque :
                mon_pickler_2 = pickle.Pickler(fichier_sim_sans_doublons_avec_manque)
                mon_pickler_2.dump(dico_manque) 


def diff_graphe(graphe1, graphe2):
    nb_diff = 0
    for noeud, data in graphe1.nodes(data=True) :
        if noeud not in graphe2.nodes() :
            print("pb1")
            nb_diff += 1
             
        else : 
            for cle in data.keys() :
                if data[cle] != graphe2.nodes[noeud][cle] :
                    print("pb2")
                    nb_diff += 1
                     
    for u,v,key,data in graphe1.edges(data=True, keys=True) :
        if (u,v) not in graphe2.edges() :
            print("pb3")
            nb_diff += 1
        else :
            if key in graphe2[u][v] :
                for cle in data.keys() :
                    if data[cle] != graphe2[u][v][key][cle] :
                        print("pb7")
                        print(u,v)
                        print(cle)
                        print(data[cle])
                        print(graphe2[u][v][key][cle])
                        nb_diff += 1
            else :
                print("pb10")
                nb_diff += 1
         
    for noeud, data in graphe2.nodes(data=True) :
        if noeud not in graphe1.nodes() :
            print("pb4")
            nb_diff += 1
             
        else : 
            for cle in data.keys() :
                if data[cle] != graphe1.nodes[noeud][cle] :
                    print("pb5")
                    nb_diff += 1
                     
    for u,v,key,data in graphe2.edges(data=True,keys=True) :
        if (u,v) not in graphe1.edges() :
            print("pb6") 
            print((u,v))
            nb_diff += 1   
        else :
            if key in graphe1[u][v] :
                for cle in data.keys() :
                    if data[cle] != graphe1[u][v][key][cle] :
                        print("pb8") 
                        nb_diff += 1
            else :
                print("pb9")
                nb_diff += 1  
    
    return nb_diff               
                    
def non_can_consecutif(graphe): 
    for noeud in graphe.nodes() :
        voisin_b53 = -1
        voisin_non_cov = []
        for voisin in graphe[noeud] :  
            for edge in graphe[noeud][voisin] :
                if graphe[noeud][voisin][edge]["label"] == 'B53' :
                    voisin_b53 = voisin
                else :
                    voisin_non_cov.append(voisin) 
        
        if voisin_b53 == -1 :
            for voisin in graphe.predecessors(noeud) :  
                for edge in graphe[voisin][noeud] :
                    if graphe[voisin][noeud][edge]["label"] == 'B53' :
                        voisin_b53 = voisin
                    else :
                        voisin_non_cov.append(voisin) 
        
        if voisin_b53 in voisin_non_cov :
            return True
    return False

if __name__ == '__main__':
    
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load()
    
        with open("dico_new_060420_avec_graphes_en_plus.pickle", 'rb') as fichier_graphe :
            mon_depickler = pickle.Unpickler(fichier_graphe)
            dico_graphe = mon_depickler.load() 
            
            with open("dico_algo_heuristique_grands_graphes_taille_4.pickle", 'rb') as fichier_graphe :
                mon_depickler = pickle.Unpickler(fichier_graphe)
                dico_graphe_non_contractes = mon_depickler.load() 
                
                compter_diff = 0
                compter_tot_sup_07 = 0
                for elt in dico_graphe.keys() :
                    
                    if elt[0][0] != '6hrm' and elt[1][0] != '6hrm' :
                        nom1 = "fichier_%s_%s_taille_4.pdb"%(elt[0][0], elt[0][1])
                        nom2 = "fichier_%s_%s_taille_4.pdb"%(elt[1][0], elt[1][1])
                        if (nom1, nom2) in rmsd.keys() :
                            val_rmsd = rmsd[(nom1, nom2)]
                        else :
                            val_rmsd = rmsd[(nom2, nom1)]
                        if elt in dico_graphe_non_contractes.keys() :
                            if(dico_graphe_non_contractes[elt]["sim"]  > 0.75 or dico_graphe[elt]["sim"] > 0.75) and val_rmsd != None and val_rmsd <= 2.5:
                                compter_tot_sup_07 += 1
                            
                            if dico_graphe[elt]["sim"] > dico_graphe_non_contractes[elt]["sim"] :
                                
                                if (dico_graphe_non_contractes[elt]["sim"]  > 0.75 or dico_graphe[elt]["sim"] > 0.75) and val_rmsd != None and val_rmsd <= 2.5:
                                    print(elt)
                                    print(dico_graphe[elt]["sim"])
                                    print(dico_graphe_non_contractes[elt]["sim"])
                                    compter_diff += 1
                        else :
                            if (dico_graphe_non_contractes[(elt[1], elt[0])]["sim"]  > 0.75 or dico_graphe[elt]["sim"] > 0.75) and val_rmsd != None and val_rmsd <= 2.5:
                                compter_tot_sup_07 += 1
                            if dico_graphe[elt]["sim"] > dico_graphe_non_contractes[(elt[1], elt[0])]["sim"] :
                                
                                if (dico_graphe_non_contractes[(elt[1], elt[0])]["sim"]  > 0.75 or dico_graphe[elt]["sim"] > 0.75) and val_rmsd != None and val_rmsd <= 2.5:
                                    print(elt)
                                    print(dico_graphe[elt]["sim"])
                                    print(dico_graphe_non_contractes[(elt[1], elt[0])]["sim"])
                                    compter_diff += 1
                print(compter_diff)
                print(compter_tot_sup_07)
            exit()
            print(len(dico_graphe))
            
            compteur = 0
            for cle in dico_graphe.keys() :
                if len(dico_graphe[cle]["graphe_en_plus"]) > 1 :
                    print(len(dico_graphe[cle]["graphe_en_plus"]))
                    #print(cle)
                    #diff_graphe(dico_graphe[cle]["graphe_en_plus"][0], dico_graphe[cle]["graphe_en_plus"][1])
                    compteur += 1
                    #exit()
                if cle == (('5dm6', 3), ('4ybb', 19)) :
                    print(cle)
                    print(len(dico_graphe[cle]["graphe_en_plus"]))
            
            print(compteur)
            exit()
        
        
        with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_graphe_1 :
            mon_depickler = pickle.Unpickler(fichier_graphe_1)
            dico_graphe_1 = mon_depickler.load() 
            
            compteur = 0
            for cle in dico_graphe.keys() :
                if cle[0][0] != '6hrm' and cle[1][0] != '6hrm' :
                    if cle in dico_graphe_1.keys() :
                        if dico_graphe[cle]["sim"] != dico_graphe_1[cle]["sim"] :
                            print("raaa")
                            compteur +=1
                    else :
                        if dico_graphe[cle]["sim"] != dico_graphe_1[(cle[1], cle[0])]["sim"] :
                            print("raaa")
                            compteur +=1
            print(compteur)
            exit()
    
    with open("Graphs/1mms.pickle", "rb") as fichier_graphe :
            mon_depickler = pickle.Unpickler(fichier_graphe)
            graphe1 = mon_depickler.load()
            
    #           
            with open("Nouvelles_donnees/fichier_1mms_1_4.pickle", "rb") as fichier_extension :
                mon_depickler = pickle.Unpickler(fichier_extension)
                extension1 = mon_depickler.load()
                
                print(graphe1[(extension1.nodes[12]["num_ch"], extension1.nodes[12]["position"][0])])
                print(graphe1[(extension1.nodes[12]["num_ch"], extension1.nodes[12]["position"][0]+1)])
                print(graphe1[(extension1.nodes[12]["num_ch"], extension1.nodes[12]["position"][0]+2)])
                
    exit()
    
    types_arn = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16S_arnm", "arnt_16S"]

    
    liste_tout = []
    with open("resolutions.pickle", 'rb') as fichier_pickle :
                mon_depickler = pickle.Unpickler(fichier_pickle)
                resolutions = mon_depickler.load()
    for elt in types_arn :      
                with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
                        mon_depickler = pickle.Unpickler(fichier_sortie)
                        liste_a_garder_noms = mon_depickler.load()    
                                    
                        for element in liste_a_garder_noms :
                            if element == ('4w2f', 16) :
                                print(element)
                                print(elt)
                            if resolutions[element[0]] <= 3.0 :
                                if element in liste_tout :
                                    print(element)
                                    print(elt)
                                if element not in liste_tout :
                                    liste_tout.append(element)
                
    print(liste_tout)  
    
    liste_chgt = []
    for nom in liste_tout :
        with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(nom[0], nom[1]), "rb") as fichier_extension :
            mon_depickler = pickle.Unpickler(fichier_extension)
            extension1 = mon_depickler.load()
        
        with open("Nouvelles_donnees/fichier_%s_%s_4.pickle"%(nom[0], nom[1]), "rb") as fichier_extension :
            mon_depickler = pickle.Unpickler(fichier_extension)
            extension2 = mon_depickler.load() 
            
            print(nom)
            nb_diff = diff_graphe(extension1, extension2)
            if nb_diff > 0 : 
                liste_chgt.append(nom)
            
            #exit()
            
    print(len(liste_chgt))
    print(liste_chgt)
    #exit()
            
            
    
#     liste_change = []
#     with open("groupes_homologues_60_2.pickle", 'rb') as fichier_pickle :
#         mon_depickler = pickle.Unpickler(fichier_pickle)
#         groupes_homologues = mon_depickler.load()
#          
#         with open("/media/coline/Maxtor/dico_new_100320_sim_par_branche_0.65.pickle", 'rb') as fichier_graphe :    
#             mon_depickler = pickle.Unpickler(fichier_graphe)
#             dico_graphe_1 = mon_depickler.load() 
#              
#             with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_graphe_2 :    
#                 mon_depickler_2 = pickle.Unpickler(fichier_graphe_2)
#                 dico_graphe_2 = mon_depickler_2.load() 
#                  
#                 compter_homol = 0
#                 compte = 0
#                 for cle in dico_graphe_1.keys() :
#                     if dico_graphe_1[cle]["sim"] != dico_graphe_2[cle]["sim"] :
#                         if cle[0][0] == '6hrm' or cle[1][0] == '6hrm' :
#                             compte += 1
#                         for groupe in groupes_homologues :
#                             if cle[0] in groupe and cle[1] in groupe :
#                                 compter_homol += 1
#                         liste_change.append(cle)
#                 print(len(liste_change))
#                 print(liste_change)
#                 print(compter_homol)
#                 print(compte)
#     exit()
        

    
    with open("all_aminor.pickle", 'rb') as fichier_all_aminor :
        mon_depickler = pickle.Unpickler(fichier_all_aminor)
        all_aminor = mon_depickler.load()
        
        for cle in all_aminor.keys() :
            for graphe in all_aminor[cle] :
                if abs(int(list(graphe.nodes())[0][1]) - int(list(graphe.nodes())[1][1])) < 20 :
                    print(cle)
                    print(graphe.nodes())
        #exit()
#      
    with open("Graphs/5ool.pickle", "rb") as fichier_graphe :
        mon_depickler = pickle.Unpickler(fichier_graphe)
        graphe1 = mon_depickler.load()
        
        for noeud in graphe1.nodes() :
            if noeud[0] == 'AV' :
                print(noeud)
                
        with open("test_pour_dessin_noeuds.csv", "w", newline="") as fichier_noeuds :
            csvwriter1 = csv.writer(fichier_noeuds)
            csvwriter1.writerow(["id", "label"])
            liste_noeuds = []
            for noeud, data in graphe1.nodes(data=True) :
                if noeud[0] == 'A' and noeud[1] >= 456 and noeud[1] <= 485 :
                    liste_noeuds.append(noeud)
                    csvwriter1.writerow([noeud, data["nt"]])
                
        with open("test_pour_dessin.csv", "w", newline="") as fichier_aretes :
            csvwriter2 = csv.writer(fichier_aretes)
            csvwriter2.writerow(["source", "target", "label"])
            
            for u,v,data in graphe1.edges(data=True) :
                if u in liste_noeuds or v in liste_noeuds  :
                    if data["label"] != 'B53' :
                        print("rapala")
                    csvwriter2.writerow([u,v,data["label"]])
        
                
#         print(graphe1.nodes())
#         print(graphe1.edges())
        #exit()
        
        with open("Graphs/5ngm.pickle", "rb") as fichier_graphe :
            mon_depickler = pickle.Unpickler(fichier_graphe)
            graphe1 = mon_depickler.load()
            
    #           
            with open("Nouvelles_donnees/fichier_5ngm_3_4.pickle", "rb") as fichier_extension :
                mon_depickler = pickle.Unpickler(fichier_extension)
                extension1 = mon_depickler.load()
            
            print(graphe1[(extension1.nodes[14]["num_ch"], extension1.nodes[14]["position"][0])])
            print(graphe1[(extension1.nodes[23]["num_ch"], extension1.nodes[23]["position"][0])])
            
        exit()
            
#             print(extension1.nodes.data())
#             print(extension1.edges.data())
#             for noeud, data in extension1.nodes(data=True) :
#                 if data["type"] != -1 :
#                     for i in range(data["position"][0], data["position"][1]+1) :
#                         print(graphe1.nodes[(data["num_ch"], i)])
#     exit()
#                  
#                 with open(NEW_EXTENSION_PATH_TAILLE+"groupes_identiques_23S.pickle", 'rb') as fichier_identiques :
#                     mon_depickler = pickle.Unpickler(fichier_identiques)
#                     groupes_identiques = mon_depickler.load() 
#                     
#                     for groupe in groupes_identiques :
#                         ok = False
#                         for elt in groupe :
#                             if elt[0] == '5dm6' :
#                                 ok = True
#                         if ok : 
#                             for elt in groupe :
#                                 if elt[0] == '4ioa' :
#                                     print(groupe)
#                             #print(groupe)
#                  
#                 with open("resolutions.pickle", 'rb') as fichier_pickle :
#                     mon_depickler = pickle.Unpickler(fichier_pickle)
#                     resolutions = mon_depickler.load() 
#                     
#                     print(resolutions['5dm6'])
#                     print(resolutions['4ioa'])
#                     
#                     
#                     for cle in all_aminor.keys() :
#                             compteur = 1
#                             if cle == "4ioa" or cle == '5dm6' :
#                                 print("rapala")
#                                 for graphe in all_aminor[cle] :
#                                     #if compteur == 16 :
#                                     print(graphe.nodes.data())
#                                     compteur += 1
#                              
#                 print(extension1.nodes[1]["num_ch"], extension1.nodes[1]["position"][0])
#                 print(extension1.nodes[2]["num_ch"], extension1.nodes[2]["position"][0])
#                 print(graphe1[(extension1.nodes[1]["num_ch"], extension1.nodes[1]["position"][0])])
#                 exit()
    #liste = [(('4y4o', 47), ('1vq8', 19)), (('4y4o', 47), ('1yhq', 24)), (('4y4o', 47), ('2qex', 19)), (('4y4o', 47), ('4u4r', 13)), (('4ybb', 13), ('1vqo', 8)), (('4w2f', 23), ('4u4r', 36)), (('1vqo', 8), ('5dm6', 1)), (('1vqo', 8), ('1mms', 1)), (('1vq8', 19), ('6eri', 1)), (('1yhq', 24), ('6eri', 1)), (('6eri', 1), ('2qex', 19))]
    #liste = [('3cc2', 17), ('6eri', 11), ('6ek0', 4), ('4ybb', 54), ('6hma', 8), ('1vq8', 21), ('5dm6', 9), ('5afi', 17), ('4y4o', 28), ('4u4r', 4)]
    
#     with open("dico_algo_heuristique_new_v_e8f97fe.pickle", 'rb') as fichier_graphe :    
#         mon_depickler = pickle.Unpickler(fichier_graphe)
#         dico_graphe_1 = mon_depickler.load()
#         
#     with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_graphe :    
#         mon_depickler = pickle.Unpickler(fichier_graphe)
#         dico_graphe_2 = mon_depickler.load()  
#         compteur = 0
#         compteur_test = 0
#         for elt in dico_graphe_1.keys() : 
#             if elt in dico_graphe_2.keys() :
#                 if dico_graphe_1[elt]["sim"] != dico_graphe_2[elt]["sim"] :
#                     print(elt)
#                     print(dico_graphe_1[elt]["sim"])
#                     print(dico_graphe_2[elt]["sim"])
#                     compteur += 1
                    
#                     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt[0][0]+"_"+str(elt[0][1])+"_2.pickle", 'rb') as fichier1 :
#                         mon_depickler_1 = pickle.Unpickler(fichier1)
#                         graphe1 = mon_depickler_1.load()
# 
#                                        
#                     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt[1][0]+"_"+str(elt[1][1])+"_2.pickle", 'rb') as fichier2 :
#                         mon_depickler = pickle.Unpickler(fichier2)
#                         graphe2 = mon_depickler.load()
#                     
#                     trouve_noeud_pas_bon = False
#                     for noeud in dico_graphe_1[elt]["graphe"] :
#                         if test_compatibilite(dico_graphe_1[elt]["graphe"], noeud, graphe1, graphe2) == False :
#                             print("rapaaaaalaaaa")
#                             trouve_noeud_pas_bon = True
#                     if trouve_noeud_pas_bon : 
#                         compteur_test += 1
#                 if dico_graphe_1[elt]["sim"] != dico_graphe_2[elt]["sim"] :
#                     print(elt)
#                     compteur += 1
#             else :
#                 print("rapala")
#         print(compteur)
#         print(compteur_test)
#         
#     exit()
    
#     liste = []
#     os.makedirs("/media/coline/Maxtor/Fichiers_mmcif/test", exist_ok = True)
#     with open("script_python.py", 'w') as fichier_python :
#         for fichier in os.listdir("/media/coline/Maxtor/Fichiers_mmcif/taille_4") :
#             for elt in liste :
#                 if str(elt[0]) + "_"+ str(elt[1]) in fichier :
#                     fichier_python.write("cp /media/coline/Maxtor/Fichiers_mmcif/taille_4/%s /media/coline/Maxtor/Fichiers_mmcif/test/%s\n"%(fichier, fichier))
            
    
    
    
    liste = [('4ybb', 12), ('6ek0', 10)]
    with open("/home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/Fichiers_mmcif/script_pymol.pml", 'w') as fichier_pymol :
        tabs = []
        chaines = []
        for elt in liste :
            with open("Graphs/%s.pickle"%elt[0], "rb") as fichier_graphe :
                mon_depickler = pickle.Unpickler(fichier_graphe)
                graphe1 = mon_depickler.load()
                 
    #         with open("Graphs/%s.pickle"%elt[1][0], "rb") as fichier_graphe :
    #             mon_depickler = pickle.Unpickler(fichier_graphe)
    #             graphe2 = mon_depickler.load()
             
            with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(elt[0], elt[1]), "rb") as fichier_extension :
                mon_depickler = pickle.Unpickler(fichier_extension)
                extension1 = mon_depickler.load()
                
                
                 
    #         with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(elt[1][0], elt[1][1]), "rb") as fichier_extension :
    #             mon_depickler = pickle.Unpickler(fichier_extension)
    #             extension2 = mon_depickler.load()
                 
                tab1 = []
                #tab2 = []
                for i in range(1,6) :
                    print((extension1.nodes[i]["num_ch"],extension1.nodes[i]["position"][0]))
                    print(graphe1.nodes[(extension1.nodes[i]["num_ch"],extension1.nodes[i]["position"][0])])
                    tab1.append(graphe1.nodes[(extension1.nodes[i]["num_ch"],extension1.nodes[i]["position"][0])]["fr3d"])
                    #tab2.append(graphe2.nodes[(extension2.nodes[i]["num_ch"],extension2.nodes[i]["position"][0])]["fr3d"])
                 
                print(elt)
                print(tab1)
                chaines.append(extension1.nodes[1]["num_ch"])
                tabs.append(tab1)
                
                if liste[len(chaines)-1][0].upper()+"_chaine_"+chaines[len(chaines)-1]+".pdb" not in os.listdir("/home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/Fichiers_mmcif/") :
                    fichier_pymol.write("load /home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/Fichiers_mmcif/%s\n"%(liste[len(chaines)-1][0].upper()+".cif"))
                    fichier_pymol.write("select /%s/%s\n"%(liste[len(chaines)-1][0].upper(), chaines[len(chaines)-1]))
                    fichier_pymol.write("create %s, sele\n"%(liste[len(chaines)-1][0].upper()+"_chaine_"+chaines[len(chaines)-1]))
                    fichier_pymol.write("save %s\n"%(liste[len(chaines)-1][0].upper()+"_chaine_"+chaines[len(chaines)-1]+".pdb"))
                else :
                    fichier_pymol.write("load /home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/Fichiers_mmcif/%s\n"%(liste[len(chaines)-1][0].upper()+"_chaine_"+chaines[len(chaines)-1]+".pdb"))
                #fichier_pymol.write("select sele%d, %s and (resi %d-%d or resi %d-%d or resi %d-%d)\n"%(len(chaines)+2, liste[len(chaines)-1][0].upper()+"_chaine_"+chaines[len(chaines)-1], int(tab1[1])-50, int(tab1[1]+50), int(tab1[2])-50, int(tab1[2])+50, int(tab1[4])-50, int(tab1[4])+50))
                
                
                #print(tab2)
        fichier_pymol.write("align %s, %s\n"%(liste[0][0].upper()+"_chaine_"+chaines[0], liste[1][0].upper()+"_chaine_"+chaines[1]))
        fichier_pymol.write("select sele1, %s and (resi %s-%s or resi %s-%s or resi %s) \n"%(liste[0][0].upper()+"_chaine_"+chaines[0],tabs[0][2],tabs[0][0], tabs[0][1], tabs[0][3], tabs[0][4] ))
        fichier_pymol.write("select sele2, %s and (resi %s-%s or resi %s-%s or resi %s) \n"%(liste[1][0].upper()+"_chaine_"+chaines[1],tabs[1][2],tabs[1][0], tabs[1][1], tabs[1][3], tabs[1][4] ))
        
        nb = 50       
        fichier_pymol.write("select sele3, %s and (resi %d-%d or resi %d-%d or resi %d-%d)\n"%(liste[0][0].upper()+"_chaine_"+chaines[0], int(tabs[0][1])-nb, int(tabs[0][1])+nb, int(tabs[0][2])-nb, int(tabs[0][2])+nb, int(tabs[0][4])-nb, int(tabs[0][4])+nb))
        fichier_pymol.write("create obj01, sele3\n")
        fichier_pymol.write("select sele4, %s and (resi %d-%d or resi %d-%d or resi %d-%d)\n"%(liste[1][0].upper()+"_chaine_"+chaines[1], int(tabs[1][1])-nb, int(tabs[1][1])+nb, int(tabs[1][2])-nb, int(tabs[1][2])+nb, int(tabs[1][4])-nb, int(tabs[1][4])+nb))
        fichier_pymol.write("create obj02, sele4\n")
        fichier_pymol.write("save /home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/Fichiers_mmcif/%s_%s_%s_%s.pse\n"%(liste[0][0].upper(), liste[0][1], liste[1][0].upper(), liste[1][1]))

        
#     with open("Graphs/4y4o.pickle", "rb") as fichier_graphe :
#         mon_depickler = pickle.Unpickler(fichier_graphe)
#         graphe = mon_depickler.load()
#          
# #         for j in range(1, 23) :
# #             if "fichier_1vq8_%d_2.pickle"%j in os.listdir("Nouvelles_donnees") :
#         with open("Nouvelles_donnees/fichier_4y4o_47_2.pickle", "rb") as fichier_extension :
#                     mon_depickler = pickle.Unpickler(fichier_extension)
#                     extension = mon_depickler.load()
# #                     print("fichier_1vq8_%d_2.pickle"%j)
# #                     print(graphe.nodes[(extension.nodes[1]["num_ch"],extension.nodes[1]["position"][0])])
#                     for i in range(1,6) :
#                         print((extension.nodes[i]["num_ch"],extension.nodes[i]["position"][0]))
#                         print(graphe.nodes[(extension.nodes[i]["num_ch"],extension.nodes[i]["position"][0])])
# 
#             
#     print("\n")
#     with open("Graphs/4u4r.pickle", "rb") as fichier_graphe :
#         mon_depickler = pickle.Unpickler(fichier_graphe)
#         graphe = mon_depickler.load()
#          
#          
# #         for j in range(1, 54) :
# #             if "fichier_4y4o_%d_2.pickle"%j in os.listdir("Nouvelles_donnees") :
#         with open("Nouvelles_donnees/fichier_4u4r_13_2.pickle", "rb") as fichier_extension :
#                     mon_depickler = pickle.Unpickler(fichier_extension)
#                     extension = mon_depickler.load()
# #                     print("fichier_4y4o_%d_2.pickle"%j)
# #                     print(graphe.nodes[(extension.nodes[1]["num_ch"],extension.nodes[1]["position"][0])])
#                     
#                     for i in range(1,6) :
#                         print((extension.nodes[i]["num_ch"],extension.nodes[i]["position"][0]))
#                         print(graphe.nodes[(extension.nodes[i]["num_ch"],extension.nodes[i]["position"][0])])

            
        
#     liste = []
#     tps1 = time.time()
#     compteur_idem = 0
#     compteur = 0
#     compteur_sans_non_can_consec = 0
#     with open("dico_algo_heuristique.pickle", "rb") as fichier_dico_graphe_algo_exact :
#         mon_depickler = pickle.Unpickler(fichier_dico_graphe_algo_exact)
#         dico_graphe_sim_d = mon_depickler.load()
#         
#         print(dico_graphe_sim_d[(('3ccm', 21), ('3t1y', 8))]["graphe"].nodes())
#         print(dico_graphe_sim_d[(('3ccm', 21), ('3t1y', 8))]["graphe"].edges())


        
#     with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_sans_doublons.pickle", "rb") as fichier_dico_graphe_algo_exact_1 :
#         mon_depickler_2 = pickle.Unpickler(fichier_dico_graphe_algo_exact_1)
#         dico_graphe_sim_p = mon_depickler_2.load()
#             
#     types_arn = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16S_arnm", "arnt_16S"]
#              
#     with open("resolutions.pickle", 'rb') as fichier_pickle :
#         mon_depickler = pickle.Unpickler(fichier_pickle)
#         resolutions = mon_depickler.load()
#                            
#         proportion_inf_3 = 0
#         total = 0
#                            
#         liste_a_garder = []
#                            
#         for typ in types_arn :
#             with open("Nouvelles_donnees/liste_representant_%s.pickle"%typ, 'rb') as fichier_repr :
#                 mon_depickler = pickle.Unpickler(fichier_repr)
#                 liste_representant = mon_depickler.load() 
#                                    
#                 #print(resolutions)
#                 for elt in liste_representant :
#                     if resolutions[elt[0]] <= 3.0 :
#                         if elt == ('4w2f', 16) :
#                             print(typ)
#                         if elt not in liste_a_garder :
#                             liste_a_garder.append(elt)
#                   
# #     with open("occ_multi_chaine.pickle", 'rb') as fichier_multi_chaines :
# #         mon_depickler = pickle.Unpickler(fichier_multi_chaines)
# #         liste_plusieurs_chaines = mon_depickler.load() 
#              
#         print(len(liste_a_garder))
#              
#         for i in range(len(liste_a_garder)) :
#             #if liste_a_garder[i] == ('3g78', 2) :
#                 for j in range(i+1, len(liste_a_garder)) :
#                     
#                        
#                         if (liste_a_garder[i], liste_a_garder[j]) in dico_graphe_sim_p.keys() :
#                             sim1 = dico_graphe_sim_p[(liste_a_garder[i], liste_a_garder[j])]["sim"]
#                             dico_graphe1 = dico_graphe_sim_p[(liste_a_garder[i], liste_a_garder[j])]["graphe"]
#                                    
#                         elif (liste_a_garder[j], liste_a_garder[i]) in dico_graphe_sim_p.keys() :
#                             sim1 = dico_graphe_sim_p[(liste_a_garder[j], liste_a_garder[i])]["sim"]
#                             dico_graphe1 = dico_graphe_sim_p[(liste_a_garder[j], liste_a_garder[i])]["graphe"]
#                         else :
#                             print("dommage")
#                             
#                         if (liste_a_garder[i], liste_a_garder[j]) in dico_graphe_sim_d.keys() :
#                             sim2 = dico_graphe_sim_d[(liste_a_garder[i], liste_a_garder[j])]["sim"]
#                             dico_graphe2 = dico_graphe_sim_d[(liste_a_garder[i], liste_a_garder[j])]["graphe"]
#                                    
#                         elif (liste_a_garder[j], liste_a_garder[i]) in dico_graphe_sim_d.keys() :
#                             sim2 = dico_graphe_sim_d[(liste_a_garder[j], liste_a_garder[i])]["sim"]
#                             dico_graphe2 = dico_graphe_sim_d[(liste_a_garder[j], liste_a_garder[i])]["graphe"]
#                         else :
#                             print("dommage") 
#                         
#                         if sim1 != sim2 :
#                             compteur += 1
#                             print(liste_a_garder[i], liste_a_garder[j])
#                             diff_graphe(dico_graphe1, dico_graphe2)
#                             if not (non_can_consecutif(dico_graphe1) and non_can_consecutif(dico_graphe2)) :
#                                 print("petit rat")
#                                 compteur_sans_non_can_consec += 1
#                             
#                             
#     print(compteur)
#     print(compteur_sans_non_can_consec)
#     
#     
#     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_4ybb_6_2.pickle", 'rb') as fichier_graphe1 :
#         mon_depickler_g1 = pickle.Unpickler(fichier_graphe1)
#         graphe1 = mon_depickler_g1.load()
#     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_4u3u_2_2.pickle"%(), 'rb') as fichier_graphe2 :
#         mon_depickler_g2 = pickle.Unpickler(fichier_graphe2)
#         graphe2 = mon_depickler_g2.load()
#          
#         print(graphe1.nodes.data())
#         print(graphe2.nodes.data())
#         print(graphe1.edges.data())
#         print(graphe2.edges.data())
# #         
#     with open("Graphs/4ybb.pickle", 'rb') as fichier_graphe :
#         mon_depickler = pickle.Unpickler(fichier_graphe)
#         graphe = mon_depickler.load()
#         #print(graphe.nodes())
#         print(graphe[('DA', 1383)])
#         
#     with open("/media/coline/Maxtor/dico_new.pickle", 'rb') as fichier_graphe :
#                     mon_depickler = pickle.Unpickler(fichier_graphe)
#                     dico_graphe = mon_depickler.load() 
#                     
#                     print(dico_graphe[(('5fk2', 1), ('5fk4', 1))]["graphe"].nodes.data())
#                     #print(dico_graphe[(('4u4r', 12), ('4u4r', 36))]["graphe"].edges.data())
#         
#         
#     types_arn = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16s_arnm", "arnt_16s"]
#     liste_tout = []
#     for elt in types_arn :
#             with open("resolutions.pickle", 'rb') as fichier_pickle :
#                 mon_depickler = pickle.Unpickler(fichier_pickle)
#                 resolutions = mon_depickler.load()
#                    
#                 with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
#                         mon_depickler = pickle.Unpickler(fichier_sortie)
#                         liste_a_garder_noms = mon_depickler.load()    
#                              
#                         for element in liste_a_garder_noms :
#                             if resolutions[element[0]] <= 3 :
#                                 liste_tout.append(element)
#          
#     print(liste_tout)  
#     print(len(liste_tout))
#     
#     liste_chgt = []
#     for elt in liste_tout :
#         with open("Nouvelles_donnees/fichier_%s_%d.pickle"%(elt[0],elt[1]), "rb") as fichier_1 :
#                 mon_depickler = pickle.Unpickler(fichier_1)
#                 graphe1 = mon_depickler.load()
#                 
#                 with open("Nouvelles_donnees/fichier_%s_%d_2.pickle"%(elt[0],elt[1]), "rb") as fichier_2 :
#                     mon_depickler = pickle.Unpickler(fichier_2)
#                     graphe2 = mon_depickler.load()
#                     
#                     for noeud, data in graphe1.nodes(data=True) :
#                         if noeud not in graphe2.nodes() :
#                             print(elt)
#                             print("pb1")
#                             if elt == ('4y4o', 23) :
#                                 print(graphe1.nodes.data())
#                                 print(graphe2.nodes.data())
#                                 print(graphe1.edges.data())
#                                 print(graphe2.edges.data())
#                             if elt not in liste_chgt :
#                                 liste_chgt.append(elt)
#                             
#                         else : 
#                             for cle in data.keys() :
#                                 if data[cle] != graphe2.nodes[noeud][cle] :
#                                     print(elt)
#                                     print(data[cle])
#                                     print(graphe2.nodes[noeud][cle])
#                                     print("pb2")
#                                     if elt not in liste_chgt :
#                                         liste_chgt.append(elt)
#                                     
#                     for u,v,key,data in graphe1.edges(data=True, keys=True) :
#                         if (u,v) not in graphe2.edges() :
#                             print(elt)
#                             print("pb3")
#                             if elt not in liste_chgt :
#                                 liste_chgt.append(elt)
#                         else :
#                             if key in graphe2[u][v] :
#                                 for cle in data.keys() :
#                                     if data[cle] != graphe2[u][v][key][cle] :
#                                         print(elt)
#                                         print("pb7")
#                                         if elt not in liste_chgt :
#                                             liste_chgt.append(elt)
#                             else :
#                                 print(elt)
#                                 print("pb10")
#                                 if elt not in liste_chgt :
#                                     liste_chgt.append(elt)
#                         
#                     for noeud, data in graphe2.nodes(data=True) :
#                         if noeud not in graphe1.nodes() :
#                             print(elt)
#                             print("pb4")
#                             if elt not in liste_chgt :
#                                 liste_chgt.append(elt)
#                             
#                         else : 
#                             for cle in data.keys() :
#                                 if data[cle] != graphe1.nodes[noeud][cle] :
#                                     print(elt)
#                                     print("pb5")
#                                     if elt not in liste_chgt :
#                                         liste_chgt.append(elt)
#                                     
#                     for u,v,key,data in graphe2.edges(data=True,keys=True) :
#                         if (u,v) not in graphe1.edges() :
#                             print(elt)
#                             print("pb6") 
#                             if elt not in liste_chgt :
#                                 liste_chgt.append(elt)     
#                         else :
#                             if key in graphe1[u][v] :
#                                 for cle in data.keys() :
#                                     if data[cle] != graphe1[u][v][key][cle] :
#                                         print(elt)
#                                         print("pb8") 
#                                         if elt not in liste_chgt :
#                                             liste_chgt.append(elt)
#                             else :
#                                 print(elt)
#                                 print("pb9")
#                                 if elt not in liste_chgt :
#                                     liste_chgt.append(elt)
#         
#     print(liste_chgt)
#     print(len(liste_chgt))
    
    
    
#     with open("Graphs/5ibb.pickle", "rb") as fichier_graphe_tot :
#         mon_depickler = pickle.Unpickler(fichier_graphe_tot)
#         graphe_tot = mon_depickler.load()
#          
#         print(graphe_tot[('14',)])
    
#     liste_pbs = []
#     for elt in liste_tout :
#         with open("Nouvelles_donnees/fichier_%s_%d.pickle"%(elt[0], elt[1]), "rb") as fichier_graphe :
#             mon_depickler = pickle.Unpickler(fichier_graphe)
#             graphe = mon_depickler.load()
#             
#             
#             for noeud in graphe.nodes() :
#                 compteur_b53 = 0
#                 for voisin in graphe[noeud] :
#                     for edge in graphe[noeud][voisin] :
#                         if graphe[noeud][voisin][edge]["label"] == "B53" :
#                             compteur_b53 += 1
# #                 for voisin in graphe.predecessors(noeud) :
# #                     for edge in graphe[voisin][noeud] :
# #                         if graphe[voisin][noeud][edge]["label"] == "B53" :
# #                             compteur_b53 += 1
#                 if compteur_b53 > 1 :
#                     print(elt)
#                     liste_pbs.append(elt)
#                     
#                     
#     for elt in liste_pbs :
#         with open("Graphs/%s.pickle"%elt[0], "rb") as fichier_graphe_tot :
#             mon_depickler = pickle.Unpickler(fichier_graphe_tot)
#             graphe_tot = mon_depickler.load()
#             print(elt)
#             
#             with open("Nouvelles_donnees/fichier_%s_%d.pickle"%(elt[0], elt[1]), "rb") as fichier_graphe :
#                 mon_depickler = pickle.Unpickler(fichier_graphe)
#                 graphe = mon_depickler.load()
#                 
#                 #print(graphe_tot.nodes.data())
#                 
#                 print(graphe[4])
#                 print(graphe[15])
#                 
#                 
# #                 for i in range(1,5) :
# #                     print(graphe.nodes[i]["position"])
# #                     for j in range(5) :
# #                         if i == 1 or i == 4 :
# #                                 if graphe.nodes[i]["position"][0]+j <= graphe_tot.number_of_nodes() :
# #                                     #print(graphe.nodes[i])
# #                                     print(graphe_tot.nodes[(graphe.nodes[i]["num_ch"], graphe.nodes[i]["position"][0]+j)])
# #                         else :
# #                                 if graphe.nodes[i]["position"][0]-j > 0 :
# #                                     print(graphe_tot.nodes[(graphe.nodes[i]["num_ch"], graphe.nodes[i]["position"][0]-j)])
#         break            
#     with open("/media/coline/Maxtor/dico_graphe_sim_rassembles_3.pickle", 'rb') as fichier_entree :
#         mon_depickler = pickle.Unpickler(fichier_entree)
#         dico_graphe_sim = mon_depickler.load()
#         
#         types_arn = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16s_arnm", "arnt_16s"]
#         liste_tout = []
#         for elt in types_arn :
#             with open("resolutions.pickle", 'rb') as fichier_pickle :
#                 mon_depickler = pickle.Unpickler(fichier_pickle)
#                 resolutions = mon_depickler.load()
#                  
#                 with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
#                         mon_depickler = pickle.Unpickler(fichier_sortie)
#                         liste_a_garder_noms = mon_depickler.load()    
#                            
#                         for element in liste_a_garder_noms :
#                             if resolutions[element[0]] <= 3 :
#                                 liste_tout.append(element)
#        
#         print(liste_tout)  
#         compteur = 0
#         for i in range(len(liste_tout)) :
#             for j in range(i+1, len(liste_tout)) :
#                 if (liste_tout[i], liste_tout[j]) not in dico_graphe_sim.keys() and (liste_tout[j], liste_tout[i]) not in dico_graphe_sim.keys() :
#                     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s.pickle"%(liste_tout[i][0], liste_tout[i][1]), 'rb') as fichier_graphe1 :
#                         mon_depickler_g1 = pickle.Unpickler(fichier_graphe1)
#                         graphe1 = mon_depickler_g1.load()
#                     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s.pickle"%(liste_tout[j][0], liste_tout[j][1]), 'rb') as fichier_graphe2 :
#                         mon_depickler_g2 = pickle.Unpickler(fichier_graphe2)
#                         graphe2 = mon_depickler_g2.load()
#                              
#                         graphe_commun_max,sim = comparaison(graphe1, graphe2, "petit rat") 
#                     
#                         dico_graphe_sim.update({(liste_tout[i], liste_tout[j]) : {"graphe" : graphe_commun_max, "sim" : sim}})
#                         
#                 compteur += 1
#                 print(compteur)
#         for fic in os.listdir("/media/coline/Maxtor") :
#             if "dico_sim_manque" in fic and '[' in fic :
#                 print(fic)
#                 with open("/media/coline/Maxtor/"+fic, "rb") as fichier_tout :
#                     mon_depickler = pickle.Unpickler(fichier_tout)
#                     dico_tout = mon_depickler.load()
#                     print(len(dico_tout))
#                     compteur = 1
#                     for cle in dico_tout.keys() :
#                         print(cle)
#                         print((cle[0].split("_")[0], int(cle[0].split("_")[1])))
#                         print((cle[1].split("_")[0], int(cle[1].split("_")[1])))
#                         
#                         with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s.pickle"%(cle[0].split("_")[0], int(cle[0].split("_")[1])), 'rb') as fichier_graphe1 :
#                             mon_depickler_g1 = pickle.Unpickler(fichier_graphe1)
#                             graphe1 = mon_depickler_g1.load()
#                         with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s.pickle"%(cle[1].split("_")[0], int(cle[1].split("_")[1])), 'rb') as fichier_graphe2 :
#                             mon_depickler_g2 = pickle.Unpickler(fichier_graphe2)
#                             graphe2 = mon_depickler_g2.load()
#                              
#                             graphe_commun_max,sim = comparaison(graphe1, graphe2, "petit rat") 
#                             
#                             if sim != dico_tout[cle]["sim"] : 
#                                 print("bizarre")
#                                 exit()
#                         
#                         if (cle[0].split("_")[0], int(cle[0].split("_")[1])) in liste_tout and (cle[1].split("_")[0], int(cle[1].split("_")[1])) in liste_tout :
#                             if cle in dico_graphe_sim.keys() and ((isinstance(dico_tout[cle], dict) and isinstance(dico_graphe_sim[cle], dict) and dico_graphe_sim[cle]["sim"] != dico_tout[cle]["sim"]) or \
#                             (not isinstance(dico_tout[cle], dict) and isinstance(dico_graphe_sim[cle], dict) and dico_graphe_sim[cle]["sim"] != dico_tout[cle]) or \
#                             (isinstance(dico_tout[cle], dict) and not isinstance(dico_graphe_sim[cle], dict) and dico_graphe_sim[cle] != dico_tout[cle]["sim"])or \
#                             (not isinstance(dico_tout[cle], dict) and not isinstance(dico_graphe_sim[cle], dict) and dico_graphe_sim[cle] != dico_tout[cle])):
#                                 print(fic)
#                                 print(cle)
#                                 print(dico_graphe_sim[cle])
#                                 print(dico_tout[cle])
#                                 print("bizarre")
#                                 exit()
#                             dico_graphe_sim.update({cle : dico_tout[cle]})
#                             print("petit rat")
#                         print(len(dico_tout))
#                         print(compteur)
#                         compteur += 1
#         print(liste_tout)
#         print(len(dico_graphe_sim))
        
#         with open("/media/coline/Maxtor/dico_graphe_sim_rassembles_4.pickle", 'wb') as fichier_sortie :
#             mon_pickler = pickle.Pickler(fichier_sortie)
#             mon_pickler.dump(dico_graphe_sim)
                
                
    
#     print("petit rat")
#     with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_avec_manque_['23S', '18S', '16S', 'Ribozyme', 'Riboswitch', 'SRP', '28S', '25S', 'Intron', 'arnt_16s_arnm', 'arnt_16s']_new_noms_res_3A.pickle", "rb") as fichier_tout :
#         mon_depickler = pickle.Unpickler(fichier_tout)
#         graphe_tout = mon_depickler.load()
#         print(graphe_tout.number_of_nodes())
#         print("petit rat")
# #         a_enlever = []
# #         for u,v,data in graphe_tout.edges(data=True) :
# #             if isinstance(data["sim"], dict) :
# #                 if data["sim"]["sim"] < 0.6 :
# #                     a_enlever.append((u,v))
# #             else :
# #                 if data["sim"] < 0.6 :
# #                     a_enlever.append((u,v))
# #                 
# #         for elt in a_enlever :
# #             graphe_tout.remove_edge(elt[0], elt[1])
#         
#         with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_avec_manque_Ribozyme_new_noms_res_3a.pickle", "rb") as fichier_ribozyme :
#             mon_depickler = pickle.Unpickler(fichier_ribozyme)
#             graphe_ribozyme = mon_depickler.load()
#             
#             
#             new_dico = {}
#             for u,v,data in graphe_ribozyme.edges(data=True) :
#                 print(u)
#                 print(v)
#                 print(data["sim"])
#                 print(graphe_tout.edges[u,v]["sim"])
#                 if data["sim"] != graphe_tout.edges[u,v]["sim"] :
#                     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s.pickle"%(u[0], u[1]), 'rb') as fichier_graphe1 :
#                         mon_depickler_g1 = pickle.Unpickler(fichier_graphe1)
#                         graphe1 = mon_depickler_g1.load()
#                     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s.pickle"%(v[0], v[1]), 'rb') as fichier_graphe2 :
#                         mon_depickler_g2 = pickle.Unpickler(fichier_graphe2)
#                         graphe2 = mon_depickler_g2.load()
#                         
#                         graphe_commun_max,sim = comparaison(graphe1, graphe2, "petit rat") 
#                         new_dico.update({(u,v) : sim})
#                         print(sim)
#             
#             for u,v,data in graphe_ribozyme.edges(data=True) :
#                 print(u)
#                 print(v)
#                 print(data["sim"])
#                 print(graphe_tout.edges[u,v]["sim"])
#                 if (u,v) in new_dico.keys() :
#                     print(new_dico[(u,v)])
#             
# #             a_enlever = []
# #             for u,v,data in graphe_ribozyme.edges(data=True) :
# #                 if isinstance(data["sim"], dict) :
# #                     if data["sim"]["sim"] < 0.6 :
# #                         a_enlever.append((u,v))
# #                 else :
# #                     if data["sim"] < 0.6 :
# #                         a_enlever.append((u,v))
# #                     
# #             for elt in a_enlever :
# #                 graphe_ribozyme.remove_edge(elt[0], elt[1])
#         
#             for u,v,data in graphe_ribozyme.edges(data=True) :
#                 if (u,v) not in graphe_tout.edges() :
#                     print("ramou")
#     G = nx.Graph()
#     for i in range(1,9) :
#         G.add_node(i)
#     
#     G.add_edge(1,2, sim=0.80)
#     G.add_edge(1,3, sim=0.87)
#     G.add_edge(2,3, sim=0.93)
#     G.add_edge(1,4, sim=0.70)
#     G.add_edge(2,4, sim=0.67)
#     G.add_edge(3,4, sim=0.63)
#     G.add_edge(1,5, sim=0.60)
#     G.add_edge(2,5, sim=0.73)
#     G.add_edge(3,5, sim=0.74)
#     G.add_edge(4,6, sim=0.68)
#     G.add_edge(4,7, sim=0.64)
#     G.add_edge(4,8, sim=0.62)
#     G.add_edge(5,6, sim=0.72)
#     G.add_edge(5,7, sim=0.70)
#     G.add_edge(5,8, sim=0.65)
#     G.add_edge(9,1, sim=0.63)
#     G.add_edge(6,7, sim=0.88)
#     G.add_edge(6,8, sim=0.91)
#     G.add_edge(7,8, sim=0.93)
#     
#     
#     pos = nx.spring_layout(G)
#     edge_labels=dict([((u,v,),d["sim"])for u,v,d in G.edges(data=True)])
# #     for i in range(9) :
# #         for j in range(i+1, 9) :
# #             G.add_edge(i,j)
#     
#     nx.draw_networkx_nodes(G, pos, node_color='black')
#     nx.draw_networkx_edge_labels(G, pos, edge_labels= edge_labels, label_pos=0.7, font_size=6)
#     nx.draw_networkx_edges(G, pos)
#     plt.axis('off')
#     plt.savefig("/home/coline/Documents/Extensions/poster_macim/exemple_clustering.svg", format='svg', transparent=True)
#     plt.show()
    
#     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_5dm6_8.pickle", "rb") as fichier_extension :
#         mon_depickler_ext = pickle.Unpickler(fichier_extension)
#         extension3 = mon_depickler_ext.load()
#         
#         with open("Graphs/5dm6.pickle", 'rb') as fichier_graphe :
#             mon_depickler = pickle.Unpickler(fichier_graphe)
#             graphe = mon_depickler.load()
#         
#             for noeud, data in extension3.nodes(data=True) :
#                 print(noeud)
#                 print(data)
#                 if data["type"] != -1 :
#                     for i in range(data["position"][0], data["position"][1]+1) :
#                         print(graphe.nodes[(data["num_ch"], i)])
#                     
#             
#             for u,v,data in extension3.edges(data=True) :
#                 print(u)
#                 print(v)
#                 print(data)
#             print(graphe.nodes[('X', 473)])   
#             print(graphe[('X', 694)])
#             print(graphe[('X', 695)])
#             print(graphe.nodes[('X', 810)])
#             print(graphe.nodes[('X', 811)])
#             print(graphe[('X', 700)])
#             print(graphe.nodes[('X', 699)])
#             print(graphe.nodes[('X', 702)])
#             print(graphe.nodes[('X', 792)])
            #print(graphe.nodes[('X', 471)])
    #liste = ["SRP", "Intron", "Ribozyme", "28S", "18S", "25S", "16S", "23S"]
    #liste = ["arnt_16s_arnm"]
    #for elt in liste :
    #    script_a_faire(elt)
#     with open("groupes_16S_homologues.pickle", 'rb') as fichier_homologues :
#         mon_depickler_1 = pickle.Unpickler(fichier_homologues)
#         groupes_homologues = mon_depickler_1.load()
#         
#         for groupe in groupes_homologues : 
#             print(len(groupe))
    
    
#     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_5gaf_4.pickle", 'rb') as fichier :
#             mon_depickler = pickle.Unpickler(fichier)
#             graphe = mon_depickler.load()
#              
#             print(graphe.nodes.data())
# #             
#             with open("Graphs/5j91.pickle", 'rb') as fichier_graphe :
#                 mon_depickler = pickle.Unpickler(fichier_graphe)
#                 graphe = mon_depickler.load()
#                 for noeud,data in graphe.nodes(data=True) :
#                     if noeud[1] < 600 :
#                         print(str(noeud)+ " "+str(data))
    #test_combien_reste()
    
    # Define input file
#     in_file = "seq1.fa"
#      
#     # Define output file
#     out_file = "aligned.aln"
#      
#     # Get the command for Clustal Omega
#     clustalw_cline = ClustalwCommandline(infile=in_file, outfile=out_file)
#      
#     # Print the executable command
#     print(clustalw_cline)
#     stdout, stderr = clustalw_cline()
#     print(stdout)
#     print(type(stdout))
#     nb_seq = 3
#     lignes = stdout.split("\n")
#     print(lignes)
#     ligne = lignes[0]
#     compteur = 0
#     while compteur < len(lignes) and "Sequences (1:%d) Aligned."%nb_seq not in ligne :
#         ligne = lignes[compteur]
#         compteur += 1
#     if compteur != len(lignes) :   
#         print(ligne) 
#         print(ligne.split(" "))
#         print(ligne.split(" ")[5])
        
#     print(ligne)
#     print(compteur)
            
    #print(stderr)
    
    #print(align.get_alignment_length())
#     with open("groupes_23S_homologues.pickle", "rb") as fichier_homologues :
#         mon_depickler_1 = pickle.Unpickler(fichier_homologues)
#         groupes_homologues = mon_depickler_1.load()
#             
#         with open("groupes_16S_homologues.pickle", "rb") as fichier_homologues_2 : 
#             mon_depickler_2 = pickle.Unpickler(fichier_homologues_2)
#             groupes_homologues_2 = mon_depickler_2.load()
#                 
#             liste_manque = []
#             for i in range(len(groupes_homologues)) : 
#                 for j in range(len(groupes_homologues_2)) :
#                     print(i)
#                     print(j)
#                     if "dico_sim_new_algo_23S_16S_homologues_groupe_%s_groupe_%s_part_1.pickle"%(i,j) not in os.listdir("/media/coline/Maxtor/Resultats/23S-16S") :#or "dico_sim_new_algo_23S_homologues_groupe_%s_groupe_%s_part_1.pickle"%(i,j) in os.listdir("/media/coline/Maxtor/Resultats_dans_dico_sim") :
#                         liste_manque.append(("23S", "16S",i,j))
#                             
#             print(liste_manque)
#             print(len(liste_manque))
            
    
    
#     liste_manque =   [(0, 92), (1, 6), (1, 7), (1, 8), (1, 12), (1, 25), (1, 52), (3, 7), (3, 10), (3, 17), (3, 31), (3, 35), (3, 88), (3, 92), (4, 5), (4, 17), (4, 35), (4, 36), (4, 88), (5, 8), (5, 38), (5, 41), (5, 47), (5, 52), (5, 74), (5, 88), (6, 7), (6, 12), (6, 36), (6, 47), (6, 52), (6, 74), (6, 88), (6, 92), (7, 8), (7, 28), (7, 47), (7, 52), (7, 88), (7, 92), (8, 9), (8, 31), (8, 36), (8, 52), (8, 88), (8, 92), (9, 88), (10, 92), (12, 35), (17, 28), (34, 93), (37, 88), (40, 92), (41, 92), (47, 92), (52, 90), (74, 92), (88, 92), (92, 93)]
#     liste_a_faire = [(14,88), (14, 48), (30, 80), (15, 103), (13, 65), (11, 16), (34, 78), (31, 74), (28, 74), (2, 92), (16, 48), (0, 1), (23, 86), (18, 104), (11, 101), (0, 84), (11, 97), (24, 85), (30, 63), (14, 38), (1, 5), (34, 74), (18, 52), (5, 6), (34, 81), (0, 50), (18, 66), (18, 48), (30, 52), (31, 55), (30, 69), (35, 68), (16, 94), (11, 35), (35, 74), (0, 9), (3, 6), (25, 52), (26, 49), (34, 44), (25, 54), (25, 36), (11, 18), (13, 45), (16, 40), (25, 39), (14, 46), (35, 84), (0, 2), (23, 104), (35, 87), (30, 54), (23, 84), (11, 77), (16, 30), (0, 70), (19, 28), (0, 88), (19, 36), (19, 35), (31, 73), (14, 50), (30, 86), (31, 72), (0, 69), (11, 45), (0, 29), (13, 30), (11, 68), (5, 7), (26, 45), (24, 98), (19, 26), (16, 65), (25, 33), (31, 57), (30, 65), (19, 41), (15, 74), (0, 54), (18, 84), (13, 24), (31, 51), (13, 73), (23, 92), (6, 89), (26, 43), (11, 57), (26, 27), (11, 46), (15, 81), (18, 105), (14, 42), (0, 16), (0, 28), (11, 51), (14, 15), (11, 61), (11, 76), (30, 62), (32, 86), (16, 98), (14, 27), (13, 39), (24, 66), (0, 30), (11, 15), (30, 84), (31, 58), (34, 48), (31, 56), (18, 58), (28, 86), (32, 92), (0, 21), (19, 32), (13, 32), (10, 103), (30, 88), (4, 92), (0, 33), (16, 67), (8, 90), (30, 82), (32, 90), (0, 32), (6, 93), (18, 51), (24, 49), (13, 27), (13, 21), (16, 89), (24, 39), (19, 31), (25, 63), (19, 66), (30, 53), (13, 72), (13, 85), (16, 71), (16, 100), (31, 50), (29, 105), (11, 24), (13, 78), (1, 4), (10, 102), (31, 88), (11, 17), (35, 88), (11, 44), (34, 75), (13, 83), (23, 105), (13, 23), (31, 60), (19, 47), (18, 70), (24, 48), (0, 53), (0, 40), (25, 65), (16, 83), (18, 44), (15, 84), (16, 34), (30, 85), (30, 57), (25, 72), (34, 52), (19, 65), (0, 58), (15, 66), (16, 63), (32, 89), (18, 50), (13, 31), (18, 29), (2, 96), (30, 50), (11, 32), (31, 53), (25, 68), (13, 75), (4, 48), (34, 46), (25, 43), (19, 72), (25, 59), (15, 96), (24, 61), (2, 97), (30, 67), (0, 15), (31, 83), (30, 37), (27, 51), (30, 48), (17, 26), (11, 83), (24, 82), (31, 76), (13, 97), (24, 105), (16, 38), (25, 48), (11, 95), (11, 55), (19, 34), (34, 56), (11, 40), (32, 88), (0, 51), (19, 60), (13, 86), (13, 82), (13, 70), (18, 102), (11, 64), (25, 40), (15, 104), (24, 89), (11, 60), (18, 75), (19, 69), (8, 91), (15, 106), (30, 77), (16, 61), (18, 99), (18, 96), (13, 105), (35, 89), (14, 31), (19, 77), (25, 37), (14, 35), (25, 50), (27, 52), (32, 83), (6, 91), (15, 71), (15, 76), (25, 27), (18, 90), (17, 25), (24, 72), (14, 40), (2, 103), (0, 65), (24, 93), (31, 66), (18, 35), (29, 99), (24, 30), (11, 92), (14, 21), (11, 90), (16, 101), (15, 59), (13, 96), (34, 65), (26, 28), (2, 98), (23, 102), (35, 66), (14, 28), (16, 103), (35, 81), (13, 61), (14, 29), (7, 89), (19, 45), (19, 68), (0, 38), (35, 83), (14, 20), (18, 20), (11, 28), (16, 41), (30, 32), (17, 19), (19, 49), (11, 98), (11, 48), (10, 88), (34, 79), (11, 65), (13, 57), (14, 47), (24, 80), (0, 56), (24, 65), (18, 72), (18, 55), (16, 28), (23, 89), (25, 66), (2, 91), (13, 55), (26, 34), (4, 101), (16, 73), (0, 13), (18, 93), (13, 84), (35, 86), (17, 21), (15, 69), (19, 39), (16, 20), (30, 81), (14, 32), (24, 33), (34, 57), (30, 56), (18, 41), (24, 103), (15, 94), (4, 49), (24, 84), (16, 86), (11, 38), (13, 14), (26, 30), (30, 33), (31, 52), (13, 53), (25, 30), (34, 66), (18, 83), (16, 51), (25, 44), (13, 81), (11, 36), (28, 88), (23, 99), (26, 40), (0, 61), (19, 30), (0, 44), (4, 45), (34, 41), (13, 92), (10, 99), (24, 55), (0, 85), (11, 49), (18, 31), (19, 27), (13, 56), (15, 75), (17, 22), (11, 94), (11, 52), (13, 79), (2, 101), (18, 69), (34, 88), (0, 5), (0, 75), (15, 102), (18, 91), (25, 74), (31, 71), (34, 82), (11, 80), (19, 81), (2, 100), (18, 98), (10, 106), (23, 96), (35, 85), (24, 28), (25, 35), (19, 43), (15, 72), (19, 54), (15, 97), (14, 45), (0, 36), (14, 24), (14, 52), (31, 89), (13, 66), (19, 90), (26, 41), (26, 42), (0, 8), (30, 73), (23, 91), (13, 104), (13, 91), (19, 33), (11, 105), (19, 57), (0, 46), (31, 75), (15, 64), (24, 96), (18, 79), (16, 74), (0, 91), (31, 92), (13, 22), (16, 50), (34, 55), (7, 91), (28, 80), (18, 56), (34, 53), (13, 89), (15, 50), (24, 90), (35, 78), (25, 71), (18, 53), (30, 74), (1, 2), (13, 33), (0, 79), (25, 32), (11, 72), (30, 41), (19, 75), (15, 53), (23, 85), (24, 104), (16, 59), (0, 12), (13, 50), (16, 33), (24, 35), (30, 43), (11, 74), (13, 90), (14, 41), (13, 54), (15, 62), (11, 53), (0, 10), (16, 17), (34, 69), (15, 54), (19, 84), (18, 36), (15, 95), (23, 93), (34, 62), (13, 25), (24, 101), (3, 5), (25, 29), (18, 95), (30, 47), (30, 58), (13, 42), (35, 76), (18, 82), (34, 64), (2, 102), (25, 62), (25, 73), (10, 104), (14, 36), (18, 42), (14, 22), (24, 77), (19, 40), (23, 97), (24, 29), (30, 87), (19, 46), (16, 49), (19, 73), (13, 26), (24, 32), (0, 57), (0, 64), (19, 56), (15, 86), (15, 99), (24, 41), (19, 59), (18, 39), (11, 56), (11, 100), (18, 47), (24, 74), (11, 54), (16, 66), (18, 60), (24, 92), (18, 86), (30, 66), (13, 93), (16, 78), (13, 59), (16, 31), (25, 47), (32, 91), (24, 38), (26, 32), (16, 26), (25, 31), (15, 98), (25, 46), (35, 69), (13, 102), (13, 94), (18, 92), (14, 25), (13, 103), (24, 88), (15, 79), (13, 71), (11, 22), (18, 64), (11, 102), (0, 68), (18, 106), (0, 19), (26, 36), (18, 94), (18, 76), (28, 82), (13, 43), (16, 21), (16, 57), (0, 24), (11, 62), (11, 58), (18, 40), (34, 37), (16, 19), (24, 26), (25, 38), (16, 52), (11, 41), (25, 28), (13, 18), (25, 49), (0, 41), (28, 84), (24, 63), (34, 77), (19, 48), (15, 65), (24, 46), (13, 63), (13, 49), (24, 40), (16, 104), (16, 76), (11, 30), (13, 106), (32, 85), (12, 13), (11, 78), (13, 52), (26, 44), (30, 68), (31, 84), (26, 46), (4, 104), (0, 22), (0, 34), (11, 14), (35, 79), (28, 87), (24, 75), (16, 25), (25, 42), (15, 101), (13, 88), (4, 91), (0, 31), (35, 71), (18, 77), (16, 62), (11, 89), (13, 29), (29, 106), (0, 35), (28, 92), (34, 47), (16, 64), (0, 27), (0, 4), (30, 90), (26, 33), (0, 49), (13, 58), (25, 67), (19, 22), (26, 51), (16, 69), (16, 68), (30, 49), (13, 35), (18, 19), (30, 42), (15, 105), (30, 38), (24, 86), (13, 80), (15, 68), (31, 77), (23, 88), (18, 101), (26, 38), (30, 64), (19, 25), (11, 12), (11, 84), (25, 45), (32, 82), (11, 19), (28, 78), (4, 102), (24, 50), (13, 69), (35, 72), (24, 56), (8, 89), (19, 23), (25, 51), (18, 87), (16, 96), (18, 28), (0, 81), (24, 87), (24, 44), (25, 57), (19, 74), (24, 81), (0, 17), (16, 102), (31, 80), (24, 31), (30, 76), (11, 34), (31, 85), (15, 60), (30, 83), (30, 59), (10, 100), (30, 36), (34, 40), (24, 64), (0, 43), (7, 90), (10, 87), (0, 45), (24, 94), (18, 85), (16, 42), (32, 93), (19, 61), (16, 81), (19, 87), (18, 89), (16, 84), (34, 63), (23, 90), (24, 34), (32, 87), (34, 72), (18, 88), (14, 34), (35, 70), (18, 71), (0, 60), (16, 54), (19, 38), (19, 44), (15, 77), (25, 60), (13, 47), (13, 62), (15, 85), (13, 98), (14, 37), (28, 91), (26, 35), (18, 24), (11, 73), (16, 44), (0, 87), (16, 45), (30, 71), (31, 79), (24, 76), (11, 67), (26, 48), (16, 90), (15, 78), (14, 23), (19, 52), (16, 18), (16, 95), (16, 27), (28, 81), (0, 83), (35, 80), (19, 80), (0, 23), (24, 25), (24, 68), (24, 91), (11, 59), (11, 21), (29, 101), (11, 81), (34, 68), (15, 58), (1, 3), (23, 100), (28, 79), (18, 68), (10, 98), (0, 55), (19, 29), (24, 79), (11, 87), (19, 64), (18, 33), (13, 19), (34, 61), (19, 86), (11, 39), (15, 91), (19, 78), (16, 91), (18, 45), (16, 80), (0, 47), (34, 76), (11, 37), (14, 19), (13, 99), (0, 77), (11, 27), (28, 89), (30, 60), (4, 50), (16, 77), (19, 50), (2, 95), (14, 39), (30, 79), (30, 45), (34, 36), (11, 93), (18, 65), (23, 95), (28, 83), (24, 95), (25, 26), (34, 54), (31, 54), (14, 49), (18, 57), (15, 92), (16, 105), (18, 73), (16, 60), (35, 75), (24, 97), (13, 100), (4, 105), (15, 82), (17, 20), (2, 99), (19, 53), (24, 57), (34, 51), (28, 90), (30, 51), (15, 90), (19, 20), (11, 31), (19, 70), (19, 42), (13, 41), (13, 77), (19, 82), (34, 80), (31, 67), (24, 51), (15, 80), (25, 58), (24, 47), (16, 88), (24, 58), (15, 83), (31, 82), (25, 41), (32, 78), (25, 69), (18, 59), (34, 50), (19, 51), (10, 101), (30, 35), (24, 59), (24, 43), (15, 52), (13, 15), (13, 76), (18, 37), (31, 62), (16, 32), (19, 55), (18, 38), (0, 66), (11, 25), (19, 92), (16, 82), (24, 73), (24, 52), (14, 17), (31, 69), (19, 71), (2, 93), (26, 29), (15, 49), (13, 51), (34, 71), (34, 38), (0, 26), (18, 62), (30, 55), (15, 88), (18, 54), (14, 33), (0, 72), (31, 70), (27, 50), (0, 89), (11, 99), (34, 67), (16, 93), (31, 64), (0, 52), (4, 106), (16, 56), (26, 39), (13, 16), (18, 34), (24, 54), (11, 104), (16, 87), (11, 70), (0, 42), (30, 46), (19, 85), (13, 60), (25, 64), (34, 59), (28, 75), (19, 67), (13, 38), (0, 86), (30, 39), (0, 20), (0, 6), (11, 47), (24, 83), (0, 7), (15, 63), (2, 94), (19, 79), (16, 92), (34, 49), (0, 39), (18, 97), (16, 58), (15, 73), (13, 37), (4, 42), (25, 34), (0, 37), (16, 70), (24, 36), (11, 71), (11, 88), (16, 43), (34, 70), (18, 63), (19, 62), (28, 76), (29, 102), (14, 18), (13, 34), (14, 44), (32, 79), (18, 103), (15, 51), (31, 63), (19, 91), (0, 11), (0, 74), (4, 103), (0, 62), (18, 23), (18, 67), (17, 18), (17, 23), (11, 63), (16, 36), (13, 101), (26, 47), (13, 20), (24, 78), (15, 67), (15, 70), (11, 103), (18, 74), (25, 55), (19, 83), (31, 61), (16, 72), (35, 82), (11, 20), (16, 22), (16, 75), (16, 99), (17, 27), (11, 23), (13, 74), (29, 103), (11, 26), (0, 18), (24, 27), (30, 89), (11, 79), (14, 43), (34, 43), (6, 90), (32, 80), (24, 102), (0, 48), (18, 26), (26, 31), (0, 3), (29, 100), (24, 100), (34, 35), (13, 68), (16, 46), (4, 51), (0, 67), (13, 36), (34, 58), (24, 60), (16, 97), (31, 90), (35, 67), (11, 91), (34, 45), (30, 70), (15, 61), (16, 79), (13, 64), (31, 59), (16, 53), (4, 43), (13, 67), (24, 106), (18, 22), (11, 42), (11, 106), (13, 87), (25, 53), (18, 25), (11, 86), (26, 37), (30, 91), (34, 42), (34, 73), (13, 44), (32, 81), (24, 53), (15, 87), (31, 65), (0, 59), (19, 88), (31, 86), (24, 42), (2, 105), (23, 98), (0, 71), (3, 4), (16, 39), (4, 46), (16, 55), (15, 57), (4, 41), (19, 24), (14, 51), (18, 30), (15, 56), (0, 90), (34, 83), (30, 75), (18, 61), (19, 21), (24, 67), (23, 87), (18, 100), (30, 31), (0, 63), (19, 58), (24, 45), (34, 87), (0, 25), (11, 69), (24, 70), (34, 39), (11, 29), (11, 96), (19, 89), (13, 95), (18, 81), (19, 37), (30, 61), (34, 85), (0, 80), (15, 93), (34, 84), (23, 103), (34, 60), (19, 63), (11, 66), (0, 73), (34, 86), (15, 100), (18, 80), (13, 28), (18, 21), (15, 89), (31, 78), (17, 24), (31, 81), (0, 14), (30, 44), (24, 62), (13, 46), (0, 82), (11, 82), (14, 30), (16, 47), (25, 61), (28, 77), (29, 104), (24, 69), (25, 70), (0, 78), (30, 72), (16, 29), (11, 75), (13, 48), (16, 35), (2, 104), (28, 85), (30, 40), (18, 78), (30, 34), (0, 76), (18, 46), (32, 84), (16, 37), (23, 94), (15, 55), (16, 24), (23, 106), (18, 49), (16, 106), (4, 47), (16, 85), (14, 26), (18, 43), (25, 56), (18, 27), (13, 40), (11, 85), (35, 77), (4, 44), (30, 78), (16, 23), (11, 43), (18, 32), (31, 91), (31, 87), (11, 33), (26, 50), (24, 71), (14, 16), (13, 17), (24, 37), (31, 68), (11, 13), (10, 97), (10, 105), (21, 25), (35, 73), (23, 101), (19, 76), (2, 106), (24, 99)]
#     
#     for elt in liste_manque :
#         if elt not in liste_a_faire :
#             print(elt)
            
            
#     for fic in os.listdir("/media/coline/Maxtor/Resultats_new") :
#         if fic not in os.listdir("/media/coline/Maxtor/Resultats") :
#             print(fic)
    
#     with open("groupes_23S_homologues.pickle", "rb") as fichier_homologues :
#         mon_depickler_1 = pickle.Unpickler(fichier_homologues)
#         groupes_homologues = mon_depickler_1.load()
#     
#         
#         for elt in liste_a_faire :
#             taille = math.ceil((len(groupes_homologues[elt[0]])*len(groupes_homologues[elt[1]]))/8000)
# #             if taille > 10 :
# #                 print(elt)
# #                 print(taille)
#             vu = [False]*taille
#             for i in range(1, taille+1) :
#                 for fic in os.listdir("/media/coline/Maxtor/Resultats") :
#                     if "groupe_%d_groupe_%d_part_%d"%(elt[0],elt[1],i) in fic :
#                         vu[i-1] = True
#             if False in vu :
#                 print(elt) 
    
#     with open("graphs_2.92.pickle", 'rb') as fichier_graphes :
#         mon_depickler = pickle.Unpickler(fichier_graphes)
#         graphes = mon_depickler.load()
#         print(graphes[('4V9F', '0')][1322])

#     with open("Graphs/1fjg.pickle", 'rb') as fichier_graphe :
#         mon_depickler = pickle.Unpickler(fichier_graphe)
#         graphe = mon_depickler.load()
#         print(graphe.nodes())
    
#     compteur = 1
#     while "fichier_5j7l_%d.pickle"%compteur in os.listdir(NEW_EXTENSION_PATH_TAILLE) :
#         with open(NEW_EXTENSION_PATH_TAILLE+"fichier_5j7l_%d.pickle"%compteur, 'rb') as fichier :
#             mon_depickler = pickle.Unpickler(fichier)
#             graphe = mon_depickler.load()
            #print(graphe.nodes.data())
#             if graphe.nodes[1]["position"][0] == 1226 and graphe.nodes[2]["position"][0] == 812 :
#                 print("fichier_5j7l_%d.pickle"%compteur)
#                 print(graphe.nodes.data())
#         compteur += 1

#     compteur = 0
#     liste_a_faire = []
#     with open("groupes_23S_homologues.pickle", "rb") as fichier_homologues :
#         mon_depickler_1 = pickle.Unpickler(fichier_homologues)
#         groupes_homologues = mon_depickler_1.load()
#        
#         for i in range(0, len(groupes_homologues)) :
#             for j in range(i+1, len(groupes_homologues)) :
#                 if "dico_sim_new_algo_23S_homologues_groupe_%s_groupe_%s_part_1.pickle"%(i,j) not in os.listdir("/home/coline/Bureau/Resultats") and "dico_sim_new_algo_23S_homologues_groupe_%s_groupe_%s_part_1.pickle"%(i,j) not in os.listdir("/home/coline/Bureau/Resultats/Resultats") :
#                     print((i,j))
#                     liste_a_faire.append((i,j))
#         print(liste_a_faire)
#         print(len(liste_a_faire))
    
    
    
    
    

#         
#         with open(NEW_EXTENSION_PATH_TAILLE+"Traitement_resultats/dico_sim.pickle", 'wb') as fichier_pickle :
#             mon_pickler = pickle.Pickler(fichier_pickle)
#             mon_pickler.dump(dico_sim)
            

#     liste_manque_morceaux = []
#     for fic in os.listdir("/home/coline/Bureau/Resultats") :
#         print(fic)
#         if os.path.isfile("/home/coline/Bureau/Resultats/"+fic) :
#             print("Dernière modification: %s" % time.ctime(os.path.getmtime("/home/coline/Bureau/Resultats/"+fic)))
#             date = int(time.ctime(os.path.getmtime("/home/coline/Bureau/Resultats/"+fic)).split(" ")[2])
#             temps = time.ctime(os.path.getmtime("/home/coline/Bureau/Resultats/"+fic)).split(" ")[3]
#             temps = int(temps.split(":")[0])*60+int(temps.split(":")[1])
#             print(date)
#             print(temps)
#             
#             if ((date < 18 ) or (date == 18 and temps < 832)) and len(fic.split("_")) > 8:
#                 un = int(fic.split("_")[7])
#                 deux = int(fic.split("_")[9].split(".")[0])
#                 if (un, deux) not in liste_manque_morceaux :
#                     liste_manque_morceaux.append((un,deux))
#     
#     print(liste_manque_morceaux)    
#     print(len(liste_manque_morceaux))

#     liste_que_sim = []
#     liste_pas_bon = []
#     for fic in os.listdir("/home/coline/Bureau/Resultats/Resultats") :
#         print(fic)
#         if os.path.isfile("/home/coline/Bureau/Resultats/Resultats/"+fic) :
#             with open("/home/coline/Bureau/Resultats/Resultats/"+fic, 'rb') as fichier :
#                 try : 
#                     mon_depickler = pickle.Unpickler(fichier)
#                     dico = mon_depickler.load()
#                     #print(type(dico[list(dico.keys())[0]]))
#                     if isinstance(dico[list(dico.keys())[0]], float):
#                         liste_que_sim.append(fic)
#                 except pickle.UnpicklingError :
#                     liste_pas_bon.append(fic)
#                 except EOFError :
#                     liste_pas_bon.append(fic)
#                 #print(dico[list(dico.keys())[0]])
#     print(liste_que_sim)
#     print(len(liste_que_sim))
#      
#     print(liste_pas_bon)
#     print(len(liste_pas_bon))
    
    
#     with open(NEW_EXTENSION_PATH_TAILLE+"Resultats/dico_sim_new_algo_18S_homologues_groupe_1_groupe_6_part_1.pickle", 'rb') as fichier:
#             mon_depickler = pickle.Unpickler(fichier)
#             dico = mon_depickler.load()
#             
#             print(len(dico))

   

#                 with open("groupes_23S_homologues.pickle", 'rb') as fichier_homologues :
#                     mon_depickler_homologues = pickle.Unpickler(fichier_homologues)
#                     groupes_homologues = mon_depickler_homologues.load()
#                     
#                     compteur1 = 0
#                     compteur2 = 0
#                     compter = 0
#                     for groupe1 in groupes_homologues :
#                         if compteur1 == 4 :
#                             for groupe2 in groupes_homologues :
#                                 if compteur2 == 92 :
#                                     for elt1 in groupe1 :
#                                         for elt2 in groupe2 :
#                                             if (str(elt1[0])+"_"+str(elt1[1]), str(elt2[0])+"_"+str(elt2[1])) not in dico_sim.keys() :
#                                                 print(elt1)
#                                                 print(elt2)
#                                                 compter += 1
#                                 compteur2 += 1 
#                         compteur1 += 1
#                         
#                     print(compter)
#                     print(len(dico_sim.keys()))
#                     break
    
    #tabnanny.check("/home/coline/Bureau/new_algo_comparaison.py")
#     

#     with open("test_dico_sim_new_algo.pickle", 'rb') as fichier :
#         mon_depickler = pickle.Unpickler(fichier)
#         new_dico_sim = mon_depickler.load()
#         
#         with open(EXTENSION_PATH%4+"sim_extensions_toutes_aretes_coeff_all1_taille_4.pickle", "rb") as fichier_old :
#             mon_depickler = pickle.Unpickler(fichier_old)
#             dico_sim = mon_depickler.load()
#             
#             compteur = 0
#             print(dico_sim.keys())
#             print(len(new_dico_sim))
#             for cle in new_dico_sim.keys() :
#                 if cle in dico_sim.keys() :
#                     old_cle = cle
#                 elif (cle[1], cle[0]) in dico_sim.keys() :
#                     old_cle = (cle[1], cle[0])
#                 else :
#                     old_cle = "ramous"
#                 #print(cle)
# #                 if cle == ('2XD0_V_36_21', '5DM6_X_48_9') :
# #                     print(cle)
# #                     print(new_dico_sim[cle])
# #                     print(dico_sim[old_cle])
#                 if new_dico_sim[cle] != dico_sim[old_cle] :
#                     
#                     if new_dico_sim[cle] > dico_sim[old_cle] :
#                         print(cle)
#                         print(new_dico_sim[cle])
#                         print(dico_sim[old_cle])
#                     compteur += 1
#             print(compteur)
            
#     with open("Graphs/6d9j.pickle", 'rb') as fichier_pickle :
#         mon_depickler = pickle.Unpickler(fichier_pickle)
#         tab_aminor = mon_depickler.load()
#          
#         print(tab_aminor)
#          
#         for noeud, data in tab_aminor.nodes(data=True) :
#             if noeud[1] > 2000 :
#                 print(str(noeud)+ " " +str(data))
             
#         for u,v,data in tab_aminor.edges(data=True) :
#             if u[1] > 1000 and v[1] > 2000 :
#                 print(str((u,v)) + " " + str(data))
        
#     with open("all_aminor.pickle", 'rb') as fichier_all_aminor :
#         mon_depickler = pickle.Unpickler(fichier_all_aminor)
#         all_aminor = mon_depickler.load()
#         
#         for cle in all_aminor.keys() :
#             if cle == '6cae' :
#                 compteur = 1
#                 for graphe in all_aminor[cle] :
#                     if compteur == 26 :
#                         print(graphe.nodes.data())
#                     compteur += 1
    
    #mmcif_dict = MMCIF2Dict.MMCIF2Dict(PATH_MMCIF+'5FDU.cif')
    
    #for elt in mmcif_dict['loop_'] :
    #    print(elt)
    
#     num_chaine = '1A'
#     
#     doc = cif.read_file(PATH_MMCIF+"5FDU.cif")
#     block = doc.sole_block()
#     print(block.name)
#     cat = block.find_mmcif_category("_struct_ref_seq.") 
#     print(list(cat.tags))
#     for row in cat :
#         if row[3] == '1A' :
#             id = row[1]
#     
#     cat = block.find_mmcif_category("_entity.") 
#     
#     for row in cat :
#         if row[0] == id :
#             type_ARN = row[3]
#             
#     cat = block.find_mmcif_category("_entity_src_nat.")         
#            
#     for row in cat :
#         if row[0] == id :
#             org = row[6]
#             
#     print(type_ARN)
#     print(org)
    
    #print(mmcif_dict['_atom_site.id'])
    #distrib_sim_range_par_taille_ext(['4V88_A5_25_30', '4YAZ_R_36_25', '1FJG_A_58_23', '5FDU_1A_30_17', '5J7L_DA_25_10'], "groupe_aleatoire_gnra")
    
#     for elt in CLUSTERING_PEREZ_VERSION_NON_CAN_2[4] :
#         with open(EXTENSION_PATH_TAILLE%6+"fichier_"+elt+".pickle", 'rb') as fichier :
#             mon_depickler = pickle.Unpickler(fichier)
#             graphes = mon_depickler.load()
#             
#             print(graphes.edges.data())
    
    
    
    #problemes sommets artificiels lies a plusieurs vrais sommets (13/07/19)
#     for i in range(1,11) :
#         list_pas_bon = []
#         for fic in os.listdir(EXTENSION_PATH_TAILLE%i) :
#             if ".pickle" in fic and "couples_possibles" not in fic and "coord" not in fic and len(fic.split("_")) == 5 :
#                 with open(EXTENSION_PATH_TAILLE%i+fic, 'rb') as fichier_graphe :
#                     mon_depickler = pickle.Unpickler(fichier_graphe)
#                     graphe = mon_depickler.load()
#                      #print(graphe.nodes.data())
#                     for noeud, data in graphe.nodes(data=True) :
#                         if data["type"] == -1 and len(graphe[noeud]) != 1 :
#                             if fic not in list_pas_bon :
#                                 list_pas_bon.append(fic)
#                              #print(noeud)
#                              
#         print(i)
#         print(len(list_pas_bon))
#         print(list_pas_bon)

#     with open(EXTENSION_PATH_TAILLE%3 + "fichier_1FJG_A_271_1_3.pickle", 'rb') as fichier :
#             mon_depickler = pickle.Unpickler(fichier)
#             graphe = mon_depickler.load()
#             
#             print(graphe.nodes.data())
#             print(graphe.edges.data())



#     with open(EXTENSION_PATH_TAILLE%4+"couples_possibles_fichier_5FDU_1A_134_3_fichier_5DM6_X_134_2.pickle", 'rb') as fichier :
#         mon_depickler = pickle.Unpickler(fichier)
#         couples = mon_depickler.load()
#         print(couples)

#     with open(PATH_MMCIF+"fichiers_rmsd_taille_4_que_carbone_1.pickle", 'rb') as fichier :
#         mon_depickler = pickle.Unpickler(fichier)
#         rmsds = mon_depickler.load()
#         groupe = list(GROUPE_GNRA)
#         groupe.extend(GROUPE_GNRA_ETENDU)
#         
#         
#         for elt in CLUSTERING_PEREZ_VERSION_NON_CAN_2[11] :
#             groupe.remove(elt)
#         
#         print(len(groupe))
#         print(groupe)
#         
#         compteur = 0
#         for elt in rmsds.keys() :
#            # print(elt[0][8:len(elt)-15])
#             if elt[0][8:len(elt)-15] in groupe and elt[1][8:len(elt)-15] in groupe : 
#                 print("ramousnif")
#                 print(elt)
#                 print(rmsds[elt])
#                 compteur +=1
#                 if rmsds[elt] == None :
#                     print("gros rat")
#         print(compteur)
        
#         for cle in graphe.keys() :
#             print(graphe[cle].edges.data())

#     with open("graphs_2.92.pickle", 'rb') as fichier_tout :
#         mon_depickler_graphes = pickle.Unpickler(fichier_tout)
#         graphes = mon_depickler_graphes.load()
#         
#         print(graphes[('3JCS', '1')].nodes[1753])
    
#     test_a_enlever = []
#     for groupe in HOMOLOGUES :
#         compteur = 0
#         for elt in groupe :
#             if compteur > 0 :
#                 test_a_enlever.append(elt)
#             compteur +=1
#             
#     print(test_a_enlever)
                        
#     with open("classer.sh", 'w') as fichier_classer :
#         for fichier in os.listdir("/home/coline/eclipse-workspace/a_minor_extension_V2/recup_data") :
#             if ".py" not in fichier :
#                 fichier_classer.write("mv %s %s\n"%(fichier, "Non_classes/"+fichier))


            
        
            
                        
                        
                        
                                                                  