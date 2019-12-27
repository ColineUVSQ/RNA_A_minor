'''
Created on 6 déc. 2018

@author: coline
Dessiner les superpositions pour les extensions version matplotlib
'''
from recup_data.constantes import EXTENSION_PATH, EXTENSION_PATH_TAILLE,\
    NEW_EXTENSION_PATH_TAILLE

import pickle
import networkx as nx
import matplotlib.pyplot as plt
from networkx.classes.function import get_node_attributes
import os
import numpy as np
from collections import OrderedDict
from recup_data.new_algo_comparaison import calcul_sim_aretes_avec_coeff
#from recup_data.calcul_sim import calcul_sim_aretes_avec_coeff

liste = ['5J7L_DA_191_3', '5J7L_DA_191_4', '5FDU_1A_301_1', '5J7L_DA_301_2', '5DM6_X_334_1', '5FDU_1A_334_2', '4V9F_0_335_1', '5J7L_DA_335_2', '3JCS_1_137_4', '4V88_A5_290_1', '4V88_A6_314_2', '5J7L_DA_218_3', '4V9F_0_251_2', '1FJG_A_62_8', '5J7L_DA_137_1', '4V9F_0_118_1', '4V9F_0_62_2', '5J7L_DA_271_2', '4V9F_0_224_1', '5DM6_X_197_1', '3GX5_A_138_6', '1FJG_A_317_2', '5J5B_BA_317_1', '1FJG_A_326_1', '5DM6_X_137_3', '5J5B_BA_314_1', '4V9F_0_134_6', '4V9F_0_328_1', '4V9F_0_197_2', '4V9F_0_62_16', '5J7L_DA_282_2', '4V88_A5_137_2', '5FDU_1A_224_3', '5J7L_DA_326_2']

''' renvoie vrai si le graphe_commun passe en argument contient le sommet passe en argument
renvoie faux sinon 
le parametre num permet de savoir quel index de noeud il faut regarder (premier index pour le premier graphe ou deuxième index pour le deuxième graphe)'''

def contient_sommet(sommet, graphe_commun, num ):
    for noeud in graphe_commun.nodes() :
        if noeud[num] == sommet :
            return True
    return False

''' renvoie vrai si le graphe_commun passe en argument contient l arete passee en argument
renvoie faux sinon 
le parametre num permet de savoir quel index d'arete il faut regarder (premier index pour le premier graphe ou deuxième index pour le deuxième graphe)'''

def contient_arete(arete, graphe_commun, num):
    for edge in graphe_commun.edges() :
        if edge[0][num] == arete[0] and edge[1][num] == arete[1] or edge[1][num] == arete[0] and edge[0][num] == arete[1]   :
            return True
    return False

def draw_isomorphisme(taille_ext):
    #with open("fichier_affichage_isomorphisme_nouvelle_metrique_coeff_all1_cww_non_can_taille_%s.txt"%taille_ext, 'w') as fichier :
        
        
#         with open(EXTENSION_PATH%taille_ext+"dico_comp_complet_metrique_toutes_aretes_coeff_all1_taille_%s.pickle"%taille_ext, 'rb') as fichier_graphe :
#             mon_depickler = pickle.Unpickler(fichier_graphe)
#             dico_graphe = mon_depickler.load()
#             print(len(dico_graphe))
#             with open(EXTENSION_PATH%taille_ext+"sim_extensions_toutes_aretes_coeff_all1_taille_%s.pickle"%taille_ext, 'rb') as fichier_sim : 
#                 mon_depickler_sim = pickle.Unpickler(fichier_sim)
#                 dico_sim = mon_depickler_sim.load()
#                
#                 dico_sim_sorted = OrderedDict(sorted(dico_sim.items(), key= lambda t: t[1], reverse=True))
#                 print(dico_sim_sorted)
                
        #for fic in os.listdir("/home/coline/Bureau/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_taille_%d/"%taille_ext) : 
            #if "graphe_comp"  in fic : 
                #with open("/home/coline/Bureau/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_taille_%d/"%taille_ext +fic, 'rb') as fichier_graphe :
#                 with open(EXTENSION_PATH_TAILLE%taille_ext +"graphe_comp_test_couples_possibles_fichier_5FDU_1A_134_3_fichier_5DM6_X_134_2.pickle", 'rb') as fichier_graphe :
#                     mon_depickler = pickle.Unpickler(fichier_graphe)
#                     dico_graphe = mon_depickler.load()   
                with open("dico_algo_heuristique.pickle", 'rb') as fichier_graphe :
                    mon_depickler = pickle.Unpickler(fichier_graphe)
                    dico_graphe = mon_depickler.load() 
                    #print(dico_graphe)
                    print(len(dico_graphe))
                    compteur = 0
                    for elt in dico_graphe.keys() :
                        #print(elt)
                        if elt[0] == ('5dm6', 4)  and elt[1] == ('4ybb', 18) :
                        #if compteur < 10 :
                                with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt[0][0]+"_"+str(elt[0][1])+"_2.pickle", 'rb') as fichier_graphe1 :
                                        mon_depickler1 = pickle.Unpickler(fichier_graphe1)
                                        graphe1 = mon_depickler1.load()
                                        
                                with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt[1][0]+"_"+str(elt[1][1])+"_2.pickle", 'rb') as fichier_graphe2 :
                                        mon_depickler2 = pickle.Unpickler(fichier_graphe2)
                                        graphe2 = mon_depickler2.load()
                                print(graphe1.edges.data())
                                print(graphe2.edges.data())
                                print(calcul_sim_aretes_avec_coeff(graphe1, graphe2, dico_graphe[elt]["graphe"], "petit rat", 1, 1, 1))
                                print(dico_graphe[elt]["graphe"].edges.data())
                                #elt = ('fichier_5J7L_DA_30_15', 'fichier_5FDU_1A_30_17')
                                #print(dico_graphe[elt].edges.data())
                                enlever = False 
        #                         for l in liste : 
        #                                 if l in elt[0] or l in elt[1] :
        #                                     enlever = True
        #                         print(enlever)
                                    
                                if enlever == False :
        #                                 cle_sim_1 = elt[0].split('_')[1] + "_" + elt[0].split('_')[2] + "_" + elt[0].split('_')[3] + "_" + elt[0].split('_')[4]
        #                                 cle_sim_2 = elt[1].split('_')[1] + "_" + elt[1].split('_')[2] + "_" + elt[1].split('_')[3] + "_" + elt[1].split('_')[4]
        #     
        #                             if (cle_sim_1, cle_sim_2) in dico_sim.keys() :
        #                                     sim = round(dico_sim[(cle_sim_1, cle_sim_2)],2)
        #                             else :
        #                                     sim = round(dico_sim[(cle_sim_2, cle_sim_1)],2)
        
                                    sim = round(calcul_sim_aretes_avec_coeff(graphe1, graphe2, dico_graphe[elt]["graphe"], "petit rat", 1, 1, 1),2)
                                    
                                    if str(sim)+"__" + str(elt[0]) + "_" + str(elt[1]) + ".png" not in os.listdir(EXTENSION_PATH_TAILLE%taille_ext) :
                                        GC = dico_graphe[elt]["graphe"].copy()
                #                         print(GC.nodes.data())
                #                         
                                        compteur_elt = 0
                #                         
                                        fig, axs=plt.subplots(figsize=(10, 12), nrows=1, ncols=2)
                #                         #fig=plt.figure()
                #             
                                        columns = 2
                                        rows = 1
                #                         
                                        nom = ""
                                        #elt = ('fichier_5FDU_1A_197_3', 'fichier_5DM6_X_48_9')
                                        print(elt)
                                        
                                        
        #                                 with open("/home/coline/Bureau/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_taille_%d/"%taille_ext+elt[0]+".pickle", "rb") as fichier_1 :
        #                                     mon_depickler_1 = pickle.Unpickler(fichier_1)
        #                                     graphe1 = mon_depickler_1.load()      
        #                                     with open("/home/coline/Bureau/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_taille_%d/"%taille_ext+elt[1]+".pickle", "rb") as fichier_2 :
        #                                         mon_depickler_2 = pickle.Unpickler(fichier_2)
        #                                         graphe2 = mon_depickler_2.load() 
        #                 
        #                                         sim = calcul_sim_aretes_avec_coeff(graphe1, graphe2, dico_graphe[elt], elt, 1, 1, 1)
                                        
                                        ''' on veut afficher la similarite et le rang '''                           
        #                                 if (cle_sim_1, cle_sim_2) in dico_sim.keys() :
        #                                     etape = 1
        #                                     rang = 0
        #                                     for elt_sim in dico_sim_sorted :
        #                                         if elt_sim == (cle_sim_1, cle_sim_2) : 
        #                                             rang = etape
        #                                         etape += 1
        #                                      
        #                                 else :
        #                                     etape = 1
        #                                     rang = 0
        #                                     for elt_sim in dico_sim_sorted :
        #                                         if elt_sim == (cle_sim_2, cle_sim_1) : 
        #                                             rang = etape
        #                                         etape += 1
        
                                        
                                        #fig.suptitle('sim = %.2f, rang = %d/4005'%(round(sim,2), rang), fontsize=16)
                                        fig.suptitle('sim = %.2f'%(round(sim,2)), fontsize=16)    
                                        
                                        for element in elt :
                                                print(element)
                                                print(compteur)
                                                '''on cree un subplot par graphe de la comparaison'''
                                                fig.add_subplot(rows, columns, compteur_elt%2+1)
                                                #fig.axis('off')
                                                #axs[compteur%2].xaxis.set_visible(False)
                                                axs[compteur_elt%2].axis("off")
                                                axs[compteur_elt%2].set_title(element)
                                            
#                                             with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+element+".pickle", 'rb') as fichier_entree :
                                                print(element)
                                                if compteur_elt == 0 :  
                                                    G = graphe1
                                                else :
                                                    G = graphe2
#                                                 mon_depickler = pickle.Unpickler(fichier_entree)
#                                                 G = mon_depickler.load()
                                                print(G.nodes.data())
                                                
                                                nx.set_node_attributes(G, (33,33), "coordonnees")
                                                G.nodes[1]["coordonnees"] = (0.0,0.5)
                                                G.nodes[2]["coordonnees"] = (2.0,0.5)
                                                G.nodes[3]["coordonnees"] = (0.0,0.0)
                                                G.nodes[4]["coordonnees"] = (2.0,0.0)
                                                G.nodes[5]["coordonnees"] = (3.0,0.5)
                                                
        #                                         fichier.write(str(element)+"\n") 
        #                                         fichier.write(str(G.number_of_nodes())+"\n") 
                                                #print(G.nodes())
                                                
                                                nodes_list = [u for u,d in G.nodes(data=True) if d["type"] != -1] 
                                                print(nodes_list)
                                                ordre_noeuds = [1,2,3,4,5]
                                                
                                                '''on ordonne les noeuds par chaine'''
                                                chaines = [[1]]
                                                for i in range(1,5) :
                                                    compteur = i
                                                    if i != 1 : chaines.append([i])
                                                    liaison_B53 = True
                                                    while liaison_B53 :
                                                        liaison_B53 = False
                                                        temp = compteur
                                
                                                        for voisin in G.successors(compteur) :
                                                            for arc in G[compteur][voisin] :
                                                                if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and G[compteur][voisin][arc]["label"] == 'B53' :
                                                                    liaison_B53 = True
                                                                    temp = voisin
                                                                    chaines[len(chaines)-1].append(voisin)
                                                                     
                                                        for voisin in G.predecessors(compteur) :
                                                            for arc in G[voisin][compteur] :
                                                                if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and G[voisin][compteur][arc]["label"] == 'B53' :
                                                                    liaison_B53 = True
                                                                    temp = voisin
                                                                    chaines[len(chaines)-1].append(voisin)
                                                        compteur = temp
                                                
                                                for i in range(4) :
                                                    for elt in chaines[i] :
                                                        if elt not in ordre_noeuds :
                                                            ordre_noeuds.append(elt)
                                                
                                                for noeud in ordre_noeuds :
                                                    #voisins = G[noeud]
                                                    coordonnees_noeud = G.nodes[noeud]["coordonnees"]
                                                    for pred in G.predecessors(noeud) :
                                                        if G.nodes[pred]["coordonnees"] == (33,33) :
                                                            coordonnees = []
                                                            for node in G.nodes() :
                                                                coordonnees.append(G.nodes[node]["coordonnees"])
                                                            for edge in G[pred][noeud] :
                                                                if G[pred][noeud][edge]["label"] == "B53" :
                                                                    if (coordonnees_noeud[0], coordonnees_noeud[1]-0.5) not in coordonnees :
                                                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.5)
                                                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.5) not in coordonnees :
                                                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.5)
                                                                    elif (coordonnees_noeud[0]-0.5, coordonnees_noeud[1]) not in coordonnees :
                                                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.5, coordonnees_noeud[1])
                                                                    elif (coordonnees_noeud[0]+0.5, coordonnees_noeud[1]) not in coordonnees :
                                                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.5, coordonnees_noeud[1])
                                                                    elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                                                    elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                                                    else :
        #                                                                 fichier.write("probleme\n")
                                                                        print("probleme")
                                                                else :
                                                                    if (coordonnees_noeud[0]-0.75, coordonnees_noeud[1]) not in coordonnees :
                                                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.75, coordonnees_noeud[1])
                                                                    elif (coordonnees_noeud[0]+0.75, coordonnees_noeud[1]) not in coordonnees :
                                                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.75, coordonnees_noeud[1])
                                                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]-0.75) not in coordonnees :
                                                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.75)
                                                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.75) not in coordonnees :
                                                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.75)
                                                                    elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                                                    elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                                                    else : 
        #                                                                 fichier.write("probleme\n")
                                                                        print("probleme") 
                                                    for succ in G.successors(noeud) :
                                                        if G.nodes[succ]["coordonnees"] == (33,33) :
                                                            coordonnees = []
                                                            for node in G.nodes() :
                                                                coordonnees.append(G.nodes[node]["coordonnees"])
                                                            for edge in G[noeud][succ] :
                                                                if G[noeud][succ][edge]["label"] == "B53" :
                                                                    if (coordonnees_noeud[0], coordonnees_noeud[1]-0.5) not in coordonnees :
                                                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.5)
                                                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.5) not in coordonnees :
                                                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.5)
                                                                    elif (coordonnees_noeud[0]-0.5, coordonnees_noeud[1]) not in coordonnees :
                                                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.5, coordonnees_noeud[1])
                                                                    elif (coordonnees_noeud[0]+0.5, coordonnees_noeud[1]) not in coordonnees :
                                                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.5, coordonnees_noeud[1])
                                                                    elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                                                    elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                                                    else :
        #                                                                 fichier.write("probleme\n")
                                                                        print("probleme")
                                                                else :
                                                                    if (coordonnees_noeud[0]-0.75, coordonnees_noeud[1]) not in coordonnees :
                                                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.75, coordonnees_noeud[1])
                                                                    elif (coordonnees_noeud[0]+0.75, coordonnees_noeud[1]) not in coordonnees :
                                                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.75, coordonnees_noeud[1])
                                                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]-0.75) not in coordonnees :
                                                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.75)
                                                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.75) not in coordonnees :
                                                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.75)
                                                                    elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                                                    elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                                                    else : 
        #                                                                 fichier.write("probleme\n")
                                                                        print("probleme") 
                                                                    
                                                
        #                                         for noeud in G.nodes() :
        #                                             fichier.write(str(noeud) + " " + str(G.nodes[noeud])+"\n")
        #                                         for u,v,edata in G.edges(data=True) :
        #                                             fichier.write(str(u)+ "'" +str(v) + " " + str(edata["label"])+"\n")
                                                
                                                #fichier.write(str(G.nodes.data())+"\n")
                                                #fichier.write(str(G.edges.data())+"\n")                
                                                
                                                #print(pos)
                                                
                                                #plt.figure(figsize =(5,12))
                                                           
                                                red_edges = [(1,2),(2,1),(1,3),(3,1),(1,5),(5,1),(2,4),(4,2),(3,4),(4,3),(2,5),(5,2)]
                                                green_edges = []
                                                blue_edges = []
                                                black_edges = []
                                                weights = []
                                                #black_edges = [edge for edge in G.edges() if edge not in red_edges]
                                                edges_list = [(u,v) for u,v,data in G.edges(data=True) if data["label"] != '0']#if data["long_range"] != None]
                                                #edges_list = [(u,v) for u,v,data in G.edges(data=True) if data["label"] != '0' and contient_arete((u,v), GC, compteur%2)]
                                                
                                                
                                                for u,v,edata, in G.edges(data=True) :
                                                    if (u,v) in edges_list :
                                                        if (u,v) not in red_edges :
                                #                             if contient_arete((u,v), GC, compteur%2) :
                                #                                 orange_edges.append((u,v))
                                                            if edata["label"] == "B53" :
                                                                green_edges.append((u,v))
                                                            elif edata["label"] == "CWW" :
                                                                blue_edges.append((u,v)) 
                                                            else :
                                                                black_edges.append((u,v))
                                                        
                                                        ''' on mettra en plus gras les aretes en commun '''
                                                        if contient_arete((u,v), GC, compteur_elt%2) :
                                                            weights.append(3)
                                                        else :
                                                            weights.append(1)
                                            
                                                edge_labels=dict([((u,v,),d["label"])for u,v,d in G.edges(data=True) if (u,v) in edges_list])
        
                                                #print(edge_labels)
                                               # node_labels=dict([(u,(d["nt"], d["type"]))for u,d in G.nodes(data=True)])## if d["type"] != None])
                                                                                           #print(node_labels)
                                                
                                                nodes_list = [u for u,d in G.nodes(data=True) if d["type"] != -1]
                                                pos = get_node_attributes(G.subgraph(nodes_list), 'coordonnees')
                                                
                                                
                                                '''on mettra en orange les sommets en commun'''
                                                orange_nodes = []
                                                pink_nodes = []
                                                for noeud in nodes_list :
                                                    if contient_sommet(noeud, GC, compteur_elt%2) :
                                                        orange_nodes.append(noeud)
                                                    else :
                                                        pink_nodes.append(noeud)
                                                node_colors = ['pink' if node in pink_nodes else 'orange' for node in nodes_list]
                                                
                                                node_labels=dict([(u, (u, d["type"], d["poids"]))for u,d in G.nodes(data=True) if d["type"] != None and d["type"] != -1 ])#and u in orange_nodes])
    
                                                #print(nodes_list)
                                                nx.draw_networkx_nodes(G, pos, nodelist=nodes_list, node_size=150, node_color=node_colors)
                                                            #nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), 
                                                            #           node_color = values, node_size = 500)
                                                nx.draw_networkx_labels(G, pos, labels = node_labels)
                                                #nx.draw_networkx_edges(G, pos, edgelist=red_edges, edge_color='r', edge_labels = edge_labels)
                                                #nx.draw_networkx_edges(G, pos, edgelist=blue_edges, edge_color='b', edge_labels = edge_labels)
                                                #nx.draw_networkx_edges(G, pos, edgelist=green_edges, edge_color='g', edge_labels = edge_labels)
                                                #nx.draw_networkx_edges(G, pos, edgelist=black_edges, edge_labels = edge_labels)
                                                
                                                
                                                edge_colors = ['black' if edge in black_edges else 'red' if edge in red_edges else 'blue' if edge in blue_edges else 'green' for edge in edges_list]
                                                
                                                
                                                #nx.draw_networkx_edge_labels(G,pos)
                                                print(len(edges_list))
                                                print(len(weights))
                                                print(weights)
                                                #print(len(edge_labels))
                                                #nx.draw_networkx_edge_labels(G,pos, edge_labels = edge_labels, font_size=6)
                                                nx.draw_networkx_edges(G,pos, edgelist=edges_list, edge_color= edge_colors, width=weights)
                                                #axs[compteur%2].set_axis_off()
                                                #plt.savefig("graphes_extension/"+element[:len(element)-7]+".png") # save as png
                                                #plt.savefig("graphes_extension/fichier_1FJG_A_48_8.png") # save as png
                                                nom = nom + "_" +element[0] +"_"+str(element[1])
                                                compteur_elt = compteur_elt+1
                                                
                                                
                                                plt.axis('off') 
                                        #plt.savefig("/home/coline/Documents/Extensions/poster_macim/comp_commun.png", format="png", transparent=True)
                 
                                        plt.show()
        #                                 if str(sim)+"_"+nom+".png" not in os.listdir(EXTENSION_PATH_TAILLE%taille_ext) :
        #                                     plt.savefig(EXTENSION_PATH_TAILLE%taille_ext+str(sim)+"_"+nom+".png") # save as png
        
                                        #plt.savefig("/home/coline/Bureau/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_taille_%d/"%taille_ext+(str(round(sim,2))) + "_"+fic[20:len(fic)-7]+".png") # save as png
                                        plt.clf()
                                        plt.close()
                                    #break
                                    #break
                            #compteur += 1
                            

if __name__ == '__main__':
    for i in range(4, 5) :
        draw_isomorphisme(i)             
    