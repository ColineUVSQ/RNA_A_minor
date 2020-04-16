'''
Created on 19 juin 2019

@author: coline
'''
import pickle
import networkx as nx
import itertools
import numpy as np
import matplotlib.pyplot as plt

import exrex

import re


from recup_data.graphe_commun_clusters import commun_cluster_clique, draw_network,\
    commun_cluster_clique_new_data, draw_new_data
from recup_data.constantes import GROUPE_ETOILE_GNRA, EXTENSION_PATH_TAILLE,\
    EXTENSION_PATH, GROUPES_TOUTES_ARETES_MAX_4_10_07,\
    CLUSTERING_PEREZ_VERSION_NON_CAN_2, GROUPE_JAUNE_VERT_BLEU,\
    NEW_EXTENSION_PATH_TAILLE
import os
from recup_data.recup_couples_V2 import recup_couples_version_graphe_moyen
from recup_data.sous_graphe_commun_max import sous_graphe_commun_max_version_graphe_moyen
from recup_data.draw_structure import draw_structure
from recup_data.new_algo_comparaison import recup_chaines
import seaborn as sns
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import MultiLabelBinarizer
import copy
import time

def recherche_graphe_commun(liste_tailles, groupe, seuil, type_sim, rep, nom):
    for i in liste_tailles :
        #digraphe_commun = commun_cluster(groupe, EXTENSION_PATH_TAILLE%i, EXTENSION_PATH%i+"dico_comp_complet_metrique_%s_taille_%s.pickle"%(type_sim, i), nom, i, seuil, "taille_max/result_k_max_4_10_toutes_aretes", False)                
        #digraphe_commun, liste_cliques = commun_cluster_clique(groupe, EXTENSION_PATH%i+"dico_comp_complet_metrique_%s_taille_%s.pickle"%(type_sim, i), EXTENSION_PATH_TAILLE%i)
        print(groupe)
        digraphe_commun, liste_cliques = commun_cluster_clique_new_data(groupe, "/media/coline/Maxtor/dico_new_100320_sim_par_branche_0.65.pickle", NEW_EXTENSION_PATH_TAILLE)
        print("rasmousnif")
        for noeud, data in digraphe_commun.nodes(data=True) :
            print(noeud, data)
        
        #draw(digraphe_commun,rep, nom, i, True) 
        #draw_new_data(digraphe_commun, rep, nom, i, True)
        #exit()
        return recup_sous_graphe_global_graphe_commun_new_data(digraphe_commun, groupe)
        #return digraphe_commun
        
def recherche_similaires(taille_ext, rep_digraphe_commun, seuil, nom):
    liste_a_faire  = []
    for fic in os.listdir(EXTENSION_PATH_TAILLE%taille_ext) :
        if ".pickle" in fic and "graphe_comp" not in fic and "couples_possibles" not in fic and len(fic.split("_")) == 5 :
            liste_a_faire.append(fic)
    print(len(liste_a_faire))
    print(liste_a_faire)   
    for elt in liste_a_faire : 
        fichier_couples_possibles = recup_couples_version_graphe_moyen(elt, taille_ext, rep_digraphe_commun, seuil, nom)
        print(fichier_couples_possibles)
        sous_graphe_commun_max_version_graphe_moyen(fichier_couples_possibles, taille_ext, rep_digraphe_commun, seuil, nom)
  
def idem_graphes_communs(digraphe1, digraphe2):
    for noeud, data in digraphe1.nodes(data=True) :
        if noeud not in digraphe2.nodes() :
            print("pb1")
            return False
        for attr in data.keys() :
            if attr not in ['espacement_motif', 'position'] and data[attr] != digraphe2.nodes[noeud][attr] :
                print(noeud)
                print(data)
                print(digraphe2.nodes[noeud])
                print(attr)
                print("pb2")
                return False
    for u,v,e, data in digraphe1.edges(keys = True, data=True) :
        if (u,v,e) not in digraphe2.edges(keys=True) :
            print("pb3")
            return False
        for attr in data.keys() :
            if data[attr] != digraphe2[u][v][e][attr] :
                print("pb4")
                return False
            
    for noeud, data in digraphe2.nodes(data=True) :
        if noeud not in digraphe1.nodes() :
            print("pb5")
            return False
        for attr in data.keys() :
            if attr not in ['espacement_motif', 'position'] and data[attr] != digraphe1.nodes[noeud][attr] :
                print(noeud)
                print(data)
                print(digraphe2.nodes[noeud])
                print(attr)
                print("pb6")
                return False
    for u,v,e, data in digraphe2.edges(keys = True, data=True) :
        if (u,v,e) not in digraphe1.edges(keys=True) :
            print("pb7")
            return False
        for attr in data.keys() :
            if data[attr] != digraphe1[u][v][e][attr] :
                print("pb8")
                return False
            
    return True

def idem_graphes_communs_version_clique(digraphe1, digraphe2):
    deja_trouves = []
    for noeud1, data1 in digraphe1.nodes(data=True) :
        ok = False
#         print(noeud1)
        for noeud2, data2 in digraphe2.nodes(data=True) :
#             print(noeud2)
#             if noeud1 == 9 and noeud2 == 10 :
#                 print(data1["type"])
#                 print(data2["type"])
#                 print(data1["poids"])
#                 print(data2["poids"])
#                 print(data1["chaine"])
#                 print(data2["chaine"])
#                 print(data1["nb_non_can"])
#                 print(data2["nb_non_can"])
#                 print(deja_trouves)
            if data1["type"] == data2["type"] and data1["poids"] == data2["poids"] and data1["chaine"] == data2["chaine"] and noeud2 not in deja_trouves:
                deja_trouves.append(noeud2)
                ok = True
                break
        if not ok :
            print("probleme1")
            print(noeud1)
            print(deja_trouves)
            return False
    
    noeuds_a_trouves = []    
    for u2, v2 in digraphe2.edges() : 
        noeuds_a_trouves.append((u2,v2))
    
    for u1, v1, data1 in digraphe1.edges(data=True) :
        ok = False
        for u2, v2, data2 in digraphe2.edges(data=True) :
            if data1["label"] == data2["label"] and (data1["long_range"] == data2["long_range"] or (data1["long_range"] == None and data2["long_range"] == None)) and (u2,v2) in noeuds_a_trouves :
#                 if (u2,v2) == (4,22) :
#                     print("ramou")
#                     print((u1,v1))
#                     print(deja_trouves)
                noeuds_a_trouves.remove((u2,v2))
                ok = True
                break
        if not ok :
            print("probleme2")
            print((u1,v1))
            print(deja_trouves)
            return False
    
    
    deja_trouves = []
    for noeud2, data2 in digraphe2.nodes(data=True) :
        ok = False
        for noeud1, data1 in digraphe1.nodes(data=True) :
            if data1["type"] == data2["type"] and data1["poids"] == data2["poids"] and data1["chaine"] == data2["chaine"] and noeud1 not in deja_trouves:
                deja_trouves.append(noeud1)
                ok = True
                break
        if not ok :
            print("probleme3")
            print(noeud1)
            print(noeud2)
        
            return False
    
    noeuds_a_trouves = []    
    for u1, v1 in digraphe1.edges() : 
        noeuds_a_trouves.append((u1,v1)) 
    for u2, v2, data2 in digraphe2.edges(data=True) :
        ok = False
        for u1, v1, data1 in digraphe1.edges(data=True) :
            if data1["label"] == data2["label"] and (data1["long_range"] == data2["long_range"] or (data1["long_range"] == None and data2["long_range"] == None)) and (u1,v1) in noeuds_a_trouves :
                noeuds_a_trouves.remove((u1,v1))
                ok = True
                break
        if not ok :
            print("probleme4")
            print((u2,v2))
            print(deja_trouves)
            return False
    
    return True

def idem_graphes_communs_version_clique_new_data(digraphe1, digraphe2):
    deja_trouves = []
    for noeud1, data1 in digraphe1.nodes(data=True) :
        ok = False
#         print(noeud1)
        for noeud2, data2 in digraphe2.nodes(data=True) :
#             print(noeud2)
#             if noeud1 == 9 and noeud2 == 10 :
#                 print(data1["type"])
#                 print(data2["type"])
#                 print(data1["poids"])
#                 print(data2["poids"])
#                 print(data1["chaine"])
#                 print(data2["chaine"])
#                 print(data1["nb_non_can"])
#                 print(data2["nb_non_can"])
#                 print(deja_trouves)
            if data1["type"] == data2["type"] and data1["poids"] == data2["poids"] and data1["chaine"] == data2["chaine"] and noeud2 not in deja_trouves:
                deja_trouves.append(noeud2)
                ok = True
                break
        if not ok :
            print("probleme1")
            print(noeud1)
            print(deja_trouves)
            return False
    
    noeuds_a_trouves = []    
    for u2, v2 in digraphe2.edges() : 
        noeuds_a_trouves.append((u2,v2))
    
    for u1, v1, data1 in digraphe1.edges(data=True) :
        ok = False
        for u2, v2, data2 in digraphe2.edges(data=True) :
            if data1["label"] == data2["label"] and (u2,v2) in noeuds_a_trouves :
#                 if (u2,v2) == (4,22) :
#                     print("ramou")
#                     print((u1,v1))
#                     print(deja_trouves)
                noeuds_a_trouves.remove((u2,v2))
                ok = True
                break
        if not ok :
            print("probleme2")
            print((u1,v1))
            print(deja_trouves)
            return False
    
    
    deja_trouves = []
    for noeud2, data2 in digraphe2.nodes(data=True) :
        ok = False
        for noeud1, data1 in digraphe1.nodes(data=True) :
            if data1["type"] == data2["type"] and data1["poids"] == data2["poids"] and data1["chaine"] == data2["chaine"] and noeud1 not in deja_trouves:
                deja_trouves.append(noeud1)
                ok = True
                break
        if not ok :
            print("probleme3")
            print(noeud1)
            print(noeud2)
        
            return False
    
    noeuds_a_trouves = []    
    for u1, v1 in digraphe1.edges() : 
        noeuds_a_trouves.append((u1,v1)) 
    for u2, v2, data2 in digraphe2.edges(data=True) :
        ok = False
        for u1, v1, data1 in digraphe1.edges(data=True) :
            if data1["label"] == data2["label"] and (u1,v1) in noeuds_a_trouves :
                noeuds_a_trouves.remove((u1,v1))
                ok = True
                break
        if not ok :
            print("probleme4")
            print((u2,v2))
            print(deja_trouves)
            return False
    
    return True

def recup_sous_graphe_global_graphe_commun(graphe_commun, groupe):
    with open("graphs_2.92.pickle", 'rb') as fichier_tout :
        mon_depickler_graphes = pickle.Unpickler(fichier_tout)
        graphes = mon_depickler_graphes.load()
        
        dico_graphes = {}
        for elt in groupe :
            motif = []
            cle = (elt.split("_")[0], elt.split("_")[1])
            graphe = nx.MultiDiGraph()
            
            for noeud, data in graphe_commun.nodes(data=True) :
                if data["type"] != -1 :
                    
                    for pos in graphe_commun.nodes[data["chaine"][0]]["num_seq"] : ## retrouver la position de l element du motif
                        if pos[0] == elt :
                            pos_chaine = pos[1]
                    
                    for pos in data["num_seq"] :
                        if pos[0] == elt :
                            positions = pos[1]
                            
                    if data["type"] in [11,12,13,14,15] :
                        motif.append(positions[0])
                    if data["type"] != 1 and data["type"] != 0 :
                        graphes[cle].nodes[positions[0]].update({"chaine" : data["chaine"], "position" : [positions[0]- pos_chaine[0]], "type" : data["type"], "en_plus" : False})
                        graphe.add_node(positions[0], **graphes[cle].nodes[positions[0]])
                    elif data["type"] == 0 :
                        for i in range(positions[0], positions[0]+data["poids"]) :
                            graphes[cle].nodes[i].update({"chaine" : data["chaine"], "position" : [i- pos_chaine[0]], "type" : data["type"], "en_plus" : False})
                            graphe.add_node(i, **graphes[cle].nodes[i])
                            if i < positions[0]+data["poids"]-1 :
                                graphe.add_edge(i, i+1, label='B53', long_range=False)
                                
                        for i in range(positions[0]+data["poids"], positions[1]+1) :
                            graphes[cle].nodes[i].update({"chaine" : data["chaine"], "position" : [i- pos_chaine[0]], "type" : data["type"], "en_plus" : True})
                            graphe.add_node(i, **graphes[cle].nodes[i])
                            if i < positions[0]+data["poids"]-1 :
                                graphe.add_edge(i, i+1, label='B53', long_range=False)
                    else :
                        for i in range(positions[0], positions[0]+data["poids"]) :
                            if i not in graphe.nodes() :
                                graphes[cle].nodes[i].update({"chaine" : data["chaine"], "position" : [i-pos_chaine[0]], "type" : data["type"], "en_plus" : False})
                                graphe.add_node(i, **graphes[cle].nodes[i])
                            if i >  positions[0]:
                                    graphe.add_edge(i-1, i, label='B53', long_range=False)
                                
                            for voisin in graphes[cle][i] :
                                if graphes[cle].edges[i,voisin]["label"] == 'CWW' :
                                    if voisin not in graphe.nodes() :
                                        graphes[cle].nodes[voisin].update({"chaine" : data["chaine"], "position" : [i-pos_chaine[0]], "type" : None, "en_plus" : False})
                                        graphe.add_node(voisin, **graphes[cle].nodes[voisin])
                                    deja_mis = False
                                    if (i,voisin) in graphe.edges() :
                                        for edge in graphe[i][voisin] :
                                            if graphe[i][voisin][edge]["label"] == 'CWW' :
                                                deja_mis = True
                                    if (i,voisin) not in graphe.edges() or not deja_mis :
                                            graphe.add_edge(i, voisin, label='CWW', long_range=False)
                                    
                                    deja_mis = False
                                    if (voisin, i) in graphe.edges() :
                                        for edge in graphe[voisin][i] :
                                            if graphe[voisin][i][edge]["label"] == 'CWW' :
                                                deja_mis = True
                                    if (voisin, i) not in graphe.edges() or not deja_mis :
                                            graphe.add_edge(voisin, i, label='CWW', long_range=False)
                        for i in range(positions[0]+data["poids"], positions[1]+1) :
                            if i not in graphe.nodes() :
                                graphes[cle].nodes[i].update({"chaine" : data["chaine"], "position" : [i-pos_chaine[0]], "type" : data["type"], "en_plus" : True})
                                graphe.add_node(i, **graphes[cle].nodes[i])
                            if i >  positions[0]:
                                graphe.add_edge(i-1, i, label='B53', long_range=False)
                                
                            for voisin in graphes[cle][i] :
                                if graphes[cle].edges[i,voisin]["label"] == 'CWW' :
                                    if voisin not in graphe.nodes() :
                                        graphes[cle].nodes[voisin].update({"chaine" : data["chaine"], "position" : [i-pos_chaine[0]], "type" : None, "en_plus" : True})
                                        graphe.add_node(voisin, **graphes[cle].nodes[voisin])
                                    
                                    deja_mis = False
                                    if (i,voisin) in graphe.edges() :
                                        for edge in graphe[i][voisin] :
                                            if graphe[i][voisin][edge]["label"] == 'CWW' :
                                                deja_mis = True
                                    if (i,voisin) not in graphe.edges() or not deja_mis :
                                            graphe.add_edge(i, voisin, label='CWW', long_range=False)
                                    
                                    deja_mis = False
                                    if (voisin, i) in graphe.edges() :
                                        for edge in graphe[voisin][i] :
                                            if graphe[voisin][i][edge]["label"] == 'CWW' :
                                                deja_mis = True
                                    if (voisin, i) not in graphe.edges() or not deja_mis :
                                            graphe.add_edge(voisin, i, label='CWW', long_range=False)
            
            for u,v,data in graphe_commun.edges(data=True) :
                if graphe_commun.nodes[u]["type"] not in [0, 1, -1] and graphe_commun.nodes[v]["type"] not in [0, 1, -1] :
                    for pos in graphe_commun.nodes[u]["num_seq"] :
                        if pos[0] == elt :
                            positions_1 = pos[1]
                            
                    for pos in graphe_commun.nodes[v]["num_seq"] :
                        if pos[0] == elt :
                            positions_2 = pos[1]
                    
                    graphe.add_edge(positions_1[0], positions_2[0], **data)
                    
            for noeud in graphe.nodes() :
                if noeud+1 in graphe.nodes() :
                    if (noeud,noeud+1) not in graphe.edges() :
                        graphe.add_edge(noeud, noeud+1, label='B53', long_range=False)
                    else :
                        deja_mis = False
                        for edge in graphe[noeud][noeud+1] :
                            if graphe[noeud][noeud+1][edge]["label"] == 'B53' :
                                deja_mis = True
                        if not deja_mis :
                            graphe.add_edge(noeud, noeud+1, label='B53', long_range=False)
                
            
            
            dico_graphes.update({elt : (graphe, motif)})
    
    for cle in dico_graphes.keys() :
        print(dico_graphes[cle][0].nodes.data())
    return dico_graphes                                        



def recup_sous_graphe_global_graphe_commun_new_data(graphe_commun, groupe):
#     with open("graphs_2.92.pickle", 'rb') as fichier_tout :
#         mon_depickler_graphes = pickle.Unpickler(fichier_tout)
#         graphes = mon_depickler_graphes.load()
    liste_graphes = []
    for elt in groupe :
        with open("Graphs/%s.pickle"%elt[0], 'rb') as fichier_tout :
            mon_depickler_graphes = pickle.Unpickler(fichier_tout)
            graphe = mon_depickler_graphes.load()
            liste_graphes.append((elt, graphe))
    
    dico_graphes = {}
    ## le premier graphe de dico graphe sera la version etendue (decontraction des noeuds contractes) du graphe commun mais ne correspondra a aucun vrai graphe
    ## avec l'ajout pour chaque sommet des nucleotides possibles
    graphe_commun_global = nx.MultiDiGraph()
    #first_compteur = max(list(graphe_commun.nodes())) + 1
    ### Traitement des noeuds
    for noeud, data in graphe_commun.nodes(data=True) :
        print("noeud")
        print(noeud)
        print(data["chaine"])
        print(data["type"])
        if data["type"] not in [0,1,-1] : ## cas facile ou il n'y a pas de decontraction a faire, on recherche juste les nucleotides correspondant au sommet dans chaque graphe du cluster
            liste_nts = {}
            for g in liste_graphes :
                for pos in data["num_seq"] :
                    if pos[0] == g[0] :
                        if g[1].nodes[(pos[2], pos[1][0])]["nt"] not in liste_nts.keys() :
                            liste_nts.update({g[1].nodes[(pos[2], pos[1][0])]["nt"] : 1})
                        else :
                            liste_nts[g[1].nodes[(pos[2], pos[1][0])]["nt"]] += 1
            #data_new = copy.deepcopy(data)
            data_new = {}
            for cle in data.keys() :
                if isinstance(data[cle], list) :
                    data_new.update({cle : list(data[cle])})
                else :
                    data_new.update({cle : data[cle]})
            data_new.update({"liste_nts" : copy.deepcopy(liste_nts)})
            
            graphe_commun_global.add_node(noeud, **data_new)
        elif data["type"] == 0 or data["type"] == 1 : ## cas ou il faut peut-etre faire une decontraction des sommets
            #data_new = copy.deepcopy(data)
            data_new = {}
            for cle in data.keys() :
                if isinstance(data[cle], list) :
                    data_new.update({cle : list(data[cle])})
                else :
                    data_new.update({cle : data[cle]})
            
            data_new["poids"] = 1
            graphe_commun_global.add_node(noeud, **data_new)
            if data["poids"] > 1 :
                ## on ajoute autant de sommets que le poids du sommet de depart et on les lie par des liaisons B53
                liste_compteurs = []
                compteur = max(max(list(graphe_commun.nodes())), max(list(graphe_commun_global.nodes())))+1
                first_compteur = compteur
                
                for i in range(1, data["poids"]) :
                    print(compteur)
                    data_new_new = {}
                    for cle in data_new.keys() :
                        if isinstance(data_new[cle], list) :
                            data_new_new.update({cle : list(data_new[cle])})
                        else :
                            data_new_new.update({cle : data_new[cle]})
                    
                    data_new_new["poids"] = 1
                    graphe_commun_global.add_node(compteur, **data_new_new)
                    liste_compteurs.append(compteur)
                    
                    if i == 1 :
                        if 1 in data["chaine"] or 4 in data["chaine"] :
                            graphe_commun_global.nodes[compteur]["position"][0] += i
                            graphe_commun_global.add_edge(noeud, compteur, label="B53", near=False)
                        else :
                            graphe_commun_global.nodes[compteur]["position"][0] -= i
                            graphe_commun_global.add_edge(compteur, noeud, label="B53", near=False)
                    else :
                        if 1 in data["chaine"] or 4 in data["chaine"] :
                            graphe_commun_global.add_edge(compteur-1, compteur, label="B53", near=False)
                            graphe_commun_global.nodes[compteur]["position"][0] += i
                        else :
                            graphe_commun_global.add_edge(compteur, compteur-1, label="B53", near=False)
                            graphe_commun_global.nodes[compteur]["position"][0] -= i
                    compteur += 1
                
                for noeud_pos, data_pos in graphe_commun.nodes(data=True) :
                    if data_pos["type"] not in [None, -1] and len([x for x in data_pos["chaine"] if x in data["chaine"]]) > 0 and noeud_pos not in liste_compteurs :
                        if (1 in data["chaine"] or 4 in data["chaine"]) and data_pos["position"][0] > data["position"][0] :
                            data_pos["position"][0] += data["poids"]
                        elif (2 in data["chaine"] or 3 in data["chaine"]) and data_pos["position"][0] < data["position"][0] :
                            data_pos["position"][0] -= data["poids"]
                            
                    
            else :
                compteur = max(max(list(graphe_commun.nodes())), max(list(graphe_commun_global.nodes())))+1
                first_compteur = compteur
            
            ## on recherche les nts correspondant a chaque sommet, en autorisant les superpositions
            compteur = noeud
            for i in range(0, data["poids"]) :
                if (1 in data["chaine"] or 4 in data["chaine"]) :
                    liste_nts = {}
                    for g in liste_graphes :
                        for pos in data["num_seq"] :
                            if pos[0] == g[0] :
                                for p in range(i+pos[1][0], pos[1][1]-data["poids"]+2+i) :
                                    if g[1].nodes[(pos[2], p)]["nt"] not in liste_nts.keys() :
                                        liste_nts.update({g[1].nodes[(pos[2], p)]["nt"] : 1})
                                    else :
                                        liste_nts[g[1].nodes[(pos[2], p)]["nt"]] += 1
                    graphe_commun_global.nodes[compteur].update({"liste_nts": copy.deepcopy(liste_nts)})
                else :
                    liste_nts = {}
                    for g in liste_graphes :
                        for pos in data["num_seq"] :
                            if pos[0] == g[0] :
                                for p in np.arange(pos[1][1]-i, pos[1][0]+data["poids"]-2-i, -1) :
                                    if g[1].nodes[(pos[2], p)]["nt"] not in liste_nts.keys() :
                                        liste_nts.update({g[1].nodes[(pos[2], p)]["nt"] : 1})
                                    else :
                                        liste_nts[g[1].nodes[(pos[2], p)]["nt"]] += 1
                            
                    graphe_commun_global.nodes[compteur].update({"liste_nts": copy.deepcopy(liste_nts)})
                if compteur == noeud :
                    compteur = first_compteur
                else :
                    compteur += 1
        elif data["type"] == -1 : ## cas des sommets artificiels (on ajoutera ceux lies aux sommets de type 1 uniquement)
            for voisin in graphe_commun[noeud] : ## normalement il n'y en a qu'un seul
                if graphe_commun.nodes[voisin]["type"] == 1 :
                    #data_new = copy.deepcopy(data)
                    data_new = {}
                    for cle in data.keys() :
                        if isinstance(data[cle], list) :
                            data_new.update({cle : list(data[cle])})
                        else :
                            data_new.update({cle : data[cle]})
                    data_new["poids"] = 1
                    graphe_commun_global.add_node(noeud, **data_new)
                    if data["poids"] > 1 :
                        ## idem qu'au-dessus, on ajoute autant de sommets que le poids du sommet de depart
                        compteur = max(max(list(graphe_commun.nodes())), max(list(graphe_commun_global.nodes())))+1
                        first_compteur = compteur
                        for i in range(1, data["poids"]) :
                            data_new["poids"] = 1
                            graphe_commun_global.add_node(compteur, **data_new)
                            if i == 1 :
                                if 1 in data["chaine"] or 4 in data["chaine"] :
                                    graphe_commun_global.add_edge(noeud, compteur, label="B53", near=False)
                                else :
                                    graphe_commun_global.add_edge(compteur, noeud, label="B53", near=False)
                            else :
                                if 1 in data["chaine"] or 4 in data["chaine"] :
                                    graphe_commun_global.add_edge(compteur-1, compteur, label="B53", near=False)
                                else :
                                    graphe_commun_global.add_edge(compteur, compteur-1, label="B53", near=False)
                            compteur += 1
                    else :
                        compteur = max(max(list(graphe_commun.nodes())), max(list(graphe_commun_global.nodes())))+1
                        first_compteur = compteur
                    
                    ## on recherche les nts correspondant en se servant des voisins (car les numeros des sommets des graphes globaux ne sont pas stockes dans l'extension)            
                    compteur = noeud
                    for i in range(0, data["poids"]) :
                        if (1 in data["chaine"] or data["chaine"] == 4) :
                            liste_nts = {}
                            for g in liste_graphes :
                                for pos in graphe_commun.nodes[voisin]["num_seq"] :
                                    if pos[0] == g[0] :
                                        for p in range(i+pos[1][0], pos[1][1]-data["poids"]+2+i) :
                                            for voisin_voisin in g[1][(pos[2], p)] :
                                                if g[1].edges[(pos[2], p),voisin_voisin]["label"] == "CWW" : ## la aussi en principe il n'y en a qu'un
                                                        if g[1].nodes[voisin_voisin]["nt"] not in liste_nts.keys() :
                                                            liste_nts.update({g[1].nodes[voisin_voisin]["nt"] : 1})
                                                        else :
                                                            liste_nts[g[1].nodes[voisin_voisin]["nt"]] += 1
                            graphe_commun_global.nodes[compteur].update({"liste_nts": copy.deepcopy(liste_nts)})
                        else :
                            liste_nts = {}
                            for g in liste_graphes :
                                for pos in graphe_commun.nodes[voisin]["num_seq"] :
                                    if pos[0] == g[0] :
                                        for p in np.arange(pos[1][1]-i, pos[1][0]+data["poids"]-2-i, -1) :
                                            for voisin_voisin in g[1][(pos[2], p)] :
                                                #for edge in g[1][(pos[2], p)][voisin_voisin] :
                                                if g[1].edges[(pos[2], p), voisin_voisin]["label"] == "CWW" : ## la aussi en principe il n'y en a qu'un
                                                        if g[1].nodes[voisin_voisin]["nt"] not in liste_nts.keys() :
                                                            liste_nts.update({g[1].nodes[voisin_voisin]["nt"] : 1})
                                                        else :
                                                            liste_nts[g[1].nodes[voisin_voisin]["nt"]] += 1
                                    
                            graphe_commun_global.nodes[compteur].update({"liste_nts": copy.deepcopy(liste_nts)})
                        if compteur == noeud :
                            compteur = first_compteur
                        else :
                            compteur += 1
        print("rapala")
        print(graphe_commun_global.nodes())
    ## Traitement des aretes
    for u,v,data in graphe_commun.edges(data=True) :
        if u in graphe_commun_global.nodes() and v in graphe_commun_global.nodes() : ## on ne rajoute pas les aretes entre des sommets de type 0 et leurs sommets artificiels car on n'a pas mis ces sommets artificiels là
            if graphe_commun.nodes[u]["poids"] == 1 or data["label"] != '0' : ## cas facile ou il n'y a ni decontraction, ni aretes artificielles
                graphe_commun_global.add_edge(u,v,**data)
#                 if u == 27 or v == 27 :
#                     print("raaaaah")
#                     print(graphe_commun.nodes[u])
#                     print(graphe_commun.nodes[v])
#                     print(data)
#                 print(u,v)
            elif graphe_commun.nodes[u]["type"] != 0 :## cas des aretes artificielles liees a des sommets de type 1 dont il faut s'occuper
                ## on ajoute les aretes CWW entre toutes les paires de sommets necessaires
    
#                 if u not in graphe_commun_global.nodes() or v not in graphe_commun_global.nodes() :
#                     print("noeud n'existe pas, bizarre")
                graphe_commun_global.add_edge(u,v, label='CWW', near=False)
#                 print(u,v)
                compteur_u = u
                compteur_v = v
                a_ajouter = []
                print(u)
                print(graphe_commun.nodes[u]["poids"])
                for i in range(1, graphe_commun.nodes[u]["poids"]) :
                    print("ouhou")
                    if (1 in graphe_commun.nodes[u]["chaine"] or 4 in graphe_commun.nodes[u]["chaine"]) :
                        for voisin_u in graphe_commun_global[compteur_u] :
                            for edge_u in graphe_commun_global[compteur_u][voisin_u] :
                                if graphe_commun_global[compteur_u][voisin_u][edge_u]["label"] == 'B53' :
                                    for voisin_v in graphe_commun_global[compteur_v] :
                                        for edge_v in graphe_commun_global[compteur_v][voisin_v] :
                                            if graphe_commun_global[compteur_v][voisin_v][edge_v]["label"] == 'B53' :
                                                a_ajouter.append((voisin_u, voisin_v))
                    else :
                        for voisin_u in graphe_commun_global.predecessors(compteur_u) :
                            for edge_u in graphe_commun_global[voisin_u][compteur_u] :
                                if graphe_commun_global[voisin_u][compteur_u][edge_u]["label"] == 'B53' :
                                    for voisin_v in graphe_commun_global.predecessors(compteur_v) :
                                        for edge_v in graphe_commun_global[voisin_v][compteur_v] :
                                            if graphe_commun_global[voisin_v][compteur_v][edge_v]["label"] == 'B53' :
                                                a_ajouter.append((voisin_u, voisin_v))
                    compteur_u = a_ajouter[len(a_ajouter)-1][0]
                    compteur_v = a_ajouter[len(a_ajouter)-1][1]
                    
                
                print("rapoulou")
                print(a_ajouter)
                print(graphe_commun_global.edges())
                print(graphe_commun_global.nodes())
                for elt in a_ajouter :
#                     if elt[0] not in graphe_commun_global.nodes() or elt[1] not in graphe_commun_global.nodes() :
#                         print("noeud n'existe pas, bizarre 2")
                    print(elt)
                    graphe_commun_global.add_edge(elt[0], elt[1], label='CWW', near=False)
                    #graphe_commun_global.add_edge(elt[1], elt[0], label='CWW', near=False)
                
    ## on enleve les liaisons b53 entre sommets artificielles lies aux sommets de type 1 (car on n'est pas surs qu'elles existent)
    a_enlever = []
    for u,v,data in graphe_commun_global.edges(data=True) :
        if data["label"] == "B53" and graphe_commun_global.nodes[u]["type"] == -1 and graphe_commun_global.nodes[v]["type"] == -1 :
            a_enlever.append((u,v))
    
    for elt in a_enlever :
        graphe_commun_global.remove_edge(elt[0], elt[1])                                  
                                
    ## on recupere chaque graphe global aussi
    for elt in groupe :
        motif = []
#             cle = (elt.split("_")[0], elt.split("_")[1])
            
        with open("Graphs/%s.pickle"%elt[0], 'rb') as fichier_tout :
            mon_depickler_graphes = pickle.Unpickler(fichier_tout)
            graphes = mon_depickler_graphes.load()
            
            graphe = nx.MultiDiGraph()
            
            for noeud, data in graphe_commun.nodes(data=True) :
                if data["type"] != -1 :
                    
                    for pos in graphe_commun.nodes[data["chaine"][0]]["num_seq"] : ## retrouver la position de l element du motif
                        if pos[0] == elt :
                            pos_chaine = pos[1]
                            
                    
                    for pos in data["num_seq"] :
                        if pos[0] == elt :
                            positions = pos[1]
                            num_ch = pos[2]
                            
                    if data["type"] in [11,12,13,14,15] :
                        motif.append((num_ch,positions[0]))
                    if data["type"] != 1 and data["type"] != 0 :
                        graphes.nodes[(num_ch,positions[0])].update({"chaine" : data["chaine"], "position" : [positions[0]- pos_chaine[0]], "type" : data["type"], "en_plus" : False, "espacement_motif" : data["espacement_motif"], "poids" : data["poids"]})
                        graphe.add_node((num_ch,positions[0]), **graphes.nodes[(num_ch,positions[0])])
                    elif data["type"] == 0 :
                        for i in range(positions[0], positions[0]+data["poids"]) :
                            graphes.nodes[(num_ch,i)].update({"chaine" : data["chaine"], "position" : [i- pos_chaine[0]], "type" : data["type"], "en_plus" : False, "espacement_motif" : data["espacement_motif"], "poids" : data["poids"]})
                            graphe.add_node((num_ch,i), **graphes.nodes[(num_ch,i)])
                            if i < positions[0]+data["poids"]-1 :
                                graphe.add_edge((num_ch,i), (num_ch,i+1), label='B53', long_range=False)

                        for i in range(positions[0]+data["poids"], positions[1]+1) :
                            graphes.nodes[(num_ch,i)].update({"chaine" : data["chaine"], "position" : [i- pos_chaine[0]], "type" : data["type"], "en_plus" : True, "espacement_motif" : data["espacement_motif"], "poids" : data["poids"]})
                            graphe.add_node((num_ch,i), **graphes.nodes[(num_ch,i)])
                            if i < positions[0]+data["poids"]-1 :
                                graphe.add_edge((num_ch,i), (num_ch,i+1), label='B53', long_range=False)
                    else :
                        for i in range(positions[0], positions[0]+data["poids"]) :
                            if (num_ch,i) not in graphe.nodes() :
                                graphes.nodes[(num_ch,i)].update({"chaine" : data["chaine"], "position" : [i-pos_chaine[0]], "type" : data["type"], "en_plus" : False, "espacement_motif" : data["espacement_motif"], "poids" : data["poids"]})
                                graphe.add_node((num_ch,i), **graphes.nodes[(num_ch,i)])
                            if i >  positions[0]:
                                    graphe.add_edge((num_ch,i-1), (num_ch,i), label='B53', long_range=False)
                                
                            for voisin in graphes[(num_ch,i)] :
                                if graphes.edges[(num_ch,i),voisin]["label"] == 'CWW' :
                                    if voisin not in graphe.nodes() :
                                        graphes.nodes[voisin].update({"chaine" : data["chaine"], "position" : [i-pos_chaine[0]], "type" : None, "en_plus" : False, "espacement_motif" : data["espacement_motif"], "poids" : data["poids"]})
                                        graphe.add_node(voisin, **graphes.nodes[voisin])
                                    deja_mis = False
                                    if ((num_ch,i),voisin) in graphe.edges() :
                                        for edge in graphe[(num_ch,i)][voisin] :
                                            if graphe[(num_ch,i)][voisin][edge]["label"] == 'CWW' :
                                                deja_mis = True
                                    if ((num_ch,i),voisin) not in graphe.edges() or not deja_mis :
                                            graphe.add_edge((num_ch,i), voisin, label='CWW', long_range=False)
                                    
                                    deja_mis = False
                                    if (voisin, (num_ch,i)) in graphe.edges() :
                                        for edge in graphe[voisin][(num_ch,i)] :
                                            if graphe[voisin][(num_ch,i)][edge]["label"] == 'CWW' :
                                                deja_mis = True
                                    if (voisin, (num_ch,i)) not in graphe.edges() or not deja_mis :
                                            graphe.add_edge(voisin, (num_ch,i), label='CWW', long_range=False)
                        for i in range(positions[0]+data["poids"], positions[1]+1) :
                            if (num_ch,i) not in graphe.nodes() :
                                graphes.nodes[(num_ch,i)].update({"chaine" : data["chaine"], "position" : [i-pos_chaine[0]], "type" : data["type"], "en_plus" : True, "espacement_motif" : data["espacement_motif"], "poids" : data["poids"]})
                                graphe.add_node((num_ch,i), **graphes.nodes[(num_ch,i)])
                            if i >  positions[0]:
                                graphe.add_edge((num_ch,i-1), (num_ch,i), label='B53', long_range=False)
                                
                            for voisin in graphes[(num_ch,i)] :
                                if graphes.edges[(num_ch,i),voisin]["label"] == 'CWW' :
                                    if voisin not in graphe.nodes() :
                                        graphes.nodes[voisin].update({"chaine" : data["chaine"], "position" : [i-pos_chaine[0]], "type" : None, "en_plus" : True, "espacement_motif" : data["espacement_motif"], "poids" : data["poids"]})
                                        graphe.add_node(voisin, **graphes.nodes[voisin])
                                    
                                    deja_mis = False
                                    if ((num_ch,i),voisin) in graphe.edges() :
                                        for edge in graphe[(num_ch,i)][voisin] :
                                            if graphe[(num_ch,i)][voisin][edge]["label"] == 'CWW' :
                                                deja_mis = True
                                    if ((num_ch,i),voisin) not in graphe.edges() or not deja_mis :
                                            graphe.add_edge((num_ch,i), voisin, label='CWW', long_range=False)
                                    
                                    deja_mis = False
                                    if (voisin, (num_ch,i)) in graphe.edges() :
                                        for edge in graphe[voisin][(num_ch,i)] :
                                            if graphe[voisin][(num_ch,i)][edge]["label"] == 'CWW' :
                                                deja_mis = True
                                    if (voisin, (num_ch,i)) not in graphe.edges() or not deja_mis :
                                            graphe.add_edge(voisin, (num_ch,i), label='CWW', long_range=False)
            
            for u,v,data in graphe_commun.edges(data=True) :
                if graphe_commun.nodes[u]["type"] not in [0, 1, -1] and graphe_commun.nodes[v]["type"] not in [0, 1, -1] :

                    for pos in graphe_commun.nodes[u]["num_seq"] :
                        if pos[0] == elt :
                            positions_1 = pos[1]
                            num_ch1 = pos[2]
                            
                    for pos in graphe_commun.nodes[v]["num_seq"] :
                        if pos[0] == elt :
                            positions_2 = pos[1]
                            num_ch2 = pos[2]
                    
                    graphe.add_edge((num_ch1,positions_1[0]), (num_ch2,positions_2[0]), **data)
            print(graphe.nodes())      
            for noeud in graphe.nodes() :
                if (noeud[0], noeud[1]+1) in graphe.nodes() :
                    if (noeud,(noeud[0], noeud[1]+1)) not in graphe.edges() :
                        graphe.add_edge(noeud, (noeud[0], noeud[1]+1), label='B53', long_range=False)
                    else :
                        deja_mis = False
                        for edge in graphe[noeud][(noeud[0], noeud[1]+1)] :
                            if graphe[noeud][(noeud[0], noeud[1]+1)][edge]["label"] == 'B53' :
                                deja_mis = True
                        if not deja_mis :
                            graphe.add_edge(noeud, (noeud[0], noeud[1]+1), label='B53', long_range=False)
                
            
            
            dico_graphes.update({elt : (graphe, motif)})
            
#     noeuds_vus = []
#     for noeud, data in graphe_commun_global.nodes(data=True) :
#         if data["type"] not in [None, -1] and noeud not in noeuds_vus:
#             noeuds_vus.append(noeud)
#             if noeud not in [1,2,3,4,5] :
#                 a_un_pred = False
#                 for voisin in graphe_commun_global.predecessors(noeud) :
#                     for edge in graphe_commun_global[voisin][noeud] :
#                         if graphe_commun_global[voisin][noeud][edge]["label"] == 'B53' :
#                             a_un_pred = True
#             
#             if noeud in [1,2,3,4,5] or not a_un_pred :
#                 liaison_B53 = True
#                 new_noeud = noeud
#                 while liaison_B53 :
#                     noeud = new_noeud
#                     liaison_B53 = False
#                     for voisin in graphe_commun_global[noeud] :
#                         for edge in graphe_commun_global[noeud][voisin] :
#                             if graphe_commun_global[noeud][voisin][edge]["label"] == 'B53' :
#                                 liaison_B53 = True
#                                 if 1 in graphe_commun_global.nodes[voisin]["chaine"] or 4 in graphe_commun_global.nodes[voisin]["chaine"] :
#                                     graphe_commun_global.nodes[voisin]["position"] = [graphe_commun_global.nodes[noeud]["position"][0]+1]
#                                 if 2 in graphe_commun_global.nodes[voisin]["chaine"] or 3 in graphe_commun_global.nodes[voisin]["chaine"] :
#                                     graphe_commun_global.nodes[voisin]["position"] = [graphe_commun_global.nodes[noeud]["position"][0]+1]
#                                 noeuds_vus.append(voisin)
#                                 new_noeud = voisin
                                
            
                
            
    
    for cle in dico_graphes.keys() :
        print(cle)
        for noeud, data in dico_graphes[cle][0].nodes(data=True) :
            print(noeud, data)
    return graphe_commun_global, dico_graphes                                        

'''02/10/19'''
def recherche_sequence_graphe_commun_groupe_11(groupe):
    liste_sequence = [[]]
    print(liste_sequence)
    compteur = 0
    for elt in groupe :
        with open(NEW_EXTENSION_PATH_TAILLE+elt, 'rb') as fichier_graphe :
            mon_depickler = pickle.Unpickler(fichier_graphe)
            extension = mon_depickler.load()
        
            with open("Graphs/"+elt.split("_")[1]+".pickle", 'rb') as fichier_graphe_tot :
                mon_depickler = pickle.Unpickler(fichier_graphe_tot)
                graphe = mon_depickler.load()
                
                k = 0
                for i in range(extension.nodes[1]["position"][0] - 6, extension.nodes[1]["position"][0]+1) :
                    #print((extension.nodes[1]["num_ch"], i))
                    #print(graphe.nodes())
                    if k != 0 and compteur == 0 :
                        liste_sequence.append([])
                    liste_sequence[k].append(graphe.nodes[(extension.nodes[1]["num_ch"], i)]["nt"])
                    k += 1
                
                print(len(liste_sequence))
                k = 7
                for i in range(extension.nodes[1]["position"][0]+1, extension.nodes[1]["position"][0]+4) :
                    if compteur == 0 :
                        liste_sequence.append([])
                    liste_sequence[k].append(graphe.nodes[(extension.nodes[1]["num_ch"], i)]["nt"])
                    k += 1
                
                print(len(liste_sequence))   
                k = 10
                for i in range(extension.nodes[4]["position"][0], extension.nodes[4]["position"][0]+5) :
                    if compteur == 0 :
                        liste_sequence.append([])
                    liste_sequence[k].append(graphe.nodes[(extension.nodes[4]["num_ch"], i)]["nt"])
                    k += 1
                    
                    if compteur == 0 :
                        liste_sequence.append([])
                    print(len(graphe[(extension.nodes[4]["num_ch"], i)]))
                    for voisin in graphe[(extension.nodes[4]["num_ch"], i)] :
                        if graphe.edges[(extension.nodes[4]["num_ch"], i), voisin]["label"] == 'CWW' and \
                        ((graphe.nodes[(extension.nodes[4]["num_ch"], i)]["nt"] == 'C' and graphe.nodes[voisin]["nt"] == 'G') or \
                          (graphe.nodes[(extension.nodes[4]["num_ch"], i)]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'C') or \
                          (graphe.nodes[(extension.nodes[4]["num_ch"], i)]["nt"] == 'A' and graphe.nodes[voisin]["nt"] == 'U') or \
                          (graphe.nodes[(extension.nodes[4]["num_ch"], i)]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'A') or \
                          (graphe.nodes[(extension.nodes[4]["num_ch"], i)]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'U') or \
                          (graphe.nodes[(extension.nodes[4]["num_ch"], i)]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'G')):
                            #print(elt)
                            #print(extension.nodes[4]["position"][0])
                            print("ramou %d"%k)
                            #print((extension.nodes[4]["num_ch"], i), voisin)
                            #print(graphe.edges[(extension.nodes[4]["num_ch"], i), voisin]["label"])
                            liste_sequence[k].append(graphe.nodes[voisin]["nt"])
                            k += 1
                
        compteur += 1
    
    liste_sequence_frequence = []   
    for pos in liste_sequence : 
        dico = {'C' : 0, 'G' : 0, 'A' : 0, 'U' : 0}
        for elt in pos :
            dico[elt] += 1
        dico_a_garder = copy.deepcopy(dico)    
        liste_sequence_frequence.append(dico_a_garder)
        
    print(len(liste_sequence))
    print(liste_sequence)
    
    print(len(liste_sequence_frequence))
    print(liste_sequence_frequence)
    
    return liste_sequence_frequence

'''02/10/19'''   
def pourcentage_signature_seq_groupe_de_28(liste_sequence_frequence):
    liste_ok = []
    liste_100 = []
    compteur = 0
    for fic in os.listdir(NEW_EXTENSION_PATH_TAILLE) : 
        if "pickle" in fic and "fichier" in fic :
            print("ramou")
            with open(NEW_EXTENSION_PATH_TAILLE+fic, 'rb') as fichier_graphe :
                mon_depickler = pickle.Unpickler(fichier_graphe)
                extension = mon_depickler.load()
            
                with open("Graphs/"+fic.split("_")[1]+".pickle", 'rb') as fichier_graphe_tot :
                    mon_depickler = pickle.Unpickler(fichier_graphe_tot)
                    graphe = mon_depickler.load()
                    
                    ok_nts  = 0
                    k = 0
                    for i in range(extension.nodes[1]["position"][0] - 6, extension.nodes[1]["position"][0]+1) :
                        #print((extension.nodes[1]["num_ch"], i))
                        #print(graphe.nodes())
                        if (extension.nodes[1]["num_ch"], i) in graphe.nodes() :
                            if liste_sequence_frequence[k][graphe.nodes[(extension.nodes[1]["num_ch"], i)]["nt"]] > 0 :
                                ok_nts += 1
                        k += 1
                    
                    k = 7
                    for i in range(extension.nodes[1]["position"][0]+1, extension.nodes[1]["position"][0]+4) :
                        if (extension.nodes[1]["num_ch"], i) in graphe.nodes() :
                            if liste_sequence_frequence[k][graphe.nodes[(extension.nodes[1]["num_ch"], i)]["nt"]] > 0 :
                                ok_nts += 1
                        k += 1
                       
                    k = 10
                    for i in range(extension.nodes[4]["position"][0], extension.nodes[4]["position"][0]+5) :
                        #print(graphe.nodes())
                        if (extension.nodes[4]["num_ch"], i) in graphe.nodes() :
                            if liste_sequence_frequence[k][graphe.nodes[(extension.nodes[4]["num_ch"], i)]["nt"]] > 0 :
                                ok_nts += 1
                            
                            
                            
                            for voisin in graphe[(extension.nodes[4]["num_ch"], i)] :
                                if graphe.edges[(extension.nodes[4]["num_ch"], i), voisin]["label"] == 'CWW' and \
                                ((graphe.nodes[(extension.nodes[4]["num_ch"], i)]["nt"] == 'C' and graphe.nodes[voisin]["nt"] == 'G') or \
                                  (graphe.nodes[(extension.nodes[4]["num_ch"], i)]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'C') or \
                                  (graphe.nodes[(extension.nodes[4]["num_ch"], i)]["nt"] == 'A' and graphe.nodes[voisin]["nt"] == 'U') or \
                                  (graphe.nodes[(extension.nodes[4]["num_ch"], i)]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'A') or \
                                  (graphe.nodes[(extension.nodes[4]["num_ch"], i)]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'U') or \
                                  (graphe.nodes[(extension.nodes[4]["num_ch"], i)]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'G')):
                                   
                                    #print(graphe.nodes())
                                    #print(voisin)
                                    if liste_sequence_frequence[k+1][graphe.nodes[voisin]["nt"]] > 0 :
                                        ok_nts += 1
                                    
                        k += 2
                    #if ok_nts/20 > pourcentage :
                    if ok_nts/20 == 1.0 :
                        liste_100.append(fic)
                    liste_ok.append(ok_nts/20)
                        
                    compteur += 1
                    print(compteur)
    return liste_ok, liste_100            
 
def draw_seq_sous_graphe_commun(G, liste_pos_motif, rep_sauve, nom_sauve):
    print(G.nodes.data())
    
    nx.set_node_attributes(G, (33,33), "coordonnees")

    G.nodes[liste_pos_motif[0]]["coordonnees"] = (0.0,0.5)
    G.nodes[liste_pos_motif[1]]["coordonnees"] = (2.0,0.5)
    G.nodes[liste_pos_motif[2]]["coordonnees"] = (0.0,0.0)
    G.nodes[liste_pos_motif[3]]["coordonnees"] = (2.0,0.0)
    G.nodes[liste_pos_motif[4]]["coordonnees"] = (3.0,0.5)
     
#                 fichier.write(str(element)+"\n") 
#                 fichier.write(str(G.number_of_nodes())+"\n") 
    print(G.edges.data())

    nodes_list = [u for u,d in G.nodes(data=True)]#and len(G[u]) > 0] 
    print(nodes_list)

    coordonnees = []
    for noeud in nodes_list :
        #voisins = G[noeud]
        print(noeud)
        if noeud not in liste_pos_motif :
            print(G.nodes[noeud])
            if G.nodes[noeud]["type"] != None :
                    chaine = G.nodes[noeud]["chaine"][0]
                    #print(chaine)
                    if chaine == 1 or chaine == 3 :
                            G.nodes[noeud]["coordonnees"] =  (G.nodes[liste_pos_motif[chaine-1]]["coordonnees"][0], G.nodes[liste_pos_motif[chaine-1]]["coordonnees"][1] + (G.nodes[noeud]["position"][0] - G.nodes[liste_pos_motif[chaine-1]]["position"][0])/2) 
                    else :
                            G.nodes[noeud]["coordonnees"] =  (G.nodes[liste_pos_motif[chaine-1]]["coordonnees"][0], G.nodes[liste_pos_motif[chaine-1]]["coordonnees"][1] + (G.nodes[liste_pos_motif[chaine-1]]["position"][0] - G.nodes[noeud]["position"][0])/2) 
        
                    coordonnees.append(G.nodes[noeud]["coordonnees"])
            else :
                    voisin = list(G[noeud])
                    print(len(voisin))
                    print(voisin)
                    if (G.nodes[voisin[0]]["coordonnees"][0] + 0.5, G.nodes[voisin[0]]["coordonnees"][1]) not in coordonnees :
                        G.nodes[noeud]["coordonnees"] = (G.nodes[voisin[0]]["coordonnees"][0] + 0.5, G.nodes[voisin[0]]["coordonnees"][1])
                        coordonnees.append(G.nodes[noeud]["coordonnees"])
                    elif  (G.nodes[voisin[0]]["coordonnees"][0] - 0.5, G.nodes[voisin[0]]["coordonnees"][1]) not in coordonnees :
                        G.nodes[noeud]["coordonnees"] = (G.nodes[voisin[0]]["coordonnees"][0] - 0.5, G.nodes[voisin[0]]["coordonnees"][1])  
                        coordonnees.append(G.nodes[noeud]["coordonnees"])  
                    else :
                        print("probleme")
                        

    print("pos")            
    #pos = nx.get_node_attributes(G, 'coordonnees')
    pos = dict([(u,(d["coordonnees"]))for u,d in G.nodes(data=True) if u in nodes_list])

    print(pos)
    for elt in pos :
        if pos[elt][0] > 30 or pos[elt][1] > 30 :
            print(elt)
     
    fig = plt.figure(figsize =(5,12))
    courbes = []      
    for noeud in nodes_list :
        courbes_x = []
        courbes_y = []  
        coordonnees_x = []
        coordonnees_y = []
        for voisin in G[noeud] :
            if noeud == (46,40) :
                print(voisin)
                print(G.nodes[voisin]["coordonnees"][0])
                print(G.nodes[noeud]["coordonnees"][0])
            #if (noeud, voisin) in G.edges() and (voisin, noeud) in G.edges() :
#                 print(voisin)
#                 print(G.nodes[voisin]["coordonnees"][1])
#                 print(noeud)
#                 print(G.nodes[noeud]["coordonnees"][1])
            if G.nodes[voisin]["coordonnees"][0] == G.nodes[noeud]["coordonnees"][0] :
                    for elt in coordonnees_x :
                        if G.nodes[voisin]["coordonnees"][0] == elt[1] :
                            courbes_x.append((elt[0], voisin))
            if G.nodes[voisin]["coordonnees"][1] == G.nodes[noeud]["coordonnees"][1] :
                    print("ramou")
                    for elt in coordonnees_y :
                        if G.nodes[voisin]["coordonnees"][1] == elt[1] :
                            courbes_y.append((elt[0], voisin))
            coordonnees_x.append((voisin, G.nodes[voisin]["coordonnees"][0]))
            coordonnees_y.append((voisin, G.nodes[voisin]["coordonnees"][1]))
        #if noeud == (33,28) :
        print(noeud)
        print(courbes_x)
        print(courbes_y)
        if len(courbes_x) > 0 :
            for elt in courbes_x :
                if (G.nodes[noeud]["coordonnees"][1] < G.nodes[elt[0]]["coordonnees"][1] and G.nodes[noeud]["coordonnees"][1] < G.nodes[elt[1]]["coordonnees"][1]) or (G.nodes[noeud]["coordonnees"][1] > G.nodes[elt[0]]["coordonnees"][1] and G.nodes[noeud]["coordonnees"][1] > G.nodes[elt[1]]["coordonnees"][1]) :
                    if abs(G.nodes[elt[0]]["coordonnees"][1] - G.nodes[noeud]["coordonnees"][1]) > abs(G.nodes[elt[1]]["coordonnees"][1] - G.nodes[noeud]["coordonnees"][1]) :
                        if (noeud, elt[0]) not in courbes and (elt[0], noeud) not in courbes  :
                            courbes.append((noeud, elt[0]))
                        
                    else :
                        if (noeud, elt[1]) not in courbes and (elt[1], noeud) not in courbes :
                            courbes.append((noeud, elt[1]))
                        
        if len(courbes_y) > 0 :
            for elt in courbes_y :
                if (G.nodes[noeud]["coordonnees"][0] < G.nodes[elt[0]]["coordonnees"][0] and G.nodes[noeud]["coordonnees"][0] < G.nodes[elt[1]]["coordonnees"][0]) or (G.nodes[noeud]["coordonnees"][0] > G.nodes[elt[0]]["coordonnees"][0] and G.nodes[noeud]["coordonnees"][0] > G.nodes[elt[1]]["coordonnees"][0]) :
                    if abs(G.nodes[elt[0]]["coordonnees"][0] - G.nodes[noeud]["coordonnees"][0]) > abs(G.nodes[elt[1]]["coordonnees"][0] - G.nodes[noeud]["coordonnees"][0]) :
                        if (noeud, elt[0]) not in courbes and (elt[0], noeud) not in courbes  :
                            courbes.append((noeud, elt[0]))
                            
                    else :
                        if (noeud, elt[1]) not in courbes and (elt[1], noeud) not in courbes :
                            courbes.append((noeud, elt[1]))
    
    for u,v,key,data in G.edges(data=True, keys=True) :
        if u in nodes_list and v in nodes_list :
            if G.nodes[u]["type"] != None and G.nodes[v]["type"] != None and data["label"] != 'B53' :
                if (u,v) not in courbes and (v,u) not in courbes :
                    courbes.append((u,v))
                            
    print(courbes)
    green_edges = []
    blue_edges = []
    black_edges = []
    grey_edges = []
    #black_edges = [edge for edge in G.edges() if edge not in red_edges]
    for u,v,edata, in G.edges(data=True) :
            if edata["label"] == "B53" :
                green_edges.append((u,v))
            elif edata["label"] == "CWW" :
                blue_edges.append((u,v)) 
            elif edata["label"] == '0' :
                grey_edges.append((u,v))
            else :
                black_edges.append((u,v))

    #edge_labels=dict([((u,v,),d["label"])for u,v,d in G.edges(data=True) if d["label"] != 'B53' and d["label"] != 'CWW' and ((u,v) not in courbes and (v,u) not in courbes)])
    #print(edge_labels)
    node_labels=dict([(u,(d["nt"]))for u,d in G.nodes(data=True) if u in nodes_list and d["type"] != -1])## if d["type"] != None])
    #node_labels=dict([(u, (u,d["type"], d["poids"])) if d["type"] != None else (u, (u)) for u,d in G.nodes(data=True) ])
    #node_labels=dict([(u, (u)) for u,d in G.nodes(data=True) ])
    #print(node_labels)
    orange_nodes = []
    grey_nodes = []
    
    for noeud, data in G.nodes(data=True) :
        print(data)
        if data["en_plus"] :
            grey_nodes.append(noeud)
        else :
            orange_nodes.append(noeud)
    
    nodes_colors = ['orange' if noeud in orange_nodes else 'grey' for noeud in G.nodes()]
    nx.draw_networkx_nodes(G, pos, node_size=150, nodelist=nodes_list, node_color=nodes_colors)
                #nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), 
                #           node_color = values, node_size = 500)
    nx.draw_networkx_labels(G, pos, labels = node_labels, font_size = 8)
    #nx.draw_networkx_edges(G, pos, edgelist=red_edges, edge_color='r', edge_labels = edge_labels)
    #nx.draw_networkx_edges(G, pos, edgelist=blue_edges, edge_color='b', edge_labels = edge_labels)
    #nx.draw_networkx_edges(G, pos, edgelist=green_edges, edge_color='g', edge_labels = edge_labels)
    #nx.draw_networkx_edges(G, pos, edgelist=black_edges, edge_labels = edge_labels)
    
    edge_colors = ['black' if (u,v) in black_edges else 'blue' if (u,v) in blue_edges else 'green' if (u,v) in green_edges else 'grey' for u,v,d in G.edges(data=True) if (u,v) not in courbes and (v,u) not in courbes and d["label"] != '0' and u in nodes_list and v in nodes_list]
#         print(black_edges)
#         print(red_edges)
#         print(edge_colors)
    #nx.draw_networkx_edge_labels(G,pos)
    edges_list = [(u,v,)for u,v,d in G.edges(data=True) if ((u,v) not in courbes and (v,u) not in courbes) and d["label"] != '0' and u in nodes_list and v in nodes_list] 
    edgelabels = dict([((u,v,), (d["label"])) for u,v,d in G.edges(data=True) if (u,v) in edges_list])
    #print(edges_list)
    #nx.draw_networkx_edge_labels(G,pos, edge_labels=edgelabels, font_size=8)
    nx.draw_networkx_edges(G,pos, edgelist = edges_list, edge_color=edge_colors)
    
    ax=plt.gca()
    draw_network(G, courbes, pos,ax)
    print("petit rat")
    plt.axis('off')
    plt.savefig(rep_sauve+nom_sauve) # save as png
    #plt.savefig("graphes_extension/fichier_1FJG_A_48_8.png") # save as png
    print(courbes)
    print(nodes_list)
    for noeud, data, in G.nodes(data=True) :
        if noeud in nodes_list :
            print(noeud)
            print(data)
    #plt.show()
    
    with open(rep_sauve+nom_sauve+".pickle", 'wb') as fichier_ecriture :
        mon_pickler = pickle.Pickler(fichier_ecriture)
        mon_pickler.dump(G)
    plt.clf()
    plt.close()


'''13/09/19 '''
def recherche_signature_sur_seq_entieres(nom, graphe_signature, pourcentage):
    with open("graphs_2.92.pickle", 'rb') as fichier_graphes :
        mon_depickler = pickle.Unpickler(fichier_graphes)
        graphes = mon_depickler.load()
        
        
'''19/09/19
creation d'une matrice cm avec :
si x = la paire ou le nt en position x sur la sequence
et y = la paire ou le nt en position y sur la sequence
cm[x][y] = le nombre d'elements contenant a la fois x et y 
est-ce une matrice de confusion au sens machine learning ? Bof  
''' 
def matrice_confusion(dico_graphes_communs_globaux):  
    liste_sequences = signature_seq_graphe_commun_par_paire(dico_graphes_communs_globaux)
    
    new_liste_sequences = []
    
    print(liste_sequences)
    ## on formatte les sequences communes: exemple :[CG1, CG2, AU3, A5,...] la chaine 1-4 jusquau 8e element et du 9e au 12e, c'est la chaine 3
    for liste in liste_sequences :
        new_liste = []
        compteur = 1
        for elt in liste :
            if len(elt) > 1 :
                #if str(elt[0]) < str(elt[1]) :
                    new_liste.append(str(elt[0])+str(elt[1])+str(compteur))
                #else :
                #    new_liste.append(str(elt[1])+str(elt[0])+str(compteur))
            else :
                new_liste.append(str(elt)+str(compteur))
            compteur +=1
        new_liste_sequences.append(new_liste)
    print(new_liste_sequences)
    print(liste_sequences)
    print(type(liste_sequences))
    
    
    ### Obtention des indices de lignes et de colonnes possibles (pour le moment on fait la difference entre CG et GC)
    ## on considere que chaque paire (resp nt) est possible pour chaque position, on fera le tri apres, dans la matrice
    labels = []
    liste_paires_possibles = ["CG", "GC", "AU", "UA", "GU", "UG"]
    #liste_paires_possibles = ["CG", "AU","GU"]
    for i in range(1,4) :
        for elt in liste_paires_possibles  :
            labels.append(str(elt)+str(i))

    liste_nts_possibles = ["A", "G", "C", "U"]
    for i in range(4,8) :
        for elt in liste_nts_possibles :
            labels.append(str(elt)+str(i))
            
    for i in range(8,13) :        
        for elt in liste_paires_possibles :
            labels.append(str(elt)+str(i))
    
    print(labels)
    print(len(labels))
    
    ## Obtention de la matrice
    cm_new = np.zeros((len(labels),len(labels)))
    for liste in new_liste_sequences :
        for i in range(len(liste)):
            num_ligne = -1
            for k in range(len(labels)) :
                if labels[k] == liste[i] :
                    num_ligne = k
            for j in range(i, len(liste)) :
                num_colonne = -1
                for k in range(len(labels)) :
                    if labels[k] == liste[j] :
                        num_colonne = k
                print(num_ligne)
                print(num_colonne)        
                
                if num_ligne > num_colonne :
                    cm_new[num_colonne][num_ligne] += 1
                    if num_ligne != num_colonne :
                        cm_new[num_ligne][num_colonne] += 1
                else :
                    cm_new[num_ligne][num_colonne] += 1
                    if num_ligne != num_colonne :
                        cm_new[num_colonne][num_ligne] += 1
    
#     for i in range(len(cm_new)) :
#         for j in range(i) :
#             cm_new[i][j] = -1            
    
    print(cm_new)
    print(type(cm_new[0][0]))
    
    ## On enlève les lignes de 0 et les colonnes de 0
    a_enlever_paires= []
    a_enlever_ligne = []  
    compteur = 0
    #print(len(cm_new))
    for ligne in cm_new :
        
        que_des_zeros = True
        for elt in ligne :
            if elt != 0 and elt != -1 :
                que_des_zeros = False
             
        if que_des_zeros :

            a_enlever_paires.append(labels[compteur])
            a_enlever_ligne.append(compteur)
             
        compteur += 1

     
    cm_new = np.delete(cm_new, a_enlever_ligne, axis=0)
    cm_new = np.delete(cm_new, a_enlever_ligne, axis=1)
     

     
    print(len(cm_new))
    for elt in a_enlever_paires :
        labels.remove(elt)
    
    np.savetxt(NEW_EXTENSION_PATH_TAILLE+"Traitement_resultats/test_matrice.txt", cm_new)      
#     with open(NEW_EXTENSION_PATH_TAILLE+"Traitement_resultats/test_matrice.txt", 'w') as fichier_matrice :
#         for ligne in cm_new :
#             fichier_matrice.write(ligne)
    
    
    cm_new_new = np.zeros((len(labels),len(labels)))
    permutations = [9, 15,  3,  7, 18, 11, 21, 23, 26, 16, 25,  1,  8, 20, 14, 12,  6, 17, 24, 19, 13, 22,  4,  2, 10,  5] ## BEA TSP
    #permutations = [11, 23,  3, 21,  5, 26, 18,  8, 14,  6, 25, 12,  1,  7, 15, 17, 10, 16, 20, 22, 19,  4, 24, 13,  2,  9] ## BEA
    #permutations = [10, 16, 17, 20, 22, 15, 19,  1, 12,  4,  7, 24, 25, 18,  6,  2,  9, 13,  8, 14,  3, 5, 11, 21, 23, 26] ## PCA
    compteur = 0
    ##permutation de colonnes
    for elt in permutations :
        for colonne in range(len(cm_new_new)) :
            #cm_new_new[ligne][compteur], cm_new_new[ligne][elt-1] = cm_new[ligne][elt-1], cm_new[ligne][compteur]
            cm_new_new[compteur][colonne] = cm_new[elt-1][colonne]   
        compteur+=1
         
    cm_new_new_int = copy.deepcopy(cm_new_new)
    ##permutation de lignes
    compteur = 0
    for elt in permutations :
        for ligne in range(len(cm_new_new)) :
            cm_new_new[ligne][compteur] = cm_new_new_int[ligne][elt-1]
#             if elt > compteur :
#                 cm_new_new[compteur][colonne], cm_new_new[elt-1][colonne] = cm_new_new_int[elt-1][colonne], cm_new_new_int[compteur][colonne]   
#             else :
#                 cm_new_new[compteur][colonne], cm_new_new[elt-1][colonne] = cm_new_new_int[elt-1][colonne], cm_new_new_int[compteur][colonne]
        compteur+=1
     
    ##permutation des labels
    new_labels = []
    for elt in permutations :
        new_labels.append(labels[elt-1])
     
    for i in range(len(cm_new_new)) :
        for j in range(i) :
            cm_new_new[i][j] = -1 
    
    ## Visualisation du resultat avec une echelle de couleurs en fonction de la valeur
    ## symetrie de la matrice => donc on ne s'interesse qu'a la partie en haut a droite
    maxi = len(cm_new_new)
    #print(cm)
    fig, ax = plt.subplots()
    im = ax.imshow(cm_new_new, interpolation='nearest', cmap=plt.cm.get_cmap('Blues'), vmin=0, vmax=5)
    ax.figure.colorbar(im, ax=ax, pad=0.2)
    # We want to show all ticks...
    ax.set(xticks=np.arange(len(cm_new_new)),
           yticks=np.arange(len(cm_new_new)),
           # ... and label them with the respective list entries
           xticklabels=new_labels, yticklabels=new_labels,
           ylabel='position',
           xlabel='position', 
           ) 
    
    plt.title("Nombre de paires de positions identiques entre éléments du groupe \n (ordre BEA-TSP)", pad=30.0)
    ax.tick_params(axis="x", bottom=True, top=True, labelbottom=True, labeltop=True)
    # Rotate and align bottom ticklabels
    plt.setp([tick.label1 for tick in ax.xaxis.get_major_ticks()], rotation=45,
             ha="right", va="center", rotation_mode="anchor")
    # Rotate and align top ticklabels
    plt.setp([tick.label2 for tick in ax.xaxis.get_major_ticks()], rotation=45,
             ha="left", va="center",rotation_mode="anchor")
    # Rotate the tick labels and set their alignment.
#     plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
#              rotation_mode="anchor")
    
    ax.tick_params(axis="y", left=True, right=True, labelright=True, labelleft=True)
    # Loop over data dimensions and create text annotations.
    fmt = '1.0f'
    thresh = maxi / 2.
    for i in range(len(cm_new_new)):
        for j in range(len(cm_new_new[i])):
            if cm_new_new[i][j] != -1 :
                ax.text(j, i, format(cm_new_new[i][j], fmt),
                        ha="center", va="center",
                        color="white" if cm_new_new[i][j] > thresh else "black")
    fig.tight_layout()
    plt.show()
    
    return cm_new_new, new_labels
    #plt.savefig("matrice_confusion_standard.png")
    

# '''09/10/19'''
# def recherche_candidat_signature(graphe_signature):
    
'''09/10/19
Obtenir la signature de séquence avec fréquence à partir du graphe de signature'''
def obtention_signature_sequence(graphe_signature):
    sequences = []
    
    predecesseurs = []
    for noeud, data in graphe_signature.nodes(data=True) :  
        predecesseur = False
        for voisin in graphe_signature.predecessors(noeud) :
            for edge in graphe_signature[voisin][noeud] :
                if graphe_signature[voisin][noeud][edge]["label"] == 'B53' :
                    predecesseur = True
        if not predecesseur :
            predecesseurs.append(noeud)
    
           
    for noeud in predecesseurs :
        sequence = [graphe_signature.nodes[noeud]["nt"]]
        liaison_B53 = True
        compteur = noeud
        while liaison_B53 :
            liaison_B53 = False
            temp = compteur
            for voisin in graphe_signature[compteur] :
                for arc in graphe_signature[compteur][voisin] :
                    if graphe_signature[compteur][voisin][arc]["label"] == 'B53' :
                        liaison_B53 = True
                        temp = voisin
                        sequence.append(graphe_signature.nodes[voisin]["nt"])
            compteur = temp
        sequences.append(sequence)
        
    print(len(sequences))
    for elt in sequences :
        print(len(elt))
        print(elt)
    print(sequences)
    return sequences

'''09/10/19'''
def recherche_signature_dans_sequence_total(sequence_signature, sequence_total, dico_graphes_communs_globaux, confusion_matrix):
    
#     print("petit rat")
#     print(sequence_signature)
    #for signature in sequence_signature :
        signature = sequence_signature[1]
        ## Generer toutes les possibilites de sequences à chercher
#         print(signature)
        expr_reg = ''
        for pos in signature :
            positions = '['
            for elt in pos :
                positions += '%s'%elt[0]
            positions += ']'
            expr_reg += '%s'%positions
#         print(expr_reg)
        generations = list(exrex.generate(expr_reg))
#         print(len(generations))
#         print(generations)
        
        
        
        index_occurrences = []
        
        compteur = 0
        for sequence in generations :
            occurrence = sequence_total.find(sequence)
            while occurrence != -1 :
                index_occurrences.append((occurrence, compteur))
                occurrence = sequence_total.find(sequence, occurrence+1)
            compteur += 1
        
        cm = confusion_matrix[0] 
        labels = confusion_matrix[1]
        
#         print(len(index_occurrences))
#         print(index_occurrences)
        
        dico_res = {}
        
        for occ in index_occurrences :
#             print(occ)
            sequence_possible = generations[occ[1]]
#             print(sequence_possible)
            # specifique de la matrice chaine 1 et 3
#             sequence_pour_matrice = []
#             compteur_position = 1
#             for i in range(3) :
#                 sequence_pour_matrice.append(sequence_possible[i]+sequence_possible[9-i]+str(compteur_position))
#                 compteur_position += 1
#              
#             for i in range(3,7) :
#                 sequence_pour_matrice.append(sequence_possible[i]+str(compteur_position))
#                 compteur_position += 1
                
                
            # specifique de la matrice chaine 4
            sequence_pour_matrice = []
            compteur_position = 8
            for i in range(5) :
                sequence_pour_matrice.append(sequence_possible[i]+str(compteur_position))
                compteur_position += 1

#             print(sequence_pour_matrice)
             
             # version chaine 1 3
#             score_tot = 0
#             for i_seq in range(len(sequence_pour_matrice)) :
#                 for j_seq in range(i_seq+1, len(sequence_pour_matrice)) :
#                     if sequence_pour_matrice[i_seq] in labels and sequence_pour_matrice[j_seq] in labels :
#                         index1 = labels.index(sequence_pour_matrice[i_seq])
#                         index2 = labels.index(sequence_pour_matrice[j_seq])
#                          
#                         if index1 < index2 :
#                             score = cm[index1][index2]
#                         else :
#                             score = cm[index2][index1]
#                         #print(score)
#                         score_tot += score
            
            # version chaine 4
            score_tot = 0
            for i_seq in range(len(sequence_pour_matrice)) :
                for j_seq in range(i_seq+1, len(sequence_pour_matrice)) :
                    for elt in labels :
                        if sequence_pour_matrice[i_seq][0] in elt and sequence_pour_matrice[i_seq][1] in elt :
                            index1 = labels.index(elt)
                        if sequence_pour_matrice[j_seq][0] in elt and sequence_pour_matrice[j_seq][1] in elt :
                            index2 = labels.index(elt)
                         
                    if index1 < index2 :
                        score = cm[index1][index2]
                    else :
                        score = cm[index2][index1]
                        #print(score)
                    score_tot += score
#             print(score_tot)
            
            dico_res.update({(sequence_possible, occ[0]) : score_tot})
            
        
#         print(generations[1289])
        return dico_res

'''09/10/19'''
def obtention_sequence_issue_de_graphe(graphe_global):
    sequence = ""
    for noeud, data in graphe_global.nodes(data=True) :  
        sequence += data["nt"]
    return sequence

'''05/09/19'''
def recherche_signature(graphe_signature, pourcentage) :
    with open("graphs_2.92.pickle", 'rb') as fichier_graphes :
        mon_depickler = pickle.Unpickler(fichier_graphes)
        graphes = mon_depickler.load()
        scores = {}
        compter = 0
        for fic in os.listdir(EXTENSION_PATH_TAILLE%4) :
            if "pickle" in fic and "couples_possibles" not in fic and "avec_coord" not in fic and len(fic.split("_")) == 6 :
                #print(fic)
                nom = (fic.split("_")[1], fic.split("_")[2])
                with open(EXTENSION_PATH_TAILLE%4+fic, 'rb') as fichier_extension :
                    mon_depickler = pickle.Unpickler(fichier_extension)
                    extension = mon_depickler.load()
                    
                    chaines = recup_chaines(extension)
                    
                    compteur = 1226
                    ok = 0
                    ok_nts = 0
                    for elt in chaines[0] :
                        for i in range(extension.nodes[elt]["poids"]) :
                            if compteur <= 1229 :
#                                 print(i)
#                                 print(extension.nodes[elt]["position"][0]+i)
#                                 print(graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"])
#                                 print(compteur)
    #                             print(i)
    #                             print(extension.nodes[elt]["position"])
    #                             print(extension.nodes[elt]["poids"])
                                vu = False 
                                for nt in graphe_signature.nodes[compteur]["nt"] :
                                    if graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"] == nt[0] :
                                        #ok_nts += nt[1]/5
                                        ok_nts += 1
                                        vu = True
#                                 if not vu :
#                                     ok_nts *= 0.000001
#                                         print("ramou1")
                                        
#                                 if graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"] not in graphe_signature.nodes[compteur]["nt"] :
#                                     ok_nts += 1*
                                    #print("ramous1")
                                
                                
                                compteur_meme_voisin = 0
                                for voisin in graphe_signature[compteur] :
                                    for edge in graphe_signature[compteur][voisin] :
                                        for voisin2 in graphes[nom][extension.nodes[elt]["position"][0]+i] :
                                            if graphes[nom].edges[extension.nodes[elt]["position"][0]+i, voisin2]["label"] == graphe_signature[compteur][voisin][edge]["label"] and graphes[nom].edges[extension.nodes[elt]["position"][0]+i, voisin2]["long_range"] == graphe_signature[compteur][voisin][edge]["long_range"] and abs(compteur-voisin) == abs(extension.nodes[elt]["position"][0]+i - voisin2):
                                                compteur_meme_voisin += 1
                                                
                                if compteur_meme_voisin != len(graphe_signature[compteur]) :
                                    ok += abs(len(graphe_signature[compteur]) -compteur_meme_voisin) 
                                    #print("ramous2")
                                
                                compteur += 1
                
                    elt = chaines[0][len(chaines[0])-1]
                    pos = extension.nodes[elt]["position"][1]+1
                       
                    while compteur <= 1229 :
                        vu = False
                        for nt in graphe_signature.nodes[compteur]["nt"] :
                            if graphes[nom].nodes[pos]["nt"] == nt[0] :
                                        #ok_nts += nt[1]/5
                                        ok_nts += 1
                                        vu = True
#                         if not vu :
#                             ok_nts *= 0.000001
#                                         print("ramou2")
                                        
#                         if graphes[nom].nodes[pos]["nt"] not in graphe_signature.nodes[compteur]["nt"] :
#                                     ok_nts += 1
                                    #print("ramous1")
                                 
                                 
                        compteur_meme_voisin = 0
                        for voisin in graphe_signature[compteur] :
                                    for edge in graphe_signature[compteur][voisin] :
                                        for voisin2 in graphes[nom][pos] :
                                            if graphes[nom].edges[pos, voisin2]["label"] == graphe_signature[compteur][voisin][edge]["label"] and graphes[nom].edges[pos, voisin2]["long_range"] == graphe_signature[compteur][voisin][edge]["long_range"] and abs(compteur-voisin) == abs(pos - voisin2):
                                                compteur_meme_voisin += 1
                                                 
                        if compteur_meme_voisin != len(graphe_signature[compteur]) :
                                    ok += abs(len(graphe_signature[compteur]) -compteur_meme_voisin) 
                                    #print("ramous2")
                                 
                        compteur += 1
                        pos += 1
                    
                        
                    compteur = 1225
#                     if fic == "fichier_4V9F_0_30_4_3.pickle" :
#                                         #print()
#                                         print(compteur)
                    for elt in chaines[2] :
                        for i in np.arange(extension.nodes[elt]["poids"]-1, -1,-1) :
                            if compteur >= 1220 :
                                
                                vu = False
                                for nt in graphe_signature.nodes[compteur]["nt"] :
                                    if graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"] == nt[0] :
                                        #ok_nts += nt[1]/5
                                        ok_nts += 1
                                        vu = True
#                                 if not vu :
#                                     ok_nts *= 0.000001
#                                         print("ramou3")
                                        
#                                 if graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"] not in graphe_signature.nodes[compteur]["nt"] :
#                                     ok_nts += 1
                                    #print("ramous3")
                                
                                
                                compteur_meme_voisin = 0
                                for voisin in graphe_signature[compteur] :
                                    for edge in graphe_signature[compteur][voisin] :
#                                         if fic == "fichier_4V9F_0_30_4_3.pickle" :
#                                             print(graphes[nom][extension.nodes[elt]["position"][0]+1])
                                        for voisin2 in graphes[nom][extension.nodes[elt]["position"][0]+i] :
#                                             if fic == "fichier_4V9F_0_30_4_3.pickle" :
#                                                 #print(nom)
#                                                 print(extension.nodes[elt]["position"][0]+i)
#                                                 print(graphes[nom].edges[extension.nodes[elt]["position"][0]+i, voisin2]["label"])
#                                                 print(elt)
#                                                 print(i)
                                            if graphes[nom].edges[extension.nodes[elt]["position"][0]+i, voisin2]["label"] == graphe_signature[compteur][voisin][edge]["label"] and graphes[nom].edges[extension.nodes[elt]["position"][0]+i, voisin2]["long_range"] == graphe_signature[compteur][voisin][edge]["long_range"] :
                                                compteur_meme_voisin += 1
                                                
                                if compteur_meme_voisin != len(graphe_signature[compteur]) :
                                    ok += abs(len(graphe_signature[compteur]) -compteur_meme_voisin)
#                                     if fic == "fichier_4V9F_0_30_4_3.pickle" :
#                                         #print()
#                                         print(compteur)
#                                         print(compteur_meme_voisin)
#                                         print(len(graphe_signature[compteur]))
                                    #print("ramous4")  
                                
                                compteur -= 1  
                                
                    elt = chaines[2][len(chaines[2])-1]
                    pos = extension.nodes[elt]["position"][0]-1
                        
                    while compteur >= 1220 :
                        vu = False
                        for nt in graphe_signature.nodes[compteur]["nt"] :
                            if graphes[nom].nodes[pos]["nt"] == nt[0] :
                                #ok_nts += 1*(nt[1]/5)
                                ok_nts += 1
                                vu = True
#                         if not vu :
#                             ok_nts *= 0.000001
#                                 print(nt[1]/5)
#                                 print(ok_nts)
#                                 print("ramou4")
                                #ok_nts += 1
#                         if graphes[nom].nodes[pos]["nt"] not in graphe_signature.nodes[compteur]["nt"] :
#                                     ok_nts += 1
                                    #print("ramous1")
                                 
                                 
                        compteur_meme_voisin = 0
                        for voisin in graphe_signature[compteur] :
                                    for edge in graphe_signature[compteur][voisin] :
                                        for voisin2 in graphes[nom][pos] :
                                            if graphes[nom].edges[pos, voisin2]["label"] == graphe_signature[compteur][voisin][edge]["label"] and graphes[nom].edges[pos, voisin2]["long_range"] == graphe_signature[compteur][voisin][edge]["long_range"] and abs(compteur-voisin) == abs(pos - voisin2):
                                                compteur_meme_voisin += 1
                                                 
                        if compteur_meme_voisin != len(graphe_signature[compteur]) :
                                    ok += abs(len(graphe_signature[compteur]) -compteur_meme_voisin)
                                     
                                    #print("ramous2")
                                 
                        compteur -= 1
                        pos -= 1
                    
                    
                    compteur = 813            
                    for elt in chaines[3] :
                        for i in range(extension.nodes[elt]["poids"]) :
                            if compteur <= 816 :
                                #print(compteur)
#                                 print(i)
#                                 print(extension.nodes[elt]["position"][0]+i)
#                                 print(graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"])
#                                 print(compteur)
    #                             print(i)
    #                             print(extension.nodes[elt]["position"])
    #                             print(extension.nodes[elt]["poids"])
                                vu = False  
                                for nt in graphe_signature.nodes[compteur]["nt"] :
                                    if graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"] == nt[0] :
                                        #ok_nts += nt[1]/5
                                        ok_nts += 1
                                        vu = True
#                                 if not vu :
#                                     ok_nts *= 0.000001
#                                         print("ramou5")
                                       
#                                 if graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"] not in graphe_signature.nodes[compteur]["nt"] :
#                                     ok_nts += 1
                                    #print("ramous1")
                                
                                
                                voisin_can = False
                                for voisin2 in graphes[nom][extension.nodes[elt]["position"][0]+i] :
                                    
#                                     print(graphes[nom].edges[extension.nodes[elt]["position"][0]+i, voisin2]["label"])
#                                     print(graphes[nom].edges[extension.nodes[elt]["position"][0]+i, voisin2]["long_range"])
                                    if graphes[nom].edges[extension.nodes[elt]["position"][0]+i, voisin2]["label"] == 'CWW' and graphes[nom].edges[extension.nodes[elt]["position"][0]+i, voisin2]["long_range"] == False \
                                     and ((graphes[nom].nodes[voisin2]["nt"] == 'G' and graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"] == 'C') or \
                                     (graphes[nom].nodes[voisin2]["nt"] == 'C' and graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"] == 'G') or \
                                     (graphes[nom].nodes[voisin2]["nt"] == 'A' and graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"] == 'U') or \
                                     (graphes[nom].nodes[voisin2]["nt"] == 'U' and graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"] == 'A') or \
                                     (graphes[nom].nodes[voisin2]["nt"] == 'G' and graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"] == 'U') or \
                                     (graphes[nom].nodes[voisin2]["nt"] == 'U' and graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"] == 'G')) :
                                        voisin_can = True
                                                
                                if not voisin_can :
                                    ok += 1
                                    #print("ramous2")
                                
                                compteur += 1
                              
                    if fic == "fichier_1FJG_A_48_8_3.pickle" :
                        print("gros rat")
                        print(ok)
                        print(ok_nts)
                    #score = ok_nts * (ok/8)
                    score = (ok_nts+ok)/28
#                     print(ok_nts)
#                     print(score)
                    scores.update({fic : score})
                    #score = 
                    #print(scores)
                    if score >= pourcentage:
                        print("ok")
                        print(fic)
                        compter += 1
                        print(ok_nts)
                        print(ok)
        
#         maxi = max(scores.values())
#         mini = min(scores.values()) 
        print(scores)
#         print("maxi")
#         print(maxi)
#         print("mini")
#         print(mini)  
#         print(maxi-mini) 
#         for cle in scores.keys() :       
#             elt = scores[cle]
#             elt = (elt-mini)/(maxi-mini)
#             #print(elt)
#             if cle == "fichier_4V9F_0_30_4_3.pickle" :
#                 print("tout petit rat")
#                 print(elt)
#             if elt >= pourcentage:
#                 print(elt)
#                 print("ok")
#                 print(cle)
#                 compter += 1
#                 print(ok_nts)
#                 print(ok)
        #return scores             
        print(compter)
        return compter
    

'''12/09/19'''
def recherche_signature_version_new_data(graphe_signature, pourcentage, liste_tout) :
#     with open("graphs_2.92.pickle", 'rb') as fichier_graphes :
#         mon_depickler = pickle.Unpickler(fichier_graphes)
#         graphes = mon_depickler.load()
        
        nx.convert_node_labels_to_integers(graphe_signature, label_attribute="premier")
        
        num_motif = [-1,-1,-1,-1]
        for noeud, data in graphe_signature.nodes(data=True) :
            for i in range(4) :
                if data["type"] == (i+1)*10 :
                    num_motif[i] = noeud 
        
        
        compter = 0
        for elt in liste_tout :
                with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s_2.pickle"%(elt[0], elt[1]), 'rb') as fichier_extension :
                    mon_depickler = pickle.Unpickler(fichier_extension)
                    extension = mon_depickler.load()
                    
                    with open("Graphs/%s.pickle"%elt[0], 'rb') as fichier_graphe :
                        mon_depickler = pickle.Unpickler(fichier_graphe)
                        graphe = mon_depickler.load()
                    
                        chaines = recup_chaines(extension)
                        
                        compteur = 1226
                        ok = 0
                        ok_nts = 0
                        for elt in chaines[0] :
                            for i in range(extension.nodes[elt]["poids"]) :
                                if compteur <= 1229 :
    #                                 print(i)
    #                                 print(extension.nodes[elt]["position"][0]+i)
    #                                 print(graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"])
    #                                 print(compteur)
        #                             print(i)
        #                             print(extension.nodes[elt]["position"])
        #                             print(extension.nodes[elt]["poids"])
                                    if graphe.nodes[(extension.nodes[elt]["num_ch"],extension.nodes[elt]["position"][0]+i)]["nt"] not in graphe_signature.nodes[compteur]["nt"] :
                                        ok_nts += 1
                                        #print("ramous1")
                                    
                                    
                                    compteur_meme_voisin = 0
                                    for voisin in graphe_signature[compteur] :
                                        for edge in graphe_signature[compteur][voisin] :
                                            for voisin2 in graphe[(extension.nodes[elt]["num_ch"],extension.nodes[elt]["position"][0]+i)] :
                                                if graphe.edges[(extension.nodes[elt]["num_ch"], extension.nodes[elt]["position"][0]+i), voisin2]["label"] == graphe_signature[compteur][voisin][edge]["label"] and abs(compteur-voisin) == abs(extension.nodes[elt]["position"][0]+i - voisin2[1]):
                                                    compteur_meme_voisin += 1
                                                    
                                    if compteur_meme_voisin != len(graphe_signature[compteur]) :
                                        ok += abs(len(graphe_signature[compteur]) -compteur_meme_voisin) 
                                        #print("ramous2")
                                    
                                    compteur += 1
                    
                        elt = chaines[0][len(chaines[0])-1]
                        pos = extension.nodes[elt]["position"][1]+1
                        num_ch = extension.nodes[elt]["num_ch"]
                           
                        while compteur <= 1229 :
                            if graphe.nodes[(num_ch, pos)]["nt"] not in graphe_signature.nodes[compteur]["nt"] :
                                        ok_nts += 1
                                        #print("ramous1")
                                     
                                     
                            compteur_meme_voisin = 0
                            for voisin in graphe_signature[compteur] :
                                        for edge in graphe_signature[compteur][voisin] :
                                            for voisin2 in graphe[(num_ch,pos)] :
                                                if graphe.edges[(num_ch, pos), voisin2]["label"] == graphe_signature[compteur][voisin][edge]["label"] and abs(compteur-voisin) == abs(pos - voisin2[1]):
                                                    compteur_meme_voisin += 1
                                                     
                            if compteur_meme_voisin != len(graphe_signature[compteur]) :
                                        ok += abs(len(graphe_signature[compteur]) -compteur_meme_voisin) 
                                        #print("ramous2")
                                     
                            compteur += 1
                            pos += 1
                        
                            
                        compteur = 1225
    #                     if fic == "fichier_4V9F_0_30_4_3.pickle" :
    #                                         #print()
    #                                         print(compteur)
                        for elt in chaines[2] :
                            for i in np.arange(extension.nodes[elt]["poids"]-1, -1,-1) :
                                if compteur >= 1220 :
                                    if graphe.nodes[(extension.nodes[elt]["num_ch"],extension.nodes[elt]["position"][0]+i)]["nt"] not in graphe_signature.nodes[compteur]["nt"] :
                                        ok_nts += 1
                                        #print("ramous3")
                                    
                                    
                                    compteur_meme_voisin = 0
                                    for voisin in graphe_signature[compteur] :
                                        for edge in graphe_signature[compteur][voisin] :
    #                                         if fic == "fichier_4V9F_0_30_4_3.pickle" :
    #                                             print(graphes[nom][extension.nodes[elt]["position"][0]+1])
                                            for voisin2 in graphe[(extension.nodes[elt]["num_ch"], extension.nodes[elt]["position"][0]+i)] :
    #                                             if fic == "fichier_4V9F_0_30_4_3.pickle" :
    #                                                 #print(nom)
    #                                                 print(extension.nodes[elt]["position"][0]+i)
    #                                                 print(graphes[nom].edges[extension.nodes[elt]["position"][0]+i, voisin2]["label"])
    #                                                 print(elt)
    #                                                 print(i)
                                                if graphe.edges[(extension.nodes[elt]["num_ch"], extension.nodes[elt]["position"][0]+i), voisin2]["label"] == graphe_signature[compteur][voisin][edge]["label"] and abs(compteur-voisin) == abs(extension.nodes[elt]["position"][0]+i - voisin2[1]):
                                                    compteur_meme_voisin += 1
                                                    
                                    if compteur_meme_voisin != len(graphe_signature[compteur]) :
                                        ok += abs(len(graphe_signature[compteur]) -compteur_meme_voisin)
    #                                     if fic == "fichier_4V9F_0_30_4_3.pickle" :
    #                                         #print()
    #                                         print(compteur)
    #                                         print(compteur_meme_voisin)
    #                                         print(len(graphe_signature[compteur]))
                                        #print("ramous4")  
                                    
                                    compteur -= 1  
                                    
                        elt = chaines[2][len(chaines[2])-1]
                        pos = extension.nodes[elt]["position"][0]-1
                        num_ch = extension.nodes[elt]["num_ch"]
                            
                        while compteur >= 1220 :
                            if graphe.nodes[(num_ch,pos)]["nt"] not in graphe_signature.nodes[compteur]["nt"] :
                                        ok_nts += 1
                                        #print("ramous1")
                                     
                                     
                            compteur_meme_voisin = 0
                            for voisin in graphe_signature[compteur] :
                                        for edge in graphe_signature[compteur][voisin] :
                                            for voisin2 in graphe[(num_ch,pos)] :
                                                if graphe.edges[(num_ch,pos), voisin2]["label"] == graphe_signature[compteur][voisin][edge]["label"]and abs(compteur-voisin) == abs(pos - voisin2[1]):
                                                    compteur_meme_voisin += 1
                                                     
                            if compteur_meme_voisin != len(graphe_signature[compteur]) :
                                        ok += abs(len(graphe_signature[compteur]) -compteur_meme_voisin)
                                         
                                        #print("ramous2")
                                     
                            compteur -= 1
                            pos -= 1
                        
                        
                        compteur = 813            
                        for elt in chaines[3] :
                            for i in range(extension.nodes[elt]["poids"]) :
                                if compteur <= 816 :
                                    #print(compteur)
    #                                 print(i)
    #                                 print(extension.nodes[elt]["position"][0]+i)
    #                                 print(graphes[nom].nodes[extension.nodes[elt]["position"][0]+i]["nt"])
    #                                 print(compteur)
        #                             print(i)
        #                             print(extension.nodes[elt]["position"])
        #                             print(extension.nodes[elt]["poids"])
                                    if graphe.nodes[(extension.nodes[elt]["num_ch"], extension.nodes[elt]["position"][0]+i)]["nt"] not in graphe_signature.nodes[compteur]["nt"] :
                                        ok_nts += 1
                                        #print("ramous1")
                                    
                                    
                                    voisin_can = False
                                    for voisin2 in graphe[(extension.nodes[elt]["num_ch"], extension.nodes[elt]["position"][0]+i)] :
                                        
    #                                     print(graphes[nom].edges[extension.nodes[elt]["position"][0]+i, voisin2]["label"])
    #                                     print(graphes[nom].edges[extension.nodes[elt]["position"][0]+i, voisin2]["long_range"])
                                        if graphe.edges[(extension.nodes[elt]["num_ch"], extension.nodes[elt]["position"][0]+i), voisin2]["label"] == 'CWW' \
                                         and ((graphe.nodes[voisin2]["nt"] == 'G' and graphe.nodes[(extension.nodes[elt]["num_ch"], extension.nodes[elt]["position"][0]+i)]["nt"] == 'C') or \
                                         (graphe.nodes[voisin2]["nt"] == 'C' and graphe.nodes[(extension.nodes[elt]["num_ch"], extension.nodes[elt]["position"][0]+i)]["nt"] == 'G') or \
                                         (graphe.nodes[voisin2]["nt"] == 'A' and graphe.nodes[(extension.nodes[elt]["num_ch"], extension.nodes[elt]["position"][0]+i)]["nt"] == 'U') or \
                                         (graphe.nodes[voisin2]["nt"] == 'U' and graphe.nodes[(extension.nodes[elt]["num_ch"], extension.nodes[elt]["position"][0]+i)]["nt"] == 'A') or \
                                         (graphe.nodes[voisin2]["nt"] == 'G' and graphe.nodes[(extension.nodes[elt]["num_ch"], extension.nodes[elt]["position"][0]+i)]["nt"] == 'U') or \
                                         (graphe.nodes[voisin2]["nt"] == 'U' and graphe.nodes[(extension.nodes[elt]["num_ch"], extension.nodes[elt]["position"][0]+i)]["nt"] == 'G')) :
                                            voisin_can = True
                                                    
                                    if not voisin_can :
                                        ok += 1
                                        #print("ramous2")
                                    
                                    compteur += 1
                                  
                        
                        if 21 - (ok + ok_nts) >= 21*pourcentage:
                            print("ok")
                            #print(fic)
                            compter += 1
                            print(ok_nts)
                            print(ok)
                    
                        
        print(compter)
        return compter  
    
''' 19/09/19 
obtention paires d elements et boucle'''
def signature_seq_graphe_commun_par_paire(dico_graphes_communs_globaux):
    
    graphe_paires = nx.MultiDiGraph()
    cle_deb = list(dico_graphes_communs_globaux.keys())[0]
    liste_triee_noeuds = list(dico_graphes_communs_globaux[cle_deb][0].nodes())
    liste_triee_noeuds.sort()
    
    liste_triee_separee = list(liste_triee_noeuds[11:])
    liste_triee_separee.extend(list(liste_triee_noeuds[:11]))
    
    
    print(liste_triee_separee)    
        
    
    noeuds_deja_vus = []
    liste_sequences = [[], [], [], [], []]
    for noeud in liste_triee_separee :
        if noeud not in noeuds_deja_vus :
            data_ref = dico_graphes_communs_globaux[cle_deb][0].nodes[noeud]
            if 1 in data_ref["chaine"] \
            or 3 in data_ref["chaine"] \
            or 4 in data_ref["chaine"] \
            or data_ref["type"] == 12 \
            or data_ref["type"] == 15 :
                paires = False
                for voisin in dico_graphes_communs_globaux[cle_deb][0][noeud] :
                    for edge in dico_graphes_communs_globaux[cle_deb][0][noeud][voisin] :
                        if dico_graphes_communs_globaux[cle_deb][0][noeud][voisin][edge]["label"] == 'CWW' :
                            noeuds_deja_vus.append(voisin)
                            paires = True
                            liste_sequences[0].append((data_ref["nt"], dico_graphes_communs_globaux[cle_deb][0].nodes[voisin]["nt"]))
                if not paires :
                    liste_sequences[0].append(data_ref["nt"])
                
                compteur = 0
                for cle in dico_graphes_communs_globaux.keys() :
                    if cle != list(dico_graphes_communs_globaux.keys())[0] :
                        for noeud_autre, data in dico_graphes_communs_globaux[cle][0].nodes(data=True) :
                            if data["chaine"] == data_ref["chaine"] and data["position"] == data_ref["position"] and data["type"] == data_ref["type"] and data["en_plus"] == data_ref["en_plus"] :
                                    
                                    if paires :
                                            for voisin_autre in dico_graphes_communs_globaux[cle][0][noeud_autre] :
                                                for edge_autre in dico_graphes_communs_globaux[cle][0][noeud_autre][voisin_autre] :
                                                    if dico_graphes_communs_globaux[cle][0][noeud_autre][voisin_autre][edge_autre]["label"] == 'CWW' :
                                                        liste_sequences[compteur].append((data["nt"], dico_graphes_communs_globaux[cle][0].nodes[voisin_autre]["nt"]))
                                    else :
                                        liste_sequences[compteur].append(data["nt"])
                    compteur += 1 
                    
    print(liste_sequences)
    return liste_sequences
                          


def nb_voisins_non_b53(noeud, graph):
    liste_voisins = []
    compteur = 0
    for voisin in graph[noeud] :
        for edge in graph[noeud][voisin] :
            if graph[noeud][voisin][edge]["label"] != "B53" :
                compteur += 1
                liste_voisins.append((voisin, edge))
    return compteur, liste_voisins

def memes_voisins(noeud1, graph1, noeud2, graph2):
    
    nb_voisins1, liste_voisins1 = nb_voisins_non_b53(noeud1, graph1)
    nb_voisins2, liste_voisins2 = nb_voisins_non_b53(noeud2, graph2)
    
    if nb_voisins1 == nb_voisins2 :
        trouve1 = []
        trouve2 = []
        for voisin1, edge1 in liste_voisins1 :
            for voisin2, edge2 in liste_voisins2 :
                if graph1[noeud1][voisin1][edge1]["label"] == graph2[noeud2][voisin2][edge2]["label"] and voisin1 not in trouve1 and voisin2 not in trouve2 :
                    trouve1.append(voisin1)
                    trouve2.append(voisin2)
                    
        if len(trouve1) == nb_voisins1 and len(trouve2) == nb_voisins2  :
            return True
        else :
            return False
    else :
        return False

def nombre_voisins_b53_meme_type(noeud1, graph1, typ):
    compteur = 0
    liaison_b53 = True
    print("succ")
    while liaison_b53 : 
        liaison_b53 = False
        new_noeud = -1
        for voisin in graph1[noeud1] :
            for edge in graph1[noeud1][voisin] :
                if graph1[noeud1][voisin][edge]["label"] == 'B53' and graph1.nodes[voisin]["type"] == typ :
                    compteur += 1
                    liaison_b53 = True    
                    new_noeud = voisin
                    print(voisin)
        noeud1 = new_noeud           
    
    return compteur

def nombre_pred_b53_meme_type(noeud1, graph1, typ):
    compteur = 0
    liaison_b53 = True
    print("pred")
    while liaison_b53 : 
        liaison_b53 = False
        new_noeud = -1
        for voisin in graph1.predecessors(noeud1) :
            for edge in graph1[voisin][noeud1] :
                if graph1[voisin][noeud1][edge]["label"] == 'B53' and graph1.nodes[voisin]["type"] == typ :
                    compteur += 1
                    liaison_b53 = True    
                    new_noeud = voisin
                    print(voisin)
        noeud1 = new_noeud           
    
    return compteur
              
''' 05/09/19
modif le 19/09/19 pour ajouter le nombre de paires de chaque type'''
def signature_sequence_graphe_commun(dico_graphes_communs_globaux):
    graphe_signature = nx.MultiDiGraph()
    
    print(dico_graphes_communs_globaux[list(dico_graphes_communs_globaux.keys())[0]][0].edges.data())
    print(dico_graphes_communs_globaux[list(dico_graphes_communs_globaux.keys())[0]][0].nodes())
    for noeud, data in dico_graphes_communs_globaux[list(dico_graphes_communs_globaux.keys())[0]][0].nodes(data=True) :
        if 1 in data["chaine"] or 3 in data["chaine"] or 4 in data["chaine"] or 2 in data["chaine"] or data["type"] == 15 :
            if data["type"] == 1 or data["type"] == 0 or data["type"] == None :
                nb_succ = nombre_voisins_b53_meme_type(noeud, dico_graphes_communs_globaux[list(dico_graphes_communs_globaux.keys())[0]][0], data["type"])
                #if nb_succ+1 > 
#             data_new = copy.deepcopy(data)
#             data_new["nt"] = [[data_new["nt"],1]]
            graphe_signature.add_node(noeud, **data)

    cle_deb = list(dico_graphes_communs_globaux.keys())[0]           
    for u,v,data in dico_graphes_communs_globaux[cle_deb][0].edges(data=True) :
        if u in graphe_signature.nodes() and v in graphe_signature.nodes() :
            data_new = copy.deepcopy(data)
            if data["label"] != 'B53' :
                print(type(dico_graphes_communs_globaux[cle_deb][0].nodes[u]["nt"]))
                print(dico_graphes_communs_globaux[cle_deb][0].nodes[u]["nt"])
                paires_nts = {(dico_graphes_communs_globaux[cle_deb][0].nodes[u]["nt"][0][0],dico_graphes_communs_globaux[cle_deb][0].nodes[v]["nt"][0][0]) : 1 }
                for cle in dico_graphes_communs_globaux.keys() :
                    if cle != list(dico_graphes_communs_globaux.keys())[0] :
                        for u_autre, v_autre, data_autre in dico_graphes_communs_globaux[cle][0].edges(data=True) :
                            if dico_graphes_communs_globaux[cle_deb][0].nodes[u]["chaine"] == dico_graphes_communs_globaux[cle][0].nodes[u_autre]["chaine"] \
                            and dico_graphes_communs_globaux[cle_deb][0].nodes[u]["position"] == dico_graphes_communs_globaux[cle][0].nodes[u_autre]["position"] \
                            and dico_graphes_communs_globaux[cle_deb][0].nodes[u]["type"] == dico_graphes_communs_globaux[cle][0].nodes[u_autre]["type"] \
                            and dico_graphes_communs_globaux[cle_deb][0].nodes[u]["en_plus"] == dico_graphes_communs_globaux[cle][0].nodes[u_autre]["en_plus"] \
                            and dico_graphes_communs_globaux[cle_deb][0].nodes[v]["chaine"] == dico_graphes_communs_globaux[cle][0].nodes[v_autre]["chaine"] \
                            and dico_graphes_communs_globaux[cle_deb][0].nodes[v]["position"] == dico_graphes_communs_globaux[cle][0].nodes[v_autre]["position"] \
                            and dico_graphes_communs_globaux[cle_deb][0].nodes[v]["type"] == dico_graphes_communs_globaux[cle][0].nodes[v_autre]["type"] \
                            and dico_graphes_communs_globaux[cle_deb][0].nodes[v]["en_plus"] == dico_graphes_communs_globaux[cle][0].nodes[v_autre]["en_plus"] :
    
                            
                                if (dico_graphes_communs_globaux[cle][0].nodes[u_autre]["nt"][0][0],dico_graphes_communs_globaux[cle][0].nodes[v_autre]["nt"][0][0]) in paires_nts.keys() :
                                    paires_nts[(dico_graphes_communs_globaux[cle][0].nodes[u_autre]["nt"][0][0], dico_graphes_communs_globaux[cle][0].nodes[v_autre]["nt"][0][0])] += 1
                                elif (dico_graphes_communs_globaux[cle][0].nodes[v_autre]["nt"][0][0], dico_graphes_communs_globaux[cle][0].nodes[u_autre]["nt"][0][0]) in paires_nts.keys() :
                                    paires_nts[(dico_graphes_communs_globaux[cle][0].nodes[v_autre]["nt"][0][0], dico_graphes_communs_globaux[cle][0].nodes[u_autre]["nt"][0][0])] += 1
                                else :
                                    paires_nts.update({(dico_graphes_communs_globaux[cle][0].nodes[u_autre]["nt"][0][0],dico_graphes_communs_globaux[cle][0].nodes[v_autre]["nt"][0][0]) : 1})
    
                
                data_new.update({"paires_nts" : paires_nts})               
            graphe_signature.add_edge(u,v,**data_new)
    
    
            
    
    for cle in dico_graphes_communs_globaux.keys() :
        #if cle != list(dico_graphes_communs_globaux.keys())[0] :
            for noeud, data in dico_graphes_communs_globaux[cle][0].nodes(data=True) :
                entre = False
                print("roupoulou")
                for noeud_ref, data_ref in graphe_signature.nodes(data=True) :
                    
                    
                    #if data["position"][0] > 50 :
#                         print("gros rat")
#                         print(str(noeud) +  " " + str(data))
#                     if data["type"] == 11 and data_ref["type"] == 11 :
#                         print("roupoulou")
#                         print(nb_voisins_non_b53(noeud, dico_graphes_communs_globaux[cle][0]))
#                         print(nb_voisins_non_b53(noeud_ref, graphe_signature))
#                         print(memes_voisins(noeud, dico_graphes_communs_globaux[cle][0], noeud_ref, graphe_signature))
#                         print(data["chaine"])
#                         print(data_ref["chaine"])
                    if data["chaine"] == data_ref["chaine"] and data["type"] == data_ref["type"] and memes_voisins(noeud, dico_graphes_communs_globaux[cle][0], noeud_ref, graphe_signature) :
                        ok = False
                        if data["type"] == 1 or data["type"] == 0 or data["type"] == None :
                            print(nombre_voisins_b53_meme_type(noeud, dico_graphes_communs_globaux[cle][0], data["type"]))
                            print(nombre_voisins_b53_meme_type(noeud_ref, graphe_signature, data["type"]))
                            print(nombre_pred_b53_meme_type(noeud, dico_graphes_communs_globaux[cle][0], data["type"]))
                            print(nombre_pred_b53_meme_type(noeud_ref, graphe_signature, data["type"]))
                            if nombre_voisins_b53_meme_type(noeud, dico_graphes_communs_globaux[cle][0], data["type"]) >= nombre_voisins_b53_meme_type(noeud_ref, graphe_signature, data["type"]) \
                            and nombre_pred_b53_meme_type(noeud, dico_graphes_communs_globaux[cle][0], data["type"]) >= nombre_pred_b53_meme_type(noeud_ref, graphe_signature, data["type"])  :
                                ok = True
                        else :
                            ok = True
                        
                        if ok :
#                             print("ripili")
                            vu = False
                            entre = True
                            for elt in graphe_signature.nodes[noeud_ref]["nt"] :
                                if elt[0] == data["nt"] :
                                    elt[1] += 1
                                    vu = True
                                    
                            if not vu :
                                graphe_signature.nodes[noeud_ref]["nt"].append([data["nt"],1])
                               
                                        #graphe_signature.nodes[noeud_ref]["nt"].append((data["nt"],1))
                            if noeud_ref == 1192 :
                                    
                                    print("gros rat")
                                    print(data["position"])
                
                if not entre :
                    print(noeud, data)
                    print("rapala")
    print(graphe_signature.nodes.data())
    for noeud, data in graphe_signature.nodes(data=True) :
        print(str(noeud) + " " + str(data))    
        
    for u,v,data in graphe_signature.edges(data=True) :
        print(str(u)+ " "+str(v) + " "+str(data))   
        
    return graphe_signature              


'''25/09/19 
recherche du graphe commun du groupe 11 dans les extensions'''
def recherche_graphe_commun_version_rapide(fic_ext):

    
    def recherche_voisin_b53(extension, noeud):
        for voisin in extension[noeud] :
            for edge in extension[noeud][voisin] :
                if extension[noeud][voisin][edge]["label"] == 'B53' :
                    return voisin
        return -1
    
    def recherche_pred_b53(extension, noeud):
        for voisin in extension.predecessors(noeud) :
            for edge in extension[voisin][noeud] :
                if extension[voisin][noeud][edge]["label"] == 'B53' :
                    return voisin
        return -1
    
    ## liaison cww de la chaine 1 liee a la chaine 3, et liaison cov avec le noeud 1
    def recherche_cww_11(extension, noeud):
        global nb_liaisons_cov_correctes, nb_liaisons_non_cov_correctes
        ok_voisin_voisin = False
        compteur = 1
        noeud_suivant = recherche_voisin_b53(extension, noeud)
        if noeud_suivant != -1 :
            if extension.nodes[noeud_suivant]["type"] == 1 :
                for voisin_voisin in extension[noeud_suivant] :
                    for edge_voisin in extension[noeud_suivant][voisin_voisin] :
                        if extension[noeud_suivant][voisin_voisin][edge_voisin]["label"] == 'CWW' and extension.nodes[voisin_voisin]["type"] == 1 and 3 in extension.nodes[voisin_voisin]["chaine"]:
                            ok_voisin_voisin = True
                            nb_liaisons_non_cov_correctes += 1
                            if noeud == 1 :
                                nb_liaisons_cov_correctes += 1
                            if extension.nodes[noeud_suivant]["poids"] >  1 :
                                compteur = 2
        return noeud_suivant, ok_voisin_voisin, compteur
    
    ## liaisons cww de la chaine 1 toutes seules et liaisons cov associees + liaisons cww de la chaine 4 et liaisons cov associees
    def recherche_cww_12_41(extension, noeud, compteur_b53, nb_liaisons):
        global nb_liaisons_cov_correctes, nb_liaisons_non_cov_correctes
        ok_voisin_voisin = False
        noeud_suivant = recherche_voisin_b53(extension, noeud)
        if noeud_suivant != -1 :
                    #print(noeud_suivant)
                    #print(extension.nodes[voisin])
            if extension.nodes[noeud_suivant]["type"] == 1 :
               #print(extension[voisin])
                for voisin_voisin in extension[noeud_suivant] :
                    #print(extension.nodes[voisin_voisin])
                    for edge_voisin in extension[noeud_suivant][voisin_voisin] :
                        if extension[noeud_suivant][voisin_voisin][edge_voisin]["label"] == '0' and extension.nodes[voisin_voisin]["type"] == -1 :
                            #print("tout petit rat")
                            ok_voisin_voisin = True
                            nb_liaisons_non_cov_correctes += min(nb_liaisons,extension.nodes[noeud_suivant]["poids"])
                            if compteur_b53 == 1 :
                                nb_liaisons_cov_correctes += min(nb_liaisons,extension.nodes[noeud_suivant]["poids"])
                            else :
                                nb_liaisons_cov_correctes += min(nb_liaisons-1,extension.nodes[noeud_suivant]["poids"]- 1)
        return noeud_suivant, ok_voisin_voisin
    
    ## aretes artificielles de la chaine 3 et liaisons cov avec le noeud 3 et associees au noeud de type 0
    def rechercher_0_3(extension, noeud):
        global nb_liaisons_cov_correctes, nb_liaisons_non_cov_correctes
        ok_type = False
        noeud_suivant = recherche_pred_b53(extension, noeud)
        #print(noeud_suivant)
        if noeud_suivant != -1 :
            #print(extension.nodes[noeud_suivant])
            if extension.nodes[noeud_suivant]["type"] == 0 :
                #print("gros rat")
                ok_type = True
                nb_liaisons_non_cov_correctes += min(2, extension.nodes[noeud_suivant]["poids"])
                if noeud == 3 :
                    nb_liaisons_cov_correctes += min(2, extension.nodes[noeud_suivant]["poids"]) 
                else :
                    nb_liaisons_cov_correctes += min(1, extension.nodes[noeud_suivant]["poids"]-1)
        return noeud_suivant, ok_type
    
    ## liaison cww avec le noeud 4
    def rechercher_liaison_4(extension, noeud):
        global nb_liaisons_non_cov_correctes
        for voisin in extension[noeud] :
            for edge in extension[noeud][voisin] :
                if extension[noeud][voisin][edge]["label"] == 'CWW' and extension.nodes[voisin]["type"] == None :
                    nb_liaisons_non_cov_correctes += 1
    
    ## liaison cww de la chaine 3 (pas compte, juste verifie) et liaison cov entre le noeud de type 0 et celui-la
    def rechercher_cww_3(extension, noeud, compteur):    
        global nb_liaisons_cov_correctes, nb_liaisons_non_cov_correctes    
        ok_voisin_voisin = False
        noeud_suivant = recherche_pred_b53(extension, noeud)
        #print(noeud_suivant)
        if noeud_suivant != -1 :
            if extension.nodes[noeud_suivant]["type"] == 1 :
                for voisin_voisin in extension[noeud_suivant] :
                    for edge_voisin in extension[noeud_suivant][voisin_voisin] :
                        if extension[noeud_suivant][voisin_voisin][edge_voisin]["label"] == 'CWW' and extension.nodes[voisin_voisin]["type"] == 1 and 1 in extension.nodes[voisin_voisin]["chaine"]:
                            ok_voisin_voisin = True
                            if compteur == 1 :
                                nb_liaisons_cov_correctes += 1
        return noeud_suivant, ok_voisin_voisin
    
    ## arete artificielle de la chaine 2
    def rechercher_0_2(extension, noeud):
        global nb_liaisons_cov_correctes, nb_liaisons_non_cov_correctes
        ok_type = False
        noeud_suivant = recherche_pred_b53(extension, noeud)
        #print(noeud_suivant)
        if noeud_suivant != -1 :
            #print(extension.nodes[noeud_suivant])
            if extension.nodes[noeud_suivant]["type"] == 0 :
                #print("gros rat")
                ok_type = True
                nb_liaisons_non_cov_correctes += 1
        return noeud_suivant, ok_type
    
    with open(NEW_EXTENSION_PATH_TAILLE+fic_ext, 'rb') as fichier_extension :
            
            mon_depickler = pickle.Unpickler(fichier_extension)
            extension = mon_depickler.load()
             
            
            ##chaine 1 
            noeud_suivant, ok_voisin_voisin, compteur = recherche_cww_11(extension, 1)
            while noeud_suivant != -1 and not ok_voisin_voisin and noeud_suivant not in [1,2,3,4]:
                noeud_suivant, ok_voisin_voisin, compteur  = recherche_cww_11(extension, noeud_suivant)
            
            
            if noeud_suivant != -1 :
                noeud_suivant, ok_voisin_voisin = recherche_cww_12_41(extension, noeud_suivant, compteur, 2)
                while noeud_suivant != -1 and not ok_voisin_voisin and noeud_suivant not in [1,2,3,4] : ## bout de la chaine
                    compteur += 1
                    noeud_suivant, ok_voisin_voisin = recherche_cww_12_41(extension, noeud_suivant, compteur,2)
                
            if nb_liaisons_cov_correctes == 0 and nb_liaisons_non_cov_correctes == 0 :
                noeud_suivant, ok_voisin_voisin = recherche_cww_12_41(extension, 1, 0, 2)
                #print(noeud_suivant)
                #print(ok_voisin_voisin)
                while noeud_suivant != -1 and not ok_voisin_voisin and noeud_suivant not in [1,2,3,4]: ## bout de la chaine
                    compteur += 1
                    #print("petit rat")
                    noeud_suivant, ok_voisin_voisin = recherche_cww_12_41(extension, noeud_suivant, 0,2)
                    
            ## chaine 3
            noeud_suivant, ok_type = rechercher_0_3(extension, 3)
            while noeud_suivant != -1 and not ok_type and noeud_suivant not in [1,2,3,4] :
                noeud_suivant, ok_type = rechercher_0_3(extension, noeud_suivant)
                 
            if noeud_suivant != -1 :
                noeud_suivant, ok_voisin_voisin = rechercher_cww_3(extension, noeud_suivant, 1)
            while noeud_suivant != -1 and not ok_voisin_voisin and noeud_suivant not in [1,2,3,4] :
                noeud_suivant, ok_voisin_voisin = rechercher_cww_3(extension, noeud_suivant, 2)
            
            ## chaine 4 
            rechercher_liaison_4(extension, 4)
            noeud_suivant, ok_voisin_voisin = recherche_cww_12_41(extension, 4, 1, 3)
            #print(noeud_suivant)
            #print(ok_voisin_voisin)
            while noeud_suivant != -1 and not ok_voisin_voisin and noeud_suivant not in [1,2,3,4]: ## bout de la chaine
                    #print(noeud_suivant)
                    #print(ok_voisin_voisin)
                    noeud_suivant, ok_voisin_voisin = recherche_cww_12_41(extension, noeud_suivant, 2,3)
                    
            
            ## chaine 2
            noeud_suivant, ok_type = rechercher_0_2(extension, 2)
            while noeud_suivant != -1 and not ok_type and noeud_suivant not in [1,2,3,4] :
                noeud_suivant, ok_type = rechercher_0_2(extension, noeud_suivant)

    return nb_liaisons_non_cov_correctes, nb_liaisons_cov_correctes                       
                        

def pourcentage_ressemblance_structure():
    ''' Pourcentage de ressemblance au niveau structure avec le groupe 11 '''
    
    
    compteur = 0
    distrib = []#             if nb_non_cov + nb_cov >= 19 :
#                 print(fic)
#                 print(nb_non_cov)
#                 print(nb_cov)
    exact = 0
    liste_100_pourcent = []
    for fic in os.listdir(NEW_EXTENSION_PATH_TAILLE) :
        if "pickle" in fic and "fichier" in fic :
            #print(compteur)
            #print(fic)
            nb_liaisons_non_cov_correctes = 0
            nb_liaisons_cov_correctes = 0
            nb_non_cov, nb_cov = recherche_graphe_commun_version_rapide(fic)
            distrib.append((nb_non_cov + nb_cov)/19)
            if nb_non_cov == 10 and nb_cov == 9 :
                exact += 1
                print(fic)
                liste_100_pourcent.append(fic)
#             if nb_non_cov + nb_cov >= 19 :
#                 print(fic)
#                 print(nb_non_cov)
#                 print(nb_cov)
                  
            #print((nb_non_cov, nb_cov))
            compteur += 1
             
    print(liste_100_pourcent)
    #print(distrib)
    liste_nb = []
    for pourcentage in np.arange(1.0, -0.1, -0.1) :
        nb = 0
        for elt in distrib :
            if elt > pourcentage :
                nb += 1
        liste_nb.append(nb)
      
    ax = plt.gca()
    #print([round(i, 1) for i in np.arange(1.0, -0.1, -0.1)])
    ax.set_xticks(np.arange(0,11,1))
    ax.set_xticklabels([round(i, 1) for i in np.arange(1.0, -0.1, -0.1)]) 
#     ax.set_yticklabels(np.arange(5, 94, 5))
#     ax.set_yticks(np.arange(5, 94, 5))             
    plt.plot(liste_nb)
    plt.title("Similarité de structure avec le graphe commun du groupe vert-jaune \n sur les extensions")
    ax.set_xlabel("Pourcentage de ressemblance")
    ax.set_ylabel("Nombre d'éléments")
    #plt.show()
    #plt.savefig(NEW_EXTENSION_PATH_TAILLE+"Traitement_resultats/similarite_structure_groupe_vert_jaune.png")
    print(exact)
    return liste_100_pourcent
   

def pourcentage_ressemblance_sequence():
    liste_nb_ok = []
    for pourcentage in np.arange(1.0, 0.6, -0.1) :
        pourcentage = 1.0   
        for i in range(4, 5) :
                    dico_graphes = recherche_graphe_commun([i], CLUSTERING_PEREZ_VERSION_NON_CAN_2[11], 0.7, "toutes_aretes_coeff_all1", "taille_4", "clustering_perez_sim_groupe_%d"%11)
                    graphe_signature = signature_sequence_graphe_commun(dico_graphes)   
                    nb_ok = recherche_signature_version_new_data(graphe_signature, pourcentage)
                    liste_nb_ok.append(nb_ok)
    print(nb_ok)
#     sns.distplot(nb_ok)
#     plt.show()
    print(liste_nb_ok) 
    ax = plt.gca()
    print([round(i, 1) for i in np.arange(1.0, -0.1, -0.1)])
    ax.set_xticks(np.arange(0,11,1))
    ax.set_xticklabels([round(i, 1) for i in np.arange(1.0, -0.1, -0.1)]) 
    ax.set_yticklabels(np.arange(5, 94, 5))
    ax.set_yticks(np.arange(5, 94, 5))             
    plt.plot(liste_nb_ok)
    plt.title("Similarité de séquence avec la signature du groupe 11 \n version frequence")
    ax.set_xlabel("Pourcentage de ressemblance")
    ax.set_ylabel("Nombre d'éléments")
    plt.show()
    
''' 09/10/19 '''
def recherche_score_occurrences(liste_100_pourcent):
    dico_graphes = recherche_graphe_commun([4], CLUSTERING_PEREZ_VERSION_NON_CAN_2[11], 0.7, "toutes_aretes_coeff_all1", "taille_4", "clustering_perez_sim_groupe_%d"%11)
    graphe_signature = signature_sequence_graphe_commun(dico_graphes) 
    sequences_signatures = obtention_signature_sequence(graphe_signature)
    cm = matrice_confusion(dico_graphes)
     
    #matrice_confusion(dico_graphes)
    liste_vraies_aminor = []
    distrib_score_vraies_aminor = []
    distrib_scores = []
#     for elt in liste_100_pourcent :
#         with open("Graphs/%s.pickle"%elt.split("_")[1], 'rb') as fichier_tout :
            #print(elt)
#             mon_depickler_graphes = pickle.Unpickler(fichier_tout)
#             graphe = mon_depickler_graphes.load()
    with open("graphs_2.92.pickle", 'rb') as fichier_tout :  
        mon_depickler_graphes = pickle.Unpickler(fichier_tout)
        graphes = mon_depickler_graphes.load()       
        print(CLUSTERING_PEREZ_VERSION_NON_CAN_2[11])
        for elt_groupe in CLUSTERING_PEREZ_VERSION_NON_CAN_2[11] :
            sequence_tot = obtention_sequence_issue_de_graphe(graphes[(elt_groupe.split("_")[0], elt_groupe.split("_")[1])])
            print(len(sequence_tot))
            print(sequence_tot)
             
            dico_res_sequence_1 = recherche_signature_dans_sequence_total(sequences_signatures, sequence_tot, dico_graphes, cm)
            print(dico_res_sequence_1)
             
            distrib_scores.extend(dico_res_sequence_1.values())
             
            with open(EXTENSION_PATH_TAILLE%4+"fichier_"+elt_groupe+".pickle", 'rb') as fichier_extension :
                mon_depickler = pickle.Unpickler(fichier_extension)
                extension = mon_depickler.load()
                 
                vraie_a_minor = 0
                print(int(extension.nodes[2]["position"][0]))
                for cle in dico_res_sequence_1 :
                    if  int(extension.nodes[2]["position"][0]) - int(cle[1]) == 0 :# and int(extension.nodes[1]["position"][0]) - int(cle[1]) > 0  :
                        vraie_a_minor = cle
            print(vraie_a_minor)
            if vraie_a_minor != 0 :
                liste_vraies_aminor.append(vraie_a_minor)
            
                distrib_score_vraies_aminor.append(dico_res_sequence_1[vraie_a_minor])
     
            
     
            #print(min(distrib_score_vraies_aminor))
            #print(max(distrib_score_vraies_aminor))
            print(distrib_scores)
             
            plt.figure(figsize=(6,6))
            #print([dico_res_sequence_1[vraie_a_minor]])
            ax = sns.distplot(list(dico_res_sequence_1.values()), kde=False, bins=60) 
            if vraie_a_minor != 0 :
                ax2 = sns.distplot([dico_res_sequence_1[vraie_a_minor]], kde=False, bins=1, color="red") 
             
            #ax.set_xticks(np.arange(0,61,5))
            #ax.set_xticklabels([round(i, 1) for i in np.arange(1.0, -0.1, -0.1)]) 
            #ax.set_yticks(np.arange(0,11,1))
            ax.set_xlabel("Score")
            ax.set_ylabel("Nombre d'occurrences")
            plt.title("Distribution des scores des occurrences de la signature \n du groupe vert jaune sur la séquence de %s \n Nombre total d'occurrences : %d "%(elt_groupe.split("_")[0], len(dico_res_sequence_1)) )
            plt.grid()
            plt.savefig("distrib_score_occ_2e_partie_groupe_11_additif_sequence_%s.png"%elt_groupe)#elt[8:len(elt)-9])
            #plt.show()   
            #plt.clf()
     
    plt.figure(figsize=(7,8))        
    ax = sns.distplot(distrib_scores, kde=False, bins=60)  
     
   # ax.set_xticks(np.arange(0,61,5))
    #ax.set_xticklabels([round(i, 1) for i in np.arange(1.0, -0.1, -0.1)]) 
    #ax.set_yticks(np.arange(0,126,5))
    ax.set_xlabel("Score")
    ax.set_ylabel("Nombre d'occurrences")
    plt.title("Distribution des scores des occurrences de la signature \n du groupe vert jaune sur les séquences du groupe étendu (28 éléments) \n Nombre total d'occurrences : %d "%(len(distrib_scores)) )
    plt.grid()
    plt.savefig("distrib_score_occ_2e_partie_groupe_11_additif.png")
    #plt.show() 
    #plt.clf()
     
    plt.figure(figsize=(7,8))
    ax = sns.distplot(distrib_score_vraies_aminor, kde=False, bins=60)  
     
    #ax.set_xticks(np.arange(0,61,5))
    #ax.set_xticklabels([round(i, 1) for i in np.arange(1.0, -0.1, -0.1)]) 
    #ax.set_yticks(np.arange(0,11,1))
    ax.set_xlabel("Score")
    ax.set_ylabel("Nombre d'occurrences")
    plt.title("Distribution des scores des vraies occurrences d'A-minor \n du groupe vert jaune sur les séquences du groupe étendu (28 éléments) \n Nombre total d'occurrences : %d "%(len(distrib_score_vraies_aminor)) )
    plt.grid()
    plt.savefig("distrib_score_vraies_occ_2e_partie_groupe_11_additif.png")
    #plt.show()  


def diff_liste_cliques(liste_cliques1, liste_cliques2):
    liste_cliques2_temp = list(liste_cliques2)
    for clique1 in liste_cliques1 : 
        ok = 0
        cl2 = -1
        compteur_cliques2 = 0
        for clique2 in liste_cliques2_temp :
            compter = 0
            for elt1 in clique1 :
                if elt1 in clique2 :
                    compter += 1 
            if compter == len(clique1) and compter == len(clique2) :
                ok += 1
                cl2 = compteur_cliques2
            compteur_cliques2 += 1
        if ok != 1 :
            return True
        else : 
            del(liste_cliques2_temp[cl2])
    return False

if __name__ == '__main__':
    
    
    with open("Nouvelles_donnees/Graphes_communs_avril_2020/graphe_commun_cluster_55_4_avec_coord.pickle", 'rb') as fichier_sortie :
        mon_depickler = pickle.Unpickler(fichier_sortie)
        graphe_commun = mon_depickler.load()
        
        for noeud, data in graphe_commun.nodes(data=True) :
            print(noeud, data)
                    #print(noeud, data["poids"], data["liste_nts"])
                   
        for u,v,data in graphe_commun.edges(data=True) :
            print(u,v,data)
        
    #exit()
    
    types = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16S_arnm", "arnt_16S"]
    with open("/media/coline/Maxtor/clustering_perez_tot_new_data_100320_sim_par_branche_0.65.pickle", 'rb') as fichier_sortie :
        mon_depickler = pickle.Unpickler(fichier_sortie)
        clusters = mon_depickler.load()
        
        print(clusters[55])
        #exit()
        
        with open("/media/coline/Maxtor/dico_new_100320_sim_par_branche_0.65.pickle", 'rb') as fichier_graphe :    
            mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
            dico_complet = mon_depickler_graphe.load() 
            
            taille_clusters = []
            taille_clusters_vrai = []
            for cluster in clusters :
#                 print("\n")
#                 print(cluster)
                groupes_egaux = []
                for i in range(len(cluster)) :
                    for j in range(i+1, len(cluster)) :
                        if (cluster[i], cluster[j]) in dico_complet.keys() :
                            sim = dico_complet[(cluster[i], cluster[j])]["sim"]
                        else :
                            sim = dico_complet[(cluster[j], cluster[i])]["sim"]
                        
                        if sim  == 1.0 :
                                ok = False
                                for groupe in groupes_egaux :
                                    if (cluster[i] in groupe or cluster[j] in groupe) and not (cluster[i] in groupe and cluster[j] in groupe) :
                                        if cluster[i] not in groupe :
                                            groupe_temp = list(groupe)
                                            groupe_temp.append(cluster[i])
                                            groupe = groupe_temp
                                        elif cluster[j] not in groupe :
                                            #print(groupe)
                                            groupe_temp = list(groupe)
                                            groupe_temp.append(cluster[j])
                                            groupe = groupe_temp
                                            #print(groupe)
                                        ok = True
                                if not ok :
                                    groupes_egaux.append([cluster[i], cluster[j]])
                    ok_1 = False   
                    for groupe in groupes_egaux :
                        if cluster[i] in groupe :
                            ok_1 = True
                    if not ok_1 :
                        groupes_egaux.append([cluster[i]])
#                 print(groupes_egaux)
                
                for groupe in groupes_egaux :
                    for i in range(len(groupe)) :
                        for j in range(i+1, len(groupe)) :
                            if (groupe[i], groupe[j]) in dico_complet.keys() :
                                sim = dico_complet[(groupe[i], groupe[j])]["sim"]
                            else :
                                sim = dico_complet[(groupe[j], groupe[i])]["sim"]
                            
                            if sim != 1 :
                                print("ahhh")
                taille_clusters.append(len(groupes_egaux))
                taille_clusters_vrai.append(len(cluster))
                if len(groupes_egaux) == 12 :
                    print(cluster)
                    print(groupes_egaux)
            
            sns.distplot(taille_clusters, kde=False)
            print(taille_clusters)
            print(taille_clusters_vrai)
            print(len(clusters))
            print(len(taille_clusters))
            plt.show()
            #exit()          
        #print(clusters[59])
        #print(clusters[51])
        #print(clusters)
        compteur = 1
        for cluster in clusters :
            
            
            #print(cluster52 : 
            if len(cluster) > 2 and compteur == 53:
            #if ('5wfs', 19) in cluster :
                #cluster.remove(('1u9s', 1))
                print(compteur)
                print(cluster)
                print(len(cluster))
                #exit()
                graphe_commun_global, dico_graphes = recherche_graphe_commun([4], cluster, 0.6, "ramou", "Graphes_communs_avril_2020", compteur)
                draw_new_data(graphe_commun_global, "Graphes_communs_avril_2020", compteur, 4, True, False)
                #exit()
#                 for cle in dico_graphes.keys() :
#                     draw_new_data(dico_graphes[cle][0], "Graphes_communs_avril_2020", compteur, 4, True, True, motif = dico_graphes[cle][1], num_graphe=cle)
                print(graphe_commun_global.number_of_nodes())
                print(cluster)
                for noeud, data in graphe_commun_global.nodes(data=True) :
                    print(noeud, data)
                    print(noeud, data["poids"], data["liste_nts"])
                exit()
                for u,v,data in graphe_commun_global.edges(data=True) :
                    print(u,v,data)
#                 print(dico_graphes.nodes.data())
#                 print(dico_graphes.edges.data())
                print(cluster)
                
                #graphe_signature = signature_sequence_graphe_commun(dico_graphes) 
                 
#                 for noeud, data in graphe_signature.nodes(data=True) :
#                     print(noeud, data)
#                 for u,v,data in graphe_signature.edges(data=True) :
#                     print(u,v,data)  
                
#             if len(cluster) == 6 :
#                 #recherche_graphe_commun([4], cluster, 0.6, "ramou", "Graphes_communs", compteur)
#               
#                 tps1 = time.time()
#                 ## verif permutations
#                 ancien_digraphe_commun = nx.MultiDiGraph()
#                 elt_ancien = -1
#                 ancienne_liste_cliques = []
#                 liste_permut = list(itertools.permutations(cluster))
#                 compter = 0
#                 for elt in liste_permut : 
#                     print(elt)
#                     digraphe_commun, liste_cliques = commun_cluster_clique_new_data(elt, "/media/coline/Maxtor/dico_new_100320_sim_par_branche_0.65.pickle", NEW_EXTENSION_PATH_TAILLE)
#                    # liste_cliques = commun_cluster_clique_new_data(elt, "/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_612.pickle", NEW_EXTENSION_PATH_TAILLE)
#                     
#                     #fichier.write(str(elt)+ " " + str(digraphe_commun.nodes.data()) + "\n\n")
#                     #fichier.write(str(digraphe_commun.nodes.data())+"\n")
#                     if ancien_digraphe_commun.number_of_nodes() > 0 and not idem_graphes_communs_version_clique_new_data(ancien_digraphe_commun, digraphe_commun) :
#         #                 fichier.write( " " + str(elt)+ " " + str(digraphe_commun.nodes.data()) + "\n"+str(ancien_digraphe_commun.nodes.data())+" %d non ok\n"%i)
#         #                 fichier.write(str(elt)+ " " + str(digraphe_commun.edges.data()) + "\n"+str(ancien_digraphe_commun.edges.data())+" %d non ok\n"%i)
#                         print("gros rat")
#                         print(digraphe_commun.number_of_nodes())
#                         print(ancien_digraphe_commun.number_of_nodes()) 
#                         print(elt)
#                         print(elt_ancien)
#                         exit()
#                       
#                           
#                     ancien_digraphe_commun = digraphe_commun.copy()
#                     #ancienne_liste_cliques = list(liste_cliques)
#                     elt_ancien = elt
#                     compter += 1
#                     print(compter)
#                     print(len(liste_permut))
#                 print("petit rat")
#                 print(ancien_digraphe_commun.nodes.data())
#                 
#                 tps2 = time.time()
#                 print(tps2 - tps1)
                #exit()
#                 
#             
            compteur += 1
        exit()
    
#     distrib = []
#     exact = 0
#     compteur = 0
#     liste_100_pourcent = []
#     for fic in os.listdir(NEW_EXTENSION_PATH_TAILLE) :
#         if "pickle" in fic and "fichier" in fic :
#             #print(compteur)
#             #print(fic)
#             nb_liaisons_non_cov_correctes = 0
#             nb_liaisons_cov_correctes = 0
#             nb_non_cov, nb_cov = recherche_graphe_commun_version_rapide(fic)
#             distrib.append((nb_non_cov + nb_cov)/19)
#             if nb_non_cov == 10 and nb_cov == 9 :
#                 exact += 1
#                 print(fic)
#                 liste_100_pourcent.append(fic)
# #             if nb_non_cov + nb_cov >= 19 :
# #                 print(fic)
# #                 print(nb_non_cov)
# #                 print(nb_cov)
#                      
#             #print((nb_non_cov, nb_cov))
#             compteur += 1
#                 
#     print(liste_100_pourcent)
#     
#     recherche_score_occurrences(liste_100_pourcent)
#     
#     
#     liste_100_vert_jaune = ['fichier_5j8a_30.pickle', 'fichier_4ybb_30.pickle', 'fichier_4v84_33.pickle', 'fichier_5gaf_4.pickle', 'fichier_4wqu_11.pickle', 'fichier_5v8i_39.pickle', 'fichier_1yjn_3.pickle', 'fichier_4v5p_32.pickle', 'fichier_4v4q_36.pickle', 'fichier_4v6d_46.pickle', 'fichier_4v8o_6.pickle', 'fichier_4v8g_8.pickle', 'fichier_6h58_29.pickle', 'fichier_5j8b_27.pickle', 'fichier_4v52_34.pickle', 'fichier_5h5u_21.pickle', 'fichier_5uym_34.pickle', 'fichier_4l71_19.pickle', 'fichier_4v7j_31.pickle', 'fichier_4u20_25.pickle', 'fichier_4v9s_51.pickle', 'fichier_4wqy_12.pickle', 'fichier_4v7y_24.pickle', 'fichier_4lnt_59.pickle', 'fichier_3jcd_11.pickle', 'fichier_4v8b_32.pickle', 'fichier_4v5n_26.pickle', 'fichier_4wq1_9.pickle', 'fichier_1yjw_4.pickle', 'fichier_4v5m_23.pickle', 'fichier_5hau_23.pickle', 'fichier_1vy7_11.pickle', 'fichier_4v6e_43.pickle', 'fichier_5fdv_12.pickle', 'fichier_6of1_3.pickle', 'fichier_3ccs_5.pickle', 'fichier_6bok_47.pickle', 'fichier_4v5e_3.pickle', 'fichier_5ib7_55.pickle', 'fichier_4wt8_52.pickle', 'fichier_5j30_47.pickle', 'fichier_4tua_19.pickle', 'fichier_5ibb_51.pickle', 'fichier_5j88_27.pickle', 'fichier_4v7z_24.pickle', 'fichier_6c4i_16.pickle', 'fichier_4zsn_12.pickle', 'fichier_5mdv_25.pickle', 'fichier_4v4h_35.pickle', 'fichier_4v5l_6.pickle', 'fichier_1yi2_5.pickle', 'fichier_4v5l_31.pickle', 'fichier_4v7y_23.pickle', 'fichier_6hrm_14.pickle', 'fichier_4u26_42.pickle', 'fichier_1vvj_13.pickle', 'fichier_4v8a_32.pickle', 'fichier_4w29_23.pickle', 'fichier_4v5f_51.pickle', 'fichier_5czp_18.pickle', 'fichier_4v5q_33.pickle', 'fichier_4p6f_47.pickle', 'fichier_4w4g_20.pickle', 'fichier_1yit_5.pickle', 'fichier_4wsd_45.pickle', 'fichier_4v6g_39.pickle', 'fichier_6nd6_4.pickle', 'fichier_4v8g_53.pickle', 'fichier_4wzo_8.pickle', 'fichier_4wsd_49.pickle', 'fichier_5el6_51.pickle', 'fichier_4v5s_64.pickle', 'fichier_3ccv_4.pickle', 'fichier_6i7v_22.pickle', 'fichier_1vy4_7.pickle', 'fichier_4v55_29.pickle', 'fichier_5j5b_32.pickle', 'fichier_4v9a_22.pickle', 'fichier_4v8d_40.pickle', 'fichier_3cc7_5.pickle', 'fichier_4w2g_30.pickle', 'fichier_4v7w_22.pickle', 'fichier_4v6c_45.pickle', 'fichier_4v5d_14.pickle', 'fichier_5wit_30.pickle', 'fichier_1yj9_5.pickle', 'fichier_4wt8_1.pickle', 'fichier_4v5d_27.pickle', 'fichier_4v8i_9.pickle', 'fichier_5el4_12.pickle', 'fichier_5hcp_10.pickle', 'fichier_3ccq_5.pickle', 'fichier_4v6p_3.pickle', 'fichier_4wpo_55.pickle', 'fichier_4v9j_6.pickle', 'fichier_4u25_19.pickle', 'fichier_6c5l_11.pickle', 'fichier_4v8u_55.pickle', 'fichier_4v5f_56.pickle', 'fichier_4wqu_34.pickle', 'fichier_4p70_14.pickle', 'fichier_6gsj_12.pickle', 'fichier_4v54_32.pickle', 'fichier_4v5c_11.pickle', 'fichier_5wis_7.pickle', 'fichier_4wqy_32.pickle', 'fichier_6bok_21.pickle', 'fichier_6q95_2.pickle', 'fichier_4p70_21.pickle', 'fichier_3jcn_15.pickle', 'fichier_4v9m_8.pickle', 'fichier_5ib8_46.pickle', 'fichier_4w29_6.pickle', 'fichier_4xej_29.pickle', 'fichier_5vpp_38.pickle', 'fichier_3j7z_13.pickle', 'fichier_6bz7_44.pickle', 'fichier_5dfe_20.pickle', 'fichier_4v6n_11.pickle', 'fichier_4v6e_25.pickle', 'fichier_5e7k_10.pickle', 'fichier_4v7l_9.pickle', 'fichier_4v8n_13.pickle', 'fichier_4ypb_19.pickle', 'fichier_4v51_26.pickle', 'fichier_4ypb_48.pickle', 'fichier_4v9h_32.pickle', 'fichier_4lsk_23.pickle', 'fichier_1vy4_52.pickle', 'fichier_5doy_28.pickle', 'fichier_3bbx_5.pickle', 'fichier_4wsm_42.pickle', 'fichier_4v9q_20.pickle', 'fichier_4wzd_13.pickle', 'fichier_4v90_28.pickle', 'fichier_4v5r_10.pickle', 'fichier_4v9q_12.pickle', 'fichier_4wsd_38.pickle', 'fichier_5mgp_14.pickle', 'fichier_4wu1_13.pickle', 'fichier_4y4o_12.pickle', 'fichier_4w2f_7.pickle', 'fichier_4v8a_49.pickle', 'fichier_5ib7_48.pickle', 'fichier_4y4o_27.pickle', 'fichier_5gae_14.pickle', 'fichier_4v5s_38.pickle', 'fichier_6czr_9.pickle', 'fichier_4lel_19.pickle', 'fichier_4p6f_19.pickle', 'fichier_3j9z_18.pickle', 'fichier_4zer_31.pickle', 'fichier_1vvj_19.pickle', 'fichier_4v67_49.pickle', 'fichier_5el4_43.pickle', 'fichier_4wqu_54.pickle', 'fichier_5ib7_16.pickle', 'fichier_4tud_47.pickle', 'fichier_5gah_11.pickle', 'fichier_3ccl_4.pickle', 'fichier_5iqr_29.pickle', 'fichier_4v5y_33.pickle', 'fichier_4v6q_8.pickle', 'fichier_4u1v_23.pickle', 'fichier_4wqr_48.pickle', 'fichier_4v9n_48.pickle', 'fichier_4w2g_11.pickle', 'fichier_5e81_49.pickle', 'fichier_6bz8_23.pickle', 'fichier_5vp2_22.pickle', 'fichier_4wq1_50.pickle', 'fichier_6gsj_50.pickle', 'fichier_5hau_33.pickle', 'fichier_4v8e_39.pickle', 'fichier_4v9c_25.pickle', 'fichier_6bok_48.pickle', 'fichier_4v95_55.pickle', 'fichier_4v8i_52.pickle', 'fichier_4v7m_55.pickle', 'fichier_4v84_8.pickle', 'fichier_3cma_6.pickle', 'fichier_6buw_19.pickle', 'fichier_4v9b_37.pickle', 'fichier_4wt1_36.pickle', 'fichier_4v6o_5.pickle', 'fichier_6c5l_34.pickle', 'fichier_5j4b_25.pickle', 'fichier_6buw_46.pickle', 'fichier_5hau_8.pickle', 'fichier_4wzd_43.pickle', 'fichier_4v8c_31.pickle', 'fichier_4tub_17.pickle', 'fichier_5hau_37.pickle', 'fichier_6cfl_10.pickle', 'fichier_4v9n_12.pickle', 'fichier_4v97_7.pickle', 'fichier_4wf1_23.pickle', 'fichier_4v5g_11.pickle', 'fichier_6q95_3.pickle', 'fichier_4v9j_35.pickle', 'fichier_4lfz_42.pickle', 'fichier_1yhq_5.pickle', 'fichier_4v4h_36.pickle', 'fichier_5hcq_24.pickle', 'fichier_5el7_13.pickle', 'fichier_4l47_14.pickle', 'fichier_6o97_30.pickle', 'fichier_4z8c_10.pickle', 'fichier_5kcr_12.pickle', 'fichier_3ccm_5.pickle', 'fichier_4v7b_27.pickle', 'fichier_5j3c_48.pickle', 'fichier_4v7j_20.pickle', 'fichier_4v5f_34.pickle', 'fichier_4v8h_22.pickle', 'fichier_4v83_35.pickle', 'fichier_4u26_23.pickle', 'fichier_4v7k_26.pickle', 'fichier_4z8c_39.pickle', 'fichier_5fdv_40.pickle', 'fichier_4v64_35.pickle', 'fichier_5mdw_25.pickle', 'fichier_5j30_18.pickle', 'fichier_4wpo_64.pickle', 'fichier_6boh_40.pickle', 'fichier_4lnt_14.pickle', 'fichier_5w4k_29.pickle', 'fichier_4v55_30.pickle', 'fichier_5hcq_34.pickle', 'fichier_6q9a_23.pickle', 'fichier_4v8e_19.pickle', 'fichier_4wqy_53.pickle', 'fichier_4v6g_25.pickle', 'fichier_4v8h_31.pickle', 'fichier_6cfj_3.pickle', 'fichier_4v8q_8.pickle', 'fichier_5vpo_15.pickle', 'fichier_6gsj_44.pickle', 'fichier_4ypb_12.pickle', 'fichier_4v72_9.pickle', 'fichier_4y4p_33.pickle', 'fichier_3jbu_15.pickle', 'fichier_4v9j_7.pickle', 'fichier_4lt8_14.pickle', 'fichier_5e81_41.pickle', 'fichier_4v8j_11.pickle', 'fichier_4v6r_6.pickle', 'fichier_3cce_3.pickle', 'fichier_4wro_7.pickle', 'fichier_4v9j_32.pickle', 'fichier_4w2i_6.pickle', 'fichier_5el5_46.pickle', 'fichier_6enu_17.pickle', 'fichier_4v5c_53.pickle', 'fichier_4v9r_10.pickle', 'fichier_5hcq_39.pickle', 'fichier_4v7t_43.pickle', 'fichier_4v7u_45.pickle', 'fichier_4wf1_44.pickle', 'fichier_3j9y_13.pickle', 'fichier_5j8a_31.pickle', 'fichier_4v9s_10.pickle', 'fichier_5fdv_27.pickle', 'fichier_4v7k_43.pickle', 'fichier_4v97_52.pickle', 'fichier_6fkr_16.pickle', 'fichier_4v57_34.pickle', 'fichier_5nrg_7.pickle', 'fichier_4w2h_31.pickle', 'fichier_4w2f_48.pickle', 'fichier_1vy6_11.pickle', 'fichier_5kpx_13.pickle', 'fichier_5ib8_12.pickle', 'fichier_4v5k_2.pickle', 'fichier_4v7x_20.pickle', 'fichier_5fdu_40.pickle', 'fichier_4wsm_8.pickle', 'fichier_6h4n_13.pickle', 'fichier_4v7p_12.pickle', 'fichier_4v7k_14.pickle', 'fichier_4v67_1.pickle', 'fichier_5j7l_30.pickle', 'fichier_4wzo_43.pickle', 'fichier_4v8u_11.pickle', 'fichier_4v6f_39.pickle', 'fichier_4u24_44.pickle', 'fichier_4lt8_22.pickle', 'fichier_4v9n_45.pickle', 'fichier_6cae_3.pickle', 'fichier_4wro_48.pickle', 'fichier_5jc9_30.pickle', 'fichier_5el5_13.pickle', 'fichier_5aka_15.pickle', 'fichier_6fkr_26.pickle', 'fichier_4w2i_49.pickle', 'fichier_4zer_37.pickle', 'fichier_4v7w_21.pickle', 'fichier_5hcr_10.pickle', 'fichier_4wro_37.pickle', 'fichier_4v9q_15.pickle', 'fichier_5el6_44.pickle', 'fichier_4z3s_5.pickle', 'fichier_4v7l_26.pickle', 'fichier_4wsm_33.pickle', 'fichier_4v9r_30.pickle', 'fichier_4tue_19.pickle', 'fichier_4v68_7.pickle', 'fichier_5j88_28.pickle', 'fichier_4u1u_43.pickle', 'fichier_4v95_11.pickle', 'fichier_4v5k_18.pickle', 'fichier_5ibb_14.pickle', 'fichier_4wro_44.pickle', 'fichier_6b4v_37.pickle', 'fichier_4woi_32.pickle', 'fichier_3ccj_3.pickle', 'fichier_4v67_7.pickle', 'fichier_4v7p_27.pickle', 'fichier_4v8g_9.pickle', 'fichier_4v8u_61.pickle', 'fichier_4v8u_36.pickle', 'fichier_4u1v_43.pickle', 'fichier_4u25_40.pickle', 'fichier_4lnt_23.pickle', 'fichier_4v8f_38.pickle', 'fichier_3cme_6.pickle', 'fichier_6n9e_8.pickle', 'fichier_1yij_5.pickle', 'fichier_4v7v_42.pickle', 'fichier_4p6f_17.pickle', 'fichier_6gsl_41.pickle', 'fichier_4v5q_10.pickle', 'fichier_5hcq_9.pickle', 'fichier_6h58_58.pickle', 'fichier_4v9s_11.pickle', 'fichier_5j4c_29.pickle', 'fichier_4tud_19.pickle', 'fichier_4wpo_2.pickle', 'fichier_4v9s_27.pickle', 'fichier_4v5e_44.pickle', 'fichier_4wqr_38.pickle', 'fichier_5el6_55.pickle', 'fichier_3cc2_6.pickle', 'fichier_3ccr_5.pickle', 'fichier_4v6f_2.pickle', 'fichier_5f8k_20.pickle', 'fichier_4v8c_33.pickle', 'fichier_4v5q_51.pickle', 'fichier_5v8i_51.pickle', 'fichier_4wq1_41.pickle', 'fichier_4z8c_26.pickle', 'fichier_4v5n_5.pickle', 'fichier_4v52_35.pickle', 'fichier_4v7z_17.pickle', 'fichier_5hcp_24.pickle', 'fichier_5el4_50.pickle', 'fichier_4w2h_13.pickle', 'fichier_3g6e_5.pickle', 'fichier_4wt8_53.pickle', 'fichier_4v8d_35.pickle', 'fichier_4wpo_28.pickle', 'fichier_4v9f_20.pickle', 'fichier_5czp_46.pickle', 'fichier_6c5l_46.pickle', 'fichier_4wfa_8.pickle', 'fichier_4v4q_35.pickle', 'fichier_4v85_27.pickle', 'fichier_4v7m_29.pickle', 'fichier_4v7m_11.pickle', 'fichier_5j4b_10.pickle', 'fichier_4v50_32.pickle', 'fichier_4v6a_26.pickle', 'fichier_4tuc_18.pickle', 'fichier_4v9i_36.pickle', 'fichier_4v9b_22.pickle', 'fichier_4ybb_31.pickle', 'fichier_4v9b_46.pickle', 'fichier_1vy6_29.pickle', 'fichier_4z3s_35.pickle', 'fichier_6nd5_32.pickle', 'fichier_6nd5_5.pickle', 'fichier_5j4d_40.pickle', 'fichier_4zer_8.pickle', 'fichier_4v95_10.pickle', 'fichier_4yzv_22.pickle', 'fichier_5hd1_10.pickle', 'fichier_4v9m_31.pickle', 'fichier_6gsl_49.pickle', 'fichier_4v5d_20.pickle', 'fichier_6cfk_22.pickle', 'fichier_5czp_44.pickle', 'fichier_6enj_9.pickle', 'fichier_5j8b_3.pickle', 'fichier_4v8i_30.pickle', 'fichier_4v9l_6.pickle', 'fichier_4lnt_24.pickle', 'fichier_4w2h_54.pickle', 'fichier_4v9k_8.pickle', 'fichier_4v5m_5.pickle', 'fichier_4v57_33.pickle', 'fichier_5jc9_29.pickle', 'fichier_4wf9_6.pickle', 'fichier_4v8g_31.pickle', 'fichier_4v8e_36.pickle', 'fichier_2j28_14.pickle', 'fichier_4wu1_51.pickle', 'fichier_4v8d_23.pickle', 'fichier_3g4s_4.pickle', 'fichier_3i56_4.pickle', 'fichier_4u27_43.pickle', 'fichier_5ndk_11.pickle', 'fichier_4v9n_11.pickle', 'fichier_4v8n_9.pickle', 'fichier_4v5j_36.pickle', 'fichier_5el5_39.pickle', 'fichier_4v9r_54.pickle', 'fichier_6o97_7.pickle', 'fichier_4v9h_6.pickle', 'fichier_4l47_21.pickle', 'fichier_6i0y_13.pickle', 'fichier_4wr6_36.pickle', 'fichier_4zsn_48.pickle', 'fichier_4v7v_24.pickle', 'fichier_4v8i_10.pickle', 'fichier_5e7k_41.pickle', 'fichier_5it8_27.pickle', 'fichier_5j3c_50.pickle', 'fichier_5j91_32.pickle', 'fichier_4v9i_11.pickle', 'fichier_4zsn_19.pickle', 'fichier_4v9b_45.pickle', 'fichier_4v7s_45.pickle', 'fichier_4w4g_12.pickle', 'fichier_4z3s_7.pickle', 'fichier_5wit_5.pickle', 'fichier_6of1_32.pickle', 'fichier_6n9f_8.pickle', 'fichier_4v54_31.pickle', 'fichier_1vy5_47.pickle', 'fichier_4w2e_27.pickle', 'fichier_5lzd_15.pickle', 'fichier_6cfl_22.pickle', 'fichier_6cae_5.pickle', 'fichier_4v9k_30.pickle', 'fichier_4v97_8.pickle', 'fichier_5it8_28.pickle', 'fichier_3i55_6.pickle', 'fichier_5j4c_3.pickle', 'fichier_4uy8_14.pickle', 'fichier_1vy7_51.pickle', 'fichier_4v9q_13.pickle', 'fichier_4v7k_32.pickle', 'fichier_5j7l_29.pickle', 'fichier_5hcr_22.pickle', 'fichier_4v8f_35.pickle', 'fichier_4v6t_26.pickle', 'fichier_4v8n_56.pickle', 'fichier_5fdu_27.pickle', 'fichier_6nd6_28.pickle', 'fichier_5ndk_44.pickle', 'fichier_5dox_31.pickle', 'fichier_4v53_33.pickle', 'fichier_4v53_32.pickle', 'fichier_5el7_48.pickle', 'fichier_3j5l_12.pickle', 'fichier_4v6s_11.pickle', 'fichier_5dox_18.pickle', 'fichier_4v5f_9.pickle', 'fichier_6cfk_9.pickle', 'fichier_5ndk_53.pickle', 'fichier_6nd6_5.pickle', 'fichier_6gc6_3.pickle', 'fichier_4v5j_11.pickle', 'fichier_6cae_34.pickle', 'fichier_4wsd_8.pickle', 'fichier_1vy7_31.pickle', 'fichier_6boh_46.pickle', 'fichier_4v9d_31.pickle', 'fichier_4wqf_54.pickle', 'fichier_5o2r_15.pickle', 'fichier_4v95_32.pickle', 'fichier_4lsk_22.pickle', 'fichier_5gad_11.pickle', 'fichier_5uyp_25.pickle', 'fichier_4v5y_32.pickle', 'fichier_4v5g_31.pickle', 'fichier_4u20_44.pickle', 'fichier_4v5r_16.pickle', 'fichier_1vy5_6.pickle', 'fichier_4v5d_56.pickle', 'fichier_6gc0_6.pickle', 'fichier_4v8n_30.pickle', 'fichier_5el6_8.pickle', 'fichier_3ccu_4.pickle', 'fichier_5j4d_46.pickle', 'fichier_4v90_7.pickle', 'fichier_4v87_41.pickle', 'fichier_4tua_48.pickle', 'fichier_3cd6_3.pickle', 'fichier_4v5r_59.pickle', 'fichier_4v7z_15.pickle', 'fichier_6cfj_29.pickle', 'fichier_1vy4_2.pickle', 'fichier_5nco_11.pickle', 'fichier_4tue_47.pickle', 'fichier_6bok_31.pickle', 'fichier_5ndj_41.pickle', 'fichier_4tua_51.pickle', 'fichier_4v97_30.pickle', 'fichier_4wqf_33.pickle', 'fichier_4v87_38.pickle', 'fichier_5vp2_8.pickle', 'fichier_4wr6_8.pickle', 'fichier_4v5q_56.pickle', 'fichier_5ndj_51.pickle', 'fichier_4v8b_38.pickle', 'fichier_5kpw_11.pickle', 'fichier_5w4k_3.pickle', 'fichier_4lel_20.pickle', 'fichier_4v67_22.pickle', 'fichier_4wt1_9.pickle', 'fichier_5hl7_8.pickle', 'fichier_4w2h_12.pickle', 'fichier_4v7m_12.pickle', 'fichier_4v8h_46.pickle', 'fichier_4v64_34.pickle', 'fichier_5dfe_47.pickle', 'fichier_4wra_20.pickle', 'fichier_4u24_24.pickle', 'fichier_5ibb_45.pickle', 'fichier_3g71_3.pickle', 'fichier_4lsk_14.pickle', 'fichier_4wt1_46.pickle', 'fichier_4wr6_45.pickle', 'fichier_5f8k_9.pickle', 'fichier_4v7z_16.pickle', 'fichier_6of1_5.pickle', 'fichier_4tua_21.pickle', 'fichier_5j5b_33.pickle', 'fichier_4v5r_34.pickle', 'fichier_5mdz_27.pickle', 'fichier_5nwy_7.pickle', 'fichier_4p6f_44.pickle', 'fichier_5czp_20.pickle', 'fichier_4wqy_11.pickle', 'fichier_4l71_13.pickle', 'fichier_4v50_33.pickle', 'fichier_6gsj_35.pickle', 'fichier_6b4v_10.pickle', 'fichier_6i7v_23.pickle', 'fichier_6fkr_22.pickle', 'fichier_5hd1_24.pickle', 'fichier_5afi_13.pickle', 'fichier_4wqf_10.pickle', 'fichier_4v7x_21.pickle', 'fichier_4v9r_11.pickle', 'fichier_4v8x_23.pickle', 'fichier_4zer_22.pickle', 'fichier_4v8j_38.pickle', 'fichier_5doy_4.pickle', 'fichier_6o97_5.pickle', 'fichier_3jbv_18.pickle', 'fichier_6n9f_20.pickle', 'fichier_5el7_54.pickle', 'fichier_5mdy_23.pickle', 'fichier_4v9d_30.pickle', 'fichier_4v71_10.pickle', 'fichier_5vpo_46.pickle', 'fichier_5vpo_17.pickle', 'fichier_4v83_9.pickle', 'fichier_4v8x_2.pickle', 'fichier_5fdv_37.pickle', 'fichier_4zsn_20.pickle', 'fichier_5j91_33.pickle', 'fichier_4u27_21.pickle', 'fichier_4v9a_49.pickle', 'fichier_5wis_22.pickle', 'fichier_4v51_25.pickle', 'fichier_3cc4_5.pickle', 'fichier_4wu1_44.pickle', 'fichier_4www_28.pickle', 'fichier_4woi_49.pickle', 'fichier_4yzv_13.pickle', 'fichier_4v5p_9.pickle', 'fichier_4w2g_52.pickle', 'fichier_4v5s_59.pickle', 'fichier_4y4p_5.pickle', 'fichier_4v8a_37.pickle', 'fichier_4v6f_50.pickle', 'fichier_4v5s_12.pickle', 'fichier_4v7u_25.pickle', 'fichier_4v9a_41.pickle', 'fichier_5fdu_9.pickle', 'fichier_6gsk_35.pickle', 'fichier_1vy6_51.pickle', 'fichier_4lfz_13.pickle']
#     liste_100_28 = ['fichier_5j8a_30.pickle', 'fichier_4ybb_30.pickle', 'fichier_4v84_33.pickle', 'fichier_5gaf_4.pickle', 'fichier_4wqu_11.pickle', 'fichier_5v8i_39.pickle', 'fichier_4v5p_32.pickle', 'fichier_4v4q_36.pickle', 'fichier_4v6d_46.pickle', 'fichier_4v8o_6.pickle', 'fichier_4v8g_8.pickle', 'fichier_6h58_29.pickle', 'fichier_5j8b_27.pickle', 'fichier_5dc3_19.pickle', 'fichier_4v52_34.pickle', 'fichier_5dge_27.pickle', 'fichier_5h5u_21.pickle', 'fichier_5uym_34.pickle', 'fichier_4l71_19.pickle', 'fichier_4v7j_31.pickle', 'fichier_4u20_25.pickle', 'fichier_4v9s_51.pickle', 'fichier_4wqy_12.pickle', 'fichier_4v7y_24.pickle', 'fichier_4lnt_59.pickle', 'fichier_6n8j_8.pickle', 'fichier_4u52_24.pickle', 'fichier_5wf0_12.pickle', 'fichier_3jcd_11.pickle', 'fichier_5dgf_15.pickle', 'fichier_4v8b_32.pickle', 'fichier_4v5n_26.pickle', 'fichier_4v8y_2.pickle', 'fichier_4wq1_9.pickle', 'fichier_4v5m_23.pickle', 'fichier_5hau_23.pickle', 'fichier_1vy7_11.pickle', 'fichier_4v6e_43.pickle', 'fichier_5fdv_12.pickle', 'fichier_6of1_3.pickle', 'fichier_4v88_7.pickle', 'fichier_6bok_47.pickle', 'fichier_4v5e_3.pickle', 'fichier_5ib7_55.pickle', 'fichier_4u3m_16.pickle', 'fichier_4wt8_52.pickle', 'fichier_6hhq_26.pickle', 'fichier_5tgm_5.pickle', 'fichier_5fcj_27.pickle', 'fichier_5j30_47.pickle', 'fichier_4tua_19.pickle', 'fichier_5ibb_51.pickle', 'fichier_5j88_27.pickle', 'fichier_4v7z_24.pickle', 'fichier_6c4i_16.pickle', 'fichier_4zsn_12.pickle', 'fichier_5mdv_25.pickle', 'fichier_4v4h_35.pickle', 'fichier_5on6_27.pickle', 'fichier_4v5l_6.pickle', 'fichier_4v5l_31.pickle', 'fichier_4v7y_23.pickle', 'fichier_6hrm_14.pickle', 'fichier_4u26_42.pickle', 'fichier_1vvj_13.pickle', 'fichier_4v8a_32.pickle', 'fichier_4w29_23.pickle', 'fichier_4v5f_51.pickle', 'fichier_5czp_18.pickle', 'fichier_4v5q_33.pickle', 'fichier_4p6f_47.pickle', 'fichier_4w4g_20.pickle', 'fichier_4u4q_15.pickle', 'fichier_4wsd_45.pickle', 'fichier_4v6g_39.pickle', 'fichier_6nd6_4.pickle', 'fichier_4v8g_53.pickle', 'fichier_4wzo_8.pickle', 'fichier_4wsd_49.pickle', 'fichier_5el6_51.pickle', 'fichier_4v5s_64.pickle', 'fichier_6i7v_22.pickle', 'fichier_6n8n_3.pickle', 'fichier_1vy4_7.pickle', 'fichier_6hhq_20.pickle', 'fichier_4v55_29.pickle', 'fichier_5j5b_32.pickle', 'fichier_4v9a_22.pickle', 'fichier_4v8d_40.pickle', 'fichier_4w2g_30.pickle', 'fichier_4v7w_22.pickle', 'fichier_4u56_15.pickle', 'fichier_4v6c_45.pickle', 'fichier_4v5d_14.pickle', 'fichier_5wit_30.pickle', 'fichier_5fci_17.pickle', 'fichier_4wt8_1.pickle', 'fichier_4v5d_27.pickle', 'fichier_4v8i_9.pickle', 'fichier_5el4_12.pickle', 'fichier_5wdt_14.pickle', 'fichier_5hcp_10.pickle', 'fichier_4v6p_3.pickle', 'fichier_4wpo_55.pickle', 'fichier_4v9j_6.pickle', 'fichier_4v6i_2.pickle', 'fichier_4u25_19.pickle', 'fichier_6c5l_11.pickle', 'fichier_4v8u_55.pickle', 'fichier_4v5f_56.pickle', 'fichier_4wqu_34.pickle', 'fichier_4p70_14.pickle', 'fichier_4u3m_25.pickle', 'fichier_4v8z_2.pickle', 'fichier_6gsj_12.pickle', 'fichier_4v54_32.pickle', 'fichier_4v5c_11.pickle', 'fichier_5wis_7.pickle', 'fichier_4wqy_32.pickle', 'fichier_5dgv_14.pickle', 'fichier_6bok_21.pickle', 'fichier_6q95_2.pickle', 'fichier_4p70_21.pickle', 'fichier_4u3u_31.pickle', 'fichier_4u4o_13.pickle', 'fichier_3jcn_15.pickle', 'fichier_4v9m_8.pickle', 'fichier_5ib8_46.pickle', 'fichier_4w29_6.pickle', 'fichier_4xej_29.pickle', 'fichier_5vpp_38.pickle', 'fichier_3j7z_13.pickle', 'fichier_6bz7_44.pickle', 'fichier_5dfe_20.pickle', 'fichier_4v6n_11.pickle', 'fichier_4u4y_14.pickle', 'fichier_4v6e_25.pickle', 'fichier_5e7k_10.pickle', 'fichier_4v7l_9.pickle', 'fichier_4v8n_13.pickle', 'fichier_4ypb_19.pickle', 'fichier_4v51_26.pickle', 'fichier_4v4n_10.pickle', 'fichier_4u4z_16.pickle', 'fichier_4u53_20.pickle', 'fichier_4ypb_48.pickle', 'fichier_4v9h_32.pickle', 'fichier_4lsk_23.pickle', 'fichier_1vy4_52.pickle', 'fichier_5doy_28.pickle', 'fichier_4u51_14.pickle', 'fichier_5juu_9.pickle', 'fichier_3bbx_5.pickle', 'fichier_4wsm_42.pickle', 'fichier_4v9q_20.pickle', 'fichier_4wzd_13.pickle', 'fichier_5fcj_19.pickle', 'fichier_4u52_16.pickle', 'fichier_4v90_28.pickle', 'fichier_5jup_13.pickle', 'fichier_4v5r_10.pickle', 'fichier_4v9q_12.pickle', 'fichier_4wsd_38.pickle', 'fichier_6n8j_9.pickle', 'fichier_5mgp_14.pickle', 'fichier_4wu1_13.pickle', 'fichier_4y4o_12.pickle', 'fichier_4w2f_7.pickle', 'fichier_4u4n_19.pickle', 'fichier_6n8m_3.pickle', 'fichier_4v8a_49.pickle', 'fichier_5ib7_48.pickle', 'fichier_4y4o_27.pickle', 'fichier_5gae_14.pickle', 'fichier_4v5s_38.pickle', 'fichier_6czr_9.pickle', 'fichier_4lel_19.pickle', 'fichier_4p6f_19.pickle', 'fichier_5wfs_11.pickle', 'fichier_3j9z_18.pickle', 'fichier_4zer_31.pickle', 'fichier_1vvj_19.pickle', 'fichier_4v67_49.pickle', 'fichier_5fci_25.pickle', 'fichier_5el4_43.pickle', 'fichier_4wqu_54.pickle', 'fichier_5ib7_16.pickle', 'fichier_4tud_47.pickle', 'fichier_5gah_11.pickle', 'fichier_5iqr_29.pickle', 'fichier_4v5y_33.pickle', 'fichier_4v6q_8.pickle', 'fichier_4u1v_23.pickle', 'fichier_4wqr_48.pickle', 'fichier_4v9n_48.pickle', 'fichier_4w2g_11.pickle', 'fichier_5e81_49.pickle', 'fichier_6bz8_23.pickle', 'fichier_5vp2_22.pickle', 'fichier_4wq1_50.pickle', 'fichier_6gsj_50.pickle', 'fichier_5hau_33.pickle', 'fichier_4v8e_39.pickle', 'fichier_4v9c_25.pickle', 'fichier_6bok_48.pickle', 'fichier_4v95_55.pickle', 'fichier_4v8i_52.pickle', 'fichier_4v7m_55.pickle', 'fichier_4v84_8.pickle', 'fichier_6buw_19.pickle', 'fichier_4v9b_37.pickle', 'fichier_4wt1_36.pickle', 'fichier_4v6o_5.pickle', 'fichier_6c5l_34.pickle', 'fichier_5j4b_25.pickle', 'fichier_6buw_46.pickle', 'fichier_5hau_8.pickle', 'fichier_4wzd_43.pickle', 'fichier_4v8c_31.pickle', 'fichier_4tub_17.pickle', 'fichier_5hau_37.pickle', 'fichier_6cfl_10.pickle', 'fichier_4v9n_12.pickle', 'fichier_4v97_7.pickle', 'fichier_4wf1_23.pickle', 'fichier_4u4u_27.pickle', 'fichier_4v5g_11.pickle', 'fichier_6q95_3.pickle', 'fichier_5tga_21.pickle', 'fichier_4v9j_35.pickle', 'fichier_4lfz_42.pickle', 'fichier_4v4h_36.pickle', 'fichier_5hcq_24.pickle', 'fichier_5el7_13.pickle', 'fichier_4l47_14.pickle', 'fichier_6o97_30.pickle', 'fichier_4z8c_10.pickle', 'fichier_5kcr_12.pickle', 'fichier_4v7b_27.pickle', 'fichier_5j3c_48.pickle', 'fichier_4v7j_20.pickle', 'fichier_4v5f_34.pickle', 'fichier_4v8h_22.pickle', 'fichier_4v83_35.pickle', 'fichier_4u3n_17.pickle', 'fichier_4u26_23.pickle', 'fichier_4v7k_26.pickle', 'fichier_4z8c_39.pickle', 'fichier_5dge_17.pickle', 'fichier_5fdv_40.pickle', 'fichier_4v64_35.pickle', 'fichier_5mdw_25.pickle', 'fichier_4u4o_19.pickle', 'fichier_5j30_18.pickle', 'fichier_4wpo_64.pickle', 'fichier_6boh_40.pickle', 'fichier_4lnt_14.pickle', 'fichier_5dc3_11.pickle', 'fichier_5w4k_29.pickle', 'fichier_4v55_30.pickle', 'fichier_5hcq_34.pickle', 'fichier_6q9a_23.pickle', 'fichier_4v8e_19.pickle', 'fichier_4wqy_53.pickle', 'fichier_4v6g_25.pickle', 'fichier_4v8h_31.pickle', 'fichier_6cfj_3.pickle', 'fichier_4v8q_8.pickle', 'fichier_5vpo_15.pickle', 'fichier_6gsj_44.pickle', 'fichier_4u3n_24.pickle', 'fichier_4ypb_12.pickle', 'fichier_4v72_9.pickle', 'fichier_4y4p_33.pickle', 'fichier_3jbu_15.pickle', 'fichier_4v9j_7.pickle', 'fichier_4lt8_14.pickle', 'fichier_5e81_41.pickle', 'fichier_4v8j_11.pickle', 'fichier_4v6r_6.pickle', 'fichier_4wro_7.pickle', 'fichier_4v9j_32.pickle', 'fichier_4w2i_6.pickle', 'fichier_5el5_46.pickle', 'fichier_6enu_17.pickle', 'fichier_4v5c_53.pickle', 'fichier_4v9r_10.pickle', 'fichier_5hcq_39.pickle', 'fichier_4v7t_43.pickle', 'fichier_4v7u_45.pickle', 'fichier_4wf1_44.pickle', 'fichier_3j9y_13.pickle', 'fichier_4u4y_21.pickle', 'fichier_5j8a_31.pickle', 'fichier_4v9s_10.pickle', 'fichier_5fdv_27.pickle', 'fichier_4v7k_43.pickle', 'fichier_4v97_52.pickle', 'fichier_6fkr_16.pickle', 'fichier_6i7o_18.pickle', 'fichier_4v57_34.pickle', 'fichier_5t6r_4.pickle', 'fichier_4w2h_31.pickle', 'fichier_4w2f_48.pickle', 'fichier_4u4r_16.pickle', 'fichier_1vy6_11.pickle', 'fichier_5kpx_13.pickle', 'fichier_5ib8_12.pickle', 'fichier_4u6f_13.pickle', 'fichier_4v5k_2.pickle', 'fichier_4v7x_20.pickle', 'fichier_5fdu_40.pickle', 'fichier_4wsm_8.pickle', 'fichier_6h4n_13.pickle', 'fichier_4v7p_12.pickle', 'fichier_5tbw_5.pickle', 'fichier_4v7k_14.pickle', 'fichier_4v67_1.pickle', 'fichier_5j7l_30.pickle', 'fichier_4wzo_43.pickle', 'fichier_4v8u_11.pickle', 'fichier_4u6f_20.pickle', 'fichier_4v6f_39.pickle', 'fichier_4u55_17.pickle', 'fichier_4u24_44.pickle', 'fichier_4lt8_22.pickle', 'fichier_4v9n_45.pickle', 'fichier_6cae_3.pickle', 'fichier_4wro_48.pickle', 'fichier_5jc9_30.pickle', 'fichier_5lyb_26.pickle', 'fichier_5el5_13.pickle', 'fichier_5aka_15.pickle', 'fichier_6fkr_26.pickle', 'fichier_4w2i_49.pickle', 'fichier_4zer_37.pickle', 'fichier_4v7w_21.pickle', 'fichier_5hcr_10.pickle', 'fichier_4wro_37.pickle', 'fichier_4v9q_15.pickle', 'fichier_5el6_44.pickle', 'fichier_4z3s_5.pickle', 'fichier_4v7l_26.pickle', 'fichier_4wsm_33.pickle', 'fichier_4v9r_30.pickle', 'fichier_4tue_19.pickle', 'fichier_4v68_7.pickle', 'fichier_5j88_28.pickle', 'fichier_4u1u_43.pickle', 'fichier_4v95_11.pickle', 'fichier_4v5k_18.pickle', 'fichier_5ndg_12.pickle', 'fichier_5ibb_14.pickle', 'fichier_4wro_44.pickle', 'fichier_6b4v_37.pickle', 'fichier_4woi_32.pickle', 'fichier_4v67_7.pickle', 'fichier_5mei_24.pickle', 'fichier_5lyb_19.pickle', 'fichier_4v7p_27.pickle', 'fichier_4v8g_9.pickle', 'fichier_4v8u_61.pickle', 'fichier_4v8u_36.pickle', 'fichier_4u1v_43.pickle', 'fichier_4u25_40.pickle', 'fichier_4lnt_23.pickle', 'fichier_4v8f_38.pickle', 'fichier_6n9e_8.pickle', 'fichier_6n8o_4.pickle', 'fichier_4v7v_42.pickle', 'fichier_4p6f_17.pickle', 'fichier_6gsl_41.pickle', 'fichier_4v5q_10.pickle', 'fichier_5hcq_9.pickle', 'fichier_6h58_58.pickle', 'fichier_4v9s_11.pickle', 'fichier_4v7r_4.pickle', 'fichier_5j4c_29.pickle', 'fichier_4tud_19.pickle', 'fichier_4wpo_2.pickle', 'fichier_4v9s_27.pickle', 'fichier_4v5e_44.pickle', 'fichier_4wqr_38.pickle', 'fichier_5ndv_15.pickle', 'fichier_5el6_55.pickle', 'fichier_4v6f_2.pickle', 'fichier_5f8k_20.pickle', 'fichier_4v8c_33.pickle', 'fichier_4v5q_51.pickle', 'fichier_5v8i_51.pickle', 'fichier_4wq1_41.pickle', 'fichier_4z8c_26.pickle', 'fichier_4v5n_5.pickle', 'fichier_4v52_35.pickle', 'fichier_4v7z_17.pickle', 'fichier_5hcp_24.pickle', 'fichier_5el4_50.pickle', 'fichier_4w2h_13.pickle', 'fichier_4wt8_53.pickle', 'fichier_4v8d_35.pickle', 'fichier_4wpo_28.pickle', 'fichier_5czp_46.pickle', 'fichier_6c5l_46.pickle', 'fichier_4v4q_35.pickle', 'fichier_4v85_27.pickle', 'fichier_5i4l_20.pickle', 'fichier_4v7m_29.pickle', 'fichier_4v7m_11.pickle', 'fichier_5j4b_10.pickle', 'fichier_4v50_32.pickle', 'fichier_4v6a_26.pickle', 'fichier_4tuc_18.pickle', 'fichier_4v9i_36.pickle', 'fichier_4v9b_22.pickle', 'fichier_4ybb_31.pickle', 'fichier_4v9b_46.pickle', 'fichier_1vy6_29.pickle', 'fichier_4z3s_35.pickle', 'fichier_5obm_14.pickle', 'fichier_6nd5_32.pickle', 'fichier_6nd5_5.pickle', 'fichier_5j4d_40.pickle', 'fichier_4zer_8.pickle', 'fichier_4v95_10.pickle', 'fichier_4yzv_22.pickle', 'fichier_5hd1_10.pickle', 'fichier_4v9m_31.pickle', 'fichier_6gsl_49.pickle', 'fichier_4v5d_20.pickle', 'fichier_6cfk_22.pickle', 'fichier_5czp_44.pickle', 'fichier_6enj_9.pickle', 'fichier_4v6u_10.pickle', 'fichier_5j8b_3.pickle', 'fichier_4v8i_30.pickle', 'fichier_4v9l_6.pickle', 'fichier_4lnt_24.pickle', 'fichier_4w2h_54.pickle', 'fichier_4v9k_8.pickle', 'fichier_4u4q_25.pickle', 'fichier_4v5m_5.pickle', 'fichier_4v57_33.pickle', 'fichier_5jc9_29.pickle', 'fichier_4u4n_30.pickle', 'fichier_4v8g_31.pickle', 'fichier_4v8e_36.pickle', 'fichier_2j28_14.pickle', 'fichier_5tga_16.pickle', 'fichier_4wu1_51.pickle', 'fichier_4v8d_23.pickle', 'fichier_4u27_43.pickle', 'fichier_5ndk_11.pickle', 'fichier_4v9n_11.pickle', 'fichier_4v8n_9.pickle', 'fichier_4v5j_36.pickle', 'fichier_5dat_15.pickle', 'fichier_4u4z_26.pickle', 'fichier_5el5_39.pickle', 'fichier_4v9r_54.pickle', 'fichier_6o97_7.pickle', 'fichier_4v9h_6.pickle', 'fichier_5we4_15.pickle', 'fichier_4l47_21.pickle', 'fichier_6i0y_13.pickle', 'fichier_4wr6_36.pickle', 'fichier_4zsn_48.pickle', 'fichier_6cb1_3.pickle', 'fichier_4v7v_24.pickle', 'fichier_4v8i_10.pickle', 'fichier_5e7k_41.pickle', 'fichier_4u4r_26.pickle', 'fichier_5it8_27.pickle', 'fichier_4u50_18.pickle', 'fichier_5j3c_50.pickle', 'fichier_5dat_24.pickle', 'fichier_5j91_32.pickle', 'fichier_4v9i_11.pickle', 'fichier_4zsn_19.pickle', 'fichier_4v9b_45.pickle', 'fichier_4v7s_45.pickle', 'fichier_4w4g_12.pickle', 'fichier_4z3s_7.pickle', 'fichier_5wit_5.pickle', 'fichier_6of1_32.pickle', 'fichier_6n9f_8.pickle', 'fichier_4v54_31.pickle', 'fichier_1vy5_47.pickle', 'fichier_5jcs_1.pickle', 'fichier_4w2e_27.pickle', 'fichier_5lzd_15.pickle', 'fichier_5m1j_6.pickle', 'fichier_6cfl_22.pickle', 'fichier_6cae_5.pickle', 'fichier_4v9k_30.pickle', 'fichier_4v97_8.pickle', 'fichier_5it8_28.pickle', 'fichier_5mei_17.pickle', 'fichier_4u4u_18.pickle', 'fichier_4v7r_8.pickle', 'fichier_5j4c_3.pickle', 'fichier_4uy8_14.pickle', 'fichier_1vy7_51.pickle', 'fichier_4v9q_13.pickle', 'fichier_4v7k_32.pickle', 'fichier_5j7l_29.pickle', 'fichier_5hcr_22.pickle', 'fichier_4v8f_35.pickle', 'fichier_4v6t_26.pickle', 'fichier_4v8n_56.pickle', 'fichier_5fdu_27.pickle', 'fichier_6nd6_28.pickle', 'fichier_5ndk_44.pickle', 'fichier_5dox_31.pickle', 'fichier_4v53_33.pickle', 'fichier_4v53_32.pickle', 'fichier_5el7_48.pickle', 'fichier_5tbw_13.pickle', 'fichier_3j5l_12.pickle', 'fichier_4v6s_11.pickle', 'fichier_5dox_18.pickle', 'fichier_4v5f_9.pickle', 'fichier_6cfk_9.pickle', 'fichier_5ndk_53.pickle', 'fichier_6nd6_5.pickle', 'fichier_6gc6_3.pickle', 'fichier_4v5j_11.pickle', 'fichier_5dgv_23.pickle', 'fichier_6cae_34.pickle', 'fichier_4wsd_8.pickle', 'fichier_4u55_25.pickle', 'fichier_4v88_3.pickle', 'fichier_1vy7_31.pickle', 'fichier_6boh_46.pickle', 'fichier_4v9d_31.pickle', 'fichier_4wqf_54.pickle', 'fichier_4u51_23.pickle', 'fichier_5o2r_15.pickle', 'fichier_5apo_8.pickle', 'fichier_4v95_32.pickle', 'fichier_5obm_23.pickle', 'fichier_4lsk_22.pickle', 'fichier_5gad_11.pickle', 'fichier_3jct_8.pickle', 'fichier_5uyp_25.pickle', 'fichier_4v5y_32.pickle', 'fichier_5t62_4.pickle', 'fichier_4v5g_31.pickle', 'fichier_4u20_44.pickle', 'fichier_4v5r_16.pickle', 'fichier_4u53_28.pickle', 'fichier_1vy5_6.pickle', 'fichier_4v5d_56.pickle', 'fichier_6gc0_6.pickle', 'fichier_4v8n_30.pickle', 'fichier_5juo_7.pickle', 'fichier_5el6_8.pickle', 'fichier_5ndw_24.pickle', 'fichier_5j4d_46.pickle', 'fichier_4v90_7.pickle', 'fichier_4u3u_20.pickle', 'fichier_5fl8_1.pickle', 'fichier_4v87_41.pickle', 'fichier_4tua_48.pickle', 'fichier_5we6_13.pickle', 'fichier_4v5r_59.pickle', 'fichier_4v7z_15.pickle', 'fichier_6cfj_29.pickle', 'fichier_1vy4_2.pickle', 'fichier_5ndv_24.pickle', 'fichier_5nco_11.pickle', 'fichier_5on6_22.pickle', 'fichier_4tue_47.pickle', 'fichier_6n8l_9.pickle', 'fichier_6c0f_2.pickle', 'fichier_6bok_31.pickle', 'fichier_5ndj_41.pickle', 'fichier_4tua_51.pickle', 'fichier_4v97_30.pickle', 'fichier_5apn_7.pickle', 'fichier_4wqf_33.pickle', 'fichier_4v87_38.pickle', 'fichier_5vp2_8.pickle', 'fichier_4wr6_8.pickle', 'fichier_4v5q_56.pickle', 'fichier_5ndj_51.pickle', 'fichier_4v8b_38.pickle', 'fichier_5kpw_11.pickle', 'fichier_5w4k_3.pickle', 'fichier_4lel_20.pickle', 'fichier_4v67_22.pickle', 'fichier_4wt1_9.pickle', 'fichier_4w2h_12.pickle', 'fichier_4v7m_12.pickle', 'fichier_4v8h_46.pickle', 'fichier_4v64_34.pickle', 'fichier_5dfe_47.pickle', 'fichier_4wra_20.pickle', 'fichier_4u24_24.pickle', 'fichier_5ibb_45.pickle', 'fichier_4lsk_14.pickle', 'fichier_4wt1_46.pickle', 'fichier_4wr6_45.pickle', 'fichier_5f8k_9.pickle', 'fichier_4v7z_16.pickle', 'fichier_6of1_5.pickle', 'fichier_6n8o_2.pickle', 'fichier_4tua_21.pickle', 'fichier_4u50_26.pickle', 'fichier_5j5b_33.pickle', 'fichier_4v5r_34.pickle', 'fichier_5mdz_27.pickle', 'fichier_5nwy_7.pickle', 'fichier_4p6f_44.pickle', 'fichier_5czp_20.pickle', 'fichier_4wqy_11.pickle', 'fichier_4l71_13.pickle', 'fichier_4v50_33.pickle', 'fichier_4v8t_12.pickle', 'fichier_6gsj_35.pickle', 'fichier_6b4v_10.pickle', 'fichier_6i7v_23.pickle', 'fichier_6fkr_22.pickle', 'fichier_5hd1_24.pickle', 'fichier_5afi_13.pickle', 'fichier_4wqf_10.pickle', 'fichier_4v7x_21.pickle', 'fichier_4v9r_11.pickle', 'fichier_4v8x_23.pickle', 'fichier_4zer_22.pickle', 'fichier_4v8j_38.pickle', 'fichier_5doy_4.pickle', 'fichier_6o97_5.pickle', 'fichier_3jbv_18.pickle', 'fichier_5ndw_16.pickle', 'fichier_6n9f_20.pickle', 'fichier_5el7_54.pickle', 'fichier_5mdy_23.pickle', 'fichier_4v9d_30.pickle', 'fichier_4v71_10.pickle', 'fichier_5vpo_46.pickle', 'fichier_5vpo_17.pickle', 'fichier_4v83_9.pickle', 'fichier_5wfk_13.pickle', 'fichier_4v8x_2.pickle', 'fichier_5fdv_37.pickle', 'fichier_4zsn_20.pickle', 'fichier_5j91_33.pickle', 'fichier_4u27_21.pickle', 'fichier_4v9a_49.pickle', 'fichier_5wis_22.pickle', 'fichier_4v51_25.pickle', 'fichier_4wu1_44.pickle', 'fichier_6elz_7.pickle', 'fichier_4www_28.pickle', 'fichier_4woi_49.pickle', 'fichier_5tgm_8.pickle', 'fichier_4yzv_13.pickle', 'fichier_4v5p_9.pickle', 'fichier_4w2g_52.pickle', 'fichier_4v5s_59.pickle', 'fichier_5i4l_13.pickle', 'fichier_4y4p_5.pickle', 'fichier_4v8a_37.pickle', 'fichier_4v6f_50.pickle', 'fichier_6q8y_3.pickle', 'fichier_4v5s_12.pickle', 'fichier_5dgf_21.pickle', 'fichier_4v7u_25.pickle', 'fichier_4v9a_41.pickle', 'fichier_5fdu_9.pickle', 'fichier_6gsk_35.pickle', 'fichier_1vy6_51.pickle', 'fichier_4lfz_13.pickle']
#     print(len(liste_100_vert_jaune))
#     print(len(liste_100_28))
#     
#     diff1 = []
#     diff2 = []
#     idem = []
#     for elt in liste_100_vert_jaune :
#         if elt in liste_100_28 : 
#             if elt not in idem :
#                 print(elt)
#                 idem.append(elt)
#         else :
#             if elt not in diff1 :
#                 print(elt)
#                 diff1.append(elt)
#     print(idem)
#     print(len(idem))
             
#     for elt in liste_100_28 : 
#         if elt in liste_100_vert_jaune :
#             if elt not in idem :
#                 idem.append(elt)
#         else :
#             if elt not in diff2 :
#                 print(elt)
#                 diff2.append(elt)
#     print(idem)
#     print(len(idem))
#      
#     print(diff1)
#     print(len(diff1))
#     print(diff2)
#     print(len(diff2))
    
    
#     with open("groupes_25S_homologues.pickle", "rb") as fichier_homologues_1 :
#         mon_depickler_1 = pickle.Unpickler(fichier_homologues_1)
#         groupes_homologues_1 = mon_depickler_1.load()
#           
#         homologues = []
#         gr_hom = []
#         for elt in idem :
#             cle = (elt.split("_")[1], int(elt.split("_")[2][:len(elt.split("_")[2])-7]))
#             compteur = 1
#            
#             for groupe in groupes_homologues_1 :
# #                 print(groupe)
# #                 print(cle)
#                 if cle in groupe :
#                     #if compteur == 87 :
#                         if compteur not in gr_hom :
#                             gr_hom.append(compteur)
#                         print("groupe : %d"%compteur)
#                         print("in")
#                         homologues.append(cle)
#                 compteur += 1
#         print(homologues)
#         print(len(homologues))
#         print(gr_hom)
    
    #plt.clf()  
    #pourcentage_ressemblance_sequence()
    
#     with open(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes"+"digraphe_commun_0.7_clustering_perez_groupe_11_taille_4.pickle", 'rb') as fichier_pickle :
#         mon_depickler = pickle.Unpickler(fichier_pickle)
#         graphe  = mon_depickler.load()
#         print(graphe.nodes.data())
#     with open("fichier_not_idem.txt", 'w') as fichier :
#         #for j in range(len(GROUPE_JAUNE_VERT_BLEU)) :
#             #if j == 0  :# and j != 12 and len(CLUSTERING_PEREZ_VERSION_NON_CAN_2[j]) > 2 and j > 6 : 
#                 liste_permut = list(itertools.permutations(GROUPE_JAUNE_VERT_BLEU))
#                 print(liste_permut)
#                 print(len(liste_permut))
# #                     
# #                 
#                 for i in range(4,11) :
#                     ancien_digraphe_commun = nx.MultiDiGraph()
#                     elt_ancien = -1
#                     ancienne_liste_cliques = []
#                     for elt in liste_permut : 
#                               
#                         digraphe_commun, liste_cliques = commun_cluster_clique(elt, EXTENSION_PATH%i+"dico_comp_complet_metrique_toutes_aretes_coeff_all1_taille_%s.pickle"%(i), EXTENSION_PATH_TAILLE%i)
#                         #fichier.write(str(elt)+ " " + str(digraphe_commun.nodes.data()) + "\n\n")
#                         #fichier.write(str(digraphe_commun.nodes.data())+"\n")
#                         if ancien_digraphe_commun.number_of_nodes() > 0 and not idem_graphes_communs_version_clique(ancien_digraphe_commun, digraphe_commun) :
#                             fichier.write( " " + str(elt)+ " " + str(digraphe_commun.nodes.data()) + "\n"+str(ancien_digraphe_commun.nodes.data())+" %d non ok\n"%i)
#                             fichier.write(str(elt)+ " " + str(digraphe_commun.edges.data()) + "\n"+str(ancien_digraphe_commun.edges.data())+" %d non ok\n"%i)
#                             print(digraphe_commun.number_of_nodes())
#                             print(ancien_digraphe_commun.number_of_nodes()) 
#                             print(elt)
#                             print(elt_ancien)
#                             exit()
#         #                 if len(liste_cliques) != len(ancienne_liste_cliques) and len(ancienne_liste_cliques) > 0 :
#         #                     fichier.write(str(elt)+ " " + str(liste_cliques) + "\n"+str(ancienne_liste_cliques)+" %d non ok\n"%i)
#         #                     print(elt)
#         #                     print(elt_ancien)
#         #                     exit()
#                         ancien_digraphe_commun = digraphe_commun.copy()
#                         ancienne_liste_cliques = list(liste_cliques)
#                         elt_ancien = elt
#                 print(ancien_digraphe_commun.nodes.data())

#     groupe = CLUSTERING_PEREZ_VERSION_NON_CAN_2[11]
#     for elt in CLUSTERING_PEREZ_VERSION_NON_CAN_2[12] :
#         if elt not in groupe :
#             groupe.append(elt)
#             
#     groupe.extend(CLUSTERING_PEREZ_VERSION_NON_CAN_2[1])
    
#     compteur = 0
#     for elt in CLUSTERING_PEREZ_VERSION_NON_CAN_2 :
#         if compteur == 14 :
#             for i in range(4, 11) :
#                 dico_graphes = recherche_graphe_commun([i], elt, 0.7, "toutes_aretes_coeff_all1", "taille_max/result_k_max_4_10_toutes_aretes", "clustering_perez_groupe_%d"%compteur)
#                 for cle in elt :
#                     os.makedirs(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes"+"seq_graphe_commun_groupe_14_clustering_perez", exist_ok = True)
#                     draw_seq_sous_graphe_commun(dico_graphes[cle][0], dico_graphes[cle][1], EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes"+"seq_graphe_commun_groupe_14_clustering_perez/", "graphe_%s_avec_coord_taille_%d"%(cle,i))
#             #             for i in range(4,11) :
# #                 recherche_similaires(i, "taille_max/result_k_max_4_10_toutes_aretes", 0.7, "clustering_perez_groupe_%d"%compteur)
#         compteur += 1

    #recherche_graphe_commun([8], CLUSTERING_PEREZ_VERSION_NON_CAN_2[11], 0.7, "toutes_aretes_coeff_all1", "taille_max/result_k_max_4_10_toutes_aretes", "clustering_perez_groupe_%d"%11)
#     with open(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes/sous_graphe_commun_clustering_perez_groupe_11_0.7/taille_8"+"couples_possibles_fichier_1FJG_A_48_8.pickle", "rb") as fichier :
#         mon_depickler = pickle.Unpickler(fichier)
#         couples_possibles = mon_depickler.load()
#         print(couples_possibles)
        
    #digraphe_commun = commun_cluster(CLUSTERING_PEREZ_VERSION_NON_CAN_2[11], EXTENSION_PATH_TAILLE%8, EXTENSION_PATH%8+"dico_comp_complet_metrique_%s_taille_%s.pickle"%("toutes_aretes_coeff_all1", 8), "clustering_perez_groupe_%d"%11, 8, 0.7, "taille_max/result_k_max_4_10_toutes_aretes")                
    #draw(digraphe_commun,"taille_max/result_k_max_4_10_toutes_aretes", "clustering_perez_groupe_%d"%11, 8)    
    #sous_graphe_commun_max_version_graphe_moyen("couples_possibles_fichier_1FJG_A_48_8.pickle", 8, "taille_max/result_k_max_4_10_toutes_aretes", 0.7, "clustering_perez_groupe_%d"%11)
        
        
    
''' 23/08/19 '''
#     with open("fichier_comp_sim_rmsd_par_seuil_tot_normalise_taille_4.pickle", 'rb') as fichier :
#         mon_depickler = pickle.Unpickler(fichier)
#         tab_max_par_seuil = mon_depickler.load()
#         
#         for elt in tab_max_par_seuil : 
#             #print(elt)
#             if round(elt[1],1) == 0.1 : 
#                 clustering_rmsd = elt[4]
#                 
#         
#         print(clustering_rmsd)
#         
#         print(len(clustering_rmsd))
#         
#         for elt in clustering_rmsd :
#             print(elt)
        
        
#         with open("fichier_not_idem.txt", 'w') as fichier :
#             for j in range(5,6) : #len(clustering_rmsd)) :
#             #if j == 0  :# and j != 12 and len(CLUSTERING_PEREZ_VERSION_NON_CAN_2[j]) > 2 and j > 6 : 
#                 liste_permut = list(itertools.permutations(clustering_rmsd[j]))
#                 print(liste_permut)
#                 print(len(liste_permut))
# #                     
# #                 
#                 for i in range(4,11) :
#                     ancien_digraphe_commun = nx.MultiDiGraph()
#                     elt_ancien = -1
#                     ancienne_liste_cliques = []
#                     for elt in liste_permut : 
#                         print(i)
#                         print(j)
#                         print(elt)
#                         if len(elt) > 2 : 
#                             digraphe_commun, liste_cliques = commun_cluster_clique(elt, EXTENSION_PATH%i+"dico_comp_complet_metrique_toutes_aretes_coeff_all1_taille_%s.pickle"%(i), EXTENSION_PATH_TAILLE%i)
#                         #fichier.write(str(elt)+ " " + str(digraphe_commun.nodes.data()) + "\n\n")
#                         #fichier.write(str(digraphe_commun.nodes.data())+"\n")
#                         if ancien_digraphe_commun.number_of_nodes() > 0 and not idem_graphes_communs_version_clique(ancien_digraphe_commun, digraphe_commun) :
#                             fichier.write( " " + str(elt)+ " " + str(digraphe_commun.nodes.data()) + "\n"+str(ancien_digraphe_commun.nodes.data())+" %d non ok\n"%i)
#                             fichier.write(str(elt)+ " " + str(digraphe_commun.edges.data()) + "\n"+str(ancien_digraphe_commun.edges.data())+" %d non ok\n"%i)
#                             print(digraphe_commun.number_of_nodes())
#                             print(ancien_digraphe_commun.number_of_nodes()) 
#                             print(elt)
#                             print(elt_ancien)
#                             exit()
#         #                 if len(liste_cliques) != len(ancienne_liste_cliques) and len(ancienne_liste_cliques) > 0 :
#         #                     fichier.write(str(elt)+ " " + str(liste_cliques) + "\n"+str(ancienne_liste_cliques)+" %d non ok\n"%i)
#         #                     print(elt)
#         #                     print(elt_ancien)
#         #                     exit()
#                         ancien_digraphe_commun = digraphe_commun.copy()
#                         ancienne_liste_cliques = list(liste_cliques)
#                         elt_ancien = elt
#                 print(ancien_digraphe_commun.nodes.data())
                

   
#     for i in range(4, 11) :
#                 dico_graphes = recherche_graphe_commun([i], clustering_rmsd[5], 0.7, "toutes_aretes_coeff_all1", "taille_max/result_k_max_4_10_toutes_aretes", "clustering_perez_rmsd_groupe_%d"%5)
#                 for cle in clustering_rmsd[5] :
#                     os.makedirs(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes"+"seq_graphe_commun_rmsd_groupe_5_clustering_perez", exist_ok = True)
#                     draw_seq_sous_graphe_commun(dico_graphes[cle][0], dico_graphes[cle][1], EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes"+"seq_graphe_commun_rmsd_groupe_5_clustering_perez/", "graphe_%s_avec_coord_taille_%d"%(cle,i))
#             #             for i in range(4,11) :
#                 recherche_similaires(i, "taille_max/result_k_max_4_10_toutes_aretes", 0.7, "clustering_perez_rmsd_groupe_%d"%5)
    
    
''' Pourcentage de ressemblance au niveau structure avec le groupe 11 '''
    
    
#     compteur = 0
#     distrib = []#             if nb_non_cov + nb_cov >= 19 :
# #                 print(fic)
# #                 print(nb_non_cov)
# #                 print(nb_cov)
#     exact = 0
#     liste_100_pourcent = []
#     for fic in os.listdir(NEW_EXTENSION_PATH_TAILLE) :
#         if "pickle" in fic and "fichier" in fic :
#             #print(compteur)
#             #print(fic)
#             nb_liaisons_non_cov_correctes = 0
#             nb_liaisons_cov_correctes = 0
#             nb_non_cov, nb_cov = recherche_graphe_commun_version_rapide(fic)
#             distrib.append((nb_non_cov + nb_cov)/19)
#             if nb_non_cov == 10 and nb_cov == 9 :
#                 exact += 1
#                 print(fic)
#                 liste_100_pourcent.append(fic)
# #             if nb_non_cov + nb_cov >= 19 :
# #                 print(fic)
# #                 print(nb_non_cov)
# #                 print(nb_cov)
#                     
#             #print((nb_non_cov, nb_cov))
#             compteur += 1
#                
#     print(liste_100_pourcent)
#     with open("groupes_23S_homologues.pickle", "rb") as fichier_homologues_1 :
#         mon_depickler_1 = pickle.Unpickler(fichier_homologues_1)
#         groupes_homologues_1 = mon_depickler_1.load()
#          
#         homologues_23s = []
#         for elt in liste_100_pourcent :
#             cle = (elt.split("_")[1], int(elt.split("_")[2][:len(elt.split("_")[2])-7]))
#             compteur = 1
#             for groupe in groupes_homologues_1 :
# #                 print(groupe)
# #                 print(cle)
#                 if cle in groupe :
#                     if compteur == 18 :
#                         print("groupe : %d"%compteur)
#                         print("in")
#                         homologues_23s.append(cle)
#                      
#                 compteur += 1
#          
#          
#          
#         with open("groupes_23S_identiques.pickle", "rb") as fichier_identiques_1 :
#             mon_depickler_1 = pickle.Unpickler(fichier_identiques_1)
#             groupes_identiques_1 = mon_depickler_1.load()
#              
#              
#             for elt in homologues_23s :
#                 compteur1 = 1
#                 for groupes in groupes_identiques_1 :
#                     compteur2 = 1
#                     for groupe in groupes :
#                         if elt in groupe :
#                             print("groupe : %s %s"%(compteur1, compteur2))
#                             print(elt)
#                             print("in")
#                         compteur2 += 1
#                     compteur1 += 1
#                      
#     with open("groupes_16S_homologues.pickle", "rb") as fichier_homologues_1 :
#         mon_depickler_1 = pickle.Unpickler(fichier_homologues_1)
#         groupes_homologues_1 = mon_depickler_1.load()
#          
#         homologues_16s = []
#         for elt in liste_100_pourcent :
#             cle = (elt.split("_")[1], int(elt.split("_")[2][:len(elt.split("_")[2])-7]))
#             compteur = 1
#             for groupe in groupes_homologues_1 :
# #                 print(groupe)
# #                 print(cle)
#                 if cle in groupe :
#                     print("groupe : %d"%compteur)
#                     print("in")
#                     homologues_16s.append(cle)
#                      
#                 compteur += 1
#                  
#                      
#         with open("groupes_16S_identiques.pickle", "rb") as fichier_identiques_1 :
#             mon_depickler_1 = pickle.Unpickler(fichier_identiques_1)
#             groupes_identiques_1 = mon_depickler_1.load()
#               
#               
#             for elt in homologues_16s :
#                 compteur1 = 1
#                 for groupes in groupes_identiques_1 :
#                     compteur2 = 1
#                     for groupe in groupes :
#                         if elt in groupe :
#                             print("groupe : %s %s"%(compteur1, compteur2))
#                             print(elt)
#                             print("in")
#                         compteur2 += 1
#                     compteur1 += 1
#          
#         with open("groupes_25S_homologues.pickle", "rb") as fichier_homologues_1 :
#             mon_depickler_1 = pickle.Unpickler(fichier_homologues_1)
#             groupes_homologues_1 = mon_depickler_1.load()
#              
#             homologues_25s = []
#             for elt in liste_100_pourcent :
#                 cle = (elt.split("_")[1], int(elt.split("_")[2][:len(elt.split("_")[2])-7]))
#                 compteur = 1
#                 for groupe in groupes_homologues_1 :
#     #                 print(groupe)
#     #                 print(cle)
#                     if cle in groupe :
#                         print("groupe : %d"%compteur)
#                         print("in")
#                         homologues_25s.append(cle)
#                          
#                     compteur += 1
#                      
#         with open("groupes_25S_identiques.pickle", "rb") as fichier_identiques_1 :
#             mon_depickler_1 = pickle.Unpickler(fichier_identiques_1)
#             groupes_identiques_1 = mon_depickler_1.load()
#               
#               
#             for elt in homologues_25s :
#                 compteur1 = 1
#                 for groupes in groupes_identiques_1 :
#                     compteur2 = 1
#                     for groupe in groupes :
#                         if elt in groupe :
#                             print("groupe : %s %s"%(compteur1, compteur2))
#                             print(elt)
#                             print("in")
#                         compteur2 += 1
#                     compteur1 += 1
#          
#         print("homologues 23s")
#         print(len(homologues_23s))
#         print(homologues_23s)
#         print("homologues 16s")
#         print(len(homologues_16s))
#         print(homologues_16s)
#         print("homologues 25s")
#         print(len(homologues_25s))
#         print(homologues_25s)
#         
#         for elt in liste_100_pourcent :
#             cle = (elt.split("_")[1], int(elt.split("_")[2][:len(elt.split("_")[2])-7]))
#             if cle not in homologues_16s and cle not in homologues_23s and cle not in homologues_25s :
#                 print(cle)
    #print(distrib)
#     liste_nb = []
#     for pourcentage in np.arange(1.0, -0.1, -0.1) :
#         nb = 0
#         for elt in distrib :
#             if elt > pourcentage :
#                 nb += 1
#         liste_nb.append(nb)
      
    #ax = plt.gca()
    #print([round(i, 1) for i in np.arange(1.0, -0.1, -0.1)])
    #ax.set_xticks(np.arange(0,11,1))
    #ax.set_xticklabels([round(i, 1) for i in np.arange(1.0, -0.1, -0.1)]) 
#     ax.set_yticklabels(np.arange(5, 94, 5))
#     ax.set_yticks(np.arange(5, 94, 5))             
    #plt.plot(liste_nb)
    #plt.title("Similarité de structure avec le graphe commun du groupe vert-jaune \n sur les extensions")
    #ax.set_xlabel("Pourcentage de ressemblance")
    #ax.set_ylabel("Nombre d'éléments")
    #plt.show()
    #plt.savefig(NEW_EXTENSION_PATH_TAILLE+"Traitement_resultats/similarite_structure_groupe_vert_jaune.png")
    #print(exact)

#     liste_vert_jaune_new_data = ["fichier_1fjg_6.pickle", "fichier_4v9f_20.pickle", "fichier_5fdu_9.pickle", "fichier_5j7l_29.pickle", "fichier_5j5b_22.pickle"]
#     liste_frequence = recherche_sequence_graphe_commun_groupe_11(liste_vert_jaune_new_data)
#     print(liste_frequence)
#     nb_ok = [0]*11
#     print(nb_ok)
#     liste_ok, liste_100 = pourcentage_signature_seq_groupe_de_28(liste_frequence)
#     print(liste_ok)
#     compteur = 0
#     for i in np.arange(1.0, -0.1,-0.1) :
#         for elt in liste_ok :
#             if elt >= i :
#                 nb_ok[compteur] += 1
#         compteur += 1
#     print(nb_ok)
#     
#     print(liste_100)
#     print(len(liste_100))

''' Groupe des 100% vert jaune et 28 '''
#     liste_100_vert_jaune = ['fichier_5j8a_30.pickle', 'fichier_4ybb_30.pickle', 'fichier_4v84_33.pickle', 'fichier_5gaf_4.pickle', 'fichier_4wqu_11.pickle', 'fichier_5v8i_39.pickle', 'fichier_1yjn_3.pickle', 'fichier_4v5p_32.pickle', 'fichier_4v4q_36.pickle', 'fichier_4v6d_46.pickle', 'fichier_4v8o_6.pickle', 'fichier_4v8g_8.pickle', 'fichier_6h58_29.pickle', 'fichier_5j8b_27.pickle', 'fichier_4v52_34.pickle', 'fichier_5h5u_21.pickle', 'fichier_5uym_34.pickle', 'fichier_4l71_19.pickle', 'fichier_4v7j_31.pickle', 'fichier_4u20_25.pickle', 'fichier_4v9s_51.pickle', 'fichier_4wqy_12.pickle', 'fichier_4v7y_24.pickle', 'fichier_4lnt_59.pickle', 'fichier_3jcd_11.pickle', 'fichier_4v8b_32.pickle', 'fichier_4v5n_26.pickle', 'fichier_4wq1_9.pickle', 'fichier_1yjw_4.pickle', 'fichier_4v5m_23.pickle', 'fichier_5hau_23.pickle', 'fichier_1vy7_11.pickle', 'fichier_4v6e_43.pickle', 'fichier_5fdv_12.pickle', 'fichier_6of1_3.pickle', 'fichier_3ccs_5.pickle', 'fichier_6bok_47.pickle', 'fichier_4v5e_3.pickle', 'fichier_5ib7_55.pickle', 'fichier_4wt8_52.pickle', 'fichier_5j30_47.pickle', 'fichier_4tua_19.pickle', 'fichier_5ibb_51.pickle', 'fichier_5j88_27.pickle', 'fichier_4v7z_24.pickle', 'fichier_6c4i_16.pickle', 'fichier_4zsn_12.pickle', 'fichier_5mdv_25.pickle', 'fichier_4v4h_35.pickle', 'fichier_4v5l_6.pickle', 'fichier_1yi2_5.pickle', 'fichier_4v5l_31.pickle', 'fichier_4v7y_23.pickle', 'fichier_6hrm_14.pickle', 'fichier_4u26_42.pickle', 'fichier_1vvj_13.pickle', 'fichier_4v8a_32.pickle', 'fichier_4w29_23.pickle', 'fichier_4v5f_51.pickle', 'fichier_5czp_18.pickle', 'fichier_4v5q_33.pickle', 'fichier_4p6f_47.pickle', 'fichier_4w4g_20.pickle', 'fichier_1yit_5.pickle', 'fichier_4wsd_45.pickle', 'fichier_4v6g_39.pickle', 'fichier_6nd6_4.pickle', 'fichier_4v8g_53.pickle', 'fichier_4wzo_8.pickle', 'fichier_4wsd_49.pickle', 'fichier_5el6_51.pickle', 'fichier_4v5s_64.pickle', 'fichier_3ccv_4.pickle', 'fichier_6i7v_22.pickle', 'fichier_1vy4_7.pickle', 'fichier_4v55_29.pickle', 'fichier_5j5b_32.pickle', 'fichier_4v9a_22.pickle', 'fichier_4v8d_40.pickle', 'fichier_3cc7_5.pickle', 'fichier_4w2g_30.pickle', 'fichier_4v7w_22.pickle', 'fichier_4v6c_45.pickle', 'fichier_4v5d_14.pickle', 'fichier_5wit_30.pickle', 'fichier_1yj9_5.pickle', 'fichier_4wt8_1.pickle', 'fichier_4v5d_27.pickle', 'fichier_4v8i_9.pickle', 'fichier_5el4_12.pickle', 'fichier_5hcp_10.pickle', 'fichier_3ccq_5.pickle', 'fichier_4v6p_3.pickle', 'fichier_4wpo_55.pickle', 'fichier_4v9j_6.pickle', 'fichier_4u25_19.pickle', 'fichier_6c5l_11.pickle', 'fichier_4v8u_55.pickle', 'fichier_4v5f_56.pickle', 'fichier_4wqu_34.pickle', 'fichier_4p70_14.pickle', 'fichier_6gsj_12.pickle', 'fichier_4v54_32.pickle', 'fichier_4v5c_11.pickle', 'fichier_5wis_7.pickle', 'fichier_4wqy_32.pickle', 'fichier_6bok_21.pickle', 'fichier_6q95_2.pickle', 'fichier_4p70_21.pickle', 'fichier_3jcn_15.pickle', 'fichier_4v9m_8.pickle', 'fichier_5ib8_46.pickle', 'fichier_4w29_6.pickle', 'fichier_4xej_29.pickle', 'fichier_5vpp_38.pickle', 'fichier_3j7z_13.pickle', 'fichier_6bz7_44.pickle', 'fichier_5dfe_20.pickle', 'fichier_4v6n_11.pickle', 'fichier_4v6e_25.pickle', 'fichier_5e7k_10.pickle', 'fichier_4v7l_9.pickle', 'fichier_4v8n_13.pickle', 'fichier_4ypb_19.pickle', 'fichier_4v51_26.pickle', 'fichier_4ypb_48.pickle', 'fichier_4v9h_32.pickle', 'fichier_4lsk_23.pickle', 'fichier_1vy4_52.pickle', 'fichier_5doy_28.pickle', 'fichier_3bbx_5.pickle', 'fichier_4wsm_42.pickle', 'fichier_4v9q_20.pickle', 'fichier_4wzd_13.pickle', 'fichier_4v90_28.pickle', 'fichier_4v5r_10.pickle', 'fichier_4v9q_12.pickle', 'fichier_4wsd_38.pickle', 'fichier_5mgp_14.pickle', 'fichier_4wu1_13.pickle', 'fichier_4y4o_12.pickle', 'fichier_4w2f_7.pickle', 'fichier_4v8a_49.pickle', 'fichier_5ib7_48.pickle', 'fichier_4y4o_27.pickle', 'fichier_5gae_14.pickle', 'fichier_4v5s_38.pickle', 'fichier_6czr_9.pickle', 'fichier_4lel_19.pickle', 'fichier_4p6f_19.pickle', 'fichier_3j9z_18.pickle', 'fichier_4zer_31.pickle', 'fichier_1vvj_19.pickle', 'fichier_4v67_49.pickle', 'fichier_5el4_43.pickle', 'fichier_4wqu_54.pickle', 'fichier_5ib7_16.pickle', 'fichier_4tud_47.pickle', 'fichier_5gah_11.pickle', 'fichier_3ccl_4.pickle', 'fichier_5iqr_29.pickle', 'fichier_4v5y_33.pickle', 'fichier_4v6q_8.pickle', 'fichier_4u1v_23.pickle', 'fichier_4wqr_48.pickle', 'fichier_4v9n_48.pickle', 'fichier_4w2g_11.pickle', 'fichier_5e81_49.pickle', 'fichier_6bz8_23.pickle', 'fichier_5vp2_22.pickle', 'fichier_4wq1_50.pickle', 'fichier_6gsj_50.pickle', 'fichier_5hau_33.pickle', 'fichier_4v8e_39.pickle', 'fichier_4v9c_25.pickle', 'fichier_6bok_48.pickle', 'fichier_4v95_55.pickle', 'fichier_4v8i_52.pickle', 'fichier_4v7m_55.pickle', 'fichier_4v84_8.pickle', 'fichier_3cma_6.pickle', 'fichier_6buw_19.pickle', 'fichier_4v9b_37.pickle', 'fichier_4wt1_36.pickle', 'fichier_4v6o_5.pickle', 'fichier_6c5l_34.pickle', 'fichier_5j4b_25.pickle', 'fichier_6buw_46.pickle', 'fichier_5hau_8.pickle', 'fichier_4wzd_43.pickle', 'fichier_4v8c_31.pickle', 'fichier_4tub_17.pickle', 'fichier_5hau_37.pickle', 'fichier_6cfl_10.pickle', 'fichier_4v9n_12.pickle', 'fichier_4v97_7.pickle', 'fichier_4wf1_23.pickle', 'fichier_4v5g_11.pickle', 'fichier_6q95_3.pickle', 'fichier_4v9j_35.pickle', 'fichier_4lfz_42.pickle', 'fichier_1yhq_5.pickle', 'fichier_4v4h_36.pickle', 'fichier_5hcq_24.pickle', 'fichier_5el7_13.pickle', 'fichier_4l47_14.pickle', 'fichier_6o97_30.pickle', 'fichier_4z8c_10.pickle', 'fichier_5kcr_12.pickle', 'fichier_3ccm_5.pickle', 'fichier_4v7b_27.pickle', 'fichier_5j3c_48.pickle', 'fichier_4v7j_20.pickle', 'fichier_4v5f_34.pickle', 'fichier_4v8h_22.pickle', 'fichier_4v83_35.pickle', 'fichier_4u26_23.pickle', 'fichier_4v7k_26.pickle', 'fichier_4z8c_39.pickle', 'fichier_5fdv_40.pickle', 'fichier_4v64_35.pickle', 'fichier_5mdw_25.pickle', 'fichier_5j30_18.pickle', 'fichier_4wpo_64.pickle', 'fichier_6boh_40.pickle', 'fichier_4lnt_14.pickle', 'fichier_5w4k_29.pickle', 'fichier_4v55_30.pickle', 'fichier_5hcq_34.pickle', 'fichier_6q9a_23.pickle', 'fichier_4v8e_19.pickle', 'fichier_4wqy_53.pickle', 'fichier_4v6g_25.pickle', 'fichier_4v8h_31.pickle', 'fichier_6cfj_3.pickle', 'fichier_4v8q_8.pickle', 'fichier_5vpo_15.pickle', 'fichier_6gsj_44.pickle', 'fichier_4ypb_12.pickle', 'fichier_4v72_9.pickle', 'fichier_4y4p_33.pickle', 'fichier_3jbu_15.pickle', 'fichier_4v9j_7.pickle', 'fichier_4lt8_14.pickle', 'fichier_5e81_41.pickle', 'fichier_4v8j_11.pickle', 'fichier_4v6r_6.pickle', 'fichier_3cce_3.pickle', 'fichier_4wro_7.pickle', 'fichier_4v9j_32.pickle', 'fichier_4w2i_6.pickle', 'fichier_5el5_46.pickle', 'fichier_6enu_17.pickle', 'fichier_4v5c_53.pickle', 'fichier_4v9r_10.pickle', 'fichier_5hcq_39.pickle', 'fichier_4v7t_43.pickle', 'fichier_4v7u_45.pickle', 'fichier_4wf1_44.pickle', 'fichier_3j9y_13.pickle', 'fichier_5j8a_31.pickle', 'fichier_4v9s_10.pickle', 'fichier_5fdv_27.pickle', 'fichier_4v7k_43.pickle', 'fichier_4v97_52.pickle', 'fichier_6fkr_16.pickle', 'fichier_4v57_34.pickle', 'fichier_5nrg_7.pickle', 'fichier_4w2h_31.pickle', 'fichier_4w2f_48.pickle', 'fichier_1vy6_11.pickle', 'fichier_5kpx_13.pickle', 'fichier_5ib8_12.pickle', 'fichier_4v5k_2.pickle', 'fichier_4v7x_20.pickle', 'fichier_5fdu_40.pickle', 'fichier_4wsm_8.pickle', 'fichier_6h4n_13.pickle', 'fichier_4v7p_12.pickle', 'fichier_4v7k_14.pickle', 'fichier_4v67_1.pickle', 'fichier_5j7l_30.pickle', 'fichier_4wzo_43.pickle', 'fichier_4v8u_11.pickle', 'fichier_4v6f_39.pickle', 'fichier_4u24_44.pickle', 'fichier_4lt8_22.pickle', 'fichier_4v9n_45.pickle', 'fichier_6cae_3.pickle', 'fichier_4wro_48.pickle', 'fichier_5jc9_30.pickle', 'fichier_5el5_13.pickle', 'fichier_5aka_15.pickle', 'fichier_6fkr_26.pickle', 'fichier_4w2i_49.pickle', 'fichier_4zer_37.pickle', 'fichier_4v7w_21.pickle', 'fichier_5hcr_10.pickle', 'fichier_4wro_37.pickle', 'fichier_4v9q_15.pickle', 'fichier_5el6_44.pickle', 'fichier_4z3s_5.pickle', 'fichier_4v7l_26.pickle', 'fichier_4wsm_33.pickle', 'fichier_4v9r_30.pickle', 'fichier_4tue_19.pickle', 'fichier_4v68_7.pickle', 'fichier_5j88_28.pickle', 'fichier_4u1u_43.pickle', 'fichier_4v95_11.pickle', 'fichier_4v5k_18.pickle', 'fichier_5ibb_14.pickle', 'fichier_4wro_44.pickle', 'fichier_6b4v_37.pickle', 'fichier_4woi_32.pickle', 'fichier_3ccj_3.pickle', 'fichier_4v67_7.pickle', 'fichier_4v7p_27.pickle', 'fichier_4v8g_9.pickle', 'fichier_4v8u_61.pickle', 'fichier_4v8u_36.pickle', 'fichier_4u1v_43.pickle', 'fichier_4u25_40.pickle', 'fichier_4lnt_23.pickle', 'fichier_4v8f_38.pickle', 'fichier_3cme_6.pickle', 'fichier_6n9e_8.pickle', 'fichier_1yij_5.pickle', 'fichier_4v7v_42.pickle', 'fichier_4p6f_17.pickle', 'fichier_6gsl_41.pickle', 'fichier_4v5q_10.pickle', 'fichier_5hcq_9.pickle', 'fichier_6h58_58.pickle', 'fichier_4v9s_11.pickle', 'fichier_5j4c_29.pickle', 'fichier_4tud_19.pickle', 'fichier_4wpo_2.pickle', 'fichier_4v9s_27.pickle', 'fichier_4v5e_44.pickle', 'fichier_4wqr_38.pickle', 'fichier_5el6_55.pickle', 'fichier_3cc2_6.pickle', 'fichier_3ccr_5.pickle', 'fichier_4v6f_2.pickle', 'fichier_5f8k_20.pickle', 'fichier_4v8c_33.pickle', 'fichier_4v5q_51.pickle', 'fichier_5v8i_51.pickle', 'fichier_4wq1_41.pickle', 'fichier_4z8c_26.pickle', 'fichier_4v5n_5.pickle', 'fichier_4v52_35.pickle', 'fichier_4v7z_17.pickle', 'fichier_5hcp_24.pickle', 'fichier_5el4_50.pickle', 'fichier_4w2h_13.pickle', 'fichier_3g6e_5.pickle', 'fichier_4wt8_53.pickle', 'fichier_4v8d_35.pickle', 'fichier_4wpo_28.pickle', 'fichier_4v9f_20.pickle', 'fichier_5czp_46.pickle', 'fichier_6c5l_46.pickle', 'fichier_4wfa_8.pickle', 'fichier_4v4q_35.pickle', 'fichier_4v85_27.pickle', 'fichier_4v7m_29.pickle', 'fichier_4v7m_11.pickle', 'fichier_5j4b_10.pickle', 'fichier_4v50_32.pickle', 'fichier_4v6a_26.pickle', 'fichier_4tuc_18.pickle', 'fichier_4v9i_36.pickle', 'fichier_4v9b_22.pickle', 'fichier_4ybb_31.pickle', 'fichier_4v9b_46.pickle', 'fichier_1vy6_29.pickle', 'fichier_4z3s_35.pickle', 'fichier_6nd5_32.pickle', 'fichier_6nd5_5.pickle', 'fichier_5j4d_40.pickle', 'fichier_4zer_8.pickle', 'fichier_4v95_10.pickle', 'fichier_4yzv_22.pickle', 'fichier_5hd1_10.pickle', 'fichier_4v9m_31.pickle', 'fichier_6gsl_49.pickle', 'fichier_4v5d_20.pickle', 'fichier_6cfk_22.pickle', 'fichier_5czp_44.pickle', 'fichier_6enj_9.pickle', 'fichier_5j8b_3.pickle', 'fichier_4v8i_30.pickle', 'fichier_4v9l_6.pickle', 'fichier_4lnt_24.pickle', 'fichier_4w2h_54.pickle', 'fichier_4v9k_8.pickle', 'fichier_4v5m_5.pickle', 'fichier_4v57_33.pickle', 'fichier_5jc9_29.pickle', 'fichier_4wf9_6.pickle', 'fichier_4v8g_31.pickle', 'fichier_4v8e_36.pickle', 'fichier_2j28_14.pickle', 'fichier_4wu1_51.pickle', 'fichier_4v8d_23.pickle', 'fichier_3g4s_4.pickle', 'fichier_3i56_4.pickle', 'fichier_4u27_43.pickle', 'fichier_5ndk_11.pickle', 'fichier_4v9n_11.pickle', 'fichier_4v8n_9.pickle', 'fichier_4v5j_36.pickle', 'fichier_5el5_39.pickle', 'fichier_4v9r_54.pickle', 'fichier_6o97_7.pickle', 'fichier_4v9h_6.pickle', 'fichier_4l47_21.pickle', 'fichier_6i0y_13.pickle', 'fichier_4wr6_36.pickle', 'fichier_4zsn_48.pickle', 'fichier_4v7v_24.pickle', 'fichier_4v8i_10.pickle', 'fichier_5e7k_41.pickle', 'fichier_5it8_27.pickle', 'fichier_5j3c_50.pickle', 'fichier_5j91_32.pickle', 'fichier_4v9i_11.pickle', 'fichier_4zsn_19.pickle', 'fichier_4v9b_45.pickle', 'fichier_4v7s_45.pickle', 'fichier_4w4g_12.pickle', 'fichier_4z3s_7.pickle', 'fichier_5wit_5.pickle', 'fichier_6of1_32.pickle', 'fichier_6n9f_8.pickle', 'fichier_4v54_31.pickle', 'fichier_1vy5_47.pickle', 'fichier_4w2e_27.pickle', 'fichier_5lzd_15.pickle', 'fichier_6cfl_22.pickle', 'fichier_6cae_5.pickle', 'fichier_4v9k_30.pickle', 'fichier_4v97_8.pickle', 'fichier_5it8_28.pickle', 'fichier_3i55_6.pickle', 'fichier_5j4c_3.pickle', 'fichier_4uy8_14.pickle', 'fichier_1vy7_51.pickle', 'fichier_4v9q_13.pickle', 'fichier_4v7k_32.pickle', 'fichier_5j7l_29.pickle', 'fichier_5hcr_22.pickle', 'fichier_4v8f_35.pickle', 'fichier_4v6t_26.pickle', 'fichier_4v8n_56.pickle', 'fichier_5fdu_27.pickle', 'fichier_6nd6_28.pickle', 'fichier_5ndk_44.pickle', 'fichier_5dox_31.pickle', 'fichier_4v53_33.pickle', 'fichier_4v53_32.pickle', 'fichier_5el7_48.pickle', 'fichier_3j5l_12.pickle', 'fichier_4v6s_11.pickle', 'fichier_5dox_18.pickle', 'fichier_4v5f_9.pickle', 'fichier_6cfk_9.pickle', 'fichier_5ndk_53.pickle', 'fichier_6nd6_5.pickle', 'fichier_6gc6_3.pickle', 'fichier_4v5j_11.pickle', 'fichier_6cae_34.pickle', 'fichier_4wsd_8.pickle', 'fichier_1vy7_31.pickle', 'fichier_6boh_46.pickle', 'fichier_4v9d_31.pickle', 'fichier_4wqf_54.pickle', 'fichier_5o2r_15.pickle', 'fichier_4v95_32.pickle', 'fichier_4lsk_22.pickle', 'fichier_5gad_11.pickle', 'fichier_5uyp_25.pickle', 'fichier_4v5y_32.pickle', 'fichier_4v5g_31.pickle', 'fichier_4u20_44.pickle', 'fichier_4v5r_16.pickle', 'fichier_1vy5_6.pickle', 'fichier_4v5d_56.pickle', 'fichier_6gc0_6.pickle', 'fichier_4v8n_30.pickle', 'fichier_5el6_8.pickle', 'fichier_3ccu_4.pickle', 'fichier_5j4d_46.pickle', 'fichier_4v90_7.pickle', 'fichier_4v87_41.pickle', 'fichier_4tua_48.pickle', 'fichier_3cd6_3.pickle', 'fichier_4v5r_59.pickle', 'fichier_4v7z_15.pickle', 'fichier_6cfj_29.pickle', 'fichier_1vy4_2.pickle', 'fichier_5nco_11.pickle', 'fichier_4tue_47.pickle', 'fichier_6bok_31.pickle', 'fichier_5ndj_41.pickle', 'fichier_4tua_51.pickle', 'fichier_4v97_30.pickle', 'fichier_4wqf_33.pickle', 'fichier_4v87_38.pickle', 'fichier_5vp2_8.pickle', 'fichier_4wr6_8.pickle', 'fichier_4v5q_56.pickle', 'fichier_5ndj_51.pickle', 'fichier_4v8b_38.pickle', 'fichier_5kpw_11.pickle', 'fichier_5w4k_3.pickle', 'fichier_4lel_20.pickle', 'fichier_4v67_22.pickle', 'fichier_4wt1_9.pickle', 'fichier_5hl7_8.pickle', 'fichier_4w2h_12.pickle', 'fichier_4v7m_12.pickle', 'fichier_4v8h_46.pickle', 'fichier_4v64_34.pickle', 'fichier_5dfe_47.pickle', 'fichier_4wra_20.pickle', 'fichier_4u24_24.pickle', 'fichier_5ibb_45.pickle', 'fichier_3g71_3.pickle', 'fichier_4lsk_14.pickle', 'fichier_4wt1_46.pickle', 'fichier_4wr6_45.pickle', 'fichier_5f8k_9.pickle', 'fichier_4v7z_16.pickle', 'fichier_6of1_5.pickle', 'fichier_4tua_21.pickle', 'fichier_5j5b_33.pickle', 'fichier_4v5r_34.pickle', 'fichier_5mdz_27.pickle', 'fichier_5nwy_7.pickle', 'fichier_4p6f_44.pickle', 'fichier_5czp_20.pickle', 'fichier_4wqy_11.pickle', 'fichier_4l71_13.pickle', 'fichier_4v50_33.pickle', 'fichier_6gsj_35.pickle', 'fichier_6b4v_10.pickle', 'fichier_6i7v_23.pickle', 'fichier_6fkr_22.pickle', 'fichier_5hd1_24.pickle', 'fichier_5afi_13.pickle', 'fichier_4wqf_10.pickle', 'fichier_4v7x_21.pickle', 'fichier_4v9r_11.pickle', 'fichier_4v8x_23.pickle', 'fichier_4zer_22.pickle', 'fichier_4v8j_38.pickle', 'fichier_5doy_4.pickle', 'fichier_6o97_5.pickle', 'fichier_3jbv_18.pickle', 'fichier_6n9f_20.pickle', 'fichier_5el7_54.pickle', 'fichier_5mdy_23.pickle', 'fichier_4v9d_30.pickle', 'fichier_4v71_10.pickle', 'fichier_5vpo_46.pickle', 'fichier_5vpo_17.pickle', 'fichier_4v83_9.pickle', 'fichier_4v8x_2.pickle', 'fichier_5fdv_37.pickle', 'fichier_4zsn_20.pickle', 'fichier_5j91_33.pickle', 'fichier_4u27_21.pickle', 'fichier_4v9a_49.pickle', 'fichier_5wis_22.pickle', 'fichier_4v51_25.pickle', 'fichier_3cc4_5.pickle', 'fichier_4wu1_44.pickle', 'fichier_4www_28.pickle', 'fichier_4woi_49.pickle', 'fichier_4yzv_13.pickle', 'fichier_4v5p_9.pickle', 'fichier_4w2g_52.pickle', 'fichier_4v5s_59.pickle', 'fichier_4y4p_5.pickle', 'fichier_4v8a_37.pickle', 'fichier_4v6f_50.pickle', 'fichier_4v5s_12.pickle', 'fichier_4v7u_25.pickle', 'fichier_4v9a_41.pickle', 'fichier_5fdu_9.pickle', 'fichier_6gsk_35.pickle', 'fichier_1vy6_51.pickle', 'fichier_4lfz_13.pickle']
#     liste_100_28 = ['fichier_5j8a_30.pickle', 'fichier_4ybb_30.pickle', 'fichier_4v84_33.pickle', 'fichier_5gaf_4.pickle', 'fichier_4wqu_11.pickle', 'fichier_5v8i_39.pickle', 'fichier_4v5p_32.pickle', 'fichier_4v4q_36.pickle', 'fichier_4v6d_46.pickle', 'fichier_4v8o_6.pickle', 'fichier_4v8g_8.pickle', 'fichier_6h58_29.pickle', 'fichier_5j8b_27.pickle', 'fichier_5dc3_19.pickle', 'fichier_4v52_34.pickle', 'fichier_5dge_27.pickle', 'fichier_5h5u_21.pickle', 'fichier_5uym_34.pickle', 'fichier_4l71_19.pickle', 'fichier_4v7j_31.pickle', 'fichier_4u20_25.pickle', 'fichier_4v9s_51.pickle', 'fichier_4wqy_12.pickle', 'fichier_4v7y_24.pickle', 'fichier_4lnt_59.pickle', 'fichier_6n8j_8.pickle', 'fichier_4u52_24.pickle', 'fichier_5wf0_12.pickle', 'fichier_3jcd_11.pickle', 'fichier_5dgf_15.pickle', 'fichier_4v8b_32.pickle', 'fichier_4v5n_26.pickle', 'fichier_4v8y_2.pickle', 'fichier_4wq1_9.pickle', 'fichier_4v5m_23.pickle', 'fichier_5hau_23.pickle', 'fichier_1vy7_11.pickle', 'fichier_4v6e_43.pickle', 'fichier_5fdv_12.pickle', 'fichier_6of1_3.pickle', 'fichier_4v88_7.pickle', 'fichier_6bok_47.pickle', 'fichier_4v5e_3.pickle', 'fichier_5ib7_55.pickle', 'fichier_4u3m_16.pickle', 'fichier_4wt8_52.pickle', 'fichier_6hhq_26.pickle', 'fichier_5tgm_5.pickle', 'fichier_5fcj_27.pickle', 'fichier_5j30_47.pickle', 'fichier_4tua_19.pickle', 'fichier_5ibb_51.pickle', 'fichier_5j88_27.pickle', 'fichier_4v7z_24.pickle', 'fichier_6c4i_16.pickle', 'fichier_4zsn_12.pickle', 'fichier_5mdv_25.pickle', 'fichier_4v4h_35.pickle', 'fichier_5on6_27.pickle', 'fichier_4v5l_6.pickle', 'fichier_4v5l_31.pickle', 'fichier_4v7y_23.pickle', 'fichier_6hrm_14.pickle', 'fichier_4u26_42.pickle', 'fichier_1vvj_13.pickle', 'fichier_4v8a_32.pickle', 'fichier_4w29_23.pickle', 'fichier_4v5f_51.pickle', 'fichier_5czp_18.pickle', 'fichier_4v5q_33.pickle', 'fichier_4p6f_47.pickle', 'fichier_4w4g_20.pickle', 'fichier_4u4q_15.pickle', 'fichier_4wsd_45.pickle', 'fichier_4v6g_39.pickle', 'fichier_6nd6_4.pickle', 'fichier_4v8g_53.pickle', 'fichier_4wzo_8.pickle', 'fichier_4wsd_49.pickle', 'fichier_5el6_51.pickle', 'fichier_4v5s_64.pickle', 'fichier_6i7v_22.pickle', 'fichier_6n8n_3.pickle', 'fichier_1vy4_7.pickle', 'fichier_6hhq_20.pickle', 'fichier_4v55_29.pickle', 'fichier_5j5b_32.pickle', 'fichier_4v9a_22.pickle', 'fichier_4v8d_40.pickle', 'fichier_4w2g_30.pickle', 'fichier_4v7w_22.pickle', 'fichier_4u56_15.pickle', 'fichier_4v6c_45.pickle', 'fichier_4v5d_14.pickle', 'fichier_5wit_30.pickle', 'fichier_5fci_17.pickle', 'fichier_4wt8_1.pickle', 'fichier_4v5d_27.pickle', 'fichier_4v8i_9.pickle', 'fichier_5el4_12.pickle', 'fichier_5wdt_14.pickle', 'fichier_5hcp_10.pickle', 'fichier_4v6p_3.pickle', 'fichier_4wpo_55.pickle', 'fichier_4v9j_6.pickle', 'fichier_4v6i_2.pickle', 'fichier_4u25_19.pickle', 'fichier_6c5l_11.pickle', 'fichier_4v8u_55.pickle', 'fichier_4v5f_56.pickle', 'fichier_4wqu_34.pickle', 'fichier_4p70_14.pickle', 'fichier_4u3m_25.pickle', 'fichier_4v8z_2.pickle', 'fichier_6gsj_12.pickle', 'fichier_4v54_32.pickle', 'fichier_4v5c_11.pickle', 'fichier_5wis_7.pickle', 'fichier_4wqy_32.pickle', 'fichier_5dgv_14.pickle', 'fichier_6bok_21.pickle', 'fichier_6q95_2.pickle', 'fichier_4p70_21.pickle', 'fichier_4u3u_31.pickle', 'fichier_4u4o_13.pickle', 'fichier_3jcn_15.pickle', 'fichier_4v9m_8.pickle', 'fichier_5ib8_46.pickle', 'fichier_4w29_6.pickle', 'fichier_4xej_29.pickle', 'fichier_5vpp_38.pickle', 'fichier_3j7z_13.pickle', 'fichier_6bz7_44.pickle', 'fichier_5dfe_20.pickle', 'fichier_4v6n_11.pickle', 'fichier_4u4y_14.pickle', 'fichier_4v6e_25.pickle', 'fichier_5e7k_10.pickle', 'fichier_4v7l_9.pickle', 'fichier_4v8n_13.pickle', 'fichier_4ypb_19.pickle', 'fichier_4v51_26.pickle', 'fichier_4v4n_10.pickle', 'fichier_4u4z_16.pickle', 'fichier_4u53_20.pickle', 'fichier_4ypb_48.pickle', 'fichier_4v9h_32.pickle', 'fichier_4lsk_23.pickle', 'fichier_1vy4_52.pickle', 'fichier_5doy_28.pickle', 'fichier_4u51_14.pickle', 'fichier_5juu_9.pickle', 'fichier_3bbx_5.pickle', 'fichier_4wsm_42.pickle', 'fichier_4v9q_20.pickle', 'fichier_4wzd_13.pickle', 'fichier_5fcj_19.pickle', 'fichier_4u52_16.pickle', 'fichier_4v90_28.pickle', 'fichier_5jup_13.pickle', 'fichier_4v5r_10.pickle', 'fichier_4v9q_12.pickle', 'fichier_4wsd_38.pickle', 'fichier_6n8j_9.pickle', 'fichier_5mgp_14.pickle', 'fichier_4wu1_13.pickle', 'fichier_4y4o_12.pickle', 'fichier_4w2f_7.pickle', 'fichier_4u4n_19.pickle', 'fichier_6n8m_3.pickle', 'fichier_4v8a_49.pickle', 'fichier_5ib7_48.pickle', 'fichier_4y4o_27.pickle', 'fichier_5gae_14.pickle', 'fichier_4v5s_38.pickle', 'fichier_6czr_9.pickle', 'fichier_4lel_19.pickle', 'fichier_4p6f_19.pickle', 'fichier_5wfs_11.pickle', 'fichier_3j9z_18.pickle', 'fichier_4zer_31.pickle', 'fichier_1vvj_19.pickle', 'fichier_4v67_49.pickle', 'fichier_5fci_25.pickle', 'fichier_5el4_43.pickle', 'fichier_4wqu_54.pickle', 'fichier_5ib7_16.pickle', 'fichier_4tud_47.pickle', 'fichier_5gah_11.pickle', 'fichier_5iqr_29.pickle', 'fichier_4v5y_33.pickle', 'fichier_4v6q_8.pickle', 'fichier_4u1v_23.pickle', 'fichier_4wqr_48.pickle', 'fichier_4v9n_48.pickle', 'fichier_4w2g_11.pickle', 'fichier_5e81_49.pickle', 'fichier_6bz8_23.pickle', 'fichier_5vp2_22.pickle', 'fichier_4wq1_50.pickle', 'fichier_6gsj_50.pickle', 'fichier_5hau_33.pickle', 'fichier_4v8e_39.pickle', 'fichier_4v9c_25.pickle', 'fichier_6bok_48.pickle', 'fichier_4v95_55.pickle', 'fichier_4v8i_52.pickle', 'fichier_4v7m_55.pickle', 'fichier_4v84_8.pickle', 'fichier_6buw_19.pickle', 'fichier_4v9b_37.pickle', 'fichier_4wt1_36.pickle', 'fichier_4v6o_5.pickle', 'fichier_6c5l_34.pickle', 'fichier_5j4b_25.pickle', 'fichier_6buw_46.pickle', 'fichier_5hau_8.pickle', 'fichier_4wzd_43.pickle', 'fichier_4v8c_31.pickle', 'fichier_4tub_17.pickle', 'fichier_5hau_37.pickle', 'fichier_6cfl_10.pickle', 'fichier_4v9n_12.pickle', 'fichier_4v97_7.pickle', 'fichier_4wf1_23.pickle', 'fichier_4u4u_27.pickle', 'fichier_4v5g_11.pickle', 'fichier_6q95_3.pickle', 'fichier_5tga_21.pickle', 'fichier_4v9j_35.pickle', 'fichier_4lfz_42.pickle', 'fichier_4v4h_36.pickle', 'fichier_5hcq_24.pickle', 'fichier_5el7_13.pickle', 'fichier_4l47_14.pickle', 'fichier_6o97_30.pickle', 'fichier_4z8c_10.pickle', 'fichier_5kcr_12.pickle', 'fichier_4v7b_27.pickle', 'fichier_5j3c_48.pickle', 'fichier_4v7j_20.pickle', 'fichier_4v5f_34.pickle', 'fichier_4v8h_22.pickle', 'fichier_4v83_35.pickle', 'fichier_4u3n_17.pickle', 'fichier_4u26_23.pickle', 'fichier_4v7k_26.pickle', 'fichier_4z8c_39.pickle', 'fichier_5dge_17.pickle', 'fichier_5fdv_40.pickle', 'fichier_4v64_35.pickle', 'fichier_5mdw_25.pickle', 'fichier_4u4o_19.pickle', 'fichier_5j30_18.pickle', 'fichier_4wpo_64.pickle', 'fichier_6boh_40.pickle', 'fichier_4lnt_14.pickle', 'fichier_5dc3_11.pickle', 'fichier_5w4k_29.pickle', 'fichier_4v55_30.pickle', 'fichier_5hcq_34.pickle', 'fichier_6q9a_23.pickle', 'fichier_4v8e_19.pickle', 'fichier_4wqy_53.pickle', 'fichier_4v6g_25.pickle', 'fichier_4v8h_31.pickle', 'fichier_6cfj_3.pickle', 'fichier_4v8q_8.pickle', 'fichier_5vpo_15.pickle', 'fichier_6gsj_44.pickle', 'fichier_4u3n_24.pickle', 'fichier_4ypb_12.pickle', 'fichier_4v72_9.pickle', 'fichier_4y4p_33.pickle', 'fichier_3jbu_15.pickle', 'fichier_4v9j_7.pickle', 'fichier_4lt8_14.pickle', 'fichier_5e81_41.pickle', 'fichier_4v8j_11.pickle', 'fichier_4v6r_6.pickle', 'fichier_4wro_7.pickle', 'fichier_4v9j_32.pickle', 'fichier_4w2i_6.pickle', 'fichier_5el5_46.pickle', 'fichier_6enu_17.pickle', 'fichier_4v5c_53.pickle', 'fichier_4v9r_10.pickle', 'fichier_5hcq_39.pickle', 'fichier_4v7t_43.pickle', 'fichier_4v7u_45.pickle', 'fichier_4wf1_44.pickle', 'fichier_3j9y_13.pickle', 'fichier_4u4y_21.pickle', 'fichier_5j8a_31.pickle', 'fichier_4v9s_10.pickle', 'fichier_5fdv_27.pickle', 'fichier_4v7k_43.pickle', 'fichier_4v97_52.pickle', 'fichier_6fkr_16.pickle', 'fichier_6i7o_18.pickle', 'fichier_4v57_34.pickle', 'fichier_5t6r_4.pickle', 'fichier_4w2h_31.pickle', 'fichier_4w2f_48.pickle', 'fichier_4u4r_16.pickle', 'fichier_1vy6_11.pickle', 'fichier_5kpx_13.pickle', 'fichier_5ib8_12.pickle', 'fichier_4u6f_13.pickle', 'fichier_4v5k_2.pickle', 'fichier_4v7x_20.pickle', 'fichier_5fdu_40.pickle', 'fichier_4wsm_8.pickle', 'fichier_6h4n_13.pickle', 'fichier_4v7p_12.pickle', 'fichier_5tbw_5.pickle', 'fichier_4v7k_14.pickle', 'fichier_4v67_1.pickle', 'fichier_5j7l_30.pickle', 'fichier_4wzo_43.pickle', 'fichier_4v8u_11.pickle', 'fichier_4u6f_20.pickle', 'fichier_4v6f_39.pickle', 'fichier_4u55_17.pickle', 'fichier_4u24_44.pickle', 'fichier_4lt8_22.pickle', 'fichier_4v9n_45.pickle', 'fichier_6cae_3.pickle', 'fichier_4wro_48.pickle', 'fichier_5jc9_30.pickle', 'fichier_5lyb_26.pickle', 'fichier_5el5_13.pickle', 'fichier_5aka_15.pickle', 'fichier_6fkr_26.pickle', 'fichier_4w2i_49.pickle', 'fichier_4zer_37.pickle', 'fichier_4v7w_21.pickle', 'fichier_5hcr_10.pickle', 'fichier_4wro_37.pickle', 'fichier_4v9q_15.pickle', 'fichier_5el6_44.pickle', 'fichier_4z3s_5.pickle', 'fichier_4v7l_26.pickle', 'fichier_4wsm_33.pickle', 'fichier_4v9r_30.pickle', 'fichier_4tue_19.pickle', 'fichier_4v68_7.pickle', 'fichier_5j88_28.pickle', 'fichier_4u1u_43.pickle', 'fichier_4v95_11.pickle', 'fichier_4v5k_18.pickle', 'fichier_5ndg_12.pickle', 'fichier_5ibb_14.pickle', 'fichier_4wro_44.pickle', 'fichier_6b4v_37.pickle', 'fichier_4woi_32.pickle', 'fichier_4v67_7.pickle', 'fichier_5mei_24.pickle', 'fichier_5lyb_19.pickle', 'fichier_4v7p_27.pickle', 'fichier_4v8g_9.pickle', 'fichier_4v8u_61.pickle', 'fichier_4v8u_36.pickle', 'fichier_4u1v_43.pickle', 'fichier_4u25_40.pickle', 'fichier_4lnt_23.pickle', 'fichier_4v8f_38.pickle', 'fichier_6n9e_8.pickle', 'fichier_6n8o_4.pickle', 'fichier_4v7v_42.pickle', 'fichier_4p6f_17.pickle', 'fichier_6gsl_41.pickle', 'fichier_4v5q_10.pickle', 'fichier_5hcq_9.pickle', 'fichier_6h58_58.pickle', 'fichier_4v9s_11.pickle', 'fichier_4v7r_4.pickle', 'fichier_5j4c_29.pickle', 'fichier_4tud_19.pickle', 'fichier_4wpo_2.pickle', 'fichier_4v9s_27.pickle', 'fichier_4v5e_44.pickle', 'fichier_4wqr_38.pickle', 'fichier_5ndv_15.pickle', 'fichier_5el6_55.pickle', 'fichier_4v6f_2.pickle', 'fichier_5f8k_20.pickle', 'fichier_4v8c_33.pickle', 'fichier_4v5q_51.pickle', 'fichier_5v8i_51.pickle', 'fichier_4wq1_41.pickle', 'fichier_4z8c_26.pickle', 'fichier_4v5n_5.pickle', 'fichier_4v52_35.pickle', 'fichier_4v7z_17.pickle', 'fichier_5hcp_24.pickle', 'fichier_5el4_50.pickle', 'fichier_4w2h_13.pickle', 'fichier_4wt8_53.pickle', 'fichier_4v8d_35.pickle', 'fichier_4wpo_28.pickle', 'fichier_5czp_46.pickle', 'fichier_6c5l_46.pickle', 'fichier_4v4q_35.pickle', 'fichier_4v85_27.pickle', 'fichier_5i4l_20.pickle', 'fichier_4v7m_29.pickle', 'fichier_4v7m_11.pickle', 'fichier_5j4b_10.pickle', 'fichier_4v50_32.pickle', 'fichier_4v6a_26.pickle', 'fichier_4tuc_18.pickle', 'fichier_4v9i_36.pickle', 'fichier_4v9b_22.pickle', 'fichier_4ybb_31.pickle', 'fichier_4v9b_46.pickle', 'fichier_1vy6_29.pickle', 'fichier_4z3s_35.pickle', 'fichier_5obm_14.pickle', 'fichier_6nd5_32.pickle', 'fichier_6nd5_5.pickle', 'fichier_5j4d_40.pickle', 'fichier_4zer_8.pickle', 'fichier_4v95_10.pickle', 'fichier_4yzv_22.pickle', 'fichier_5hd1_10.pickle', 'fichier_4v9m_31.pickle', 'fichier_6gsl_49.pickle', 'fichier_4v5d_20.pickle', 'fichier_6cfk_22.pickle', 'fichier_5czp_44.pickle', 'fichier_6enj_9.pickle', 'fichier_4v6u_10.pickle', 'fichier_5j8b_3.pickle', 'fichier_4v8i_30.pickle', 'fichier_4v9l_6.pickle', 'fichier_4lnt_24.pickle', 'fichier_4w2h_54.pickle', 'fichier_4v9k_8.pickle', 'fichier_4u4q_25.pickle', 'fichier_4v5m_5.pickle', 'fichier_4v57_33.pickle', 'fichier_5jc9_29.pickle', 'fichier_4u4n_30.pickle', 'fichier_4v8g_31.pickle', 'fichier_4v8e_36.pickle', 'fichier_2j28_14.pickle', 'fichier_5tga_16.pickle', 'fichier_4wu1_51.pickle', 'fichier_4v8d_23.pickle', 'fichier_4u27_43.pickle', 'fichier_5ndk_11.pickle', 'fichier_4v9n_11.pickle', 'fichier_4v8n_9.pickle', 'fichier_4v5j_36.pickle', 'fichier_5dat_15.pickle', 'fichier_4u4z_26.pickle', 'fichier_5el5_39.pickle', 'fichier_4v9r_54.pickle', 'fichier_6o97_7.pickle', 'fichier_4v9h_6.pickle', 'fichier_5we4_15.pickle', 'fichier_4l47_21.pickle', 'fichier_6i0y_13.pickle', 'fichier_4wr6_36.pickle', 'fichier_4zsn_48.pickle', 'fichier_6cb1_3.pickle', 'fichier_4v7v_24.pickle', 'fichier_4v8i_10.pickle', 'fichier_5e7k_41.pickle', 'fichier_4u4r_26.pickle', 'fichier_5it8_27.pickle', 'fichier_4u50_18.pickle', 'fichier_5j3c_50.pickle', 'fichier_5dat_24.pickle', 'fichier_5j91_32.pickle', 'fichier_4v9i_11.pickle', 'fichier_4zsn_19.pickle', 'fichier_4v9b_45.pickle', 'fichier_4v7s_45.pickle', 'fichier_4w4g_12.pickle', 'fichier_4z3s_7.pickle', 'fichier_5wit_5.pickle', 'fichier_6of1_32.pickle', 'fichier_6n9f_8.pickle', 'fichier_4v54_31.pickle', 'fichier_1vy5_47.pickle', 'fichier_5jcs_1.pickle', 'fichier_4w2e_27.pickle', 'fichier_5lzd_15.pickle', 'fichier_5m1j_6.pickle', 'fichier_6cfl_22.pickle', 'fichier_6cae_5.pickle', 'fichier_4v9k_30.pickle', 'fichier_4v97_8.pickle', 'fichier_5it8_28.pickle', 'fichier_5mei_17.pickle', 'fichier_4u4u_18.pickle', 'fichier_4v7r_8.pickle', 'fichier_5j4c_3.pickle', 'fichier_4uy8_14.pickle', 'fichier_1vy7_51.pickle', 'fichier_4v9q_13.pickle', 'fichier_4v7k_32.pickle', 'fichier_5j7l_29.pickle', 'fichier_5hcr_22.pickle', 'fichier_4v8f_35.pickle', 'fichier_4v6t_26.pickle', 'fichier_4v8n_56.pickle', 'fichier_5fdu_27.pickle', 'fichier_6nd6_28.pickle', 'fichier_5ndk_44.pickle', 'fichier_5dox_31.pickle', 'fichier_4v53_33.pickle', 'fichier_4v53_32.pickle', 'fichier_5el7_48.pickle', 'fichier_5tbw_13.pickle', 'fichier_3j5l_12.pickle', 'fichier_4v6s_11.pickle', 'fichier_5dox_18.pickle', 'fichier_4v5f_9.pickle', 'fichier_6cfk_9.pickle', 'fichier_5ndk_53.pickle', 'fichier_6nd6_5.pickle', 'fichier_6gc6_3.pickle', 'fichier_4v5j_11.pickle', 'fichier_5dgv_23.pickle', 'fichier_6cae_34.pickle', 'fichier_4wsd_8.pickle', 'fichier_4u55_25.pickle', 'fichier_4v88_3.pickle', 'fichier_1vy7_31.pickle', 'fichier_6boh_46.pickle', 'fichier_4v9d_31.pickle', 'fichier_4wqf_54.pickle', 'fichier_4u51_23.pickle', 'fichier_5o2r_15.pickle', 'fichier_5apo_8.pickle', 'fichier_4v95_32.pickle', 'fichier_5obm_23.pickle', 'fichier_4lsk_22.pickle', 'fichier_5gad_11.pickle', 'fichier_3jct_8.pickle', 'fichier_5uyp_25.pickle', 'fichier_4v5y_32.pickle', 'fichier_5t62_4.pickle', 'fichier_4v5g_31.pickle', 'fichier_4u20_44.pickle', 'fichier_4v5r_16.pickle', 'fichier_4u53_28.pickle', 'fichier_1vy5_6.pickle', 'fichier_4v5d_56.pickle', 'fichier_6gc0_6.pickle', 'fichier_4v8n_30.pickle', 'fichier_5juo_7.pickle', 'fichier_5el6_8.pickle', 'fichier_5ndw_24.pickle', 'fichier_5j4d_46.pickle', 'fichier_4v90_7.pickle', 'fichier_4u3u_20.pickle', 'fichier_5fl8_1.pickle', 'fichier_4v87_41.pickle', 'fichier_4tua_48.pickle', 'fichier_5we6_13.pickle', 'fichier_4v5r_59.pickle', 'fichier_4v7z_15.pickle', 'fichier_6cfj_29.pickle', 'fichier_1vy4_2.pickle', 'fichier_5ndv_24.pickle', 'fichier_5nco_11.pickle', 'fichier_5on6_22.pickle', 'fichier_4tue_47.pickle', 'fichier_6n8l_9.pickle', 'fichier_6c0f_2.pickle', 'fichier_6bok_31.pickle', 'fichier_5ndj_41.pickle', 'fichier_4tua_51.pickle', 'fichier_4v97_30.pickle', 'fichier_5apn_7.pickle', 'fichier_4wqf_33.pickle', 'fichier_4v87_38.pickle', 'fichier_5vp2_8.pickle', 'fichier_4wr6_8.pickle', 'fichier_4v5q_56.pickle', 'fichier_5ndj_51.pickle', 'fichier_4v8b_38.pickle', 'fichier_5kpw_11.pickle', 'fichier_5w4k_3.pickle', 'fichier_4lel_20.pickle', 'fichier_4v67_22.pickle', 'fichier_4wt1_9.pickle', 'fichier_4w2h_12.pickle', 'fichier_4v7m_12.pickle', 'fichier_4v8h_46.pickle', 'fichier_4v64_34.pickle', 'fichier_5dfe_47.pickle', 'fichier_4wra_20.pickle', 'fichier_4u24_24.pickle', 'fichier_5ibb_45.pickle', 'fichier_4lsk_14.pickle', 'fichier_4wt1_46.pickle', 'fichier_4wr6_45.pickle', 'fichier_5f8k_9.pickle', 'fichier_4v7z_16.pickle', 'fichier_6of1_5.pickle', 'fichier_6n8o_2.pickle', 'fichier_4tua_21.pickle', 'fichier_4u50_26.pickle', 'fichier_5j5b_33.pickle', 'fichier_4v5r_34.pickle', 'fichier_5mdz_27.pickle', 'fichier_5nwy_7.pickle', 'fichier_4p6f_44.pickle', 'fichier_5czp_20.pickle', 'fichier_4wqy_11.pickle', 'fichier_4l71_13.pickle', 'fichier_4v50_33.pickle', 'fichier_4v8t_12.pickle', 'fichier_6gsj_35.pickle', 'fichier_6b4v_10.pickle', 'fichier_6i7v_23.pickle', 'fichier_6fkr_22.pickle', 'fichier_5hd1_24.pickle', 'fichier_5afi_13.pickle', 'fichier_4wqf_10.pickle', 'fichier_4v7x_21.pickle', 'fichier_4v9r_11.pickle', 'fichier_4v8x_23.pickle', 'fichier_4zer_22.pickle', 'fichier_4v8j_38.pickle', 'fichier_5doy_4.pickle', 'fichier_6o97_5.pickle', 'fichier_3jbv_18.pickle', 'fichier_5ndw_16.pickle', 'fichier_6n9f_20.pickle', 'fichier_5el7_54.pickle', 'fichier_5mdy_23.pickle', 'fichier_4v9d_30.pickle', 'fichier_4v71_10.pickle', 'fichier_5vpo_46.pickle', 'fichier_5vpo_17.pickle', 'fichier_4v83_9.pickle', 'fichier_5wfk_13.pickle', 'fichier_4v8x_2.pickle', 'fichier_5fdv_37.pickle', 'fichier_4zsn_20.pickle', 'fichier_5j91_33.pickle', 'fichier_4u27_21.pickle', 'fichier_4v9a_49.pickle', 'fichier_5wis_22.pickle', 'fichier_4v51_25.pickle', 'fichier_4wu1_44.pickle', 'fichier_6elz_7.pickle', 'fichier_4www_28.pickle', 'fichier_4woi_49.pickle', 'fichier_5tgm_8.pickle', 'fichier_4yzv_13.pickle', 'fichier_4v5p_9.pickle', 'fichier_4w2g_52.pickle', 'fichier_4v5s_59.pickle', 'fichier_5i4l_13.pickle', 'fichier_4y4p_5.pickle', 'fichier_4v8a_37.pickle', 'fichier_4v6f_50.pickle', 'fichier_6q8y_3.pickle', 'fichier_4v5s_12.pickle', 'fichier_5dgf_21.pickle', 'fichier_4v7u_25.pickle', 'fichier_4v9a_41.pickle', 'fichier_5fdu_9.pickle', 'fichier_6gsk_35.pickle', 'fichier_1vy6_51.pickle', 'fichier_4lfz_13.pickle']
#     print(len(liste_100_vert_jaune))
#     print(len(liste_100_28))
#     diff1 = []
#     diff2 = []
#     idem = []
#     for elt in liste_100_vert_jaune :
#         if elt in liste_100_28 : 
#             if elt not in idem :
#                 print(elt)
#                 idem.append(elt)
#         else :
#             if elt not in diff1 :
#                 print(elt)
#                 diff1.append(elt)
#     print(idem)
#     print(len(idem))
#             
#     for elt in liste_100_28 : 
#         if elt in liste_100_vert_jaune :
#             if elt not in idem :
#                 idem.append(elt)
#         else :
#             if elt not in diff2 :
#                 print(elt)
#                 diff2.append(elt)
#     print(idem)
#     print(len(idem))
#     
#     print(diff1)
#     print(len(diff1))
#     print(diff2)
#     print(len(diff2))
    
    #nb_ok = [598, 969, 3124, 4959, 11244, 17357, 20902, 21281, 21306, 21309, 21309] ## avec signature des 5 éléments
#     nb_ok = [674, 1135, 2925, 3839, 10014, 13905, 19224, 20864, 21294, 21309, 21309] ## avec signature des 28 éléments
#     
#     ax = plt.gca()
#     ax.set_xticks(np.arange(0,11,1))
#     ax.set_xticklabels([round(i, 1) for i in np.arange(1.0, -0.1, -0.1)]) 
#     ax.set_yticks(np.arange(0,25000, 2000))
#     ax.set_xlabel("Pourcentage de ressemblance")
#     ax.set_ylabel("Nombre d'éléments")
#     for xy in zip(np.arange(0,11,1), nb_ok) :                                     # <--
#         ax.annotate('%s' % xy[1], xy=xy, textcoords='data')
#     plt.title("Similarité de séquence avec le groupe 11")
#     plt.plot(nb_ok)
#     plt.grid()
#     plt.show()
    
    

    
    
#     dico_graphes = recherche_graphe_commun([4], CLUSTERING_PEREZ_VERSION_NON_CAN_2[11], 0.7, "toutes_aretes_coeff_all1", "taille_4", "clustering_perez_sim_groupe_%d"%11)





    
        