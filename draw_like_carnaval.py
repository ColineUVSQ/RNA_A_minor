'''
Created on 16 juin 2019

@author: coline

Dessiner les extensions version CaRNAVal
'''

import pickle
import operator
import itertools
import subprocess
import networkx as nx
import numpy as np
import os
from recup_data.constantes import EXTENSION_PATH_TAILLE, EXTENSION_PATH,\
    NEW_EXTENSION_PATH_TAILLE
from collections import OrderedDict
from recup_data.draw_extension import draw_extension
from recup_data.draw_isomorphism_structure import draw_isomorphism_struct

'''issu du code fourni par Vladimir utilise dans CaRNAval'''
def neighbors(s1,s2):

    """

        Check whether two strands, identified by their strand orders or by their x-ordinates, are adjacent.

    """

    return abs(abs(s1)-abs(s2))==1

'''je ne sais plus a quoi ça sert'''
def motif_commun_to_2D(graphe_commun):
    liste_noeuds_faciles = []
    for noeud, data in graphe_commun.nodes(data=True) :
        if data["type"] != 1 and data["type"] != 0 :
            liste_noeuds_faciles.append(noeud)
    
    graphe_commun_2D = graphe_commun.subgraph(liste_noeuds_faciles)
    
    compteur = graphe_commun_2D.number_of_nodes()
    for noeud, data in graphe_commun.nodes(data=True) :
        can_interne = False
        b53 = -1
        if data["type"] == 1 :
            for voisin in graphe_commun[noeud] :
                for edge in graphe_commun[noeud][voisin] :
                    if graphe_commun[noeud][voisin][edge]["label"] == 'CWW' :
                        if graphe_commun.nodes[voisin]["type"] != -1 :
                            can_interne = True
                    elif graphe_commun[noeud][voisin][edge]["label"] == 'B53' :
                        b53 = ("succ", voisin)
                        
            for pred in graphe_commun.predecessors(noeud) :
                for edge in graphe_commun[pred][noeud] :
                    if graphe_commun[pred][noeud][edge]["label"] == 'B53' :
                        b53 = ("pred", pred)
                               
            if not can_interne :
                graphe_commun_2D.add_node(compteur, **data)
                graphe_commun_2D.add_node(compteur+1, **data)
                graphe_commun_2D.add_edge(compteur, compteur+1, label='CWW', long_range=False)
                if b53 != -1 :
                    if b53[0] == "pred" :
                        graphe_commun_2D.add_edge(compteur, b53[1], label='B53', long_range=False)
                    else :
                        graphe_commun_2D.add_edge(b53[1], compteur, label='B53', long_range=False)
                compteur += 2
                
                for _ in range(1,data["poids"]) :
                    graphe_commun_2D.add_node(compteur, **data)
                    graphe_commun_2D.add_node(compteur+1, **data)
                    graphe_commun_2D.add_edge(compteur, compteur+1, label='CWW', long_range=False)
                    graphe_commun_2D.add_edge(compteur-1, compteur, label='B53', long_range=False)
                    compteur += 2
                    
            else :
                for voisin in graphe_commun[noeud] :
                    for edge in graphe_commun[noeud][voisin] :
                        if graphe_commun[noeud][voisin][edge]["label"] == 'CWW' :
                            if graphe_commun.nodes[noeud]["chaine"][0] < graphe_commun.nodes[voisin]["chaine"][0] :
                                graphe_commun_2D.add_node(compteur, **data)
                                graphe_commun_2D.add_node(compteur+1, **data)
                                graphe_commun_2D.add_edge(compteur, compteur+1, label='CWW', long_range=False)
                                
                                if b53 != -1 :
                                    if b53[0] == "pred" :
                                        graphe_commun_2D.add_edge(compteur, b53[1], label='B53', long_range=False)
                                    else :
                                        graphe_commun_2D.add_edge(b53[1], compteur, label='B53', long_range=False)
                                        
                                for voisin_2 in graphe_commun[voisin] :
                                    for edge in graphe_commun[voisin][voisin_2] :
                                        if graphe_commun[voisin][voisin_2][edge]["label"] == 'B53' :   
                                            b53 = ("succ", voisin_2)  
                                              
                                for pred in graphe_commun.predecessors(voisin) :
                                    for edge in graphe_commun[pred][voisin] :
                                        if graphe_commun[pred][voisin][edge]["label"] == 'B53' :
                                            b53 = ("pred", pred)
                                            
                                if b53 != -1 :
                                    if b53[0] == "pred" :
                                        graphe_commun_2D.add_edge(compteur+1, b53[1], label='B53', long_range=False)
                                    else :
                                        graphe_commun_2D.add_edge(b53[1], compteur+1, label='B53', long_range=False)
                                
                                compteur += 2
                                
                                for _ in range(1,data["poids"]) :
                                    graphe_commun_2D.add_node(compteur, **data)
                                    graphe_commun_2D.add_node(compteur+1, **data)
                                    graphe_commun_2D.add_edge(compteur, compteur+1, label='CWW', long_range=False)
                                    graphe_commun_2D.add_edge(compteur-2, compteur, label='B53', long_range=False)
                                    graphe_commun_2D.add_edge(compteur-1, compteur+1, label='B53', long_range=False)
                                    compteur +=2
                                
                                
                
        if data["type"] == 0 :
            for voisin in graphe_commun[noeud] :
                for edge in graphe_commun[noeud][voisin] :
                    if graphe_commun[noeud][voisin][edge]["label"] == 'B53' :
                        b53 = ("succ", voisin)
                        
            for pred in graphe_commun.predecessors(noeud) :
                for edge in graphe_commun[pred][noeud] :
                    if graphe_commun[pred][noeud][edge]["label"] == 'B53' :
                        b53 = ("pred", pred)
            
            graphe_commun_2D.add_node(compteur, **data)                 
            if b53 != -1 :
                if b53[0] == "pred" :
                    graphe_commun_2D.add_edge(compteur, b53[1], label='B53', long_range=False)
                else :
                    graphe_commun_2D.add_edge(b53[1], compteur, label='B53', long_range=False)  
                compteur += 1
                         
            for _ in range(1,data["poids"]) :
                    graphe_commun_2D.add_node(compteur, **data)
                    graphe_commun_2D.add_edge(compteur-1, compteur, label='B53', long_range=False)
                    compteur += 1                   
    
                    
    return graphe_commun_2D                   
    
            
'''Transformation du graphe commun a deux graphes d'extension en un graphe global
(decontraction des sommets de type 1 et de type 0, recuperation des aretes non covalentes associes aux sommets de type 1 etc ...) 
ne marche pas pour tous les cas'''
def isomorphisme_extension_to_global(graphe_commun, graphe_ext_1, graphe_ext_2, graphe_global_1, graphe_global_2):
    liste_motif = [(1,1), (2,2), (3,3), (4,4), (5,5)]
    liste_aretes_motif = [((1,1),(2,2)), ((2,2),(1,1)), ((1,1), (5,5)), ((5,5), (1,1)), ((3,3),(4,4)), ((4,4), (3,3)), ((2,2), (5,5)), ((5,5), (2,2)), ((3,3), (1,1)), ((2,2), (4,4))]
    
    graphe_global_commun_1 = nx.MultiDiGraph()
    graphe_global_commun_2 = nx.MultiDiGraph()
    for noeud in graphe_commun.nodes() :
        nb_voisins = 0
        for voisin in graphe_commun[noeud] :
            for edge in graphe_commun[noeud][voisin] :
                if graphe_commun[noeud][voisin][edge]["label"] != 'B53' :
                    nb_voisins += 1
        if graphe_ext_1.nodes[noeud[0]]["type"] != -1 and nb_voisins > 0 :
            pos_1 = graphe_ext_1.nodes[noeud[0]]["position"]
            pos_2 = graphe_ext_2.nodes[noeud[1]]["position"]
            
            mini_poids = min(graphe_ext_1.nodes[noeud[0]]["poids"], graphe_ext_2.nodes[noeud[1]]["poids"])
            
            ## verfie si le noeud n est pas implique dans une liaison 1-1 interne au graphe
            can_interne = False
            if graphe_ext_1.nodes[noeud[0]]["type"] == 1 :
                for voisin in graphe_commun[noeud] :
                    for edge in graphe_commun[noeud][voisin] : 
                        if graphe_commun[noeud][voisin][edge]["label"] == 'CWW' and graphe_ext_1.nodes[voisin[0]]["label"] == 1 :
                            can_interne = True
#             print(pos_1)
#             print(pos_2)       
#             print(graphe_ext_1.nodes[noeud[0]]["chaine"])         
#             print(can_interne)
            if not can_interne :
                for i in range(pos_1[0], min(pos_1[0]+mini_poids, pos_1[1]+1)) :
                    if noeud in liste_motif : 
                        graphe_global_1.nodes[i].update({'motif' : True})
                    else :
                        graphe_global_1.nodes[i].update({'motif' : False})
                    graphe_global_commun_1.add_node(i, **graphe_global_1.nodes[i])
            else :
                if 1 in graphe_ext_1.nodes[noeud[0]]["chaine"] or 4 in graphe_ext_1.nodes[noeud[0]]["chaine"] :
                    for i in range(pos_1[0], min(pos_1[0]+mini_poids, pos_1[1]+1)) :
                        if noeud in liste_motif : 
                            graphe_global_1.nodes[i].update({'motif' : True})
                        else :
                            graphe_global_1.nodes[i].update({'motif' : False})
                        graphe_global_commun_1.add_node(i, **graphe_global_1.nodes[i]) 
                else :
                    print("petit rat")
                    for i in np.arange(pos_1[1], min(pos_1[1]-mini_poids, pos_1[0]), -1) :
                        print("tout petit rat")
                        print(i)
                        if noeud in liste_motif : 
                            graphe_global_1.nodes[i].update({'motif' : True})
                        else :
                            graphe_global_1.nodes[i].update({'motif' : False})
                        graphe_global_commun_1.add_node(i, **graphe_global_1.nodes[i])
            
            can_interne = False
            if graphe_ext_2.nodes[noeud[1]]["label"] == 1 :
                for voisin in graphe_commun[noeud] :
                    for edge in graphe_commun[noeud][voisin] : 
                        if graphe_commun[noeud][voisin][edge]["label"] == 'CWW' and graphe_ext_2.nodes[voisin[1]]["label"] == 1 :
                            can_interne = True
#             print(pos_1)
#             print(pos_2) 
#             print(graphe_ext_2.nodes[noeud[1]]["chaine"])
#             print(can_interne)
            
            if not can_interne :
                for i in range(pos_2[0], min(pos_2[0]+mini_poids, pos_2[1]+1)) :
                    if noeud in liste_motif : 
                        graphe_global_2.nodes[i].update({'motif' : True})
                    else :
                        graphe_global_2.nodes[i].update({'motif' : False})
                    graphe_global_commun_2.add_node(i, **graphe_global_2.nodes[i])
            else :
                if 1 in graphe_ext_2.nodes[noeud[1]]["chaine"] or 4 in graphe_ext_2.nodes[noeud[1]]["chaine"] :
                    for i in range(pos_2[0], min(pos_2[0]+mini_poids, pos_2[1]+1)) :
                        if noeud in liste_motif : 
                            graphe_global_2.nodes[i].update({'motif' : True})
                        else :
                            graphe_global_2.nodes[i].update({'motif' : False})
                        graphe_global_commun_2.add_node(i, **graphe_global_2.nodes[i]) 
                else :
                    print("petit rat")
                    for i in np.arange(pos_2[1], min(pos_2[1]-mini_poids, pos_2[0]), -1) :
                        print("tout petit rat")
                        print(i)
                        if noeud in liste_motif : 
                            graphe_global_2.nodes[i].update({'motif' : True})
                        else :
                            graphe_global_2.nodes[i].update({'motif' : False})
                        graphe_global_commun_2.add_node(i, **graphe_global_2.nodes[i])
    print("gros gros rat")
    print(graphe_global_commun_1.nodes.data())            
    for u,v,data in graphe_commun.edges(data=True) :
#         print((u,v))
        if data["type"] != '0' :
            pos_1_1 = graphe_ext_1.nodes[u[0]]["position"]
            pos_1_2 = graphe_ext_1.nodes[v[0]]["position"]
            
            pos_2_1 = graphe_ext_2.nodes[u[1]]["position"]
            pos_2_2 = graphe_ext_2.nodes[v[1]]["position"]
            
            for i in range(pos_1_1[0], pos_1_1[1]+1) :
                for voisin in graphe_global_1[i] :
                    for edge in graphe_global_1[i][voisin] :
                        if voisin >= pos_1_2[0] and voisin <= pos_1_2[1] and graphe_global_1[i][voisin][edge]["label"][:3] == data["type"][:3] :
                            
                            if (i,voisin) not in graphe_global_commun_1 and i in graphe_global_commun_1.nodes() and voisin in graphe_global_commun_1.nodes():
                                if (u,v) in liste_aretes_motif :
                                    graphe_global_1[i][voisin][edge].update({"motif" : True})
                                else :
                                    graphe_global_1[i][voisin][edge].update({"motif" : False})
                                graphe_global_commun_1.add_edge(i,voisin, **graphe_global_1[i][voisin][edge])
                            if (voisin, i) not in graphe_global_commun_1 and graphe_global_1[i][voisin][edge]["label"] != 'B53' and i in graphe_global_commun_1.nodes() and voisin in graphe_global_commun_1.nodes() :
                                if (v,u) in liste_aretes_motif :
                                    graphe_global_1[voisin][i][edge].update({"motif" : True})
                                else :
                                    graphe_global_1[voisin][i][edge].update({"motif" : False})
                                graphe_global_commun_1.add_edge(voisin, i, **graphe_global_1[voisin][i][edge])
#                             print("petit rat")
#                             print(i)
#                             print(voisin)
            
            
            for i in range(pos_2_1[0], pos_2_1[1]+1) :
                for voisin in graphe_global_2[i] :
                    for edge in graphe_global_2[i][voisin] :
                        if voisin >= pos_2_2[0] and voisin <= pos_2_2[1] and graphe_global_2[i][voisin][edge]["label"][:3] == data["type"][:3] :
                            if (i,voisin) not in graphe_global_commun_2 and i in graphe_global_commun_2.nodes() and voisin in graphe_global_commun_2.nodes():
                                if (u,v) in liste_aretes_motif :
                                    graphe_global_2[i][voisin][edge].update({"motif" : True})
                                else :
                                    graphe_global_2[i][voisin][edge].update({"motif" : False})
                                graphe_global_commun_2.add_edge(i,voisin, **graphe_global_2[i][voisin][edge])
                            if (voisin, i) not in graphe_global_commun_2 and graphe_global_2[i][voisin][edge]["label"] != 'B53' and i in graphe_global_commun_2.nodes() and voisin in graphe_global_commun_2.nodes() :
                                if (v,u) in liste_aretes_motif :
                                    graphe_global_2[voisin][i][edge].update({"motif" : True})
                                else :
                                    graphe_global_2[voisin][i][edge].update({"motif" : False})                                
                                graphe_global_commun_2.add_edge(voisin, i, **graphe_global_2[voisin][i][edge])
#                             print("petit rat")
#                             print(i)
#                             print(voisin)
            
            for i in range(pos_1_1[0], pos_1_1[1]+1) :
                for pred_1 in graphe_global_1.predecessors(i) :
                    for edge_1 in graphe_global_1[pred_1][i] :
                        for j in range(pos_2_1[0], pos_2_1[1]+1) :
                            for pred_2 in graphe_global_2.predecessors(j) :
                                for edge_2 in graphe_global_2[pred_2][j] :
                                    if pred_1 >= pos_1_2[0] and pred_1 <= pos_1_2[1] and graphe_global_1[pred_1][i][edge_1]["label"] == data["type"] and data["type"] == 'B53' and pred_2 >= pos_2_2[0] and pred_2 <= pos_2_2[1] and graphe_global_2[pred_2][j][edge_2]["label"] == data["type"] and data["type"] == 'B53' :
                                        if (pred_1, i) not in graphe_global_commun_1.edges() and pred_1 in graphe_global_commun_1.nodes() and i in graphe_global_commun_1.nodes() and (pred_2, j) not in graphe_global_commun_2.edges() and pred_2 in graphe_global_commun_2.nodes() and j in graphe_global_commun_2.nodes() :
                                            if (u,v) in liste_aretes_motif or (v,u) in liste_aretes_motif :
                                                graphe_global_1[pred_1][i][edge_1].update({"motif" : True})
                                                graphe_global_2[pred_2][j][edge_2].update({"motif" : True})
                                            else :
                                                graphe_global_1[pred_1][i][edge_1].update({"motif" : False})
                                                graphe_global_2[pred_2][j][edge_2].update({"motif" : False})
                                            graphe_global_commun_1.add_edge(pred_1, i, **graphe_global_1[pred_1][i][edge_1])
                                            graphe_global_commun_2.add_edge(pred_2, j, **graphe_global_2[pred_2][j][edge_2])
                
            if graphe_ext_1.nodes[u[0]]["type"] == 0 or graphe_ext_1.nodes[v[0]]["type"] == 0 :
                    if graphe_ext_1.nodes[u[0]]["type"] == 0 :
                        pos_1 = graphe_ext_1.nodes[u[0]]["position"]
                        pos_2 = graphe_ext_2.nodes[u[1]]["position"]
                        mini_poids = min(graphe_ext_1.nodes[u[0]]["poids"], graphe_ext_2.nodes[u[1]]["poids"])
                    else :
                        pos_1 = graphe_ext_1.nodes[v[0]]["position"]
                        pos_2 = graphe_ext_2.nodes[v[1]]["position"]
                        mini_poids = min(graphe_ext_1.nodes[v[0]]["poids"], graphe_ext_2.nodes[v[1]]["poids"])
                    compteur = 0
                    for i in range(pos_1[0], min(pos_1[0]+mini_poids, pos_1[1]+1)) :
                        if compteur > 0 :
                            graphe_global_commun_1.add_edge(i-1, i, label='B53', long_range=False, motif=False)
                        compteur += 1
                    
                    compteur = 0
                    for i in range(pos_2[0], min(pos_2[0]+mini_poids, pos_2[1]+1)) :
                        if compteur > 0 :
                            graphe_global_commun_2.add_edge(i-1, i, label='B53', long_range=False, motif=False)
                        compteur += 1
                          
        elif graphe_ext_1.nodes[u[0]]["type"] == 1 or graphe_ext_1.nodes[v[0]]["type"] == 1 :
            if graphe_ext_1.nodes[u[0]]["type"] == 1 :
                pos_1 = graphe_ext_1.nodes[u[0]]["position"]
                pos_2 = graphe_ext_2.nodes[u[1]]["position"]
                mini_poids = min(graphe_ext_1.nodes[u[0]]["poids"], graphe_ext_2.nodes[u[1]]["poids"])
            else :
                pos_1 = graphe_ext_1.nodes[v[0]]["position"]
                pos_2 = graphe_ext_2.nodes[v[1]]["position"]
                mini_poids = min(graphe_ext_1.nodes[v[0]]["poids"], graphe_ext_2.nodes[v[1]]["poids"])
            
            
            compteur = 0  
            for i in range(pos_1[0], min(pos_1[0]+mini_poids, pos_1[1]+1)) :
                if compteur > 0 :
                    graphe_global_commun_1.add_edge(i-1, i, label='B53', long_range=False, motif=False)
                for voisin in graphe_global_1[i] :
                    for edge in graphe_global_1[i][voisin] :
                        if graphe_global_1[i][voisin][edge]["label"] == 'CWW' :
                            graphe_global_1.nodes[voisin].update({'motif' : False})
                            graphe_global_commun_1.add_node(voisin, **graphe_global_1.nodes[voisin])
                            if (i,voisin) not in graphe_global_commun_1.edges() and i in graphe_global_commun_1.nodes() and voisin in graphe_global_commun_1.nodes():
                                graphe_global_1[i][voisin][edge].update({'motif' : False})
                                graphe_global_commun_1.add_edge(i,voisin, **graphe_global_1[i][voisin][edge])
                            if (voisin, i) not in graphe_global_commun_1.edges() and i in graphe_global_commun_1.nodes() and voisin in graphe_global_commun_1.nodes() :
                                graphe_global_1[voisin][i][edge].update({'motif' : False})
                                graphe_global_commun_1.add_edge(voisin, i, **graphe_global_1[voisin][i][edge])
                                
                            if min(pos_1[0]+mini_poids, pos_1[1]+1) - pos_1[0] - 1 > 0 and compteur > 0:
#                                 print("ramou")
#                                 print(voisin)
#                                 print(pos_1)
#                                 print(mini_poids)
                                graphe_global_commun_1.add_edge(voisin, voisin+1, label='B53', long_range=False, motif=False)
#                                 print(graphe_global_commun_1.edges.data())
                compteur += 1
            
            compteur = 0                  
            for i in range(pos_2[0], min(pos_2[0]+mini_poids, pos_2[1]+1)) :
                if compteur > 0 :
                    graphe_global_commun_2.add_edge(i-1, i, label='B53', long_range=False, motif=False)
                for voisin in graphe_global_2[i] :
                    for edge in graphe_global_2[i][voisin] :
                        if graphe_global_2[i][voisin][edge]["label"] == 'CWW' :
                            graphe_global_2.nodes[voisin].update({'motif' : False})
                            graphe_global_commun_2.add_node(voisin, **graphe_global_2.nodes[voisin])
                            if (i,voisin) not in graphe_global_commun_2.edges() and i in graphe_global_commun_2.nodes() and voisin in graphe_global_commun_2.nodes():
                                graphe_global_2[i][voisin][edge].update({'motif' : False})
                                graphe_global_commun_2.add_edge(i,voisin, **graphe_global_2[i][voisin][edge])
                            if (voisin, i) not in graphe_global_commun_2.edges() and i in graphe_global_commun_2.nodes() and voisin in graphe_global_commun_2.nodes():
                                graphe_global_2[voisin][i][edge].update({'motif' : False})
                                graphe_global_commun_2.add_edge(voisin, i, **graphe_global_2[voisin][i][edge])
                                
                            if min(pos_2[0]+mini_poids, pos_2[1]+1) - pos_2[0] - 1 > 0 and compteur > 0 :
#                                 print("ramou")
#                                 print(voisin)
                                graphe_global_commun_2.add_edge(voisin, voisin+1, label='B53', long_range=False, motif=False)
                compteur += 1
    
    ##liens B53 motif qui manquent
    pos_2_1 = graphe_ext_1.nodes[2]["position"]
    pos_4_1 = graphe_ext_1.nodes[4]["position"]
    
    pos_2_2 = graphe_ext_2.nodes[2]["position"]
    pos_4_2 = graphe_ext_2.nodes[4]["position"]
    
    for voisin_1 in graphe_global_1[pos_2_1[0]] : 
#         print(voisin)
        for edge_1 in graphe_global_1[pos_2_1[0]][voisin_1] :
                if  graphe_global_1[pos_2_1[0]][voisin_1][edge_1]["label"] == 'CWW' :
                    for voisin_2 in graphe_global_2[pos_2_2[0]] : 
                        for edge_2 in graphe_global_2[pos_2_2[0]][voisin_2] :
                            if  graphe_global_2[pos_2_2[0]][voisin_2][edge_2]["label"] == 'CWW' :
                                if (voisin_2-1,voisin_2) not in graphe_global_commun_2.edges() and voisin_2 in graphe_global_commun_2.nodes() and voisin_2-1 in graphe_global_commun_2.nodes() and (voisin_1-1,voisin_1) not in graphe_global_commun_1.edges() and voisin_1 in graphe_global_commun_1.nodes() and voisin_1-1 in graphe_global_commun_1.nodes():
                                    graphe_global_commun_1.add_edge(voisin_1-1,voisin_1, label='B53', long_range=False, motif=True)
                                    graphe_global_commun_2.add_edge(voisin_2-1,voisin_2, label='B53', long_range=False, motif=True)     
                        
    for voisin_1 in graphe_global_1[pos_2_1[0]-1] : 
#         print(voisin)
        for edge_1 in graphe_global_1[pos_2_1[0]-1][voisin_1] :
                if  graphe_global_1[pos_2_1[0]-1][voisin_1][edge_1]["label"] == 'CWW' :
                    for voisin_2 in graphe_global_2[pos_2_2[0]-1] : 
                        for edge_2 in graphe_global_2[pos_2_2[0]-1][voisin_2] :
                            if  graphe_global_2[pos_2_2[0]-1][voisin_2][edge_2]["label"] == 'CWW' :
                                if (voisin_2-1,voisin_2) not in graphe_global_commun_2.edges() and voisin_2 in graphe_global_commun_2.nodes() and voisin_2-1 in graphe_global_commun_2.nodes() and (voisin_1-1,voisin_1) not in graphe_global_commun_1.edges() and voisin_1 in graphe_global_commun_1.nodes() and voisin_1-1 in graphe_global_commun_1.nodes():
                                    graphe_global_commun_1.add_edge(voisin_1-1,voisin_1, label='B53', long_range=False, motif=False)
                                    graphe_global_commun_2.add_edge(voisin_2-1,voisin_2, label='B53', long_range=False, motif=False)  
    
    for voisin_1 in graphe_global_1[pos_4_1[0]+1] : 
#         print(voisin)
        for edge_1 in graphe_global_1[pos_4_1[0]+1][voisin_1] :
                if  graphe_global_1[pos_4_1[0]+1][voisin_1][edge_1]["label"] == 'CWW' :
                    for voisin_2 in graphe_global_2[pos_4_2[0]+1] : 
                        for edge_2 in graphe_global_2[pos_4_2[0]+1][voisin_2] :
                            if  graphe_global_2[pos_4_2[0]+1][voisin_2][edge_2]["label"] == 'CWW' :
                                if (voisin_2,voisin_2+1) not in graphe_global_commun_2.edges() and voisin_2 in graphe_global_commun_2.nodes() and voisin_2+1 in graphe_global_commun_2.nodes() and (voisin_1,voisin_1+1) not in graphe_global_commun_1.edges() and voisin_1 in graphe_global_commun_1.nodes() and voisin_1+1 in graphe_global_commun_1.nodes():
                                    graphe_global_commun_1.add_edge(voisin_1,voisin_1+1, label='B53', long_range=False, motif=False)
                                    graphe_global_commun_2.add_edge(voisin_2,voisin_2+1, label='B53', long_range=False, motif=False)  
    
         
    
            
#     if (pos_1[0], pos_3[0]) not in graphe_global_commun_1.edges() :
#         graphe_global_commun_1.add_edge(pos_3[0], pos_1[0], label='B53', long_range=False)
#         
#     pos_1 = graphe_ext_2.nodes[1]["position"]
#     pos_3 = graphe_ext_2.nodes[3]["position"]
#     if (pos_1[0], pos_3[0]) not in graphe_global_commun_2.edges() :
#         graphe_global_commun_2.add_edge(pos_3[0], pos_1[0], label='B53', long_range=False)
    
#     for pred in graphe_commun.predecessors((2,2)) :
#         pos_pred = graphe_ext_1.nodes[pred[0]]["position"] 
#         pos_2 = graphe_ext_1.nodes[2]["position"]
#         
#         for i in range(pos_pred[0], )
                           
    print(graphe_global_commun_1.edges.data())  
    print(graphe_global_commun_2.nodes.data()) 
    
    return graphe_global_commun_1, graphe_global_commun_2


'''03/03/20 Transformation du graphe commun a deux graphes d'extension en un graphe global
(decontraction des sommets de type 1 et de type 0, recuperation des aretes non covalentes associes aux sommets de type 1 etc ...) 
ne marche pas pour tous les cas'''
def isomorphisme_extension_to_global_version_new_data(graphe_commun, graphe_ext_1, graphe_ext_2, graphe_global_1, graphe_global_2):
    liste_motif = [(1,1), (2,2), (3,3), (4,4), (5,5)]
    liste_aretes_motif = [((1,1),(2,2)), ((2,2),(1,1)), ((1,1), (5,5)), ((5,5), (1,1)), ((3,3),(4,4)), ((4,4), (3,3)), ((2,2), (5,5)), ((5,5), (2,2)), ((3,3), (1,1)), ((2,2), (4,4))]
    
    graphe_global_commun_1 = nx.MultiDiGraph()
    graphe_global_commun_2 = nx.MultiDiGraph()
    for noeud in graphe_commun.nodes() :
        nb_voisins = 0
        for voisin in graphe_commun[noeud] :
            for edge in graphe_commun[noeud][voisin] :
                if graphe_commun[noeud][voisin][edge]["label"] != 'B53' :
                    nb_voisins += 1
        if graphe_ext_1.nodes[noeud[0]]["type"] != -1 and nb_voisins > 0 :
            pos_1 = graphe_ext_1.nodes[noeud[0]]["position"]
            pos_2 = graphe_ext_2.nodes[noeud[1]]["position"]
            
            mini_poids = min(graphe_ext_1.nodes[noeud[0]]["poids"], graphe_ext_2.nodes[noeud[1]]["poids"])
            
            ## verfie si le noeud n est pas implique dans une liaison 1-1 interne au graphe
            can_interne = False
            if graphe_ext_1.nodes[noeud[0]]["type"] == 1 :
                for voisin in graphe_commun[noeud] :
                    for edge in graphe_commun[noeud][voisin] : 
                        if graphe_commun[noeud][voisin][edge]["label"] == 'CWW' and graphe_ext_1.nodes[voisin[0]]["type"] == 1 :
                            can_interne = True
#             print(pos_1)
#             print(pos_2)       
#             print(graphe_ext_1.nodes[noeud[0]]["chaine"])         
#             print(can_interne)
            if not can_interne :
                for i in range(pos_1[0], min(pos_1[0]+mini_poids, pos_1[1]+1)) :
                    if noeud in liste_motif : 
                        graphe_global_1.nodes[(graphe_ext_1.nodes[noeud[0]]["num_ch"],i)].update({'motif' : True})
                    else :
                        graphe_global_1.nodes[(graphe_ext_1.nodes[noeud[0]]["num_ch"],i)].update({'motif' : False})
                    graphe_global_commun_1.add_node((graphe_ext_1.nodes[noeud[0]]["num_ch"],i), **graphe_global_1.nodes[(graphe_ext_1.nodes[noeud[0]]["num_ch"],i)])
            else :
                if 1 in graphe_ext_1.nodes[noeud[0]]["chaine"] or 4 in graphe_ext_1.nodes[noeud[0]]["chaine"] :
                    for i in range(pos_1[0], min(pos_1[0]+mini_poids, pos_1[1]+1)) :
                        if noeud in liste_motif : 
                            graphe_global_1.nodes[(graphe_ext_1.nodes[noeud[0]]["num_ch"],i)].update({'motif' : True})
                        else :
                            graphe_global_1.nodes[(graphe_ext_1.nodes[noeud[0]]["num_ch"],i)].update({'motif' : False})
                        graphe_global_commun_1.add_node((graphe_ext_1.nodes[noeud[0]]["num_ch"],i), **graphe_global_1.nodes[(graphe_ext_1.nodes[noeud[0]]["num_ch"],i)]) 
                else :
                    print("petit rat")
                    for i in np.arange(pos_1[1], min(pos_1[1]-mini_poids, pos_1[0]), -1) :
                        print("tout petit rat")
                        print(i)
                        if noeud in liste_motif : 
                            graphe_global_1.nodes[(graphe_ext_1.nodes[noeud[0]]["num_ch"],i)].update({'motif' : True})
                        else :
                            graphe_global_1.nodes[(graphe_ext_1.nodes[noeud[0]]["num_ch"],i)].update({'motif' : False})
                        graphe_global_commun_1.add_node((graphe_ext_1.nodes[noeud[0]]["num_ch"],i), **graphe_global_1.nodes[(graphe_ext_1.nodes[noeud[0]]["num_ch"],i)])
            
            can_interne = False
            if graphe_ext_2.nodes[noeud[1]]["type"] == 1 :
                for voisin in graphe_commun[noeud] :
                    for edge in graphe_commun[noeud][voisin] : 
                        if graphe_commun[noeud][voisin][edge]["label"] == 'CWW' and graphe_ext_2.nodes[voisin[1]]["type"] == 1 :
                            can_interne = True
#             print(pos_1)
#             print(pos_2) 
#             print(graphe_ext_2.nodes[noeud[1]]["chaine"])
#             print(can_interne)
            
            if not can_interne :
                for i in range(pos_2[0], min(pos_2[0]+mini_poids, pos_2[1]+1)) :
                    if noeud in liste_motif : 
                        graphe_global_2.nodes[(graphe_ext_2.nodes[noeud[1]]["num_ch"],i)].update({'motif' : True})
                    else :
                        graphe_global_2.nodes[(graphe_ext_2.nodes[noeud[1]]["num_ch"],i)].update({'motif' : False})
                    graphe_global_commun_2.add_node((graphe_ext_2.nodes[noeud[1]]["num_ch"],i), **graphe_global_2.nodes[(graphe_ext_2.nodes[noeud[1]]["num_ch"],i)])
            else :
                if 1 in graphe_ext_2.nodes[noeud[1]]["chaine"] or 4 in graphe_ext_2.nodes[noeud[1]]["chaine"] :
                    for i in range(pos_2[0], min(pos_2[0]+mini_poids, pos_2[1]+1)) :
                        if noeud in liste_motif : 
                            graphe_global_2.nodes[(graphe_ext_2.nodes[noeud[1]]["num_ch"],i)].update({'motif' : True})
                        else :
                            graphe_global_2.nodes[(graphe_ext_2.nodes[noeud[1]]["num_ch"],i)].update({'motif' : False})
                        graphe_global_commun_2.add_node((graphe_ext_2.nodes[noeud[1]]["num_ch"],i), **graphe_global_2.nodes[(graphe_ext_2.nodes[noeud[1]]["num_ch"],i)]) 
                else :
                    print("petit rat")
                    for i in np.arange(pos_2[1], min(pos_2[1]-mini_poids, pos_2[0]), -1) :
                        print("tout petit rat")
                        print(i)
                        if noeud in liste_motif : 
                            graphe_global_2.nodes[(graphe_ext_2.nodes[noeud[1]]["num_ch"],i)].update({'motif' : True})
                        else :
                            graphe_global_2.nodes[(graphe_ext_2.nodes[noeud[1]]["num_ch"],i)].update({'motif' : False})
                        graphe_global_commun_2.add_node((graphe_ext_2.nodes[noeud[1]]["num_ch"],i), **graphe_global_2.nodes[(graphe_ext_2.nodes[noeud[1]]["num_ch"],i)])
    print("gros gros rat")
    print(graphe_global_commun_1.nodes.data())            
    for u,v,data in graphe_commun.edges(data=True) :
#         print((u,v))
        if data["label"] != '0' :
            pos_1_1 = graphe_ext_1.nodes[u[0]]["position"]
            pos_1_2 = graphe_ext_1.nodes[v[0]]["position"]
            
            pos_2_1 = graphe_ext_2.nodes[u[1]]["position"]
            pos_2_2 = graphe_ext_2.nodes[v[1]]["position"]
            
            for i in range(pos_1_1[0], pos_1_1[1]+1) :
                num_noeud = (graphe_ext_1.nodes[u[0]]["num_ch"], i)
                for voisin in graphe_global_1[num_noeud] :
                    for edge in graphe_global_1[num_noeud][voisin] :
                        if voisin[0] == graphe_ext_1.nodes[u[0]]["num_ch"] and voisin[1] >= pos_1_2[0] and voisin[1] <= pos_1_2[1] and graphe_global_1[num_noeud][voisin][edge]["label"][:3] == data["label"][:3] :
                            
                            if (num_noeud,voisin) not in graphe_global_commun_1 and num_noeud in graphe_global_commun_1.nodes() and voisin in graphe_global_commun_1.nodes():
                                if (u,v) in liste_aretes_motif :
                                    graphe_global_1[num_noeud][voisin][edge].update({"motif" : True})
                                else :
                                    graphe_global_1[num_noeud][voisin][edge].update({"motif" : False})
                                graphe_global_commun_1.add_edge(num_noeud,voisin, **graphe_global_1[num_noeud][voisin][edge])
                            if (voisin, num_noeud) not in graphe_global_commun_1 and graphe_global_1[num_noeud][voisin][edge]["label"] != 'B53' and num_noeud in graphe_global_commun_1.nodes() and voisin in graphe_global_commun_1.nodes() :
                                if (v,u) in liste_aretes_motif :
                                    graphe_global_1[voisin][num_noeud][edge].update({"motif" : True})
                                else :
                                    graphe_global_1[voisin][num_noeud][edge].update({"motif" : False})
                                graphe_global_commun_1.add_edge(voisin, num_noeud, **graphe_global_1[voisin][num_noeud][edge])
#                             print("petit rat")
#                             print(i)
#                             print(voisin)
            
            
            for i in range(pos_2_1[0], pos_2_1[1]+1) :
                num_noeud = (graphe_ext_2.nodes[u[1]]["num_ch"], i)
                for voisin in graphe_global_2[num_noeud] :
                    for edge in graphe_global_2[num_noeud][voisin] :
                        if voisin[0] == graphe_ext_2.nodes[u[1]]["num_ch"] and voisin[1] >= pos_2_2[0] and voisin[1] <= pos_2_2[1] and graphe_global_2[num_noeud][voisin][edge]["label"][:3] == data["label"][:3] :
                            if (num_noeud,voisin) not in graphe_global_commun_2 and num_noeud in graphe_global_commun_2.nodes() and voisin in graphe_global_commun_2.nodes():
                                if (u,v) in liste_aretes_motif :
                                    graphe_global_2[num_noeud][voisin][edge].update({"motif" : True})
                                else :
                                    graphe_global_2[num_noeud][voisin][edge].update({"motif" : False})
                                graphe_global_commun_2.add_edge(num_noeud,voisin, **graphe_global_2[num_noeud][voisin][edge])
                            if (voisin, num_noeud) not in graphe_global_commun_2 and graphe_global_2[num_noeud][voisin][edge]["label"] != 'B53' and num_noeud in graphe_global_commun_2.nodes() and voisin in graphe_global_commun_2.nodes() :
                                if (v,u) in liste_aretes_motif :
                                    graphe_global_2[voisin][num_noeud][edge].update({"motif" : True})
                                else :
                                    graphe_global_2[voisin][num_noeud][edge].update({"motif" : False})                                
                                graphe_global_commun_2.add_edge(voisin, num_noeud, **graphe_global_2[voisin][num_noeud][edge])
#                             print("petit rat")
#                             print(i)
#                             print(voisin)
            
            for i in range(pos_1_1[0], pos_1_1[1]+1) :
                num_noeud = (graphe_ext_1.nodes[u[0]]["num_ch"], i)
                for pred_1 in graphe_global_1.predecessors(num_noeud) :
                    for edge_1 in graphe_global_1[pred_1][num_noeud] :
                        for j in range(pos_2_1[0], pos_2_1[1]+1) :
                            num_noeud_2 = (graphe_ext_2.nodes[u[1]]["num_ch"], j)
                            for pred_2 in graphe_global_2.predecessors(num_noeud_2) :
                                for edge_2 in graphe_global_2[pred_2][num_noeud_2] :
                                    if pred_1 == graphe_ext_1.nodes[u[0]]["num_ch"] and pred_1[1] >= pos_1_2[0] and pred_1[1] <= pos_1_2[1] and graphe_global_1[pred_1][num_noeud][edge_1]["label"] == data["label"] and data["label"] == 'B53' and pred_2 >= pos_2_2[0] and pred_2 <= pos_2_2[1] and graphe_global_2[pred_2][num_noeud_2][edge_2]["label"] == data["label"] and data["label"] == 'B53' :
                                        if (pred_1, num_noeud) not in graphe_global_commun_1.edges() and pred_1 in graphe_global_commun_1.nodes() and num_noeud in graphe_global_commun_1.nodes() and (pred_2, num_noeud_2) not in graphe_global_commun_2.edges() and pred_2 in graphe_global_commun_2.nodes() and num_noeud_2 in graphe_global_commun_2.nodes() :
                                            if (u,v) in liste_aretes_motif or (v,u) in liste_aretes_motif :
                                                graphe_global_1[pred_1][num_noeud][edge_1].update({"motif" : True})
                                                graphe_global_2[pred_2][num_noeud_2][edge_2].update({"motif" : True})
                                            else :
                                                graphe_global_1[pred_1][num_noeud][edge_1].update({"motif" : False})
                                                graphe_global_2[pred_2][num_noeud_2][edge_2].update({"motif" : False})
                                            graphe_global_commun_1.add_edge(pred_1, num_noeud, **graphe_global_1[pred_1][num_noeud][edge_1])
                                            graphe_global_commun_2.add_edge(pred_2, num_noeud_2, **graphe_global_2[pred_2][num_noeud_2][edge_2])
                
            if graphe_ext_1.nodes[u[0]]["type"] == 0 or graphe_ext_1.nodes[v[0]]["type"] == 0 :
                    if graphe_ext_1.nodes[u[0]]["type"] == 0 :
                        pos_1 = graphe_ext_1.nodes[u[0]]["position"]
                        pos_2 = graphe_ext_2.nodes[u[1]]["position"]
                        ch_noeud_1 = graphe_ext_1.nodes[u[0]]["num_ch"]
                        ch_noeud_2 = graphe_ext_2.nodes[u[1]]["num_ch"]
                        mini_poids = min(graphe_ext_1.nodes[u[0]]["poids"], graphe_ext_2.nodes[u[1]]["poids"])
                    else :
                        pos_1 = graphe_ext_1.nodes[v[0]]["position"]
                        pos_2 = graphe_ext_2.nodes[v[1]]["position"]
                        ch_noeud_1 = graphe_ext_1.nodes[v[0]]["num_ch"]
                        ch_noeud_2 = graphe_ext_2.nodes[v[1]]["num_ch"]
                        mini_poids = min(graphe_ext_1.nodes[v[0]]["poids"], graphe_ext_2.nodes[v[1]]["poids"])
                    compteur = 0
                    for i in range(pos_1[0], min(pos_1[0]+mini_poids, pos_1[1]+1)) :
                        if compteur > 0 :
                            graphe_global_commun_1.add_edge((ch_noeud_1, i-1), (ch_noeud_1, i), label='B53', motif=False)
                        compteur += 1
                    
                    compteur = 0
                    for i in range(pos_2[0], min(pos_2[0]+mini_poids, pos_2[1]+1)) :
                        if compteur > 0 :
                            graphe_global_commun_2.add_edge((ch_noeud_2, i-1), (ch_noeud_2, i), label='B53', motif=False)
                        compteur += 1
                          
        elif graphe_ext_1.nodes[u[0]]["type"] == 1 or graphe_ext_1.nodes[v[0]]["type"] == 1 :
            if graphe_ext_1.nodes[u[0]]["type"] == 1 :
                pos_1 = graphe_ext_1.nodes[u[0]]["position"]
                pos_2 = graphe_ext_2.nodes[u[1]]["position"]
                ch_noeud_1 = graphe_ext_1.nodes[u[0]]["num_ch"]
                ch_noeud_2 = graphe_ext_2.nodes[u[1]]["num_ch"]
                mini_poids = min(graphe_ext_1.nodes[u[0]]["poids"], graphe_ext_2.nodes[u[1]]["poids"])
            else :
                pos_1 = graphe_ext_1.nodes[v[0]]["position"]
                pos_2 = graphe_ext_2.nodes[v[1]]["position"]
                ch_noeud_1 = graphe_ext_1.nodes[v[0]]["num_ch"]
                ch_noeud_2 = graphe_ext_2.nodes[v[1]]["num_ch"]
                mini_poids = min(graphe_ext_1.nodes[v[0]]["poids"], graphe_ext_2.nodes[v[1]]["poids"])
            
            
            compteur = 0  
            for i in range(pos_1[0], min(pos_1[0]+mini_poids, pos_1[1]+1)) :
                if compteur > 0 :
                    graphe_global_commun_1.add_edge((ch_noeud_1, i-1), (ch_noeud_1, i), label='B53', motif=False)
                for voisin in graphe_global_1[(ch_noeud_1, i)] :
                    for edge in graphe_global_1[(ch_noeud_1, i)][voisin] :
                        if graphe_global_1[(ch_noeud_1, i)][voisin][edge]["label"] == 'CWW' :
                            graphe_global_1.nodes[voisin].update({'motif' : False})
                            graphe_global_commun_1.add_node(voisin, **graphe_global_1.nodes[voisin])
                            if ((ch_noeud_1, i),voisin) not in graphe_global_commun_1.edges() and (ch_noeud_1, i) in graphe_global_commun_1.nodes() and voisin in graphe_global_commun_1.nodes():
                                graphe_global_1[(ch_noeud_1, i)][voisin][edge].update({'motif' : False})
                                graphe_global_commun_1.add_edge((ch_noeud_1, i),voisin, **graphe_global_1[(ch_noeud_1, i)][voisin][edge])
                            if (voisin,(ch_noeud_1, i)) not in graphe_global_commun_1.edges() and (ch_noeud_1, i) in graphe_global_commun_1.nodes() and voisin in graphe_global_commun_1.nodes() :
                                graphe_global_1[voisin][(ch_noeud_1, i)][edge].update({'motif' : False})
                                graphe_global_commun_1.add_edge(voisin, (ch_noeud_1, i), **graphe_global_1[voisin][(ch_noeud_1, i)][edge])
                                
                            if min(pos_1[0]+mini_poids, pos_1[1]+1) - pos_1[0] - 1 > 0 and compteur > 0:
#                                 print("ramou")
#                                 print(voisin)
#                                 print(pos_1)
#                                 print(mini_poids)
                                graphe_global_commun_1.add_edge(voisin, (voisin[0], voisin[1]+1), label='B53', motif=False)
#                                 print(graphe_global_commun_1.edges.data())
                compteur += 1
            
            compteur = 0                  
            for i in range(pos_2[0], min(pos_2[0]+mini_poids, pos_2[1]+1)) :
                if compteur > 0 :
                    graphe_global_commun_2.add_edge((ch_noeud_2, i-1), (ch_noeud_2, i), label='B53', motif=False)
                for voisin in graphe_global_2[(ch_noeud_2, i)] :
                    for edge in graphe_global_2[(ch_noeud_2, i)][voisin] :
                        if graphe_global_2[(ch_noeud_2, i)][voisin][edge]["label"] == 'CWW' :
                            graphe_global_2.nodes[voisin].update({'motif' : False})
                            graphe_global_commun_2.add_node(voisin, **graphe_global_2.nodes[voisin])
                            if ((ch_noeud_2, i),voisin) not in graphe_global_commun_2.edges() and (ch_noeud_2, i) in graphe_global_commun_2.nodes() and voisin in graphe_global_commun_2.nodes():
                                graphe_global_2[(ch_noeud_2, i)][voisin][edge].update({'motif' : False})
                                graphe_global_commun_2.add_edge((ch_noeud_2, i),voisin, **graphe_global_2[(ch_noeud_2, i)][voisin][edge])
                            if (voisin, (ch_noeud_2, i)) not in graphe_global_commun_2.edges() and (ch_noeud_2, i) in graphe_global_commun_2.nodes() and voisin in graphe_global_commun_2.nodes():
                                graphe_global_2[voisin][(ch_noeud_2, i)][edge].update({'motif' : False})
                                graphe_global_commun_2.add_edge(voisin, (ch_noeud_2, i), **graphe_global_2[voisin][(ch_noeud_2, i)][edge])
                                
                            if min(pos_2[0]+mini_poids, pos_2[1]+1) - pos_2[0] - 1 > 0 and compteur > 0 :
#                                 print("ramou")
#                                 print(voisin)
                                graphe_global_commun_2.add_edge(voisin, (voisin[0], voisin[1]+1), label='B53', motif=False)
                compteur += 1
    
    ##liens B53 motif qui manquent
    pos_2_1 = (graphe_ext_1.nodes[2]["num_ch"], graphe_ext_1.nodes[2]["position"][0])
    pos_4_1 = (graphe_ext_1.nodes[4]["num_ch"], graphe_ext_1.nodes[4]["position"][0])
    
    pos_2_2 = (graphe_ext_2.nodes[2]["num_ch"], graphe_ext_2.nodes[2]["position"][0])
    pos_4_2 = (graphe_ext_2.nodes[4]["num_ch"], graphe_ext_2.nodes[4]["position"][0])
    
    for voisin_1 in graphe_global_1[pos_2_1] : 
#         print(voisin)
        for edge_1 in graphe_global_1[pos_2_1][voisin_1] :
                if  graphe_global_1[pos_2_1][voisin_1][edge_1]["label"] == 'CWW' :
                    for voisin_2 in graphe_global_2[pos_2_2] : 
                        for edge_2 in graphe_global_2[pos_2_2][voisin_2] :
                            if  graphe_global_2[pos_2_2][voisin_2][edge_2]["label"] == 'CWW' :
                                if ((voisin_2[0], voisin_2[1]-1),voisin_2) not in graphe_global_commun_2.edges() and voisin_2 in graphe_global_commun_2.nodes() and (voisin_2[0], voisin_2[1]-1) in graphe_global_commun_2.nodes() and ((voisin_1[0], voisin_1[1]-1),voisin_1) not in graphe_global_commun_1.edges() and voisin_1 in graphe_global_commun_1.nodes() and (voisin_1[0], voisin_1[1]-1) in graphe_global_commun_1.nodes():
                                    graphe_global_commun_1.add_edge((voisin_1[0], voisin_1[1]-1),voisin_1, label='B53', motif=True)
                                    graphe_global_commun_2.add_edge((voisin_2[0], voisin_2[1]-1),voisin_2, label='B53', motif=True)     
                        
    for voisin_1 in graphe_global_1[(pos_2_1[0], pos_2_1[1]-1)] : 
#         print(voisin)
        for edge_1 in graphe_global_1[(pos_2_1[0], pos_2_1[1]-1)][voisin_1] :
                if  graphe_global_1[(pos_2_1[0], pos_2_1[1]-1)][voisin_1][edge_1]["label"] == 'CWW' :
                    for voisin_2 in graphe_global_2[(pos_2_2[0], pos_2_2[1]-1)] : 
                        for edge_2 in graphe_global_2[(pos_2_2[0], pos_2_2[1]-1)][voisin_2] :
                            if  graphe_global_2[(pos_2_2[0], pos_2_2[1]-1)][voisin_2][edge_2]["label"] == 'CWW' :
                                if ((voisin_2[0], voisin_2[1]-1),voisin_2) not in graphe_global_commun_2.edges() and voisin_2 in graphe_global_commun_2.nodes() and (voisin_2[0], voisin_2[1]-1) in graphe_global_commun_2.nodes() and ((voisin_1[0], voisin_1[1]-1),voisin_1) not in graphe_global_commun_1.edges() and voisin_1 in graphe_global_commun_1.nodes() and (voisin_1[0], voisin_1[1]-1) in graphe_global_commun_1.nodes():
                                    graphe_global_commun_1.add_edge((voisin_1[0], voisin_1[1]-1),voisin_1, label='B53', motif=False)
                                    graphe_global_commun_2.add_edge((voisin_2[0], voisin_2[1]-1),voisin_2, label='B53', motif=False)  
    
    for voisin_1 in graphe_global_1[(pos_4_1[0], pos_4_1[1]+1)] : 
#         print(voisin)
        for edge_1 in graphe_global_1[(pos_4_1[0], pos_4_1[1]+1)][voisin_1] :
                if  graphe_global_1[(pos_4_1[0], pos_4_1[1]+1)][voisin_1][edge_1]["label"] == 'CWW' :
                    for voisin_2 in graphe_global_2[(pos_4_2[0], pos_4_2[1]+1)] : 
                        for edge_2 in graphe_global_2[(pos_4_2[0], pos_4_2[1]+1)][voisin_2] :
                            if  graphe_global_2[(pos_4_2[0], pos_4_2[1]+1)][voisin_2][edge_2]["label"] == 'CWW' :
                                if (voisin_2,(voisin_2[0], voisin_2[1]+1)) not in graphe_global_commun_2.edges() and voisin_2 in graphe_global_commun_2.nodes() and (voisin_2[0], voisin_2[1]+1) in graphe_global_commun_2.nodes() and (voisin_1,(voisin_1[0], voisin_1[1]+1)) not in graphe_global_commun_1.edges() and voisin_1 in graphe_global_commun_1.nodes() and (voisin_1[0], voisin_1[1]+1) in graphe_global_commun_1.nodes():
                                    graphe_global_commun_1.add_edge(voisin_1,(voisin_1[0], voisin_1[1]+1), label='B53', motif=False)
                                    graphe_global_commun_2.add_edge(voisin_2,(voisin_2[0], voisin_2[1]+1), label='B53', motif=False)  
    
         
    
            
#     if (pos_1[0], pos_3[0]) not in graphe_global_commun_1.edges() :
#         graphe_global_commun_1.add_edge(pos_3[0], pos_1[0], label='B53', long_range=False)
#         
#     pos_1 = graphe_ext_2.nodes[1]["position"]
#     pos_3 = graphe_ext_2.nodes[3]["position"]
#     if (pos_1[0], pos_3[0]) not in graphe_global_commun_2.edges() :
#         graphe_global_commun_2.add_edge(pos_3[0], pos_1[0], label='B53', long_range=False)
    
#     for pred in graphe_commun.predecessors((2,2)) :
#         pos_pred = graphe_ext_1.nodes[pred[0]]["position"] 
#         pos_2 = graphe_ext_1.nodes[2]["position"]
#         
#         for i in range(pos_pred[0], )
                           
    print(graphe_global_commun_1.edges.data())  
    print(graphe_global_commun_2.nodes.data()) 
    
    return graphe_global_commun_1, graphe_global_commun_2


'''renvoie le rang associe a une valeur de sim'''
def rang_sim(dico_sim, sim):
    dico_sim_sorted = OrderedDict(sorted(dico_sim.items(), key= lambda t: t[1], reverse=True))
    
    etape = 1
    rang = 0
    for elt in dico_sim_sorted.values() :
        if elt == sim :
            rang = etape
        etape += 1
        
    return rang

'''Preparation de la creation de la figure tikz de superposition de deux graphes dont le graphe commun est passe en parametre'''
def draw_iso_ext_to_global(graphe_commun, cle, rep, i, sim):
        #print(sim)
        with open(rep+cle[0]+".pickle", 'rb') as fichier_1 :
            mon_depickler_1 = pickle.Unpickler(fichier_1)
            graphe_ext_1 = mon_depickler_1.load()
            
        with open(rep+cle[1]+".pickle", 'rb') as fichier_2 :
            mon_depickler_2 = pickle.Unpickler(fichier_2)
            graphe_ext_2 = mon_depickler_2.load()
        
        with open("Graphes_globaux/graphes_sequence/taille_%d/graphe_global_"%i+cle[0][8:]+"_avec_coord.pickle", 'rb') as fichier_pickle_1 :
            mon_depickler_coord_1 = pickle.Unpickler(fichier_pickle_1)
            graphe_global_1 = mon_depickler_coord_1.load()
            
        with open("Graphes_globaux/graphes_sequence/taille_%d/graphe_global_"%i+cle[1][8:]+"_avec_coord.pickle", 'rb') as fichier_pickle_2 :
            mon_depickler_coord_2 = pickle.Unpickler(fichier_pickle_2)
            graphe_global_2 = mon_depickler_coord_2.load()
        cle_0 = cle[0].split("_")[1] + "-" + cle[0].split("_")[2] + "-" +cle[0].split("_")[3] + "-"+ cle[0].split("_")[4]
        cle_1 = cle[1].split("_")[1] + "-" + cle[1].split("_")[2] + "-" +cle[1].split("_")[3] + "-"+ cle[1].split("_")[4]
        
        #if cle_0+"_"+cle_1+".pdf" not in os.listdir("Graphes_globaux/graphes_sequence/taille_%d/"%i) :#and (cle_0 != '3JCS-2-30-56' or cle_1 != '5DM6-X-328-2') :
        #if cle_0 == '5J5B-BA-294-2' and cle_1 == '3JCS-1-25-46' :
        #if cle_0 == '5J7L-DA-272-2' and cle_1 == '5J7L-DA-197-4' :
        
        #if cle_0 == '4V9F-0-25-56' and cle_1 == '2XD0-V-36-21' :
        #if cle_0 == '1FJG-A-48-8' and cle_1 == '1U9S-A-58-11' :
        print(cle_0)
        print(cle_1)
        graphe_global_commun_1, graphe_global_commun_2 = isomorphisme_extension_to_global(graphe_commun, graphe_ext_1, graphe_ext_2, graphe_global_1, graphe_global_2)
        createtikz_comp(graphe_global_1, graphe_global_2, graphe_global_commun_1, graphe_global_commun_2, cle_0+"_"+cle_1, "Graphes_globaux/graphes_sequence/taille_%d/"%i, i, cle_0, cle_1, sim[0], sim[1])
                    #print(cle)


''' modifie a partir de la fonction createTikz fournie par Vladimir utilisee dans CaRNAval
creation d'une figure tikz stockee dans un fichier pdf a partir des graphes passe en parametre 
(version superposition de deux graphes)
encore quelques soucis'''
def createtikz_comp(graphe_1, graphe_2, graphe_commun_1, graphe_commun_2, suffix, rep, i, cle_1, cle_2, sim, rang):

    tikzfile = rep+f'{suffix}.pgf'

    fout = open(tikzfile,"w")
    
    
    #fout.write("\\vspace*{-0.5cm}\n")
    fout.write("\\hspace*{-5cm}\n")
    fout.write("\\begin{minipage}{0.4\\textwidth}\n")
    
    fout.write("\\begin{center}\n")
    fout.write("\\begin{tikzpicture}[transform canvas={scale=0.5}]\n")

    compteur = 1
    maxi_coordonnees_y = 0
    mini_coordonnees_y = 0
    for n, data in graphe_1.nodes(data=True):
            
            
        
        #if data["type"] != -1 :
            x = data["coordonnees"][0]*1.25*3
            y = data["coordonnees"][1]*2*3
            
            if maxi_coordonnees_y < y :
                maxi_coordonnees_y = y
                
            if mini_coordonnees_y > y :
                mini_coordonnees_y = y
            
            color = "black"
            if n in graphe_commun_1.nodes() :
                ##pour noeuds "fictifs" qui correspondent aux noeuds des liaisons non can au milieu des helices can
                if 'motif' not in graphe_commun_1.nodes[n].keys() :
                    graphe_commun_1.nodes[n].update({'motif' : False})
                if graphe_commun_1.nodes[n]["motif"] :
                    color = "gray"
                else :
                    color = "orange"
    
            fout.write("  \\draw (\XScale{%s},\YScale{%s}) node[nt, label={[black]below right: %s}, label={[%s]center : %s}] (n-%s){} ;\n"%(x,y,n,color,data["nt"],n))
            compteur += 1
    
    fout.write("  \\draw (\XScale{%s},\YScale{%s}) node[label={[black]above: %s}] (title){} ;\n"%(0,maxi_coordonnees_y+2,cle_1))
    fout.write("\\draw (\XScale{%s},\YScale{%s}) node[label={[black]below: sim=%1.2f, rang=%d/4005}] (sim){} ;\n"%(0,mini_coordonnees_y-2.5,sim, rang))
    
    ##ajout des liaisons b53 qui manquent (qu on ne met pas pour les extensions)
    for n, data in graphe_1.nodes(data=True) :
        if n+1 in graphe_1.nodes() and (n,n+1) not in graphe_1.edges() :
            graphe_1.add_edge(n,n+1,label='B53', long_range=False)
            #if n in graphe_commun_1.nodes() and n+1 in graphe_commun_1.nodes() :
            #    graphe_commun_1.add_edge(n,n+1,label='B53', long_range=False)
        
        
    
    first = True
    print("gros rat")
    for e in graphe_1.edges(data=True):
            print(e)
            print(graphe_commun_1.edges.data())
            
            n1,n2 = e[0],e[1]

            attr1 = e[2]
        
        #if attr1["label"] != '0' :
        
            (x1,y1) = graphe_1.nodes[n1]["coordonnees"]
    
            (x2,y2) = graphe_1.nodes[n2]["coordonnees"]
    
            style = ""
            
            if (n1<n2 and not((graphe_1.nodes[n1]["motif"] == 1 and graphe_1.nodes[n2]["motif"] == 2) or \
                 (graphe_1.nodes[n1]["motif"] == 3 and graphe_1.nodes[n2]["motif"] == 4) or \
                 (graphe_1.nodes[n1]["motif"] == 1 and graphe_1.nodes[n2]["motif"] == 5)) or \
                 (graphe_1.nodes[n1]["motif"] == 2 and graphe_1.nodes[n2]["motif"] == 1) or \
                 (graphe_1.nodes[n1]["motif"] == 4 and graphe_1.nodes[n2]["motif"] == 3) or \
                 (graphe_1.nodes[n1]["motif"] == 5 and graphe_1.nodes[n2]["motif"] == 1))and \
                 attr1["label"] != '0' :
    
                if attr1["label"] != "B53":
    
                    if attr1['long_range']:
    
                        style += ",lr"
    
                    else:
    
                        style += ",bp"
    
                if not neighbors(x1,x2) and attr1["label"] != "B53":
    
                    if attr1["label"] == 'TSS' and first :
                        style +=  ",bend right=10"
                        first = False
                    elif abs(x1-x2) >= 1 or abs(y1-y2) >= 1: 
                        style += ",bend left=10"
                
                if (n1,n2) in graphe_commun_1.edges() :
                    motif = False
                    for edge in graphe_commun_1[n1][n2] :
                        print(e)
                        if graphe_commun_1[n1][n2][edge]["motif"] :
                            motif = True
                    if motif :
                        style += ",color=lightgray"
                    else :
                        style += ",color=orange"
                print(style)
                fout.write("  \\draw (n-%s) edge[%s] (n-%s);\n"%(n1,attr1["label"]+style,n2))
                

    fout.write("\\end{tikzpicture}\n")
    fout.write("\\end{center}\n")
    fout.write("\\end{minipage}\n")
    fout.write("\\hspace{0.5cm}\n")
    
    fout.write("\\begin{minipage}{0.4\\textwidth}\n")
    fout.write("\\begin{center}\n")
    fout.write("\\begin{tikzpicture}[transform canvas={scale=0.5}]\n")

    compteur = 1
    maxi_coordonnees_y = 0
    for n, data in graphe_2.nodes(data=True):
            
        #if data["type"] != -1 :
            x = data["coordonnees"][0]*1.25*3
            y = data["coordonnees"][1]*2*3
            
            if maxi_coordonnees_y < y :
                maxi_coordonnees_y = y
            
            color = "black"
            if n in graphe_commun_2.nodes() :
                print(n)
                print(graphe_commun_2.nodes[n])
                ##pour noeuds "fictifs" qui correspondent aux noeuds des liaisons non can au milieu des helices can
                if 'motif' not in graphe_commun_2.nodes[n].keys() :
                    graphe_commun_2.nodes[n].update({'motif' : False})
                    
                if graphe_commun_2.nodes[n]["motif"] :
                    color = "gray"
                else :
                    color = "orange"
    
            fout.write("  \\draw (\XScale{%s},\YScale{%s}) node[nt, label={[black]below right: %s}, label={[%s]center :%s}] (n-%s){} ;\n"%(x,y,n,color,data["nt"],n))
            compteur += 1
    fout.write("  \\draw (\XScale{%s},\YScale{%s}) node[label={[black]above:%s}] (title){} ;\n"%(0,maxi_coordonnees_y+2,cle_2))
    ##ajout des liaisons b53 qui manquent (qu on ne met pas pour les extensions)
    for n, data in graphe_2.nodes(data=True) :
        if n+1 in graphe_2.nodes() and (n,n+1) not in graphe_2.edges() :
            graphe_2.add_edge(n,n+1,label='B53', long_range=False)
            #if n in graphe_commun_2.nodes() and n+1 in graphe_commun_2.nodes() :
            #    graphe_commun_2.add_edge(n,n+1,label='B53', long_range=False)
        
    
    first = True
    for e in graphe_2.edges(data=True):
            #print(e)
            n1,n2 = e[0],e[1]

            attr1 = e[2]
        
        #if attr1["label"] != '0' :
        
            (x1,y1) = graphe_2.nodes[n1]["coordonnees"]
    
            (x2,y2) = graphe_2.nodes[n2]["coordonnees"]
    
            style = ""
            
            if (n1<n2 and not((graphe_2.nodes[n1]["motif"] == 1 and graphe_2.nodes[n2]["motif"] == 2) or \
                 (graphe_2.nodes[n1]["motif"] == 3 and graphe_2.nodes[n2]["motif"] == 4) or \
                 (graphe_2.nodes[n1]["motif"] == 1 and graphe_2.nodes[n2]["motif"] == 5)) or \
                 (graphe_2.nodes[n1]["motif"] == 2 and graphe_2.nodes[n2]["motif"] == 1) or \
                 (graphe_2.nodes[n1]["motif"] == 4 and graphe_2.nodes[n2]["motif"] == 3) or \
                 (graphe_2.nodes[n1]["motif"] == 5 and graphe_2.nodes[n2]["motif"] == 1))and \
                 attr1["label"] != '0' :
    
                if attr1["label"] != "B53":
    
                    if attr1['long_range']:
    
                        style += ",lr"
    
                    else:
    
                        style += ",bp"
    
                if not neighbors(x1,x2) and attr1["label"] != "B53":
    
                    if attr1["label"] == 'TSS' and first :
                        style +=  ",bend right=10"
                        first = False
                    elif abs(x1-x2) >= 1 or abs(y1-y2) >= 1: 
                        style += ",bend left=10"
                
#                 attr_sans_motif = e[2]
#                 del(attr_sans_motif["motif"])
                
                if (n1,n2) in graphe_commun_2.edges() :
                    motif = False
                    for edge in graphe_commun_2[n1][n2] :
                        if graphe_commun_2[n1][n2][edge]["motif"] :
                            motif = True
                    if motif :
                        style += ",color=lightgray"
                    else :
                        style += ",color=orange"
                
                fout.write("  \\draw (n-%s) edge[%s] (n-%s);\n"%(n1,attr1["label"]+style,n2))
                

    fout.write("\\end{tikzpicture}\n")
    fout.write("\\end{center}\n")
    fout.write("\\end{minipage}\n")
    fout.close()

    srcfile = rep+f'{suffix}.tex'

    fout = open(srcfile,"w")


    fout.write(r"\documentclass{article}[20pt]")

    #fout.write("\\documentclass[tikz,border=10pt]{standalone}\n")

    fout.write("\\input{%s}\n"%TIKZ_HEADER)

    fout.write(r"""\geometry{paperwidth=%dcm,paperheight=%1.2fcm,margin=0cm}
            \tikzset{
  font={\fontsize{23pt}{12}\selectfont}}
            \begin{document}
                \begin{figure}[p]
                
                \centering"""%(55, 50))



    #fout.write("\\begin{document}\n")

    fout.write("\\input{%s}\n"%tikzfile)
    
    fout.write(r"\end{figure}")

    fout.write("\\end{document}\n")

    fout.write("\\input{%s}\n"%TIKZ_FOOTER)

    fout.close()

    #add_number_pgf(tikzfile)


    subprocess.call(["/usr/bin/pdflatex",'-output-directory', "Graphes_globaux/graphes_sequence/taille_%d/"%i, srcfile,"--quiet"])

    # subprocess.call(["convert", '-density' ,'300', srcfile[:-3] + 'pdf', srcfile[:-3] + 'png'])#imagemagick

    os.remove(srcfile.replace(".tex",".log"))

    os.remove(srcfile.replace(".tex",".aux"))          

'''lancement de la creation des figures tikz pour toutes les superpositions'''
def test():
    for i in range(9,10) :
        with open(EXTENSION_PATH%i+"dico_comp_complet_metrique_toutes_aretes_coeff_all1_taille_%d.pickle"%i, 'rb') as fichier_comp :
            mon_depickler = pickle.Unpickler(fichier_comp)
            dico_comp = mon_depickler.load()
            
            with open(EXTENSION_PATH%i+"sim_extensions_toutes_aretes_coeff_all1_taille_%d.pickle"%i, 'rb') as fichier_sim :
                mon_depickler = pickle.Unpickler(fichier_sim)
                dico_sim = mon_depickler.load()
                
                for cle in dico_comp.keys() :
                    cle_sim_1 = cle[0][8:]
                    cle_sim_2 = cle[1][8:]
                    if (cle_sim_1, cle_sim_2) in dico_sim.keys() :
                        rang = rang_sim(dico_sim, dico_sim[(cle_sim_1, cle_sim_2)])
                        draw_iso_ext_to_global(dico_comp[cle], cle, EXTENSION_PATH_TAILLE%i, i, (dico_sim[(cle_sim_1, cle_sim_2)],rang))
                    else :
                        rang = rang_sim(dico_sim, dico_sim[(cle_sim_2, cle_sim_1)])
                        draw_iso_ext_to_global(dico_comp[cle], cle, EXTENSION_PATH_TAILLE%i, i, (dico_sim[(cle_sim_2, cle_sim_1)],rang))

'''(non utilise) issu du code fourni par Vladimir utilise dans CaRNAval'''    
def nodeToNeighborsMap(g,coords):

    """Builds a node to neighbor mapping"""

    # print('coords',coords)

    node2neighbors = {n:[] for n in coords}

    for e in g.edges(data=True):

        n1,n2 = e[0],e[1]

        attr = e[2]

        # print('n1,n2,attr',n1,n2,attr)

        node2neighbors[n1].append((n2,attr["label"]))

        node2neighbors[n2].append((n1,attr["label"]))

    return node2neighbors

'''(non utilise) issu du code fourni par Vladimir utilise dans CaRNAval'''
def getMajority(votes):

    """

        Returns the element having max #occ. in the list.

        Returns the smallest element in case of a draw.

    """

    counts = {i:0 for i in votes}

    for i in  votes:

        counts[i]+=1

    l = [(counts[i],i) for i in counts]

    l.sort(reverse=True)

    wc,w = l[0]

    return w

'''(non utilise) issu du code fourni par Vladimir utilise dans CaRNAval'''
def computeCoords(g,strands,order):

    # print("computeCoords(g,strands,order)",strands,order)

    coords = {}

    # for s in [0,1]:#TAG

    for s in range(0,len(strands)):

        for i in range(len(strands[s])):

            o = order[s][i]

            x = abs(o)

            d = o/x

            by = 0

            if d<0:

                by = len(strands[s][i])-1

            for j in range(len(strands[s][i])):

                y = by+j*d

                n = strands[s][i][j]

                coords[n] = (x,y)

    return coords

'''(non utilise) issu du code fourni par Vladimir utilise dans CaRNAval'''
def doesCross(n1,n2,m1,m2,node2strand):

    """

        Check whether two base-pairs cross in the current strand ordering/orientations

    """

    if node2strand[n1] == node2strand[m2]:

        m1,m2 = m2,m1

    if node2strand[n1] == node2strand[m1] and node2strand[n2] == node2strand[m2] :

       return (-1*(n1-m1)*(n2-m2)*node2strand[n1]*node2strand[n2] > 0)

    return False

'''(non utilise) issu du code fourni par Vladimir utilise dans CaRNAval'''
def conflicts(orderStrands,strands,g):

    """

        Computes and returns a quality score for the current strand ordering,

        defined as a linear combination of the #canonical pairs on

        adjacent strands (+), the #non-canonical pairs on adjacent

        strands (+), and the #crossing pairs of interactions on adjacent strands.

    """

    node2strand = {}

    for s in strands.keys():#TAG

        for i in range(len(orderStrands[s])):

            # print(i)

            for n in  strands[s][i]:

                node2strand[n] = orderStrands[s][i]

    numCrossing = 0

    WCCNeighbors = 0

    NCNeighbors = 0

    for e in g.edges(data=True):

        n1,n2 = e[0],e[1]

        attr1 = e[2]

        if neighbors(node2strand[n1],node2strand[n2]) and n1<n2:

            if attr1["label"] == "CWW":

                WCCNeighbors += 1

            else:

                NCNeighbors += 1

        if n1<n2 and attr1["label"] != "B53":

            for f in g.edges(data=True):

                attr2 = e[2]

                m1,m2 = f[0],f[1]

                if e!=f and m1<m2 and attr2["label"] != "B53":

                    if doesCross(n1,n2,m1,m2,node2strand):

                        numCrossing += 1

    return (2*WCCNeighbors+2*NCNeighbors-numCrossing,WCCNeighbors,numCrossing)

'''(non utilise) issu du code fourni par Vladimir utilise dans CaRNAval'''
def buildStrands(g):

    """

        Builds the strands, i.e. the sequences of positions connected by

        backbone connections, and regroups them within a two level structure.

        First level: Left/right of long range interaction

        Second level: list of strands

        TODO (for the brave only):

        1) Distinguish between strands and hairpin loops

        2) Special treatment for isolated nucleotides

    """

    children = {n:[] for n in g.nodes()}

    nonRoots = {}

    # strands = {0:[],1:[]}#TAG

    strands={}

    print("Nodes data",g.nodes(data=True))

    print("Edges data",g.edges(data=True))

    s = sorted(set(n[1]['part_id'] for n in g.nodes(data=True)))

    sides = {n[0]:s.index(n[1]['part_id']) for n in g.nodes(data=True)}

    print("s,sides",s,sides)

    for n in g.nodes(data=True):

        n1 = n[0]

        data = n[1]

    for e in g.edges(data=True):

        n1,n2 = e[0],e[1]

        attr = e[2]

        if attr["label"] == "B53":

            children[n1].append(n2)

            nonRoots[n2] = 0

    for n in g.nodes(data=True):

        n1 = n[0]

        print("n,n1,nonRoots",n,n1,nonRoots)

        if n1 not in nonRoots:

            s = sides[n1]

            print('n1,s',n1,s)

            nlist = [n1]

            l = children[n1][:]

            nprev = n1

            while len(l)>0:

                nnext = l.pop()

                nlist.append(nnext)

                l = l + children[nnext]

                nprev = nnext

            tmp=strands.get(s,[])

            tmp.append(nlist)

            strands[s]=tmp

    #We want the keys to be consecutive so we perform a simple permutation from {keys} to {0,1,...}

    continuous_strands={}

    i=0

    for key in sorted(strands.keys()):

        continuous_strands[i]=strands[key]

        i+=1

    return continuous_strands

'''(non utilise) issu du code fourni par Vladimir utilise dans CaRNAval'''
def placeStrands(g,strands,motif):

    """

        Brute force computation of the strand ordering/orientation having maximum

        quality score.

        Warning: Complexities are in k1!*k2!*M^2, where k1,k2 are #strands on both

        sides of the long-range interaction, and M is the maximum #nucleotides on a strand.

        Use at your own risk!

    """

    # Version de Yann pour seulement deux strands

    # doit etre modifiee pour gerer un nombre arbitraire de strands

    # pas de soucis d'efficacite pour l'instant

    ls = []

    n_tot=0

    for i in range(0,len(strands)):
    
        n = len(strands[i])

        la1 = list(itertools.product([-1,1],repeat=n))

        lb1 = list(itertools.permutations(range(n_tot+1,n_tot+n+1)))

        n_tot+=n

        #zip :  returns a list of tuples, where the i-th tuple contains the i-th element from each of the argument sequences or iterables.

        l1 = [[a*b for a,b in zip(x,y)] for x,y in list(itertools.product(la1,lb1))]

        # print('l1',l1)

        ls.append(l1)

    # print(ls)

    ls = tuple(ls)

    # print(list(itertools.product(*ls)))

    l = [{i:a for i,a in enumerate(o)} for o in list(itertools.product(*ls))]

    #tout ce qui a precede vise a produire l'ensemble des possibilites

    #les lignes suivantes vont evaluer ces possibilites

    # print('l#',l)

    f = [(conflicts(orderStrands,strands,g),orderStrands) for orderStrands in l]

    # print('f#',f)

    # [((8, 2, 2), {0: [-1], 1: [-2]}), ((10, 2, 0), {0: [-1], 1: [2]}), ((10, 2, 0), {0: [1], 1: [-2]}), ((8, 2, 2), {0: [1], 1: [2]})]

    #faut dire au sort de ne considerer que le tuple et pas le dictionnaire qui va avec

    #putain d'internet fini a la pisse



    f.sort(reverse=True,key=operator.itemgetter(0))

    return f[0][1]

'''(non utilise) issu du code fourni par Vladimir utilise dans CaRNAval'''
def shiftHeight(g,strands,order):

    """

        Shifts strands uniformly on the y-axis in order to draw base-pairs

        horizontally.

        TODO: Shift outliers, stretching the backbone.

    """

    coords = computeCoords(g,strands,order)

    node2neighbors = nodeToNeighborsMap(g,coords)

    lx = {coords[n][0]:[] for n in coords}

    for n in coords:

        x,y = coords[n]

        lx[x].append(n)

    for i in sorted(lx.keys()):

        print("  Strand",i)

        votes = []

        for n in lx[i]:

            x,y = coords[n]

            for nn,t in node2neighbors[n]:

                xn,yn = coords[nn]

                if neighbors(xn,x) and xn<x:

                    if t != "B53":

                        votes.append(yn-y)

                    votes.append(yn-y)

        print("  ",votes)

        if len(votes)>0:

            votes.sort()

            delta = getMajority(votes)

            print("  delta=",delta)

            for n in lx[i]:

                x,y = coords[n]

                coords[n] = (x,y+delta)

    return coords

'''(non utilise) issu du code fourni par Vladimir utilise dans CaRNAval'''
def createTikz(i,suffix,g,coords):

    """

        Compiles motif layouts into TikZ and compiles them using PDFLatex

    """

    num2Nodes = {n[0]:n for n in g.nodes(data=True)}

    sort = sorted(g.nodes())

    print(sort)

    tikzfile = f'motif-{i+1}{suffix}.pgf'

    fout = open(tikzfile,"w")

    fout.write("\\begin{tikzpicture}\n")

    for n in coords:

        (x,y) = coords[n]

        fout.write("  \\node[nt] (n-%s) at (\XScale{%s},\YScale{%s}) {\\tiny %s};\n"%(n,x,y,sort.index(n)+1))

    for e in g.edges(data=True):

        n1,n2 = e[0],e[1]

        attr1 = e[2]

        (x1,y1) = coords[n1]

        (x2,y2) = coords[n2]

        style = ""

        if n1<n2:

            if attr1["label"] != "B53":

                if attr1['long_range']:

                    style += ",lr"

                else:

                    style += ",bp"

            if not neighbors(x1,x2) and attr1["label"] != "B53":

                style += ",bend left=10"

            fout.write("  \\draw (n-%s) edge[%s] (n-%s);\n"%(n1,attr1["label"]+style,n2))

    fout.write("\\end{tikzpicture}\n")

    fout.close()

    srcfile = f'motif-{i+1}{suffix}.tex'

    fout = open(srcfile,"w")



    fout.write(r"\documentclass{article}")

    #fout.write("\\documentclass[tikz,border=10pt]{standalone}\n")

    fout.write("\\input{%s}\n"%TIKZ_HEADER)

    fout.write(r"""\geometry{

            paperwidth=8cm,

            paperheight=8cm,

            margin=0cm

            }

            \begin{document}

                \begin{figure}[p]

                \centering""")



    #fout.write("\\begin{document}\n")

    fout.write("\\input{%s}\n"%tikzfile)

    fout.write(r"\end{figure}")

    fout.write("\\end{document}\n")

    fout.write("\\input{%s}\n"%TIKZ_FOOTER)

    fout.close()

''' modifie a partir de la fonction createTikz fournie par Vladimir utilisee dans CaRNAval
creation d'une figure tikz stockee dans un fichier pdf a partir du graphe passe en parametre 
'''
def createTikz_modif(suffix,graphe, rep, i):

    """

        Compiles motif layouts into TikZ and compiles them using PDFLatex

    """

    tikzfile = rep+f'{suffix}.pgf'

    fout = open(tikzfile,"w")

    fout.write("\\begin{tikzpicture}[inner sep = 0mm]\n")

    compteur = 1
    for n, data in graphe.nodes(data=True):
        
        #if data["type"] != -1 :
            x = data["coordonnees"][0]*1.25
            y = data["coordonnees"][1]*2
    
            #fout.write("  \\draw (\XScale{%s},\YScale{%s}) node[nt, label={[black]below right: \\tiny %s}] (n-%s){\\tiny %s} ;\n"%(x,y,n,n,data["nt"]))
            fout.write("  \\draw (\XScale{%s},\YScale{%s}) node[nt, label={[black]below right: \\tiny %s}] (n-%s){} ;\n"%(x,y,n,n))

            compteur += 1
    
    ##ajout des liaisons b53 qui manquent (qu on ne met pas pour les extensions)
    for n, data in graphe.nodes(data=True) :
        if n+1 in graphe.nodes() and (n,n+1) not in graphe.edges() :
            graphe.add_edge(n,n+1,label='B53', long_range=False)
        
    
    first = True
    for e in graphe.edges(data=True):
            print(e)
            n1,n2 = e[0],e[1]

            attr1 = e[2]
        
        #if attr1["label"] != '0' :
        
            (x1,y1) = graphe.nodes[n1]["coordonnees"]
    
            (x2,y2) = graphe.nodes[n2]["coordonnees"]
    
            style = ""
            print("petit rat")
            print(graphe.nodes[n1])
            print(graphe.nodes[n2])
            
            if (n1<n2 and not((graphe.nodes[n1]["motif"] == 1 and graphe.nodes[n2]["motif"] == 2) or \
                 (graphe.nodes[n1]["motif"] == 3 and graphe.nodes[n2]["motif"] == 4) or \
                 (graphe.nodes[n1]["motif"] == 1 and graphe.nodes[n2]["motif"] == 5)) or \
                 (graphe.nodes[n1]["motif"] == 2 and graphe.nodes[n2]["motif"] == 1) or \
                 (graphe.nodes[n1]["motif"] == 4 and graphe.nodes[n2]["motif"] == 3) or \
                 (graphe.nodes[n1]["motif"] == 5 and graphe.nodes[n2]["motif"] == 1))and \
                 attr1["label"] != '0' :
    
                if attr1["label"] != "B53":
    
                    if attr1['long_range']:
    
                        style += ",lr"
    
                    else:
    
                        style += ",bp"
    
                if not neighbors(x1,x2) and attr1["label"] != "B53":
    
                    if attr1["label"] == 'TSS' and first :
                        style +=  ",bend right=10"
                        first = False
                    elif abs(x1-x2) >= 1 or abs(y1-y2) >= 1: 
                        style += ",bend left=10"
    
                fout.write("  \\draw (n-%s) edge[%s] (n-%s);\n"%(n1,attr1["label"]+style,n2))
                

    fout.write("\\end{tikzpicture}\n")

    fout.close()

    srcfile = rep+f'{suffix}.tex'

    fout = open(srcfile,"w")



    fout.write(r"\documentclass{article}[10pt]")

    #fout.write("\\documentclass[tikz,border=10pt]{standalone}\n")

    fout.write("\\input{%s}\n"%TIKZ_HEADER)

    fout.write(r"""\geometry{paperwidth=%dcm,paperheight=%1.2fcm,margin=0cm}
            
            \begin{document}
                \begin{figure}[p]
                \caption{%s}
                
                \centering"""%(8, 8+i+1.75, suffix))



    #fout.write("\\begin{document}\n")

    fout.write("\\input{%s}\n"%tikzfile)

    fout.write(r"\end{figure}")

    fout.write("\\end{document}\n")

    fout.write("\\input{%s}\n"%TIKZ_FOOTER)

    fout.close()

    #add_number_pgf(tikzfile)


    subprocess.call(["/usr/bin/pdflatex",'-output-directory', "/media/coline/Maxtor/", srcfile,"--quiet"])

    # subprocess.call(["convert", '-density' ,'300', srcfile[:-3] + 'pdf', srcfile[:-3] + 'png'])#imagemagick

#     os.remove(srcfile.replace(".tex",".log"))
# 
#     os.remove(srcfile.replace(".tex",".aux"))

''' modifie a partir de la fonction createTikz fournie par Vladimir utilisee dans CaRNAval
creation d'une figure tikz stockee dans un fichier pdf a partir du graphe passe en parametre
(version graphe commun a un groupe) 
'''
def createTikz_modif_commun(suffix,graphe, rep, i):

    """

        Compiles motif layouts into TikZ and compiles them using PDFLatex

    """
    print(graphe.nodes.data())
    tikzfile = rep+f'{suffix}.pgf'

    fout = open(tikzfile,"w")

    fout.write("\\begin{tikzpicture}[inner sep = 0mm]\n")

    compteur = 1
    for n, data in graphe.nodes(data=True):
        
        #if data["type"] != -1 :
            x = data["coordonnees"][0]*1.25*5
            y = data["coordonnees"][1]*2*5
    
            fout.write("  \\draw (\XScale{%s},\YScale{%s}) node[nt] (n-%s){} ;\n"%(x,y,n))
            compteur += 1
    
    ##ajout des liaisons b53 qui manquent (qu on ne met pas pour les extensions)
    for n, data in graphe.nodes(data=True) :
        if n+1 in graphe.nodes() and (n,n+1) not in graphe.edges() :
            graphe.add_edge(n,n+1,label='B53', long_range=False)
        
    
    first = True
    for e in graphe.edges(data=True):
            print(e)
            n1,n2 = e[0],e[1]

            attr1 = e[2]
        
        #if attr1["label"] != '0' :
        
            (x1,y1) = graphe.nodes[n1]["coordonnees"]
    
            (x2,y2) = graphe.nodes[n2]["coordonnees"]
    
            style = ""
            print("petit rat")
            print(graphe.nodes[n1])
            print(graphe.nodes[n2])
            
            if n1<n2 and \
                 attr1["label"] != '0' :
    
                if attr1["label"] != "B53":
    
                    if attr1['long_range']:
    
                        style += ",lr"
    
                    else:
    
                        style += ",bp"
    
                if not neighbors(x1,x2) and attr1["label"] != "B53":
    
                    if attr1["label"] == 'TSS' and first :
                        style +=  ",bend right=10"
                        first = False
                    elif abs(x1-x2) >= 1 or abs(y1-y2) >= 1: 
                        style += ",bend left=10"
    
                fout.write("  \\draw (n-%s) edge[%s] (n-%s);\n"%(n1,attr1["label"]+style,n2))
                

    fout.write("\\end{tikzpicture}\n")

    fout.close()

    srcfile = rep+f'{suffix}.tex'

    fout = open(srcfile,"w")



    fout.write(r"\documentclass{article}[10pt]")

    #fout.write("\\documentclass[tikz,border=10pt]{standalone}\n")

    fout.write("\\input{%s}\n"%TIKZ_HEADER)

    fout.write(r"""\geometry{paperwidth=%dcm,paperheight=%1.2fcm,margin=0cm}
            
            \begin{document}
                \begin{figure}[p]
                \caption{%s}
                
                \centering"""%(50, 50, suffix))



    #fout.write("\\begin{document}\n")

    fout.write("\\input{%s}\n"%tikzfile)

    fout.write(r"\end{figure}")

    fout.write("\\end{document}\n")

    fout.write("\\input{%s}\n"%TIKZ_FOOTER)

    fout.close()

    #add_number_pgf(tikzfile)


    subprocess.call(["/usr/bin/pdflatex",'-output-directory', EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes", srcfile,"--quiet"])

    # subprocess.call(["convert", '-density' ,'300', srcfile[:-3] + 'pdf', srcfile[:-3] + 'png'])#imagemagick

    os.remove(srcfile.replace(".tex",".log"))

    os.remove(srcfile.replace(".tex",".aux"))


TIKZ_HEADER="Graphes_globaux/graphes_sequence/taille_%d/"%4+"header.inc.tex"

TIKZ_FOOTER="Graphes_globaux/graphes_sequence/taille_%d/"%4+"footer.inc.tex"


# def recherche_categories(graphe):
#     chaines_1 = [[1]]
#     for i in range(1,5) :
#             compteur = i
#             if i != 1 : chaines_1.append([i])
#             liaison_B53 = True
#             while liaison_B53 :
#                 liaison_B53 = False
#                 temp = compteur
#                 for voisin in graphe.successors(compteur) :
#                     for arc in graphe[compteur][voisin] :
#                         if voisin not in [1,2,3,4] and voisin not in chaines_1[len(chaines_1)-1] and graphe[compteur][voisin][arc]["label"] == 'B53' :
#                             liaison_B53 = True
#                             temp = voisin
#                             chaines_1[len(chaines_1)-1].append(voisin)
#                             
#                 for voisin in graphe.predecessors(compteur) :
#                     for arc in graphe[voisin][compteur] :
#                         if voisin not in [1,2,3,4] and voisin not in chaines_1[len(chaines_1)-1] and graphe[voisin][compteur][arc]["label"] == 'B53' :
#                             liaison_B53 = True
#                             temp = voisin
#                             chaines_1[len(chaines_1)-1].append(voisin)
#                 compteur = temp
#     lien_chaine = False            
#     for i in range(len(chaines_1)) :
#         for elt in chaines_1[i] :
#             for voisin in graphe[elt] :
#                 if graphe.nodes[voisin]
        
        
if __name__ == '__main__':
    #with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_graphe :    
    #with open("/media/coline/Maxtor/dico_new_100320_sim_par_branche_0.65.pickle", 'rb') as fichier_graphe :    
    with open("dico_algo_heuristique_new_v_33bfb11.pickle", 'rb') as fichier_graphe :    
    
        mon_depickler = pickle.Unpickler(fichier_graphe)
        dico_graphe = mon_depickler.load() 
        
        with open("grands_graphes_new_data_taille_4.pickle", 'rb') as fichier :
            mon_depickler = pickle.Unpickler(fichier)
            dico_graphe_global = mon_depickler.load()
            #liste_a_tester = [(('4y4o', 23), ('5afi', 24)),(('4v67', 7), ('4u3u', 7)),(('4ybb', 12), ('4u3u', 7)),(('4ybb', 12), ('6ek0', 7)),(('4v67', 7), ('6ek0', 7)),(('4y4o', 23), ('4u27', 3)), (('5wfs', 19), ('4u3u', 7)),(('5wfs', 19), ('6ek0', 7)),(('2zjr', 3), ('4u4r', 18))]
            liste_a_tester = [(('4ybb', 6), ('1vq8', 19))]

            for couple in liste_a_tester :
                for elt in dico_graphe.keys() :
                    if elt[0] in couple and elt[1] in couple :
                        with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt[0][0]+"_"+str(elt[0][1])+"_2.pickle", 'rb') as fichier_graphe1 :
                            mon_depickler1 = pickle.Unpickler(fichier_graphe1)
                            graphe1 = mon_depickler1.load()
                                
                        with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt[1][0]+"_"+str(elt[1][1])+"_2.pickle", 'rb') as fichier_graphe2 :
                            mon_depickler2 = pickle.Unpickler(fichier_graphe2)
                            graphe2 = mon_depickler2.load()
                            print(type(dico_graphe_global[elt[0]][0]))
                            graphe_commun_global_1, graphe_commun_global_2 = isomorphisme_extension_to_global_version_new_data(dico_graphe[elt]["graphe"], graphe1, graphe2, dico_graphe_global[elt[0]][0], dico_graphe_global[elt[1]][0])
                            print("ramou")
                            print(dico_graphe_global[elt[0]][0].edges.data())
                            print(dico_graphe_global[elt[1]][0].edges.data())
                            print("ramousnif")
                            sim = round(calcul_sim_aretes_avec_coeff(graphe1, graphe2, dico_graphe[elt]["graphe"], "petit rat", 1, 1, 1),2)
                            print(graphe_commun_global_1.nodes.data())
                            print(graphe_commun_global_1.edges.data())
                            print(graphe_commun_global_2.nodes.data())
                            print(graphe_commun_global_2.edges.data())
                            draw_isomorphism_struct(elt[0], elt[1], graphe_commun_global_1, graphe_commun_global_2, dico_graphe_global[elt[0]][0], dico_graphe_global[elt[1]][0], sim)
                            
                            #exit()
    
    exit()
    nom_fichier ="1c2w_1"
    
    
    with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_avec_coord.pickle"%nom_fichier, 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        graphe = mon_depickler.load()
        
        with open("Graphs/"+nom_fichier.split("_")[0]+".pickle", 'rb') as fichier_pickle_grand :
            mon_depickler_grand = pickle.Unpickler(fichier_pickle_grand)
            graphe_grand = mon_depickler_grand.load()
            
            liste_a_garder = []
            for noeud, data in graphe.nodes(data=True) :
                for i in range(data["position"][0], data["position"][1]+1) :
                    liste_a_garder.append(i)
            
            sous_graphe = graphe_grand.subgraph(liste_a_garder)
        
            createTikz_modif(nom_fichier, sous_graphe, "", 4)
#     
#     with open("grands_graphes_taille_2.pickle", 'rb') as fichier :
#         mon_depickler = pickle.Unpickler(fichier)
#         dico_graphes = mon_depickler.load()
#         
#         
#         for cle in dico_graphes.keys() :
#             print(cle)
#             print(dico_graphes[cle].nodes.data())
#             strands = buildStrands(dico_graphes[cle])
#             print(strands)
#             orderStrands = placeStrands(dico_graphes[cle],strands)
#             print("ramou")
#             print(orderStrands)
#             #orderStrands = {0: [1,2,3,4], 1: [5,6,7], 2: [-1,-2], 3: [-3,-4], 4: [8]}
#             coords = shiftHeight(dico_graphes[cle],strands,orderStrands)
#             createTikz(0,"test",dico_graphes[cle],coords)
#             break
    
    #subprocess.call(["/usr/bin/pdflatex",'-output-directory', "Graphes_globaux/graphes_sequence/taille_%d/"%8, '/home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/Graphes_globaux/graphes_sequence/taille_8/5J7L-DA-50-21-taille-8.tex',"--quiet"])
#     
#     for i in range(8, 11) :
#         for fic in os.listdir("Graphes_globaux/graphes_sequence/taille_%d"%i) :
#             if "pickle" in fic :
#                 with open("Graphes_globaux/graphes_sequence/taille_%d/"%i+fic, 'rb') as fichier_pickle :
#                     mon_depickler_coord = pickle.Unpickler(fichier_pickle)
#                     graphe = mon_depickler_coord.load()
#                     cle = fic.split("_")[2] + "-" +   fic.split("_")[3] + "-"+ fic.split("_")[4] + "-" + fic.split("_")[5]
#                     createTikz_modif(cle+"-taille-%d"%i, graphe, "Graphes_globaux/graphes_sequence/taille_%d/"%i, i)
    
    #test()
    
#     with open(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes"+"graphe_commun_toutes_aretes_max_0.7_groupe_clustering_perez_12_taille_8_avec_coord.pickle", 'rb') as fichier_pickle :
#         mon_depickler_coord = pickle.Unpickler(fichier_pickle)
#         graphe = mon_depickler_coord.load()
#         #cle = fic.split("_")[2] + "-" +   fic.split("_")[3] + "-"+ fic.split("_")[4] + "-" + fic.split("_")[5]
#         createTikz_modif_commun("test", graphe, EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes", 8)
    
    