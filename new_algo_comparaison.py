'''
Created on 29 août 2019

@author: coline

Nouvel algo pour comparer les graphes d'extension deux à deux (Methode Branch and Cut)
(version toutes donnees PDB)
'''
import recup_data.clustering_perez
from recup_data.extension_new import obtenir_extension_un_elt
from recup_data.clustering_perez import algo_principal
from recup_data import graphe_commun_clusters
liste_pbs = [('6az1', 5), ('6ek0', 9), ('4wqf', 14), ('6ek0', 1), ('6ek0', 8), ('6ek0', 9), ('6az1', 4), ('5wdt', 11), ('6gyv', 1), ('4y4o', 18), ('5dm6', 6), ('4ybb', 27), ('5ibb', 34), ('1vqo', 21), ('5afi', 15), ('6h4n', 23), ('4v9d', 41), ('1k8a', 8), ('4u4r', 28), ('4u3u', 5), ('6ek0', 11), ('2nz4', 3), ('4ena', 1), ('3t1y', 2)]
import networkx as nx
import os
import pickle
import json
import time
import copy
import matplotlib.pyplot as plt
from recup_data.constantes import EXTENSION_PATH, EXTENSION_PATH_TAILLE,\
    NEW_EXTENSION_PATH_TAILLE, deux_chaines, PATH_MMCIF
import multiprocessing
import csv
import seaborn as sns
import numpy as np

''' Partie liee a l'algo'''

def chaines_reliees(graphe):
    ''' renvoie les chaines du graphe qui sont reliees soit par liaison covalente soit par liaison hydrogene '''
    paires_chaines = []
    for noeud,data in graphe.nodes(data=True) :
        if data["type"] != None and data["type"] != '0' and noeud not in [1,2,3,4,5] :
            if len(data["chaine"]) >  1 :
                if data["chaine"] not in paires_chaines :
                    paires_chaines.append([e for e in data["chaine"]if e != 5])
            
            for voisin in graphe[noeud] :
                for edge in graphe[noeud][voisin] :
                    if graphe[noeud][voisin][edge]["label"] != 'B53':
                        for elt1 in data["chaine"] :
                            for elt2 in graphe.nodes[voisin]["chaine"] :
                                if [elt1, elt2] not in paires_chaines and [elt2, elt1] not in paires_chaines and elt1 != elt2 and elt1 != 5 and elt2 != 5 :
                                    paires_chaines.append([elt1, elt2])
                    
    return paires_chaines

'''idem que dans calcul_sim '''
def calcul_aretes_avec_coeff(graphe, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
    
    for u,v,data in graphe.edges(data=True) :
        if data["label"] != 'B53' :
            if data["label"] == '0' :
                if coeffa == 1 :
                    somme_aretes += graphe.nodes[u]["poids"]
                    
            elif graphe.nodes[u]["type"] == 1 and graphe.nodes[v]["type"] == 1 :
                if coeffc == 1 :
                    somme_aretes += graphe.nodes[u]["poids"]
                    
            else :
                if coeffn == 1 :
                    somme_aretes += graphe.nodes[u]["poids"]
                    
    somme_aretes = somme_aretes/2 - 4
#     print(somme_aretes)
    
    return somme_aretes

'''idem que dans calcul_sim '''
def calcul_aretes_communes_avec_coeff(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
    for u,v,data in graphe_commun.edges(data=True) :
        if data["label"] != 'B53' :
            if data["label"] == '0' :
                if coeffa == 1 :
                    print(u,v)
                    print(graphe1.nodes[u[0]]["poids"])
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
#                     print(u,v)
#                     print(min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]))
#                     print(somme_aretes)
                    #print((u,v))
            elif graphe1.nodes[u[0]]["type"] == 1 and graphe1.nodes[v[0]]["type"] == 1 :
                if coeffc == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
#                     print(u,v)
#                     print(min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]))
#                     print(somme_aretes)
                    #print((u,v))
            else :
                if coeffn == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
#                     print(u,v)
#                     print(min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) )
#                     print(somme_aretes)
                    #print((u,v))
#     print("rapoulou")
#     print(somme_aretes/2)     
             
    somme_aretes = somme_aretes/2 - 4
#     print(somme_aretes)  
    return somme_aretes
       
'''idem que dans calcul_sim '''
def calcul_sim_aretes_avec_coeff(graphe1, graphe2, graphe_commun, cle, coeffc, coeffa, coeffn):
    aretes_1 = calcul_aretes_avec_coeff(graphe1, cle, coeffc, coeffa, coeffn)
    aretes_2 = calcul_aretes_avec_coeff(graphe2, cle, coeffc, coeffa, coeffn)
    aretes_commun = calcul_aretes_communes_avec_coeff(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn)
#     print(aretes_1)
#     print(aretes_2)
#     print(aretes_commun)
#     print(aretes_1)
#     print(aretes_2)
#     print(aretes_commun)
    
    return aretes_commun/max(aretes_1, aretes_2)

''' calcul de la similarite entre les deux graphes en prenant le nombre d'aretes en commun divise par le nombre d'aretes dans le 2e graphe '''
def calcul_sim_aretes_avec_coeff_graphe_moyen(graphe1, graphe2, graphe_commun, cle, coeffc, coeffa, coeffn):
    aretes_1 = calcul_aretes_avec_coeff(graphe1, cle, coeffc, coeffa, coeffn)
    aretes_2 = calcul_aretes_avec_coeff(graphe2, cle, coeffc, coeffa, coeffn)
    aretes_commun = calcul_aretes_communes_avec_coeff(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn)
#     print(aretes_1)
#     print(aretes_2)
#     print(aretes_commun)
#     print(aretes_1)
#     print(aretes_2)
#     print(aretes_commun)
    
    return aretes_commun/aretes_2

# def distance_entre_noeuds_meme_chaine(graphe, noeud1, noeud2):
# #     print(noeud1)
# #     print(noeud2)
#     chaine = -1
#     for elt in graphe.nodes[noeud1]["chaine"] :
#         if elt in graphe.nodes[noeud2]["chaine"] :
#             chaine = elt
#             
#     if chaine != -1 :
#         if graphe.nodes[noeud1]["position"][0] < graphe.nodes[noeud2]["position"][0] :
#             mini = noeud1
#             maxi = noeud2
#         else :
#             mini = noeud2
#             maxi = noeud1
#             
#         compteur = mini
#         dist = 0
#         temp = compteur
#         liaison_B53 = True
#         while temp != maxi and liaison_B53:
# #             print("ramou")
#             liaison_B53 = False
#             #if chaine in [2,3] :
#             for voisin in graphe.successors(compteur) :
#                     for arc in graphe[compteur][voisin] :
#                         if voisin not in [1,2,3,4] and graphe[compteur][voisin][arc]["label"] == 'B53' :
# #                             print("voisin")
# #                             print(voisin)
#                             temp = voisin
#                             dist += 1
#                             liaison_B53 = True
#                             
# #             if chaine in [1,4] :
# #                 for voisin in graphe.predecessors(compteur) :
# #                     for arc in graphe[voisin][compteur] :
# #                         if voisin not in [1,2,3,4] and graphe[voisin][compteur][arc]["label"] == 'B53' :
# #                             temp = voisin
# #                             dist += 1
# #                             liaison_B53 = True
# #             print(temp)           
#             compteur = temp
#             
#         return dist
#     return None

def test_compatibilite(graphe_commun, noeud, graphe1, graphe2):
    ''' renvoie vrai si ajouter le noeud au graphe_commun ne provoquera pas l'apparition d'une incompatibilité 
    de séquence dans le graphe_commun (deux superpositions de noeuds qui sont dans un ordre dans le 1er graphe
    et dans l'autre ordre dans le 2e graphe)'''
    
    for i in range(len(list(graphe_commun.nodes()))) :
            noeud2 = list(graphe_commun.nodes())[i]
            meme_chaine = False
            
            
            
            modele = graphe1.nodes[1]["position"][0]
            
            for elt in graphe1.nodes[noeud[0]]["chaine"] :
                if elt in graphe1.nodes[noeud2[0]]["chaine"] :
                    meme_chaine = True
                    
            for elt in graphe2.nodes[noeud[1]]["chaine"] :
                if elt in graphe2.nodes[noeud2[1]]["chaine"] :
                    meme_chaine = True
#             if noeud == (24,27) :
#                 print(noeud)
#                 print(graphe1.nodes[noeud[0]]["type"])
#                 print(graphe1.nodes.data())
#                 print(graphe2.nodes.data())
                   
#             if noeud == (1033,1034) :
#                 print(graphe1.nodes[noeud[0]]["chaine"])
            if noeud != noeud2 and meme_chaine and graphe1.nodes[noeud[0]]["type"] not in [None,-1] and graphe1.nodes[noeud2[0]]["type"] not in [None,-1] and noeud[0] not in [1,2,3,4,5] and noeud2[0] not in [1,2,3,4,5]:
                
#                 print("test compa")
#                 print(graphe1.nodes[noeud[0]]["position"])
#                 print(graphe1.nodes[noeud2[0]]["position"])
#                 
#                 print(graphe1.nodes[noeud[0]]["position"])
#                 print(graphe1.nodes[noeud2[0]]["position"])
                if (graphe1.nodes[noeud[0]]["position"][0] < graphe1.nodes[noeud2[0]]["position"][0] and graphe2.nodes[noeud[1]]["position"][0] > graphe2.nodes[noeud2[1]]["position"][0]) or (graphe1.nodes[noeud[0]]["position"][0] > graphe1.nodes[noeud2[0]]["position"][0] and graphe2.nodes[noeud[1]]["position"][0] < graphe2.nodes[noeud2[1]]["position"][0]) :
                    
#                     print("test compa")
#                     print(noeud)
#                     print(chaine_noeud11)
#                     print(chaine_noeud21)
#                     print(distance_entre_noeuds_meme_chaine(graphe1, chaine_noeud11, noeud[0]))
#                     print(distance_entre_noeuds_meme_chaine(graphe1, chaine_noeud21, noeud2[0]))
                    #if (distance_entre_noeuds_meme_chaine(graphe1, chaine_noeud11, noeud[0]) < distance_entre_noeuds_meme_chaine(graphe1, chaine_noeud21, noeud2[0]) and distance_entre_noeuds_meme_chaine(graphe2, chaine_noeud12, noeud[1]) > distance_entre_noeuds_meme_chaine(graphe2, chaine_noeud22, noeud2[1])) or  (distance_entre_noeuds_meme_chaine(graphe1, chaine_noeud11, noeud[0]) > distance_entre_noeuds_meme_chaine(graphe1, chaine_noeud21, noeud2[0]) and distance_entre_noeuds_meme_chaine(graphe2, chaine_noeud12, noeud[1]) < distance_entre_noeuds_meme_chaine(graphe2, chaine_noeud22, noeud2[1])) :
#                     if noeud == (17,18) and noeud2 == (14,15) :
#                             print("tout petit rat")
#                             print(noeud)
                    #if (abs(noeud[0]-modele) < abs(noeud2[0]-modele) and abs(noeud[1] - modele) > abs(noeud2[1]-modele)) or (abs(noeud[0]-modele) > abs(noeud2[0]-modele) and abs(noeud[1] - modele) < abs(noeud2[1]-modele)) :
                        
                    return False
    
    return True

def recup_chaines(graphe):
    ''' renvoie les numeros des sommets appartenant a chaque chaine (du motif vers le bout de la chaine) '''
    chaines = [[1]]
    for i in range(1,5) :
            compteur = i
            if i != 1 : chaines.append([i])
            liaison_B53 = True
            while liaison_B53 :
                liaison_B53 = False
                temp = compteur
                for voisin in graphe.successors(compteur) :
                    for arc in graphe[compteur][voisin] :
                        if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and graphe[compteur][voisin][arc]["label"] == 'B53' :
                            liaison_B53 = True
                            temp = voisin
                            chaines[len(chaines)-1].append(voisin)
                            
                for voisin in graphe.predecessors(compteur) :
                    for arc in graphe[voisin][compteur] :
                        if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and graphe[voisin][compteur][arc]["label"] == 'B53' :
                            liaison_B53 = True
                            temp = voisin
                            chaines[len(chaines)-1].append(voisin)
                compteur = temp
    return chaines

def liste_meme_type(chaines, graphe, chaines_comp, graphe_comp):
    '''renvoie une liste de dictionnaires dans lequel sont indiques pour chaque sommet du graphe (cle du dictionnaire)
    les sommets (de graphe_comp) de la meme chaine qui sont de meme type (et aussi -1 pour pouvoir considerer le cas ou le sommet n'est pas superpose) '''
    #groupe_types = [[1,2], [2,3], [None, -1]]
    #groupe_types = [[2,3]]
    groupe_types = [[1,2], [None, -1]]
    dico_chaines = [{}]
    for i in range(len(chaines)) :
        if i != 0 :
            dico_chaines.append({})
        for elt in chaines[i] :
            for elt_comp in chaines_comp[i] :
                if graphe.nodes[elt]["type"] == graphe_comp.nodes[elt_comp]["type"] :
                #if graphe.nodes[elt]["type"] == graphe_comp.nodes[elt_comp]["type"] or (graphe.nodes[elt]["type"] in groupe_types[0] and graphe_comp.nodes[elt_comp]["type"] in groupe_types[0]) or (graphe.nodes[elt]["type"] in groupe_types[1] and graphe_comp.nodes[elt_comp]["type"] in groupe_types[1]) or (graphe.nodes[elt]["type"] in groupe_types[2] and graphe_comp.nodes[elt_comp]["type"] in groupe_types[2]) :
                #if graphe.nodes[elt]["type"] == graphe_comp.nodes[elt_comp]["type"] or (graphe.nodes[elt]["type"] in groupe_types[0] and graphe_comp.nodes[elt_comp]["type"] in groupe_types[0])  :
                #if graphe.nodes[elt]["type"] == graphe_comp.nodes[elt_comp]["type"] or (graphe.nodes[elt]["type"] in groupe_types[0] and graphe_comp.nodes[elt_comp]["type"] in groupe_types[0]) or (graphe.nodes[elt]["type"] in groupe_types[1] and graphe_comp.nodes[elt_comp]["type"] in groupe_types[1]) :
                #if graphe.nodes[elt]["type"] == graphe_comp.nodes[elt_comp]["type"] or (graphe.nodes[elt]["type"] == 1 and graphe_comp.nodes[elt_comp]["type"] == 3 and arete_CWWn(graphe_comp, elt_comp)) \
                   # or (graphe.nodes[elt]["type"] == 3 and graphe_comp.nodes[elt_comp]["type"] == 1 and arete_CWWn(graphe, elt)) \
                    #or (graphe.nodes[elt]["type"] == 2 and graphe_comp.nodes[elt_comp]["type"] == 3 and arete_CWWn(graphe_comp, elt_comp))  \
                    #or (graphe.nodes[elt]["type"] == 3 and graphe_comp.nodes[elt_comp]["type"] == 2 and arete_CWWn(graphe, elt)) :

                    if elt not in dico_chaines[len(dico_chaines)-1].keys() :
                        dico_chaines[len(dico_chaines)-1].update({elt : [elt_comp]})
                    else :
                        dico_chaines[len(dico_chaines)-1][elt].append(elt_comp)
    
    for i in range(len(chaines)) :
        chaine = chaines[i]
        for elt in chaine :
            if elt not in dico_chaines[i].keys() :
                dico_chaines[i].update({elt : [-1]})
            else :
                dico_chaines[i][elt].append(-1)
            
        
    return dico_chaines        

def dans_graphe(graphe, couple_a_chercher):
    ''' renvoie vrai si l'un des noeuds du couple_a_chercher est deja present dans une paire dans le graphe et pas l'autre
    faux sinon '''
    for noeud in graphe.nodes() :
        #print(noeud)
        if (noeud[0] == couple_a_chercher[0] and noeud[1] != couple_a_chercher[1]) or (noeud[1] == couple_a_chercher[1] and noeud[0] != couple_a_chercher[0])  :
            return True
    return False


''' 21/03/20 test pour accepter superposition CWW/CWWn '''
def arete_CWWn(graphe, noeud):
    for voisin in graphe[noeud] :
        for edge in graphe[noeud][voisin] :
            if graphe[noeud][voisin][edge]["label"] == 'CWWn' :
                return True
    return False

def meme_type_meme_chaine(noeud1, noeud2, graphe1, graphe2):
    ''' renvoie vrai si les 2 noeuds (issus respectivement de graphe1 et graphe2) sont de meme type et appartiennent a la meme chaine '''
    
    #groupe_types = [[1,2], [2,3], [None, -1]]
    #groupe_types = [[2,3]]
    groupe_types = [[1,2], [None, -1]]
    meme_chaine = False
    for elt in graphe1.nodes[noeud1]["chaine"] :
        if elt in graphe2.nodes[noeud2]["chaine"] :
            meme_chaine = True
    if (graphe1.nodes[noeud1]["type"] == graphe2.nodes[noeud2]["type"] and meme_chaine) or (noeud1 in [1,2,3,4,5] and noeud2 in [1,2,3,4,5]):
    #if ((graphe1.nodes[noeud1]["type"] == graphe2.nodes[noeud2]["type"] or (graphe1.nodes[noeud1]["type"] in groupe_types[0] and graphe2.nodes[noeud2]["type"] in groupe_types[0]) or (graphe1.nodes[noeud1]["type"] in groupe_types[1] and graphe2.nodes[noeud2]["type"] in groupe_types[1]) or (graphe1.nodes[noeud1]["type"] in groupe_types[2] and graphe2.nodes[noeud2]["type"] in groupe_types[2])) and meme_chaine) or (noeud1 in [1,2,3,4,5] and noeud2 in [1,2,3,4,5]) :
    #if ((graphe1.nodes[noeud1]["type"] == graphe2.nodes[noeud2]["type"] or (graphe1.nodes[noeud1]["type"] in groupe_types[0] and graphe2.nodes[noeud2]["type"] in groupe_types[0])) and meme_chaine) or (noeud1 in [1,2,3,4,5] and noeud2 in [1,2,3,4,5]) :
    #if ((graphe1.nodes[noeud1]["type"] == graphe2.nodes[noeud2]["type"] or (graphe1.nodes[noeud1]["type"] in groupe_types[0] and graphe2.nodes[noeud2]["type"] in groupe_types[0]) or (graphe1.nodes[noeud1]["type"] in groupe_types[1] and graphe2.nodes[noeud2]["type"] in groupe_types[1])) and meme_chaine) or (noeud1 in [1,2,3,4,5] and noeud2 in [1,2,3,4,5]) :
    #if ((graphe1.nodes[noeud1]["type"] == graphe2.nodes[noeud2]["type"] or (graphe1.nodes[noeud1]["type"] == 1 and graphe2.nodes[noeud2]["type"] == 3 and arete_CWWn(graphe2, noeud2)) \
    #    or (graphe1.nodes[noeud1]["type"] == 3 and graphe2.nodes[noeud2]["type"] == 1 and arete_CWWn(graphe1, noeud1)) \
    #    or (graphe1.nodes[noeud1]["type"] == 2 and graphe2.nodes[noeud2]["type"] == 3 and arete_CWWn(graphe2, noeud2))  \
    #    or (graphe1.nodes[noeud1]["type"] == 3 and graphe2.nodes[noeud2]["type"] == 2 and arete_CWWn(graphe1, noeud1)) \
    #    or (graphe1.nodes[noeud1]["type"] == None and graphe2.nodes[noeud2]["type"] == -1 and arete_CWWn(graphe1, noeud1)) \
    #    or (graphe1.nodes[noeud1]["type"] == -1 and graphe2.nodes[noeud2]["type"] == None and arete_CWWn(graphe2, noeud2))) and meme_chaine) or (noeud1 in [1,2,3,4,5] and noeud2 in [1,2,3,4,5]):
        return True
    else :
        return False

def algo_principal(graphe1, graphe2, chaines_1, dico_chaines_1, graphe_commun_temp, cle):
    ''' 
    while pile non vide :
        x,y = depiler
        ajouter la paire x,y au graphe commun
        if num_chaine(x et y) == 4 and x == dernier element de la derniere chaine :
            on arrive en fin de tour => on regarde si le graphe courant est le meilleur vu jusque la
            if val_comp > max_val_comp :
                max_val_comp = val_comp
                graphe_commun_max = graphe_courant
                
            
        
        else : on est encore en phase d'empiler des trucs
            if x+1 < len(chaines_1[num_chaine(x)]) :
                for super_possible in dico_chaines_1[num_chaine(x+1)][x+1] :
                    if compatibles(super_possible, x+1) and alone(super_possible) and alone(x+1)
                        empiler la paire  
            else :
                on change de chaine 
                for super_possible in dico_chaines_1[num_chaine(x)+1][0] :
                    if compatibles(super_possible, chaines_1[num_chaine(x)+1][0]) and alone(super_possible) and alone(chaines_1[num_chaine(x)+1][0])
                        empiler la paire 
       
             '''  
    
    ### Initialisation
    pile = []
#     print(chaines_1[0][0])
    for super_possible in dico_chaines_1[0][chaines_1[0][0]] :
        pile.append((chaines_1[0][0], super_possible, 0, 0, 0)) 
    
    
    graphe_commun_max = graphe_commun_temp.copy()
    sim_max = 0.0
    
    liste_etapes = []

    while len(pile) > 0 :
#         print("pile courante")
#         print(pile)
#         print("graphe courant")
#         print(graphe_commun_temp.nodes.data())
#         print(graphe_commun_temp.edges.data())
        
        
        
        x,y,num_chaine,num_x,etape = pile.pop()
        #if x == 15 and y == 12 :
        #print(x,y)
#         print(graphe_commun_temp.nodes())
        #print(test_compatibilite(graphe_commun_temp, (x,y), graphe1, graphe2))
#         print(etape)
        #print(graphe_commun_temp.edges.data())
            
        for i in range(etape, len(liste_etapes)) :
    #             print(liste_etapes[i])
#                 print("gros rat")
#                 print(etape)
#                 print(len(liste_etapes))
                if "noeuds" in liste_etapes[i].keys() :
                    graphe_commun_temp.remove_nodes_from(liste_etapes[i]["noeuds"])
                if "aretes" in liste_etapes[i].keys() :
                    graphe_commun_temp.remove_edges_from(liste_etapes[i]["aretes"])
        
#         if x == 27 and y == 27 : 
#             print(graphe_commun_temp.nodes.data())
#             print(dans_graphe(graphe_commun_temp, (x,y)))
#             print(test_compatibilite(graphe_commun_temp, (x,y), graphe1, graphe2))
#             print(liste_etapes)
#             print(etape)
#         if x == 17 and y == 18 :
#             print("gros tas")
#             print(graphe_commun_temp.nodes.data())
#             print(test_compatibilite(graphe_commun_temp, (x,y), graphe1, graphe2))
        if y == -1 or (not dans_graphe(graphe_commun_temp, (x,y)) and test_compatibilite(graphe_commun_temp, (x,y), graphe1, graphe2)) :
            
#             if x == 15 and y == 12 :
#                 print(x,y)
#                 print(graphe_commun_temp.nodes())
#             if x == 6 and y == 6 :
#                 print(x,y)
#                 print(graphe_commun_temp.nodes())
            
            #print(liste_etapes)
    #         print(etape)
            
#             for i in range(etape, len(liste_etapes)) :
#     #             print(liste_etapes[i])
#                 if "noeuds" in liste_etapes[i].keys() :
#                     graphe_commun_temp.remove_nodes_from(liste_etapes[i]["noeuds"])
#                 if "aretes" in liste_etapes[i].keys() :
#                     graphe_commun_temp.remove_edges_from(liste_etapes[i]["aretes"])
            
            if y != -1 :
                
                ## ajout de la paire (x,y) au graphe commun
                taille_liste_etapes = len(liste_etapes)
                if not dans_graphe(graphe_commun_temp, (x,y)) and (x,y) not in graphe_commun_temp.nodes() :
                    graphe_commun_temp.add_node((x,y))
                    liste_etapes.append({"noeuds" :[(x,y)]})
                
                ## recherche de voisins non cov de x et y a ajouter aussi
                
                voisins_1 = graphe1[x]
                voisins_2 = graphe2[y]
                
                ## on construit un dico des labels identiques au prealable pour parer au probleme de plusieurs aretes de meme type incidentes au meme sommet
                dico_labels = {}
                for voisin_1 in voisins_1 :
                    for voisin_2 in voisins_2 :
                            for edge_1 in graphe1[x][voisin_1] :
                                for edge_2 in graphe2[y][voisin_2] :
                                    if graphe1[x][voisin_1][edge_1]["label"] != 'B53' and graphe2[y][voisin_2][edge_2]["label"] != 'B53' :
                                        label1 = graphe1[x][voisin_1][edge_1]["label"]
                                        label2 = graphe2[y][voisin_2][edge_2]["label"]
                                        if label1 == label2 :
                                            if (voisin_1, edge_1) not in dico_labels.keys() :
                                                dico_labels.update({(voisin_1, edge_1) : [(voisin_2, edge_2)]})
                                            else :
                                                dico_labels[(voisin_1, edge_1)].append((voisin_2, edge_2))
#                 print("roupoulou")
#                 print(dico_labels)
                
                ## on choisit la paire ou le nombre de voisins est le plus proche                        
                for cle in dico_labels.keys() :
                    if len(dico_labels[cle]) > 1 :
                        diff_voisins = 5
                        meilleur_voisin = -1
                        compteur_v = 0
                        for elt in dico_labels[cle] :
                            if diff_voisins > abs(len(graphe1[cle[0]])-len(graphe2[elt[0]])) :
                                diff_voisins = abs(len(graphe1[cle[0]])-len(graphe2[elt[0]]))
                                meilleur_voisin = compteur_v
                            compteur_v += 1
                        dico_labels[cle] = [dico_labels[cle][meilleur_voisin]]
#                 
#                 print("rapala")
#                 print(dico_labels)
                
                
                voisins = {}
                #for voisin1 in graphe1[x] :
                    
#                     for edge1 in graphe1[x][voisin1] :
#                         if graphe1[x][voisin1][edge1]["label"] != 'B53' :
#                             label1 = graphe1[x][voisin1][edge1]["label"]
#                             #long_range1 = graphe1[x][voisin1][edge1]["long_range"]
#                             #if graphe1[x][voisin1][edge1]["near"] != True :
#                             for voisin2 in graphe2[y] :
#                                     for edge2 in graphe2[y][voisin2] :
#                                         if graphe2[y][voisin2][edge2]["label"] != 'B53' :
                for cle in dico_labels.keys() :
                                            voisin1 = cle[0]
                                            edge1 = cle[1]
                                            voisin2 = dico_labels[cle][0][0]
                                            edge2 = dico_labels[cle][0][1]
                                            label1 = graphe1[x][voisin1][edge1]["label"]
                                            label2 = graphe2[y][voisin2][edge2]["label"]
                                            
#                                             print((x,voisin1), (y,voisin2))
#                                             print(label1, label2)
                                            
    #                                         print("rapoulou")
    #                                         print(voisin1, voisin2)
    #                                         print(label1)
    #                                         print(graphe2[y][voisin2][edge2]["label"])
    #                                         print(label1 == graphe2[y][voisin2][edge2]["label"] )
    #                                         print("ramou")
    #                                         print((voisin1, voisin2))
    #                                         print(label1)
    #                                         print(graphe2[y][voisin2][edge2]["label"])
                                            #if label1 == graphe2[y][voisin2][edge2]["label"] or (label1 in ['CWW', '0'] and  graphe2[y][voisin2][edge2]["label"] in ['CWW', '0']): #and long_range1 == graphe2[y][voisin2][edge2]["long_range"] : ## meme label des aretes
                                            #if graphe2[y][voisin2][edge2]["near"] != True and label1 == graphe2[y][voisin2][edge2]["label"] : #and long_range1 == graphe2[y][voisin2][edge2]["long_range"] : ## meme label des aretes
                                            #if label1 == graphe2[y][voisin2][edge2]["label"] or (label1 in ['0', 'CWW', 'CWWn'] and graphe2[y][voisin2][edge2]["label"] in ['0', 'CWW', 'CWWn']):    
                                            #if label1 == graphe2[y][voisin2][edge2]["label"] : 
                                            
                                                
                                                #print(voisin1, voisin2)
                                            if label1 != 'B53' and label2 != 'B53' :
                                                
                                                if meme_type_meme_chaine(voisin1, voisin2, graphe1, graphe2)  :
                                                   
                                                    if not dans_graphe(graphe_commun_temp, (voisin1, voisin2)) and test_compatibilite(graphe_commun_temp, (voisin1, voisin2), graphe1, graphe2) : ## verif elements de la paire des voisins ne sont pas deja dans ue autre paire
                                                             
            #                                                 print(liste_etapes)
                                                        #if graphe1[x][voisin1][edge1]["near"] == False and graphe2[y][voisin2][edge2]["near"] == False :
        
                                                            if (voisin1, voisin2) not in graphe_commun_temp.nodes() :
                                                                if taille_liste_etapes == len(liste_etapes) :
                                                                    liste_etapes.append({"noeuds" :[(voisin1, voisin2)]})
                                                                else :
                                                                    if "noeuds" not in liste_etapes[len(liste_etapes)-1].keys() :
                                                                        liste_etapes[len(liste_etapes)-1].update({"noeuds" :[(voisin1, voisin2)]})
                                                                    else :
                                                                        liste_etapes[len(liste_etapes)-1]["noeuds"].append((voisin1, voisin2))
                                                                
                                                                graphe_commun_temp.add_node((voisin1, voisin2))
                                                                
                                                            if ((x,y), (voisin1, voisin2)) in graphe_commun_temp.edges() :
                                                                deja_vu = False
                                                                for edge in graphe_commun_temp[(x,y)][(voisin1,voisin2)] :
                                                                    if graphe_commun_temp[(x,y)][(voisin1,voisin2)][edge]["label"] == label1 :
                                                                        deja_vu = True
                                     
                                                            if ((x,y), (voisin1, voisin2)) not in graphe_commun_temp.edges() or not deja_vu:
#                                                                 print(((x,y), (voisin1, voisin2)))
                                                                graphe_commun_temp.add_edge((x,y), (voisin1, voisin2), label=label1, motif=False)
                                                                if len(label1) == 3 :
                                                                    graphe_commun_temp.add_edge((voisin1, voisin2), (x,y), label=label1[0]+label1[2]+label1[1], motif=False)
                                                                elif len(label1) == 4 :
                                                                    graphe_commun_temp.add_edge((voisin1, voisin2), (x,y), label=label1[0]+label1[2]+label1[1]+label1[3:], motif=False)
                                                                else : ## aretes artificielles
                                                                    graphe_commun_temp.add_edge((voisin1, voisin2), (x,y), label=label1, motif=False)
                                                                    
                                                                if taille_liste_etapes == len(liste_etapes) :
                                                                    liste_etapes.append({"aretes" :[((x,y),(voisin1, voisin2)), ((voisin1, voisin2), (x,y))]})
                                                                else :
                                                                    if "aretes" not in liste_etapes[len(liste_etapes)-1].keys() :
                                                                        liste_etapes[len(liste_etapes)-1].update({"aretes" :[((x,y),(voisin1, voisin2)), ((voisin1, voisin2), (x,y))]})
                                                                    else :
                                                                        liste_etapes[len(liste_etapes)-1]["aretes"].append(((x,y),(voisin1, voisin2)))
                                                                        liste_etapes[len(liste_etapes)-1]["aretes"].append(((voisin1, voisin2), (x,y))) 
                    
            ## si c'est la fin d'un tour on calcule la sim et on compare au max
            
            if num_chaine == len(chaines_1) -1 and x == chaines_1[len(chaines_1)-1][len(chaines_1[len(chaines_1)-1])-1] :
                #sim = calcul_sim_aretes_avec_coeff(graphe1, graphe2, graphe_commun_temp, cle, 1, 1, 1) 
                
                #print(graphe_commun_temp.nodes.data())
                sim = calcul_aretes_communes_avec_coeff(graphe_commun_temp, graphe1, graphe2, cle, 1, 1, 1)
#                 print(sim)
#                 print(sim_max)
                if sim > sim_max :
                    sim_max = sim
                    graphe_commun_max = graphe_commun_temp.copy() 
#                     print("ripili")
#                     print(graphe_commun_max.nodes.data())  
#                     print(graphe_commun_max.edges.data())
#                     print(sim)
#                     print("petit rat")
#                     print(graphe_commun_max.nodes.data())
    #             print("sim courante")
    #             print(sim)
    #             print("graphe commun courant")
    #             print(graphe_commun_temp.nodes.data())
    #             if (8,6) in graphe_commun_temp.nodes() :
    #                 break
    #             print("graphe commun max courant")
    #             print(graphe_commun_max.nodes.data())
    #             
    #             for i in range(etape, len(liste_etapes)) :
    #                 graphe_commun_temp.remove_nodes_from(liste_etapes[i]["noeuds"])
    #                 graphe_commun_temp.remove_edges_from(liste_etapes[i]["aretes"])
    #             print("verif")
    #             print(graphe_commun_temp.edges.data()) 
            
            ## ce n'est pas la fin d'un tour : on continue de descendre dans l'arbre des possibles
            else :
                if num_x+1 < len(chaines_1[num_chaine]) : 
                    ## on n'a pas fini la chaine courante
                    #print(num_x)
                    suivant = chaines_1[num_chaine][num_x+1]
                    for super_possible in dico_chaines_1[num_chaine][suivant] :
#                         if suivant == 13 and super_possible == 13 :
#                             print(test_compatibilite(graphe_commun_temp, (suivant, super_possible), graphe1, graphe2))
                        if super_possible == -1 or (test_compatibilite(graphe_commun_temp, (suivant, super_possible), graphe1, graphe2) and not dans_graphe(graphe_commun_temp, (suivant, super_possible)) ):#and (suivant, super_possible) not in graphe_commun_temp.nodes() ):
                            pile.append((suivant, super_possible, num_chaine, num_x+1, len(liste_etapes))) 
#                             print("empiler")
#                             print(suivant, super_possible) 
                            
                else :
                    ## on change de chaine
                    suivant = chaines_1[num_chaine+1][0]
    #                 print(suivant)
                    for super_possible in dico_chaines_1[num_chaine+1][suivant] :
                        if super_possible == -1 or suivant in [1,2,3,4] or (test_compatibilite(graphe_commun_temp, (suivant, super_possible), graphe1, graphe2) and not dans_graphe(graphe_commun_temp, (suivant, super_possible)) ):#and (suivant, super_possible) not in graphe_commun_temp.nodes() ):
                            pile.append((suivant, super_possible, num_chaine+1, 0, len(liste_etapes)))
                            #print("empiler")
                            #print(suivant, super_possible) 
            
    return graphe_commun_max
                        
#                 print("pile courante")
#                 print(pile)



def comparaison(graphe1, graphe2, cle):
    ''' execute la recherche de graphe commun max entre les deux graphes '''
    chaines_1 = recup_chaines(graphe1)
    chaines_2 = recup_chaines(graphe2)
    
    
    
    print(chaines_1)
    print(chaines_2)
    
    dico_chaines_1 = liste_meme_type(chaines_1, graphe1, chaines_2, graphe2)
    #dico_chaines_2 = liste_meme_type(chaines_2, graphe2, chaines_1, graphe1)
    print(dico_chaines_1)
    
    chaines_reliees_1 = chaines_reliees(graphe1)
    chaines_reliees_2 = chaines_reliees(graphe2)
    print(chaines_reliees_1)
    print(chaines_reliees_2)
    
    chaines_reliees_tot = []
    chaines_reliees_tot.extend(list(chaines_reliees_1))
    print(chaines_reliees_1)
    chaines_reliees_tot.extend(list(chaines_reliees_2))

    petit_dico = {}
    for i in range(1,5) :
        for elt in chaines_reliees_tot :
            if i in elt :
                if i not in petit_dico.keys() :
                    petit_dico.update({i : [e for e in elt if e != i]})
                else :
                    for e in elt :
                        if e != i :
                            if e not in petit_dico[i] :
                                petit_dico[i].append(e)
                            
    chaines_a_mettre_ensemble = [[1]]
    dico_vu = {1:True, 2:False, 3:False, 4:False}
    petite_pile = [(1, 0)]
    while False in dico_vu.values() :
        if len(petite_pile) > 0 :
            cle, num_groupe = petite_pile.pop()
            if cle in petit_dico.keys() :
                for elt in petit_dico[cle] :
                    if not dico_vu[elt]:
                        petite_pile.append((elt, num_groupe))
                        dico_vu[elt] = True
                        chaines_a_mettre_ensemble[num_groupe].append(elt)
                
        else :
   
            for cle in dico_vu.keys() :
                if not dico_vu[cle] :
                    num_groupe += 1
                    petite_pile.append((cle, num_groupe))
                    dico_vu[cle] = True
                    chaines_a_mettre_ensemble.append([])
                    chaines_a_mettre_ensemble[num_groupe].append(cle)
                    break
                    
        
        
            
            
#             #if elt > cle :
#                 vu = False
#                 for elt2 in chaines_a_mettre_ensemble :
#                     for e in elt2 :
#                         if e == cle :
#                             if elt not in elt2 :
#                                 elt2.append(elt)
#                             if dico_vu[elt] == False :
#                                 dico_vu[elt] = True
#                             vu = True
#                 if not vu :
#                     if dico_vu[cle] == False :
#                         dico_vu[cle] = True
#                     if dico_vu[elt] == False :
#                         dico_vu[elt] = True
#                     chaines_a_mettre_ensemble.append([cle, elt])
#                 temp = 
#     
#     for cle in dico_vu.keys() :
#         if not dico_vu[cle] :
#             chaines_a_mettre_ensemble.append([cle])
                
    print("ramou")
    print(chaines_reliees_1)
    print(chaines_reliees_2)
    print(chaines_reliees_tot)
    print(petit_dico)
    print(chaines_a_mettre_ensemble)
        
                
    
    graphe_commun_temp = nx.MultiDiGraph()
    
    for i in range(1,6) :
        graphe_commun_temp.add_node((i,i))
    graphe_commun_temp.add_edge((1,1),(2,2), label="CSS", motif=True)
    graphe_commun_temp.add_edge((2,2),(1,1), label="CSS", motif=True)
    graphe_commun_temp.add_edge((1,1),(5,5), label="TSS", motif=True)
    graphe_commun_temp.add_edge((5,5),(1,1), label="TSS", motif=True)
    graphe_commun_temp.add_edge((2,2),(5,5), label="CWW", motif=True)
    graphe_commun_temp.add_edge((5,5),(2,2), label="CWW", motif=True)
    graphe_commun_temp.add_edge((3,3),(4,4), label="CSS", motif=True)
    graphe_commun_temp.add_edge((4,4),(3,3), label="CSS", motif=True)
           
    #print(chaines_reliees_tot)
    for elt in chaines_a_mettre_ensemble :
        chaines = []
        dico_chaines = []
        for chaine in elt :
            print(chaine)
            chaines.append(chaines_1[chaine-1])
            dico_chaines.append(dico_chaines_1[chaine-1])
        graphe_commun_temp = algo_principal(graphe1, graphe2, chaines, dico_chaines, graphe_commun_temp, cle)
        
    graphe_commun_max = graphe_commun_temp.copy()
    sim_max = calcul_sim_aretes_avec_coeff(graphe1, graphe2, graphe_commun_max, cle, 1,1,1)
    
    print(graphe_commun_max.edges.data())
        
    ## recherche s'il faut ajouter une liaison cov => a la fin     
    for noeud in graphe_commun_max.nodes() :
        num_voisin1_b53 = -1#             print("sim courante")
#             print(sim)
#             print("graphe commun courant")
#             print(graphe_commun_temp.nodes.data())
        for voisin in graphe1[noeud[0]] :
            for edge in graphe1[noeud[0]][voisin] :
                if graphe1[noeud[0]][voisin][edge]["label"] == 'B53' :
                    num_voisin1_b53 = voisin
                    
        num_voisin2_b53 = -1
        for voisin in graphe2[noeud[1]] :
            for edge in graphe2[noeud[1]][voisin] :
                if graphe2[noeud[1]][voisin][edge]["label"] == 'B53' :
                    num_voisin2_b53 = voisin
        
        if num_voisin1_b53 != -1 and num_voisin2_b53 != -1 and (num_voisin1_b53, num_voisin2_b53) in graphe_commun_max.nodes() :
            graphe_commun_max.add_edge(noeud, (num_voisin1_b53, num_voisin2_b53), label='B53', motif=False)
            
    if sim_max >= 0.65 :
        new_sim_max = verif_sim_branche(graphe_commun_max, sim_max, chaines_1, chaines_2, graphe1, graphe2, 1,1,1, 0.4)
     
        print(sim_max)
        if new_sim_max != sim_max :
            return graphe_commun_max, new_sim_max, True
        else :
            return graphe_commun_max, new_sim_max, False
    return graphe_commun_max, sim_max, False

def verif_sim_branche(graphe_commun_max, sim_max, chaines_1, chaines_2, graphe1, graphe2, coeffc, coeffa, coeffn, seuil):
    
    
    compteur = 0
    compter_branche = 0
    for chaine in chaines_1 :
        poids_aretes_commun = 0
        poids_aretes_1 = 0
        poids_aretes_2 = 0
        for u,v, data in graphe_commun_max.edges(data=True) :
#             print(u,v)
#             print(data["motif"])
            if ((u[0] in chaine and u[1] in chaines_2[compteur]) or (v[0] in chaine and v[1] in chaines_2[compteur])) and data["motif"] == False :
                #print(u,v)
                if data["label"] != 'B53' :
                    if data["label"] == '0' :
                        if coeffa == 1 :
                            poids_aretes_commun += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
                    elif graphe1.nodes[u[0]]["type"] == 1 and graphe1.nodes[v[0]]["type"] == 1 :
                        if coeffc == 1 :
                            poids_aretes_commun += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
                    else :
                        if coeffn == 1 :
                            poids_aretes_commun += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"])
        
        for u,v, data in graphe1.edges(data=True) :
            if u in chaine or v in chaine :
                if (u,v,data["label"]) not in [(1,2,'CSS'), (2,1,'CSS'), (1,5,'TSS'), (5,1,'TSS'), (2,5,'CWW'), (5,2, 'CWW'), (3,4,'CSS'), (4,3,'CSS')] :
                    if data["label"] != 'B53' :
                        if data["label"] == '0' :
                            if coeffa == 1 :
                                poids_aretes_1 += graphe1.nodes[u]["poids"]
                                
                        elif graphe1.nodes[u]["type"] == 1 and graphe1.nodes[v]["type"] == 1 :
                            if coeffc == 1 :
                                poids_aretes_1 += graphe1.nodes[u]["poids"]
                                
                        else :
                            if coeffn == 1 :
                                poids_aretes_1 += graphe1.nodes[u]["poids"]
        
        for u,v, data in graphe2.edges(data=True) :
            if u in chaines_2[compteur] or v in chaines_2[compteur] :
                if (u,v,data["label"]) not in [(1,2,'CSS'), (2,1,'CSS'), (1,5,'TSS'), (5,1,'TSS'), (2,5,'CWW'), (5,2, 'CWW'), (3,4,'CSS'), (4,3,'CSS')] :
                    if data["label"] != 'B53' :
                        if data["label"] == '0' :
                            if coeffa == 1 :
                                poids_aretes_2 += graphe2.nodes[u]["poids"]
                                
                        elif graphe2.nodes[u]["type"] == 1 and graphe2.nodes[v]["type"] == 1 :
                            if coeffc == 1 :
                                poids_aretes_2 += graphe2.nodes[u]["poids"]
                                
                        else :
                            if coeffn == 1 :
                                poids_aretes_2 += graphe2.nodes[u]["poids"]
        print(poids_aretes_commun)
        print(poids_aretes_1)
        print(poids_aretes_2)                       
        sim = (poids_aretes_commun)/max(poids_aretes_1, poids_aretes_2)
        print(sim)
        if sim <= seuil  :
            compter_branche += 1
        

        compteur += 1
    print("gros tas")
    print(compter_branche)
    sim_max = (1-(compter_branche/4))*sim_max
    return sim_max


                

''' Obtention de resultats '''

''' 16/03/20 '''
def genere_graphe_seuil_sim(graphe, seuil) :
    a_enlever = []
    for u,v,data in graphe.edges(data=True) :
        if data["sim"] < seuil :
            a_enlever.append((u,v))
    for elt in a_enlever :
        graphe.remove_edge(elt[0], elt[1])

''' 16/03/20 '''
def genere_graphe_seuil_rmsd(graphe, seuil) :
    a_enlever = []
    for u,v,data in graphe.edges(data=True) :
        if data["rmsd"] == None or data["rmsd"] > seuil :
            a_enlever.append((u,v))
    for elt in a_enlever :
        graphe.remove_edge(elt[0], elt[1])
        

def nb_liaison_near(graphe):
    nb_near = 0
    for u,v, data in graphe.edges(data=True) :   
        if data["near"] == True and (u,v) not in [(1,2), (2,1), (3,4), (4,3), (1,5), (5,1)] :
            nb_near += 1
            
    return int(nb_near/2)

def enum_graphe_liaison_near(graphe, nom):
    liste_liaison_near = []
    liste_aretes = []
    for u,v, key, data in graphe.edges(keys=True, data=True) :
        
        if  ((isinstance(data["near"], bool) and data["near"] == True) or (isinstance(data["near"], list) and True in data["near"])) and (u,v) not in liste_aretes and (u,v) not in [(1,2), (2,1), (3,4), (4,3), (1,5), (5,1)]:
            #print(u,v)
#             print(u,v)
#             print(data["label"])
#             print(data["near"])
            if (v,u) in graphe.edges() :
                for edge in graphe[v][u] :
                    if len(data["label"]) == 3 : 
                        label_inv = data["label"][0] + data["label"][2] + data["label"][1]
                    elif len(data["label"]) == 4 :
                        label_inv = data["label"][0] + data["label"][2] + data["label"][1] + data["label"][3]
                    else :
                        print("probleme")
                    if graphe[v][u][edge]["label"] == label_inv and graphe[v][u][edge]["near"] == data["near"] :
                        liste_liaison_near.append(((u,v,key), (v,u,edge))) 
                liste_aretes.append((u,v))
                liste_aretes.append((v,u))
            else :
                print("bizarre")
                
            if isinstance(data["near"], list) and True in data["near"] :
                with open("fichier_near_true_false.txt", 'a') as fichier_txt :
                    fichier_txt.write(str(nom) + " " + str(u) + "," + str(v) + " " + str(data["label"]) + " " + str(data["near"]) + "\n")
    
#     nb = nb_liaison_near(graphe) 
#     if nb != len(liste_liaison_near) : 
#         print("ripili")
#         print(nb)
#         print(len(liste_liaison_near))
#         print(liste_liaison_near)
#         print(liste_aretes)
             
    return liste_liaison_near


''' issu de https://python.jpvweb.com/python/mesrecettespython/doku.php?id=parties_ensemble&s[]=ensemble&s[]=parties&s[]=un 
donne la liste des ensemble des parties d'un ensemble (en utilisant notation binaire)'''    
def partiesliste(seq):
    p = []
    i, imax = 0, 2**len(seq)-1
    while i <= imax:
        s = []
        j, jmax = 0, len(seq)-1
        while j <= jmax:
            if (i>>j)&1 == 1:
                s.append(seq[j])
            j += 1
        p.append(s)
        i += 1 
    return p  
            


def obtention_comparaison_avec_liaison_near(liste_num_ARN):
    ### Recup des occurrences non redondantes et avec une bonne resolution ###
    liste_tout = []
    for elt in liste_num_ARN :
        with open("resolutions.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            resolutions = mon_depickler.load()
             
            with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
                    mon_depickler = pickle.Unpickler(fichier_sortie)
                    liste_a_garder_noms = mon_depickler.load()    
                       
                    for element in liste_a_garder_noms :
                        if resolutions[element[0]] <= 3 and element not in liste_pbs :
                            liste_tout.append((elt, element))
    
    parties_liste_taille_1_4 = []
    liste = []
    for i in range(1,5) :
        liste.append(i)
        parties_liste_taille_1_4.append(partiesliste(liste))
    print(parties_liste_taille_1_4)
    dico_graphes = {}
    compter = 0
    for i in range(len(liste_tout)) :
        with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+liste_tout[i][1][0]+"_"+str(liste_tout[i][1][1])+"_2.pickle", 'rb') as fichier1 :
            mon_depickler_1 = pickle.Unpickler(fichier1)
            graphe1 = mon_depickler_1.load()

            liste_near = enum_graphe_liaison_near(graphe1)
            #nb_near = nb_liaison_near(graphe1)
            print(liste_tout[i])
            print(len(liste_near))
            print(liste_near)
        
            if len(liste_near) != 0 :
                if len(liste_near) == 3 :
                    print("rapoulou")
                    compter += 1
                liste_temp_graphes = []
                print("premier")
                print(graphe1.edges.data())
                print("modifie")
                for comb in parties_liste_taille_1_4[len(liste_near)-1] :
                    graph_copy = graphe1.copy()
                    for elt in comb :
                        graph_copy.remove_edge(liste_near[elt-1][0][0], liste_near[elt-1][0][1], key=liste_near[elt-1][0][2])
                        graph_copy.remove_edge(liste_near[elt-1][1][0], liste_near[elt-1][1][1], key=liste_near[elt-1][1][2])
                    print(graph_copy.edges.data())
                    liste_temp_graphes.append(graph_copy)
                dico_graphes.update({liste_tout[i][1] : liste_temp_graphes})
            else :
                dico_graphes.update({liste_tout[i][1] : [graphe1]})
            print("apres")
            print(graphe1.edges.data())
            print(dico_graphes[liste_tout[i][1]])
    print(compter)
    #exit()
    dico_new = {}
    compteur = 0
    #liste_faux_neg = [(('4ybb', 12), ('4u4r', 21)),(('4w2g', 52), ('6az1', 3)),(('4ybb', 12), ('6az1', 3)),(('2zjr', 3), ('6az1', 3)),(('4v67', 7), ('6az1', 3)),(('5wfs', 19), ('6az1', 3)),(('4u4r', 18), ('4y4o', 11)),(('4ybb', 1), ('4y4o', 22)),(('4ybb', 12), ('4y4o', 22)),(('4y4o', 25), ('4y4o', 22)),(('2zjr', 3), ('4y4o', 22)),(('4v67', 7), ('4y4o', 22)),(('6hma', 1), ('4y4o', 22)),(('5wfs', 19), ('4y4o', 22)),(('5dm6', 2), ('5nwy', 17)),(('1vq8', 16), ('3mum', 1)),(('6hma', 1), ('6ek0', 7)),(('6hma', 14), ('4faw', 2)),(('1vq8', 18), ('4faw', 2)),(('5ngm', 15), ('4faw', 2)),(('6eri', 17), ('4faw', 2)),(('6eri', 17), ('4y1n', 1)),(('5dm6', 4), ('4faw', 2)),(('4w2g', 52), ('4y4o', 22)),(('6hma', 8), ('4faw', 2)),(('6eri', 12), ('4faw', 2)),(('6qul', 4), ('4y4o', 38)),(('4ybb', 22), ('4faw', 2)),(('4ybb', 54), ('4woi', 62)),(('4ybb', 54), ('5nwy', 12)),(('4w2f', 37), ('4faw', 2)),(('4y4o', 43), ('4faw', 2)),(('5afi', 17), ('4woi', 62)),(('5afi', 17), ('5nwy', 12)),(('4ybb', 30), ('1u9s', 1)),(('1vq8', 16), ('3mur', 1)),(('6eri', 16), ('3mum', 1)),(('3cc7', 17), ('3mum', 1)),(('3ccm', 16), ('3mum', 1)),(('4ybb', 12), ('5ngm', 4)),(('4w2g', 52), ('4faw', 1)),(('4y4o', 25), ('6ek0', 7)),(('5wfs', 19), ('5ngm', 4)),(('4ybb', 12), ('4ybb', 7)),(('5wfs', 19), ('4ybb', 7)),(('4y4o', 38), ('4faw', 1)),(('4ybb', 21), ('4faw', 1)),(('5afi', 17), ('4faw', 2)),(('4y4o', 58), ('4faw', 2)),(('4y4o', 58), ('4y1n', 1))]
    liste_faux_neg = [(('3cc7', 17), ('3mum', 1))]
    for i in range(len(liste_tout)) :
        for j in range(i+1, len(liste_tout)) :
            #print(compteur)
#             print(liste_tout[i])
            if (liste_tout[i][1], liste_tout[j][1]) in liste_faux_neg or (liste_tout[j][1], liste_tout[i][1]) in liste_faux_neg  :
                
                print(compteur)
#                 print(liste_tout[i])
#                 print(liste_tout[j])
                graphe_commun_max = nx.MultiDiGraph()
                sim_max = 0.0
                for graphe1 in dico_graphes[liste_tout[i][1]] :
                    for graphe2 in dico_graphes[liste_tout[j][1]] :
                        print(graphe2.edges.data())
                        print("nombre d'aretes")
                        print(len([(u,v) for u,v, data in graphe2.edges(data=True) if data["label"] != 'B53']))
                        graphe_commun, sim = comparaison(graphe1, graphe2, "petit rat") 
                        print(sim)
                        if sim > sim_max :
                            sim_max = sim
                            graphe_commun_max = graphe_commun.copy()
                print(len(dico_graphes[liste_tout[j][1]]))
                print(sim_max)
                dico_new.update({(liste_tout[i][1], liste_tout[j][1]) : {"graphe" : graphe_commun_max, "sim" : sim_max} })
                compteur += 1
                
                
    with open("fichier_csv_liaison_near.csv", "w") as fichier_csv :
        csvwriter = csv.writer(fichier_csv)
        csvwriter.writerow(["Paire", "Sim"])
        for cle in dico_new.keys() :
            csvwriter.writerow([cle, dico_new[cle]["sim"]])
            print(cle)
            print(dico_new[cle]["sim"])
    #print(dico_new)
            
def obtention_differents_graphes_avec_liaison_near_v2(liste_num_ARN) :
        ### Recup des occurrences non redondantes et avec une bonne resolution ###
    liste_tout = []
    for elt in liste_num_ARN :
        with open("resolutions.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            resolutions = mon_depickler.load()
             
            with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
                    mon_depickler = pickle.Unpickler(fichier_sortie)
                    liste_a_garder_noms = mon_depickler.load()    
                       
                    for element in liste_a_garder_noms :
                        if resolutions[element[0]] <= 3 and element[0] != '6hrm' :
                            liste_tout.append((elt, element))
    
    parties_liste_taille_1_9 = []
    liste = []
    for i in range(1,10) :
        liste.append(i)
        parties_liste_taille_1_9.append(partiesliste(liste))
    print(parties_liste_taille_1_9)
    dico_graphes = {}
    compter = 0
    compteur = 0
    for i in range(len(liste_tout)) :
        #if liste_tout[i][1] == ('5ngm', 3) :
            with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+liste_tout[i][1][0]+"_"+str(liste_tout[i][1][1])+"_5.pickle", 'rb') as fichier1 :
                mon_depickler_1 = pickle.Unpickler(fichier1)
                graphe1 = mon_depickler_1.load()
                
                with open("Graphs/%s.pickle"%liste_tout[i][1][0], 'rb') as fichier_tout :
                    mon_depickler_graphes = pickle.Unpickler(fichier_tout)
                    graphe = mon_depickler_graphes.load()
                
                    liste_near = enum_graphe_liaison_near(graphe1, liste_tout[i])
                    #nb_near = nb_liaison_near(graphe1)
                    print(liste_tout[i])
                    print(len(liste_near))
                    print(liste_near)
                    print("compteur")
                    print(compteur)
                
                    if len(liste_near) != 0 :
                        if len(liste_near) == 3 :
                            print("rapoulou")
                            compter += 1
                        liste_temp_graphes = []
                        print("premier")
                        print(graphe1.edges.data())
                        print("modifie")
                        print(len(liste_near))
                        print(len(parties_liste_taille_1_9))
                        print(liste_near)
                        for comb in parties_liste_taille_1_9[len(liste_near)-1] :
                            graphe_copy = graphe.copy()
                            for elt in comb :
                                
                                print(graphe[(graphe1.nodes[liste_near[elt-1][0][0]]["num_ch"], graphe1.nodes[liste_near[elt-1][0][0]]["position"][0])])
                                print(graphe[(graphe1.nodes[liste_near[elt-1][1][0]]["num_ch"], graphe1.nodes[liste_near[elt-1][1][0]]["position"][0])])
    
                                graphe_copy.remove_edge((graphe1.nodes[liste_near[elt-1][0][0]]["num_ch"], graphe1.nodes[liste_near[elt-1][0][0]]["position"][0]), (graphe1.nodes[liste_near[elt-1][0][1]]["num_ch"], graphe1.nodes[liste_near[elt-1][0][1]]["position"][0]))
                                graphe_copy.remove_edge((graphe1.nodes[liste_near[elt-1][1][0]]["num_ch"], graphe1.nodes[liste_near[elt-1][1][0]]["position"][0]), (graphe1.nodes[liste_near[elt-1][1][1]]["num_ch"], graphe1.nodes[liste_near[elt-1][1][1]]["position"][0]))
                            graph_copy = obtenir_extension_un_elt(liste_tout[i][1], graphe_copy, 4)
                            print(graph_copy.edges.data())
                            liste_temp_graphes.append(graph_copy)
                        dico_graphes.update({liste_tout[i][1] : liste_temp_graphes})
                    else :
                        dico_graphes.update({liste_tout[i][1] : [graphe1]})
                    print("apres")
                    print(graphe1.edges.data())
                    print(dico_graphes[liste_tout[i][1]])
                    compteur += 1
                
    with open("Resultats/dico_graphes_liaisons_near_avec_ou_sans.pickle", 'wb') as fichier_liaisons_near :
        mon_pickler = pickle.Pickler(fichier_liaisons_near)
        mon_pickler.dump(dico_graphes)
    print(compter)
    
''' 26/03/20 '''
def obtention_graphe_extension_sans_liaison_near(liste_num_ARN) :
        ### Recup des occurrences non redondantes et avec une bonne resolution ###
    liste_tout = []
    for elt in liste_num_ARN :
        with open("resolutions.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            resolutions = mon_depickler.load()
             
            with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
                    mon_depickler = pickle.Unpickler(fichier_sortie)
                    liste_a_garder_noms = mon_depickler.load()    
                       
                    for element in liste_a_garder_noms :
                        if resolutions[element[0]] <= 3 and element[0] != '6hrm' :
                            liste_tout.append((elt, element))
    
    dico_graphes = {}
    compter = 0
    for i in range(len(liste_tout)) :
        #if liste_tout[i][1] == ('3mum', 1) :
            with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+liste_tout[i][1][0]+"_"+str(liste_tout[i][1][1])+"_5.pickle", 'rb') as fichier1 :
                mon_depickler_1 = pickle.Unpickler(fichier1)
                graphe1 = mon_depickler_1.load()
                
                with open("Graphs/%s.pickle"%liste_tout[i][1][0], 'rb') as fichier_tout :
                    mon_depickler_graphes = pickle.Unpickler(fichier_tout)
                    graphe = mon_depickler_graphes.load()
                    print(type(graphe))
                
                    liste_near = enum_graphe_liaison_near(graphe1, liste_tout[i])
                    #nb_near = nb_liaison_near(graphe1)
                    print(liste_tout[i])
                    print(len(liste_near))
                    print(liste_near)
                
                    if len(liste_near) != 0 :
                        if len(liste_near) == 3 :
                            print("rapoulou")
                            compter += 1
                        print("premier")
                        print(graphe1.edges.data())
                        print("modifie")
                        print(liste_near)
                         
                        graphe_copy = graphe.copy()
                        for elt in liste_near :
                            graphe_copy.remove_edge((graphe1.nodes[elt[0][0]]["num_ch"], graphe1.nodes[elt[0][0]]["position"][0]), (graphe1.nodes[elt[0][1]]["num_ch"], graphe1.nodes[elt[0][1]]["position"][0]))
                            graphe_copy.remove_edge((graphe1.nodes[elt[1][0]]["num_ch"], graphe1.nodes[elt[1][0]]["position"][0]), (graphe1.nodes[elt[1][1]]["num_ch"], graphe1.nodes[elt[1][1]]["position"][0]))
                            if liste_tout[i][1] == ('3mum', 1) :
                                print("ahahahah")
                        graph_copy = obtenir_extension_un_elt(liste_tout[i][1], graphe_copy, 4)
                        dico_graphes.update({liste_tout[i][1] : graph_copy})
                    else :
                        dico_graphes.update({liste_tout[i][1] : graphe1})
                    print("apres")
                    print(graphe1.edges.data())
                    print(dico_graphes[liste_tout[i][1]])
                
    with open("dico_graphes_sans_liaison_near.pickle", 'wb') as fichier_liaisons_near :
        mon_pickler = pickle.Pickler(fichier_liaisons_near)
        mon_pickler.dump(dico_graphes)
    print(compter)


def obtention_comparaison_avec_liaison_near_v2(liste_num_ARN):    
    ### Recup des occurrences non redondantes et avec une bonne resolution ###
    liste_tout = []
    for elt in liste_num_ARN :
        with open("resolutions.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            resolutions = mon_depickler.load()
             
            with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
                    mon_depickler = pickle.Unpickler(fichier_sortie)
                    liste_a_garder_noms = mon_depickler.load()    
                       
                    for element in liste_a_garder_noms :
                        if resolutions[element[0]] <= 3 and element not in liste_pbs :
                            liste_tout.append((elt, element))
    
    with open("dico_graphes_liaisons_near_avec_ou_sans.pickle", 'rb') as fichier_liaisons_near :
        mon_depickler = pickle.Unpickler(fichier_liaisons_near)
        dico_graphes = mon_depickler.load()
                            
        dico_new = {}
        compteur = 0
        #liste_faux_neg = [(('4ybb', 12), ('4u4r', 21)),(('4w2g', 52), ('6az1', 3)),(('4ybb', 12), ('6az1', 3)),(('2zjr', 3), ('6az1', 3)),(('4v67', 7), ('6az1', 3)),(('5wfs', 19), ('6az1', 3)),(('4u4r', 18), ('4y4o', 11)),(('4ybb', 1), ('4y4o', 22)),(('4ybb', 12), ('4y4o', 22)),(('4y4o', 25), ('4y4o', 22)),(('2zjr', 3), ('4y4o', 22)),(('4v67', 7), ('4y4o', 22)),(('6hma', 1), ('4y4o', 22)),(('5wfs', 19), ('4y4o', 22)),(('5dm6', 2), ('5nwy', 17)),(('1vq8', 16), ('3mum', 1)),(('6hma', 1), ('6ek0', 7)),(('6hma', 14), ('4faw', 2)),(('1vq8', 18), ('4faw', 2)),(('5ngm', 15), ('4faw', 2)),(('6eri', 17), ('4faw', 2)),(('6eri', 17), ('4y1n', 1)),(('5dm6', 4), ('4faw', 2)),(('4w2g', 52), ('4y4o', 22)),(('6hma', 8), ('4faw', 2)),(('6eri', 12), ('4faw', 2)),(('6qul', 4), ('4y4o', 38)),(('4ybb', 22), ('4faw', 2)),(('4ybb', 54), ('4woi', 62)),(('4ybb', 54), ('5nwy', 12)),(('4w2f', 37), ('4faw', 2)),(('4y4o', 43), ('4faw', 2)),(('5afi', 17), ('4woi', 62)),(('5afi', 17), ('5nwy', 12)),(('4ybb', 30), ('1u9s', 1)),(('1vq8', 16), ('3mur', 1)),(('6eri', 16), ('3mum', 1)),(('3cc7', 17), ('3mum', 1)),(('3ccm', 16), ('3mum', 1)),(('4ybb', 12), ('5ngm', 4)),(('4w2g', 52), ('4faw', 1)),(('4y4o', 25), ('6ek0', 7)),(('5wfs', 19), ('5ngm', 4)),(('4ybb', 12), ('4ybb', 7)),(('5wfs', 19), ('4ybb', 7)),(('4y4o', 38), ('4faw', 1)),(('4ybb', 21), ('4faw', 1)),(('5afi', 17), ('4faw', 2)),(('4y4o', 58), ('4faw', 2)),(('4y4o', 58), ('4y1n', 1))]
        #liste_faux_neg = [(('3cc7', 17), ('3mum', 1))]
        for i in range(len(liste_tout)) :
            for j in range(i+1, len(liste_tout)) :
                #print(compteur)
    #             print(liste_tout[i])
                #if (liste_tout[i][1], liste_tout[j][1]) in liste_faux_neg or (liste_tout[j][1], liste_tout[i][1]) in liste_faux_neg  :
                    
                    print(compteur)
    #                 print(liste_tout[i])
    #                 print(liste_tout[j])
                    graphe_commun_max = nx.MultiDiGraph()
                    sim_max = 0.0
                    compteur_1 = 0
                    compteur_2 = 0
                    compteurs = (0,0)
                    for graphe1 in dico_graphes[liste_tout[i][1]] :
                        compteur_2 = 0
                        for graphe2 in dico_graphes[liste_tout[j][1]] :
                            #print(graphe1.edges.data())
                            print(graphe2.edges.data())
                            print("nombre d'aretes")
                            print(len([(u,v) for u,v, data in graphe2.edges(data=True) if data["label"] != 'B53']))
                            
                            graphe_commun, sim, diff = comparaison(graphe1, graphe2, "petit rat") 
                            print(sim)
                            if sim > sim_max :
                                sim_max = sim
                                graphe_commun_max = graphe_commun.copy()
                                compteurs = (compteur_1, compteur_2)
                            compteur_2 += 1
                        compteur_1 += 1
                    print(len(dico_graphes[liste_tout[j][1]]))
                    print(sim_max)
                    dico_new.update({(liste_tout[i][1], liste_tout[j][1]) : {"graphe" : graphe_commun_max, "sim" : sim_max, "nums_graphe" : compteurs} })
                    compteur += 1
                    
                    
        with open("fichier_csv_liaison_near.csv", "w") as fichier_csv :
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["Paire", "Sim"])
            for cle in dico_new.keys() :
                csvwriter.writerow([cle, dico_new[cle]["sim"]])
                print(cle)
                print(dico_new[cle]["sim"])
        
        with open("/media/coline/Maxtor/dico_new_260320_sim_par_branche_0.65_avec_liaison_near_plus_infos_new.pickle", 'wb') as fichier_graphes_tot :
            mon_pickler = pickle.Pickler(fichier_graphes_tot)
            mon_pickler.dump(dico_new)        
        
''' 26/03/20
comparaison en enlevant les liaisons near du modèle (sauf celles du motif of course) '''     
def obtention_comparaison_sans_liaison_near(liste_num_ARN):    
    ### Recup des occurrences non redondantes et avec une bonne resolution ###
    liste_tout = []
    for elt in liste_num_ARN :
        with open("resolutions.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            resolutions = mon_depickler.load()
             
            with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
                    mon_depickler = pickle.Unpickler(fichier_sortie)
                    liste_a_garder_noms = mon_depickler.load()    
                       
                    for element in liste_a_garder_noms :
                        if resolutions[element[0]] <= 3 and element not in liste_pbs :
                            liste_tout.append((elt, element))
    
                            
    dico_new = {}
    compteur = 0
    
    with open("dico_graphes_sans_liaison_near.pickle", 'rb') as fichier_liaisons_near :
        mon_depickler = pickle.Unpickler(fichier_liaisons_near)
        dico_graphes = mon_depickler.load()
    #liste_faux_neg = [(('4ybb', 12), ('4u4r', 21)),(('4w2g', 52), ('6az1', 3)),(('4ybb', 12), ('6az1', 3)),(('2zjr', 3), ('6az1', 3)),(('4v67', 7), ('6az1', 3)),(('5wfs', 19), ('6az1', 3)),(('4u4r', 18), ('4y4o', 11)),(('4ybb', 1), ('4y4o', 22)),(('4ybb', 12), ('4y4o', 22)),(('4y4o', 25), ('4y4o', 22)),(('2zjr', 3), ('4y4o', 22)),(('4v67', 7), ('4y4o', 22)),(('6hma', 1), ('4y4o', 22)),(('5wfs', 19), ('4y4o', 22)),(('5dm6', 2), ('5nwy', 17)),(('1vq8', 16), ('3mum', 1)),(('6hma', 1), ('6ek0', 7)),(('6hma', 14), ('4faw', 2)),(('1vq8', 18), ('4faw', 2)),(('5ngm', 15), ('4faw', 2)),(('6eri', 17), ('4faw', 2)),(('6eri', 17), ('4y1n', 1)),(('5dm6', 4), ('4faw', 2)),(('4w2g', 52), ('4y4o', 22)),(('6hma', 8), ('4faw', 2)),(('6eri', 12), ('4faw', 2)),(('6qul', 4), ('4y4o', 38)),(('4ybb', 22), ('4faw', 2)),(('4ybb', 54), ('4woi', 62)),(('4ybb', 54), ('5nwy', 12)),(('4w2f', 37), ('4faw', 2)),(('4y4o', 43), ('4faw', 2)),(('5afi', 17), ('4woi', 62)),(('5afi', 17), ('5nwy', 12)),(('4ybb', 30), ('1u9s', 1)),(('1vq8', 16), ('3mur', 1)),(('6eri', 16), ('3mum', 1)),(('3cc7', 17), ('3mum', 1)),(('3ccm', 16), ('3mum', 1)),(('4ybb', 12), ('5ngm', 4)),(('4w2g', 52), ('4faw', 1)),(('4y4o', 25), ('6ek0', 7)),(('5wfs', 19), ('5ngm', 4)),(('4ybb', 12), ('4ybb', 7)),(('5wfs', 19), ('4ybb', 7)),(('4y4o', 38), ('4faw', 1)),(('4ybb', 21), ('4faw', 1)),(('5afi', 17), ('4faw', 2)),(('4y4o', 58), ('4faw', 2)),(('4y4o', 58), ('4y1n', 1))]
    #liste_faux_neg = [(('3cc7', 17), ('3mum', 1))]
        for i in range(len(liste_tout)) :
            for j in range(i+1, len(liste_tout)) :
                    #print(compteur)
        #             print(liste_tout[i])
                    #if (liste_tout[i][1], liste_tout[j][1]) in liste_faux_neg or (liste_tout[j][1], liste_tout[i][1]) in liste_faux_neg  :
                        
                        print(compteur)
        #                 print(liste_tout[i])
        #                 print(liste_tout[j])
                        graphe_commun_max = nx.MultiDiGraph()
                        sim_max = 0.0
                        
                        graphe1 = dico_graphes[liste_tout[i][1]]
                        graphe2 = dico_graphes[liste_tout[j][1]]
                        
                        graphe_commun, sim, diff = comparaison(graphe1, graphe2, "petit rat") 
                        print(sim)
                        if sim > sim_max :
                            sim_max = sim
                            graphe_commun_max = graphe_commun.copy()
                        print(len(dico_graphes[liste_tout[j][1]]))
                        print(sim_max)
                        dico_new.update({(liste_tout[i][1], liste_tout[j][1]) : {"graphe" : graphe_commun_max, "sim" : sim_max} })
                        compteur += 1
                    
                    
        with open("fichier_csv_liaison_near.csv", "w") as fichier_csv :
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["Paire", "Sim"])
            for cle in dico_new.keys() :
                csvwriter.writerow([cle, dico_new[cle]["sim"]])
                print(cle)
                print(dico_new[cle]["sim"])
        
        with open("/media/coline/Maxtor/dico_new_260320_sim_par_branche_0.65_sans_liaison_near_changer_modele.pickle", 'wb') as fichier_graphes_tot :
            mon_pickler = pickle.Pickler(fichier_graphes_tot)
            mon_pickler.dump(dico_new)        


def creation_graphe_complet(liste_num_ARN):
    ''' creation du graphe complet pondere par les sim '''
    
    ### Recup des occurrences non redondantes et avec une bonne resolution ###
    liste_tout = []
    for elt in liste_num_ARN :
        with open("resolutions.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            resolutions = mon_depickler.load()
             
            with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
                    mon_depickler = pickle.Unpickler(fichier_sortie)
                    liste_a_garder_noms = mon_depickler.load()    
                       
                    for element in liste_a_garder_noms :
                        if resolutions[element[0]] <= 3 and element[0] != '6hrm' : #and element not in liste_pbs :
                            liste_tout.append((elt, element))
   
    #print(liste_tout)  
    
    
    ### Creation du graphe complet qui va contenir les occurrences et leur sim ###
    graphe_complet = nx.Graph()
                
    compteur = 1
    for elt in liste_tout :
        graphe_complet.add_node(elt[1], type=elt[0], nom=elt[1])
        compteur += 1
    #print(graphe_complet.nodes.data())
    print(graphe_complet.number_of_nodes())
 
    
    #with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_graphe_sim : ## Fichier contenant les graphes communs max et les sim pour chaque paire
    with open("/media/coline/Maxtor/dico_new_100320_sim_par_branche_0.65.pickle", 'rb') as fichier_graphe_sim :    

            mon_depickler_2 = pickle.Unpickler(fichier_graphe_sim)
            dico_graphe_sim = mon_depickler_2.load()
            
            for cle in dico_graphe_sim.keys() :
                if cle[0][0] != '6hrm' and cle[1][0] != '6hrm' : # and cle[0] not in liste_pbs and cle[1] not in liste_pbs :
                    graphe_complet.add_edge(cle[0], cle[1], sim=dico_graphe_sim[cle]["sim"])
                    
            print(list(graphe_complet.edges.data())[0])
            #print(graphe_complet.nodes.data())
            print(graphe_complet.number_of_nodes())
            print(graphe_complet.number_of_edges())
            for noeud, data in graphe_complet.nodes(data=True) :
                print(noeud, data)
                
            if ('6ek0', 8) not in graphe_complet.nodes()  :
                print("ahhhh")
            #exit()
            
    return graphe_complet

def creation_graphe_complet_rmsd(liste_num_ARN):
    ''' creation du graphe complet pondere par les rmsd '''
    
    ### Recup des occurrences non redondantes et avec une bonne resolution ###
    liste_tout = []
    for elt in liste_num_ARN :
        with open("resolutions.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            resolutions = mon_depickler.load()
             
            with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
                    mon_depickler = pickle.Unpickler(fichier_sortie)
                    liste_a_garder_noms = mon_depickler.load()    
                       
                    for element in liste_a_garder_noms :
                        if resolutions[element[0]] <= 3 and element[0] != '6hrm' : # and element not in liste_pbs :
                            liste_tout.append((elt, element))
   
    #print(liste_tout)  
    
    
    ### Creation du graphe complet qui va contenir les occurrences et leur sim ###
    graphe_complet = nx.Graph()
                
    compteur = 1
    for elt in liste_tout :
        graphe_complet.add_node(elt[1], type=elt[0], nom=elt[1])
        compteur += 1
    #print(graphe_complet.nodes.data())
    print(graphe_complet.number_of_nodes())
 
    
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load() 
            
        for cle in rmsd.keys() :
            cle_0 = (cle[0].split("_")[1], int(cle[0].split("_")[2]))
            cle_1 = (cle[1].split("_")[1], int(cle[1].split("_")[2]))
            #print(cle_0)
            #print(cle_1)
            if cle_0[0] != '6hrm' and cle_1[0] != '6hrm' :# and cle_0 not in liste_pbs and cle_1 not in liste_pbs :
                graphe_complet.add_edge(cle_0, cle_1, rmsd=rmsd[cle])
        print(len(rmsd.keys()))        
        print(list(graphe_complet.edges.data())[0])
        #print(graphe_complet.nodes.data())
        print(graphe_complet.number_of_nodes())
        print(graphe_complet.number_of_edges())
            
    return graphe_complet

    

def distrib_rmsd_clusters(liste_num_ARN, clusters, graphe_sim):
    
    liste_verte = [(('1vqo', 12), ('6nd6', 7)),(('4y4o', 40), ('4u4r', 8)),(('6eri', 9), ('4u4r', 8)), (('6hma', 15), ('1vq8', 15)), (('5dm6', 2), ('6eri', 1)), (('4u27', 5), ('1vqp', 19)),(('4u27', 5), ('3ccr', 4)),(('4u27', 5), ('3cc7', 4)),(('5dm6', 7), ('1vq8', 11)),(('5dm6', 7), ('3ccr', 12)),(('5dm6', 7), ('3ccl', 11)),(('5dm6', 7), ('3ccs', 12)),(('5dm6', 7), ('3cc7', 12)),(('5dm6', 7), ('3ccq', 12)),(('5dm6', 7), ('3cce', 9)),(('5dm6', 7), ('3ccu', 11)),(('4ybb', 8), ('1vq8', 11)),(('4ybb', 8), ('3ccr', 12)),(('4ybb', 8), ('3ccl', 11)),(('4ybb', 8), ('3ccs', 12)),(('4ybb', 8), ('3cc7', 12)),(('4ybb', 8), ('3ccq', 12)),(('4ybb', 8), ('3cce', 9)),(('4ybb', 8), ('3ccu', 11)),(('3cc2', 2), ('6hma', 1)),(('5dm6', 13), ('1vqp', 19)),(('5dm6', 13), ('3ccr', 4)),(('5dm6', 13), ('3cc7', 4)), (('6hma', 9), ('4u4r', 8)), (('4y4o', 25), ('1vq8', 5)),(('4y4o', 25), ('3cc2', 2)),(('6hma', 15), ('4u4r', 16)), (('1vqp', 19), ('6eri', 9)), (('3ccr', 4), ('6eri', 9)),(('6eri', 9), ('3cc7', 4)),(('5afi', 15), ('4u4r', 28)),(('6hma', 9), ('1vqp', 19)),(('6hma', 9), ('3ccr', 4)),(('6hma', 9), ('3cc7', 4)),(('1vq8', 11), ('6hma', 6)),(('3ccr', 12), ('6hma', 6)),(('6hma', 6), ('3ccl', 11)),(('6hma', 6), ('3ccs', 12)),(('6hma', 6), ('3cc7', 12)),(('6hma', 6), ('3ccq', 12)),(('6hma', 6), ('3cce', 9)),(('6hma', 6), ('3ccu', 11)),(('4y4o', 40), ('1vqp', 19)),(('4y4o', 40), ('3ccr', 4)),(('4y4o', 40), ('3cc7', 4)),(('4y4o', 3), ('3ccl', 11)),(('4y4o', 3), ('3cce', 9)),(('4v51', 22), ('3ccl', 11)),(('4v51', 22), ('3cce', 9)),(('5e81', 12), ('3ccl', 11)),(('5e81', 12), ('3cce', 9)),(('3ccl', 11), ('4v90', 13)),(('4v90', 13), ('3cce', 9)),(('1vqo', 12), ('4u4r', 15)),(('4ybb', 1), ('1vq8', 5)),(('4y4o', 3), ('3ccr', 12)),(('4y4o', 3), ('3ccs', 12)),(('4y4o', 3), ('3cc7', 12)),(('4y4o', 3), ('3ccu', 11)),(('4v51', 22), ('3ccr', 12)),(('4v51', 22), ('3ccs', 12)),(('4v51', 22), ('3cc7', 12)),(('4v51', 22), ('3ccu', 11)),(('5e81', 12), ('3ccr', 12)),(('5e81', 12), ('3ccs', 12)),(('5e81', 12), ('3cc7', 12)),(('5e81', 12), ('3ccu', 11)),(('3ccr', 12), ('4v90', 13)),(('4v90', 13), ('3ccs', 12)),(('4v90', 13), ('3cc7', 12)),(('4v90', 13), ('3ccu', 11)),(('4y4o', 3), ('3ccq', 12)),(('4v51', 22), ('3ccq', 12)),(('5e81', 12), ('3ccq', 12)),(('4v90', 13), ('3ccq', 12)),(('4ybb', 1), ('3cc2', 2)),(('4y4o', 3), ('1vq8', 11)),(('4v51', 22), ('1vq8', 11)),(('5e81', 12), ('1vq8', 11)),(('1vq8', 5), ('6hma', 1)),(('1vq8', 11), ('4v90', 13)), (('5dm6', 9), ('1vq8', 21)), (('6hma', 2), ('4u4r', 19)), (('4ybb', 6), ('1yhq', 24)), (('6hma', 15), ('3cc2', 6))]
#     # 16S
    liste_verte.extend([(('5nwy', 17), ('4y4o', 4)), (('6ek0', 2), ('6az1', 1)), (('4u27', 54), ('4y4o', 22)), (('4ybb', 35), ('4u3u', 22))])
    ## homologues cherches a la main avec bonne sim et rmsd
    liste_verte.extend([(('4y4o', 13), ('4u3u', 22)),(('3cc2', 1), ('4u4r', 36)),(('1vqo', 14), ('4u4r', 36)),(('6eri', 3), ('4u3u', 2)),(('4ybb', 6), ('4u3u', 2)),(('4y4o', 9), ('4u4r', 11)),(('4v67', 23), ('4u4r', 11)),(('4v67', 23), ('6ek0', 6)),(('4y4o', 9), ('6ek0', 6)),(('4y4o', 8), ('3cc2', 1)),(('3t1y', 8), ('4u4r', 11)),(('3t1y', 8), ('6ek0', 6)),(('2qex', 19), ('4u3u', 2)),(('1yhq', 24), ('4u3u', 2)),(('4y4o', 13), ('6az1', 6)),(('4y4o', 8), ('1vqo', 14)),(('5ngm', 1), ('4u4r', 11)),(('6eri', 15), ('4u4r', 11)),(('5ngm', 1), ('6ek0', 6)),(('6eri', 15), ('6ek0', 6)),(('4ybb', 9), ('4u4r', 11)),(('5j7l', 10), ('4u4r', 11)),(('4u27', 32), ('4u4r', 11)),(('5j7l', 10), ('6ek0', 6)),(('4ybb', 9), ('6ek0', 6)),(('4u27', 32), ('6ek0', 6)),(('1vq8', 19), ('4u3u', 2)),(('4ybb', 35), ('6az1', 6)),(('2vqe', 5), ('4u4r', 11)),(('2vqe', 5), ('6ek0', 6))])
        
    ## homologues trouves parmi les 23S et cie qu'on met ensemble avec notre sim  a 0.75 et pas avec la rmsd a 2
    liste_verte.extend([(('4y4o', 25),('1vq8', 5)),(('4y4o', 25),('3cc2', 2)),(('4ybb', 1),('1vq8', 5)),(('4ybb', 1),('3cc2', 2)), (('6hma', 1),('1vq8', 5)),(('6hma', 1),('3cc2', 2))])
        
    ## homologues du 05/03/20 : cherches dans les bleus totaux faux neg, faux pos et vrais pos
    liste_verte.extend([(('4u27', 54), ('6eri', 14)), (('2r8s', 1), ('1l8v', 1)), (('4y4o', 51), ('4u4r', 21)), (('4ybb', 7), ('4u4r', 21)),(('5ngm', 4), ('4u4r', 21)),(('4y4o', 11), ('4u4r', 24)), (('4y4o', 51), ('6az1', 3)), (('4ybb', 7), ('6az1', 3)), (('5ngm', 4), ('6az1', 3)), (('1vqp', 19), ('4u4r', 8)), (('4u4u', 11), ('6az1', 1)), (('4y4o', 8), ('4u4r', 36)), (('3ccr', 4), ('4u4r', 8)),(('3cc7', 4), ('4u4r', 8)), (('4y4o', 11), ('6eri', 10)), (('1gid', 1), ('1l8v', 1))])
    
    
    dico_pos = {(('5dm6', 3), ('4y4o', 2)): 169.0, (('5dm6', 3), ('4ybb', 5)): 175.0, (('5dm6', 3), ('6eri', 21)): 166.5, (('5dm6', 8), ('4ybb', 28)): 132.0, (('5dm6', 8), ('4y4o', 35)): 123.0, (('5dm6', 8), ('4v51', 12)): 118.5, (('5dm6', 8), ('4wsd', 40)): 125.5, (('5dm6', 8), ('4v67', 15)): 128.5, (('5dm6', 8), ('6eri', 22)): 110.5, (('5dm6', 4), ('4ybb', 22)): 169.0, (('5dm6', 4), ('4w2f', 37)): 247.0, (('5dm6', 4), ('4y4o', 43)): 247.0, (('5dm6', 4), ('6hma', 14)): 196.0, (('5dm6', 4), ('6eri', 12)): 193.0, (('5dm6', 7), ('4ybb', 8)): 157.0, (('5dm6', 7), ('4v51', 22)): 172.5, (('5dm6', 7), ('5e81', 12)): 184.0, (('5dm6', 7), ('6eri', 13)): 213.0, (('5dm6', 7), ('4v90', 13)): 184.0, (('5dm6', 10), ('4ybb', 18)): 120.0, (('5dm6', 10), ('4y4o', 53)): 180.0, (('5dm6', 10), ('4v51', 49)): 180.0, (('5dm6', 10), ('5wdt', 12)): 120.0, (('5dm6', 9), ('4ybb', 54)): 137.0, (('5dm6', 9), ('5afi', 17)): 137.0, (('5dm6', 9), ('6eri', 11)): 117.5, (('5dm6', 2), ('4ybb', 48)): 67.0, (('5dm6', 2), ('6qul', 16)): 67.0, (('4ybb', 6), ('6eri', 3)): 99.5, (('4u27', 5), ('5dm6', 13)): 162.5, (('4u27', 5), ('6eri', 9)): 122.5, (('4ybb', 22), ('4w2f', 37)): 168.5, (('4ybb', 22), ('4y4o', 43)): 168.5, (('4ybb', 22), ('6hma', 14)): 188.5, (('4ybb', 22), ('6eri', 12)): 162.5, (('4ybb', 14), ('5mdv', 21)): 301.0, (('4ybb', 14), ('4v51', 5)): 91.0, (('4ybb', 14), ('6eri', 23)): 91.5, (('4ybb', 14), ('6qul', 18)): 251.0, (('4ybb', 36), ('4y4o', 17)): 141.0, (('4ybb', 36), ('5dm6', 6)): 138.5, (('4ybb', 36), ('4v51', 29)): 144.5, (('4ybb', 36), ('4wsd', 4)): 150.0, (('4ybb', 36), ('4v67', 5)): 159.0, (('4ybb', 36), ('6eri', 2)): 150.5, (('4ybb', 17), ('5dm6', 14)): 187.0, (('4ybb', 17), ('4y4o', 8)): 69.0, (('4ybb', 19), ('4w2g', 53)): 238.0, (('4ybb', 19), ('4y4o', 37)): 231.0, (('4ybb', 19), ('6hma', 17)): 238.0, (('4ybb', 48), ('6qul', 16)): 256.5, (('4ybb', 20), ('4y4o', 34)): 102.5, (('4ybb', 20), ('5dm6', 5)): 113.0, (('4ybb', 20), ('4v67', 43)): 102.5, (('4ybb', 20), ('6hma', 13)): 90.0, (('4ybb', 20), ('5wdt', 5)): 274.0, (('4ybb', 8), ('4v51', 22)): 172.5, (('4ybb', 8), ('5e81', 12)): 184.0, (('4ybb', 8), ('6eri', 13)): 139.0, (('4ybb', 8), ('4v90', 13)): 184.0, (('4u26', 22), ('5j5b', 15)): 260.0, (('4ybb', 37), ('4y4o', 19)): 150.0, (('4ybb', 37), ('6h4n', 7)): 301.0, (('4ybb', 37), ('6hma', 10)): 196.0, (('4ybb', 37), ('6eri', 16)): 115.0, (('4ybb', 37), ('6qul', 4)): 289.5, (('4ybb', 10), ('6eri', 5)): 76.0, (('4ybb', 10), ('6qul', 12)): 257.0, (('4ybb', 28), ('4y4o', 35)): 150.0, (('4ybb', 28), ('4v51', 12)): 147.5, (('4ybb', 28), ('4wsd', 40)): 150.0, (('4ybb', 28), ('4v67', 15)): 150.0, (('4ybb', 28), ('6eri', 22)): 167.0, (('4ybb', 18), ('4y4o', 53)): 98.5, (('4ybb', 18), ('4v51', 49)): 98.5, (('4ybb', 18), ('5wdt', 12)): 301.0, (('4ybb', 24), ('4y4o', 50)): 120.0, (('4ybb', 24), ('5dm6', 11)): 114.0, (('4ybb', 24), ('5wfs', 6)): 301.0, (('4ybb', 24), ('6eri', 24)): 127.5, (('4ybb', 30), ('4y4o', 12)): 77.0, (('4ybb', 30), ('2zjr', 15)): 113.5, (('4ybb', 30), ('6eri', 8)): 114.0, (('4ybb', 30), ('5wfs', 11)): 283.0, (('4ybb', 54), ('5afi', 17)): 301.0, (('4ybb', 54), ('6eri', 11)): 170.0, (('4ybb', 23), ('5dm6', 12)): 238.0, (('4ybb', 23), ('4y4o', 6)): 231.0, (('4ybb', 23), ('4w2f', 4)): 238.0, (('4ybb', 23), ('6hma', 4)): 238.0, (('4ybb', 23), ('6eri', 20)): 220.0, (('4ybb', 12), ('2zjr', 3)): 97.0, (('4ybb', 12), ('4v67', 7)): 96.5, (('4ybb', 12), ('5wfs', 19)): 249.0, (('4ybb', 27), ('5afi', 15)): 301.0, (('4y4o', 23), ('1vq8', 8)): 119.5, (('4y4o', 23), ('1yhq', 9)): 119.5, (('4y4o', 23), ('2qex', 12)): 84.0, (('4y4o', 50), ('5dm6', 11)): 165.0, (('4y4o', 50), ('5ngm', 8)): 144.0, (('4y4o', 50), ('5wfs', 6)): 120.0, (('4y4o', 50), ('6eri', 24)): 146.0, (('4y4o', 40), ('6hma', 9)): 115.0, (('4y4o', 40), ('6eri', 9)): 90.0, (('4y4o', 34), ('5dm6', 5)): 147.5, (('4y4o', 34), ('4v67', 43)): 301.0, (('4y4o', 34), ('6hma', 13)): 73.5, (('4y4o', 34), ('5wdt', 5)): 96.5, (('4y4o', 12), ('2zjr', 15)): 93.5, (('4y4o', 12), ('6eri', 8)): 94.0, (('4y4o', 12), ('5wfs', 11)): 83.5, (('4y4o', 29), ('6hma', 11)): 151.0, (('4w2f', 37), ('4y4o', 43)): 303.0, (('4w2f', 37), ('6hma', 14)): 223.0, (('4w2f', 37), ('1vq8', 18)): 142.0, (('4w2f', 37), ('3ccl', 14)): 142.0, (('4w2f', 37), ('6eri', 12)): 220.0, (('4w2f', 37), ('3cce', 12)): 142.0, (('4y4o', 3), ('6hma', 6)): 175.0, (('4y4o', 3), ('6eri', 13)): 175.0, (('4y4o', 35), ('4v51', 12)): 284.0, (('4y4o', 35), ('4wsd', 40)): 301.0, (('4y4o', 35), ('4v67', 15)): 286.0, (('4y4o', 35), ('6eri', 22)): 147.5, (('4y4o', 28), ('6hma', 8)): 106.5, (('4y4o', 25), ('6hma', 1)): 175.0, (('4y4o', 17), ('5dm6', 6)): 135.0, (('4y4o', 17), ('4v51', 29)): 267.5, (('4y4o', 17), ('4wsd', 4)): 292.0, (('4y4o', 17), ('4v67', 5)): 265.0, (('4y4o', 17), ('6eri', 2)): 135.0, (('4y4o', 19), ('1vq8', 16)): 110.5, (('4y4o', 19), ('1yhq', 17)): 110.5, (('4y4o', 19), ('6h4n', 7)): 141.0, (('4y4o', 19), ('3ccr', 17)): 101.5, (('4y4o', 19), ('6hma', 10)): 199.0, (('4y4o', 19), ('3cd6', 14)): 110.5, (('4y4o', 19), ('3ccs', 16)): 101.5, (('4y4o', 19), ('6eri', 16)): 142.0, (('4y4o', 19), ('3cc7', 17)): 101.5, (('4y4o', 19), ('6qul', 4)): 156.0, (('4y4o', 19), ('3ccq', 18)): 104.0, (('4y4o', 19), ('3ccm', 16)): 110.5, (('4y4o', 19), ('3ccu', 16)): 101.5, (('4y4o', 19), ('2qex', 20)): 110.5, (('4y4o', 47), ('1vq8', 19)): 42.0, (('4y4o', 47), ('1yhq', 24)): 39.5, (('4y4o', 47), ('6eri', 1)): 47.5, (('4y4o', 47), ('2qex', 19)): 42.0, (('4y4o', 47), ('4u4r', 13)): 33.0, (('4w2g', 53), ('4y4o', 37)): 303.0, (('4w2g', 53), ('1vqo', 2)): 140.0, (('4w2g', 53), ('6hma', 17)): 274.0, (('4w2g', 53), ('3ccl', 22)): 131.0, (('4w2g', 53), ('3cd6', 19)): 140.0, (('4w2g', 53), ('3cce', 20)): 131.0, (('4w2g', 53), ('3ccm', 21)): 134.5, (('4y4o', 2), ('4ybb', 5)): 181.5, (('4y4o', 2), ('6eri', 21)): 135.0, (('4y4o', 53), ('4v51', 49)): 292.0, (('4y4o', 53), ('5wdt', 12)): 98.5, (('4ybb', 5), ('6eri', 21)): 134.5, (('5mdv', 21), ('4v51', 5)): 91.0, (('5mdv', 21), ('6eri', 23)): 91.5, (('5mdv', 21), ('6qul', 18)): 251.0, (('4u27', 3), ('5afi', 24)): 301.0, (('4u27', 3), ('6h4n', 24)): 301.0, (('4ybb', 13), ('1vqo', 8)): 30.0, (('4ybb', 13), ('5dm6', 1)): 183.0, (('4ybb', 13), ('1mms', 1)): 125.0, (('4ybb', 13), ('6eri', 4)): 166.0, (('5dm6', 12), ('4y4o', 6)): 292.0, (('5dm6', 12), ('4w2f', 4)): 292.0, (('5dm6', 12), ('6hma', 4)): 274.0, (('5dm6', 12), ('6eri', 20)): 247.0, (('5dm6', 6), ('4v51', 29)): 122.0, (('5dm6', 6), ('4wsd', 4)): 129.5, (('5dm6', 6), ('4v67', 5)): 144.0, (('5dm6', 6), ('6eri', 2)): 122.5, (('5dm6', 14), ('4y4o', 8)): 73.0, (('2zjr', 1), ('4v7l', 44)): 102.5, (('4y4o', 6), ('4w2f', 4)): 303.0, (('4y4o', 6), ('1vqo', 6)): 127.5, (('4y4o', 6), ('6hma', 4)): 267.0, (('4y4o', 6), ('6eri', 20)): 240.0, (('4y4o', 37), ('6hma', 17)): 267.0, (('4w2f', 23), ('6eri', 4)): 202.0, (('4w2f', 23), ('4u4r', 36)): 42.0, (('4y4o', 43), ('6hma', 14)): 223.0, (('4y4o', 43), ('1vq8', 18)): 135.5, (('4y4o', 43), ('3ccl', 14)): 135.5, (('4y4o', 43), ('6eri', 12)): 220.0, (('4y4o', 43), ('3cce', 12)): 135.5, (('6hma', 9), ('6eri', 9)): 154.5, (('6hma', 14), ('1vq8', 18)): 133.5, (('6hma', 14), ('3ccl', 14)): 133.5, (('6hma', 14), ('6eri', 12)): 196.0, (('6hma', 14), ('3cce', 12)): 133.5, (('6hma', 5), ('6eri', 21)): 162.0, (('4w2f', 4), ('6hma', 4)): 274.0, (('4w2f', 4), ('6eri', 20)): 247.0, (('4v51', 12), ('4wsd', 40)): 278.5, (('4v51', 12), ('4v67', 15)): 283.0, (('4v51', 12), ('6eri', 22)): 139.5, (('4v51', 29), ('4wsd', 4)): 267.5, (('4v51', 29), ('4v67', 5)): 275.0, (('4v51', 29), ('6eri', 2)): 142.0, (('4v51', 22), ('5e81', 12)): 289.5, (('4v51', 22), ('6eri', 13)): 154.5, (('4v51', 22), ('4v90', 13)): 280.5, (('4v51', 49), ('5wdt', 12)): 98.5, (('2zjr', 3), ('4v67', 7)): 98.5, (('2zjr', 3), ('5wfs', 19)): 78.5, (('5e81', 12), ('6eri', 13)): 166.0, (('5e81', 12), ('4v90', 13)): 292.0, (('4wsd', 4), ('4v67', 5)): 274.0, (('4wsd', 4), ('6eri', 2)): 138.5, (('4wsd', 40), ('4v67', 15)): 284.0, (('4wsd', 40), ('6eri', 22)): 150.0, (('5ibb', 34), ('6h4n', 23)): 133.0, (('5ibb', 34), ('4v9d', 41)): 133.0, (('4v67', 15), ('6eri', 22)): 146.5, (('1vq8', 1), ('3cc2', 20)): 274.0, (('1vq8', 1), ('2qex', 18)): 295.0, (('1vq8', 3), ('3cc2', 21)): 265.0, (('1vq8', 3), ('3cpw', 3)): 301.0, (('1vq8', 4), ('3cc2', 23)): 301.0, (('1vq8', 5), ('3cc2', 2)): 265.0, (('1vqo', 8), ('5dm6', 1)): 52.5, (('1vqo', 8), ('1mms', 1)): 26.5, (('1vq8', 8), ('1yhq', 9)): 301.0, (('1vq8', 8), ('2qex', 12)): 200.0, (('1vqo', 10), ('1yhq', 6)): 274.0, (('1vq8', 10), ('3cc2', 4)): 292.0, (('1vq8', 10), ('2qex', 24)): 295.0, (('1vq8', 11), ('3ccr', 12)): 301.0, (('1vq8', 11), ('3ccl', 11)): 301.0, (('1vq8', 11), ('3ccs', 12)): 301.0, (('1vq8', 11), ('3cc7', 12)): 301.0, (('1vq8', 11), ('3ccq', 12)): 301.0, (('1vq8', 11), ('3cce', 9)): 301.0, (('1vq8', 11), ('3ccu', 11)): 301.0, (('1vqo', 14), ('3cc2', 1)): 301.0, (('1vq8', 12), ('3cc2', 19)): 265.0, (('1vq8', 12), ('3cpw', 14)): 301.0, (('1vq8', 13), ('3cc2', 14)): 292.0, (('1vq8', 13), ('4v9f', 8)): 283.0, (('1vq8', 13), ('2qex', 8)): 201.0, (('1vq8', 15), ('3cc2', 6)): 265.0, (('1vq8', 16), ('1yhq', 17)): 301.0, (('1vq8', 16), ('3ccr', 17)): 301.0, (('1vq8', 16), ('6hma', 10)): 110.5, (('1vq8', 16), ('3cd6', 14)): 301.0, (('1vq8', 16), ('3ccs', 16)): 301.0, (('1vq8', 16), ('6eri', 16)): 103.0, (('1vq8', 16), ('3cc7', 17)): 301.0, (('1vq8', 16), ('3ccq', 18)): 301.0, (('1vq8', 16), ('3ccm', 16)): 301.0, (('1vq8', 16), ('3ccu', 16)): 301.0, (('1vq8', 16), ('2qex', 20)): 262.0, (('1vq8', 18), ('3ccl', 14)): 301.0, (('1vq8', 18), ('6eri', 12)): 113.0, (('1vq8', 18), ('3cce', 12)): 301.0, (('1vq8', 19), ('1yhq', 24)): 247.0, (('1vq8', 19), ('6eri', 1)): 32.5, (('1vq8', 19), ('2qex', 19)): 295.0, (('1vqo', 24), ('1yij', 19)): 301.0, (('1vq8', 21), ('3cc2', 17)): 274.0, (('1vq8', 22), ('3cc2', 12)): 301.0, (('5dm6', 5), ('4v67', 43)): 138.5, (('5dm6', 5), ('6hma', 13)): 119.5, (('5dm6', 5), ('5wdt', 5)): 113.0, (('5dm6', 11), ('5wfs', 6)): 114.0, (('5dm6', 11), ('6eri', 24)): 140.5, (('5dm6', 13), ('6eri', 9)): 131.0, (('5dm6', 1), ('1mms', 1)): 112.0, (('5dm6', 1), ('6eri', 4)): 204.0, (('1vy7', 20), ('5dm6', 15)): 153.0, (('1vy7', 20), ('6i7v', 27)): 157.0, (('3cc2', 4), ('2qex', 24)): 277.0, (('1yhq', 9), ('2qex', 12)): 191.0, (('3cc2', 14), ('4v9f', 8)): 292.0, (('3cc2', 14), ('2qex', 8)): 192.0, (('1yhq', 17), ('3ccr', 17)): 301.0, (('1yhq', 17), ('6hma', 10)): 110.5, (('1yhq', 17), ('3cd6', 14)): 301.0, (('1yhq', 17), ('3ccs', 16)): 301.0, (('1yhq', 17), ('6eri', 16)): 103.0, (('1yhq', 17), ('3cc7', 17)): 301.0, (('1yhq', 17), ('3ccq', 18)): 301.0, (('1yhq', 17), ('3ccm', 16)): 301.0, (('1yhq', 17), ('3ccu', 16)): 301.0, (('1yhq', 17), ('2qex', 20)): 253.0, (('3cc2', 19), ('3cpw', 14)): 256.0, (('3cc2', 20), ('2qex', 18)): 259.0, (('3cc2', 21), ('3cpw', 3)): 256.0, (('3cc2', 22), ('1vqo', 16)): 292.0, (('3cc2', 22), ('2qex', 16)): 292.0, (('1yhq', 24), ('6eri', 1)): 39.5, (('1yhq', 24), ('2qex', 19)): 247.0, (('5afi', 24), ('6h4n', 24)): 301.0, (('6h4n', 7), ('6hma', 10)): 187.0, (('6h4n', 7), ('6eri', 16)): 106.0, (('6h4n', 7), ('6qul', 4)): 280.5, (('5afi', 17), ('6eri', 11)): 170.0, (('1vqp', 19), ('3ccr', 4)): 301.0, (('1vqp', 19), ('3cc7', 4)): 301.0, (('1vqo', 2), ('6hma', 17)): 131.0, (('1vqo', 2), ('3ccl', 22)): 301.0, (('1vqo', 2), ('3cd6', 19)): 301.0, (('1vqo', 2), ('3cce', 20)): 301.0, (('1vqo', 2), ('3ccm', 21)): 301.0, (('5dm6', 15), ('6i7v', 27)): 121.0, (('4v67', 5), ('6eri', 2)): 144.0, (('4v67', 7), ('5wfs', 19)): 99.0, (('4v67', 43), ('6hma', 13)): 73.5, (('4v67', 43), ('5wdt', 5)): 96.5, (('2zjr', 15), ('6eri', 8)): 80.5, (('2zjr', 15), ('5wfs', 11)): 110.0, (('1vqo', 6), ('6hma', 4)): 127.5, (('1vqo', 6), ('6eri', 20)): 127.5, (('4v51', 5), ('6eri', 23)): 76.0, (('4v51', 5), ('6qul', 18)): 110.0, (('4v9f', 8), ('2qex', 8)): 201.0, (('6h4n', 23), ('4v9d', 41)): 301.0, (('6nd6', 7), ('5wdt', 11)): 172.5, (('3ccr', 4), ('3cc7', 4)): 292.0, (('3ccr', 12), ('3ccl', 11)): 301.0, (('3ccr', 12), ('3ccs', 12)): 292.0, (('3ccr', 12), ('3cc7', 12)): 292.0, (('3ccr', 12), ('3ccq', 12)): 301.0, (('3ccr', 12), ('3cce', 9)): 301.0, (('3ccr', 12), ('3ccu', 11)): 292.0, (('3ccr', 17), ('6hma', 10)): 108.0, (('3ccr', 17), ('3cd6', 14)): 301.0, (('3ccr', 17), ('3ccs', 16)): 292.0, (('3ccr', 17), ('6eri', 16)): 97.0, (('3ccr', 17), ('3cc7', 17)): 292.0, (('3ccr', 17), ('3ccq', 18)): 301.0, (('3ccr', 17), ('3ccm', 16)): 301.0, (('3ccr', 17), ('3ccu', 16)): 292.0, (('3ccr', 17), ('2qex', 20)): 253.0, (('6hma', 6), ('6eri', 13)): 139.0, (('6hma', 13), ('5wdt', 5)): 93.0, (('6hma', 17), ('3ccl', 22)): 122.0, (('6hma', 17), ('3cd6', 19)): 131.0, (('6hma', 17), ('3cce', 20)): 122.0, (('6hma', 17), ('3ccm', 21)): 131.0, (('6hma', 10), ('3cd6', 14)): 110.5, (('6hma', 10), ('3ccs', 16)): 101.5, (('6hma', 10), ('6eri', 16)): 148.0, (('6hma', 10), ('3cc7', 17)): 101.5, (('6hma', 10), ('6qul', 4)): 192.0, (('6hma', 10), ('3ccq', 18)): 108.0, (('6hma', 10), ('3ccm', 16)): 110.5, (('6hma', 10), ('3ccu', 16)): 101.5, (('6hma', 10), ('2qex', 20)): 110.5, (('6hma', 7), ('5ngm', 3)): 301.0, (('6hma', 4), ('6eri', 20)): 265.0, (('1mms', 1), ('6eri', 4)): 108.0, (('3ccl', 11), ('3ccs', 12)): 301.0, (('3ccl', 11), ('3cc7', 12)): 301.0, (('3ccl', 11), ('3ccq', 12)): 301.0, (('3ccl', 11), ('3cce', 9)): 301.0, (('3ccl', 11), ('3ccu', 11)): 301.0, (('3ccl', 14), ('6eri', 12)): 113.0, (('3ccl', 14), ('3cce', 12)): 301.0, (('3ccl', 22), ('3cd6', 19)): 301.0, (('3ccl', 22), ('3cce', 20)): 301.0, (('3ccl', 22), ('3ccm', 21)): 301.0, (('6eri', 1), ('2qex', 19)): 32.5, (('6eri', 8), ('5wfs', 11)): 103.0, (('6eri', 13), ('4v90', 13)): 175.0, (('3cd6', 14), ('3ccs', 16)): 301.0, (('3cd6', 14), ('6eri', 16)): 103.0, (('3cd6', 14), ('3cc7', 17)): 301.0, (('3cd6', 14), ('3ccq', 18)): 301.0, (('3cd6', 14), ('3ccm', 16)): 292.0, (('3cd6', 14), ('3ccu', 16)): 301.0, (('3cd6', 14), ('2qex', 20)): 253.0, (('3cd6', 19), ('3cce', 20)): 301.0, (('3cd6', 19), ('3ccm', 21)): 292.0, (('5wfs', 6), ('6eri', 24)): 127.5, (('3ccs', 12), ('3cc7', 12)): 292.0, (('3ccs', 12), ('3ccq', 12)): 292.0, (('3ccs', 12), ('3cce', 9)): 301.0, (('3ccs', 12), ('3ccu', 11)): 301.0, (('3ccs', 16), ('6eri', 16)): 94.0, (('3ccs', 16), ('3cc7', 17)): 292.0, (('3ccs', 16), ('3ccq', 18)): 292.0, (('3ccs', 16), ('3ccm', 16)): 301.0, (('3ccs', 16), ('3ccu', 16)): 301.0, (('3ccs', 16), ('2qex', 20)): 253.0, (('1hc8', 2), ('1y39', 2)): 148.0, (('6eri', 16), ('3cc7', 17)): 94.0, (('6eri', 16), ('6qul', 4)): 114.0, (('6eri', 16), ('3ccq', 18)): 103.0, (('6eri', 16), ('3ccm', 16)): 103.0, (('6eri', 16), ('3ccu', 16)): 94.0, (('6eri', 16), ('2qex', 20)): 103.0, (('6eri', 23), ('6qul', 18)): 94.0, (('6eri', 5), ('6qul', 12)): 86.0, (('6eri', 12), ('3cce', 12)): 113.0, (('3cc7', 12), ('3ccq', 12)): 292.0, (('3cc7', 12), ('3cce', 9)): 301.0, (('3cc7', 12), ('3ccu', 11)): 292.0, (('3cc7', 17), ('3ccq', 18)): 292.0, (('3cc7', 17), ('3ccm', 16)): 301.0, (('3cc7', 17), ('3ccu', 16)): 292.0, (('3cc7', 17), ('2qex', 20)): 253.0, (('3ccq', 12), ('3cce', 9)): 301.0, (('3ccq', 12), ('3ccu', 11)): 292.0, (('3ccq', 18), ('3ccm', 16)): 301.0, (('3ccq', 18), ('3ccu', 16)): 292.0, (('3ccq', 18), ('2qex', 20)): 253.0, (('3cce', 9), ('3ccu', 11)): 301.0, (('3cce', 20), ('3ccm', 21)): 301.0, (('3ccm', 16), ('3ccu', 16)): 301.0, (('3ccm', 16), ('2qex', 20)): 253.0, (('1vqo', 16), ('2qex', 16)): 295.0, (('3ccu', 16), ('2qex', 20)): 253.0, (('4ybb', 35), ('4y4o', 13)): 166.0, (('6i7v', 39), ('4w2f', 33)): 159.5, (('6i7v', 39), ('4wsd', 1)): 159.5, (('6i7v', 39), ('3t1y', 6)): 159.5, (('6i7v', 39), ('4woi', 43)): 235.0, (('6i7v', 39), ('4v67', 25)): 159.5, (('6i7v', 39), ('4v90', 2)): 159.5, (('6i7v', 39), ('4v8b', 31)): 158.0, (('6i7v', 39), ('4v9h', 14)): 159.5, (('6i7v', 39), ('4u26', 1)): 235.0, (('6i7v', 39), ('5e81', 1)): 159.5, (('4ybb', 46), ('4y4o', 58)): 127.0, (('4ybb', 46), ('5ngm', 15)): 124.0, (('4ybb', 2), ('4y4o', 31)): 142.0, (('4ybb', 2), ('4v8d', 26)): 137.0, (('4ybb', 2), ('4v7l', 15)): 127.0, (('4ybb', 2), ('4u27', 2)): 205.0, (('4ybb', 2), ('5f8k', 24)): 142.0, (('4ybb', 2), ('3t1y', 10)): 142.0, (('4ybb', 2), ('4v67', 6)): 142.0, (('4ybb', 2), ('4v90', 11)): 142.0, (('4ybb', 2), ('6h4n', 17)): 205.0, (('4ybb', 2), ('5ngm', 18)): 151.0, (('4ybb', 2), ('5e81', 10)): 142.0, (('4u27', 10), ('4y4o', 39)): 95.5, (('4u27', 10), ('6eri', 6)): 131.5, (('4ybb', 51), ('5ngm', 9)): 71.0, (('4ybb', 51), ('6hrm', 11)): 175.0, (('4u27', 42), ('4y4o', 26)): 148.0, (('4ybb', 4), ('4y4o', 54)): 179.0, (('4ybb', 11), ('4y4o', 11)): 187.0, (('4ybb', 7), ('4y4o', 51)): 111.0, (('4ybb', 7), ('5ngm', 4)): 166.0, (('4u27', 32), ('4y4o', 9)): 158.0, (('4u27', 32), ('4ybb', 9)): 230.0, (('4u27', 32), ('5j7l', 10)): 230.0, (('4u27', 32), ('2vqe', 5)): 158.0, (('4u27', 32), ('3t1y', 8)): 161.0, (('4u27', 32), ('4v67', 23)): 159.5, (('4u27', 32), ('5ngm', 1)): 155.5, (('5nwy', 17), ('4y4o', 4)): 39.5, (('4u27', 54), ('4y4o', 22)): 74.0, (('4y4o', 58), ('5ngm', 15)): 99.5, (('4y4o', 9), ('4ybb', 9)): 155.0, (('4y4o', 9), ('5j7l', 10)): 155.0, (('4y4o', 9), ('2vqe', 5)): 235.0, (('4y4o', 9), ('3t1y', 8)): 215.0, (('4y4o', 9), ('4v67', 23)): 235.0, (('4y4o', 9), ('5ngm', 1)): 208.0, (('4y4o', 51), ('5ngm', 4)): 136.5, (('4y4o', 31), ('4v8d', 26)): 200.0, (('4y4o', 31), ('4v7l', 15)): 190.0, (('4y4o', 31), ('4u27', 2)): 172.0, (('4y4o', 31), ('5f8k', 24)): 250.0, (('4y4o', 31), ('3t1y', 10)): 223.5, (('4y4o', 31), ('4v67', 6)): 255.0, (('4y4o', 31), ('4v90', 11)): 245.0, (('4y4o', 31), ('6h4n', 17)): 147.0, (('4y4o', 31), ('5ngm', 18)): 202.0, (('4y4o', 31), ('5e81', 10)): 240.0, (('4w2f', 33), ('4wsd', 1)): 276.0, (('4w2f', 33), ('3t1y', 6)): 253.5, (('4w2f', 33), ('4woi', 43)): 159.5, (('4w2f', 33), ('4v67', 25)): 285.0, (('4w2f', 33), ('4v90', 2)): 275.0, (('4w2f', 33), ('4v8b', 31)): 230.0, (('4w2f', 33), ('4v9h', 14)): 280.0, (('4w2f', 33), ('4u26', 1)): 159.5, (('4w2f', 33), ('5e81', 1)): 270.0, (('4y4o', 38), ('4ybb', 21)): 97.5, (('4ybb', 9), ('5j7l', 10)): 235.0, (('4ybb', 9), ('2vqe', 5)): 155.0, (('4ybb', 9), ('3t1y', 8)): 161.0, (('4ybb', 9), ('4v67', 23)): 159.5, (('4ybb', 9), ('5ngm', 1)): 155.0, (('5j7l', 10), ('2vqe', 5)): 155.0, (('5j7l', 10), ('3t1y', 8)): 161.0, (('5j7l', 10), ('4v67', 23)): 159.5, (('5j7l', 10), ('5ngm', 1)): 155.0, (('4y4o', 20), ('5nwy', 22)): 201.0, (('2vqe', 5), ('3t1y', 8)): 215.0, (('2vqe', 5), ('4v67', 23)): 240.0, (('2vqe', 5), ('5ngm', 1)): 213.0, (('4y4o', 39), ('6eri', 6)): 95.0, (('4wsd', 1), ('3t1y', 6)): 247.0, (('4wsd', 1), ('4woi', 43)): 159.5, (('4wsd', 1), ('4v67', 25)): 276.0, (('4wsd', 1), ('4v90', 2)): 266.0, (('4wsd', 1), ('4v8b', 31)): 230.0, (('4wsd', 1), ('4v9h', 14)): 271.0, (('4wsd', 1), ('4u26', 1)): 159.5, (('4wsd', 1), ('5e81', 1)): 270.0, (('4v8d', 26), ('4v7l', 15)): 190.0, (('4v8d', 26), ('4u27', 2)): 137.0, (('4v8d', 26), ('5f8k', 24)): 200.0, (('4v8d', 26), ('3t1y', 10)): 200.0, (('4v8d', 26), ('4v67', 6)): 200.0, (('4v8d', 26), ('4v90', 11)): 200.0, (('4v8d', 26), ('6h4n', 17)): 137.0, (('4v8d', 26), ('5ngm', 18)): 182.0, (('4v8d', 26), ('5e81', 10)): 200.0, (('3t1y', 8), ('4v67', 23)): 215.0, (('4v7l', 15), ('4u27', 2)): 127.0, (('4v7l', 15), ('5f8k', 24)): 190.0, (('4v7l', 15), ('3t1y', 10)): 190.0, (('4v7l', 15), ('4v67', 6)): 190.0, (('4v7l', 15), ('4v90', 11)): 190.0, (('4v7l', 15), ('6h4n', 17)): 127.0, (('4v7l', 15), ('5ngm', 18)): 172.0, (('4v7l', 15), ('5e81', 10)): 190.0, (('4u27', 2), ('5f8k', 24)): 172.0, (('4u27', 2), ('3t1y', 10)): 154.0, (('4u27', 2), ('4v67', 6)): 172.0, (('4u27', 2), ('4v90', 11)): 172.0, (('4u27', 2), ('6h4n', 17)): 210.0, (('4u27', 2), ('5ngm', 18)): 181.0, (('4u27', 2), ('5e81', 10)): 172.0, (('5f8k', 24), ('3t1y', 10)): 218.5, (('5f8k', 24), ('4v67', 6)): 250.0, (('5f8k', 24), ('4v90', 11)): 245.0, (('5f8k', 24), ('6h4n', 17)): 147.0, (('5f8k', 24), ('5ngm', 18)): 202.0, (('5f8k', 24), ('5e81', 10)): 240.0, (('3t1y', 6), ('4woi', 43)): 159.5, (('3t1y', 6), ('4v67', 25)): 253.5, (('3t1y', 6), ('4v90', 2)): 247.0, (('3t1y', 6), ('4v8b', 31)): 230.0, (('3t1y', 6), ('4v9h', 14)): 248.5, (('3t1y', 6), ('4u26', 1)): 159.5, (('3t1y', 6), ('5e81', 1)): 247.0, (('3t1y', 10), ('4v67', 6)): 223.5, (('3t1y', 10), ('4v90', 11)): 217.0, (('3t1y', 10), ('6h4n', 17)): 147.0, (('3t1y', 10), ('5e81', 10)): 217.0, (('4woi', 43), ('4v67', 25)): 159.5, (('4woi', 43), ('4v90', 2)): 159.5, (('4woi', 43), ('4v8b', 31)): 158.0, (('4woi', 43), ('4v9h', 14)): 159.5, (('4woi', 43), ('4u26', 1)): 265.0, (('4woi', 43), ('5e81', 1)): 159.5, (('4v67', 6), ('4v90', 11)): 245.0, (('4v67', 6), ('6h4n', 17)): 147.0, (('4v67', 6), ('5ngm', 18)): 202.0, (('4v67', 6), ('5e81', 10)): 240.0, (('4v67', 23), ('5ngm', 1)): 213.0, (('4v67', 25), ('4v90', 2)): 275.0, (('4v67', 25), ('4v8b', 31)): 230.0, (('4v67', 25), ('4v9h', 14)): 280.0, (('4v67', 25), ('4u26', 1)): 159.5, (('4v67', 25), ('5e81', 1)): 270.0, (('4v90', 2), ('4v8b', 31)): 230.0, (('4v90', 2), ('4v9h', 14)): 275.0, (('4v90', 2), ('4u26', 1)): 159.5, (('4v90', 2), ('5e81', 1)): 270.0, (('4v90', 11), ('6h4n', 17)): 147.0, (('4v90', 11), ('5ngm', 18)): 202.0, (('4v90', 11), ('5e81', 10)): 240.0, (('6h4n', 17), ('5ngm', 18)): 156.0, (('6h4n', 17), ('5e81', 10)): 147.0, (('4v8b', 31), ('4v9h', 14)): 230.0, (('4v8b', 31), ('4u26', 1)): 158.0, (('4v8b', 31), ('5e81', 1)): 230.0, (('4v9h', 14), ('4u26', 1)): 159.5, (('4v9h', 14), ('5e81', 1)): 270.0, (('4u26', 1), ('5e81', 1)): 159.5, (('5ngm', 9), ('6hrm', 11)): 86.0, (('5ngm', 18), ('5e81', 10)): 202.0, (('6ek0', 2), ('6az1', 1)): 57.0}

    
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load() 
    
    
    with open("groupes_homologues_60_2.pickle", 'rb') as fichier_hom :
        mon_depickler = pickle.Unpickler(fichier_hom)
        homologues = mon_depickler.load()
    
    distrib_rmsd_intra = []
    distrib_rmsd_inter = []
    distrib_rmsd_homologues_intra = []
    distrib_rmsd_homologues_inter = []
    
    for cle in rmsd.keys() :
        if rmsd[cle] != None :
            nom1 = (cle[0].split("_")[1], int(cle[0].split("_")[2]))
            nom2 = (cle[1].split("_")[1], int(cle[1].split("_")[2]))
            if nom1[0] != '6hrm' and nom2[0] != '6hrm' :
                trouve = False
                pas_homologues = True
                for groupe in homologues : 
                    if nom1 in groupe and nom2 in groupe :
                        pas_homologues = False
                for c in clusters :
                    if nom1 in c and nom2 in c and not trouve :
                        if pas_homologues and (nom1, nom2) not in liste_verte and (nom2, nom1) not in liste_verte and (nom1, nom2) not in dico_pos.keys() and (nom2, nom1) not in dico_pos.keys():
                            distrib_rmsd_intra.append(rmsd[cle])
    #                         print(nom1, nom2)
    #                         print(graphe_sim.edges[nom1, nom2]["sim"])
                        else :
                            print("rapoulou")
                            distrib_rmsd_homologues_intra.append(rmsd[cle])
                        trouve = True
                if not trouve :
                    if pas_homologues and (nom1, nom2) not in liste_verte and (nom2, nom1) not in liste_verte and (nom1, nom2) not in dico_pos.keys() and (nom2, nom1) not in dico_pos.keys():
                        distrib_rmsd_inter.append(rmsd[cle])
                        print(nom1, nom2)
                        print(rmsd[cle])
                    else :
                        print("rapili")
                        distrib_rmsd_homologues_inter.append(rmsd[cle])
    
    #print(distrib_rmsd_inter)
    for elt in distrib_rmsd_inter : 
        if elt > 20 :
            print(elt)
    print(len(distrib_rmsd_intra))
    print(min(distrib_rmsd_intra))
    print(max(distrib_rmsd_intra))
    fig, ax = plt.subplots()
    axs2 = ax.twinx()
    sns.distplot(distrib_rmsd_inter, kde=False, ax=ax)
    sns.distplot(distrib_rmsd_intra, kde=False, ax=axs2, color="orange")
    ax.set_xlabel("RMSD")
    ax.set_ylabel("Nombre de paires")
    #sns.distplot(distrib_rmsd_homologues_inter, kde=False, norm_hist= True)
    #sns.distplot(distrib_rmsd_homologues_intra, kde=False, norm_hist=True)
    plt.show()


def creation_fichier_gephi_seuil_rmsd(liste_num_ARN, seuil_rmsd):
        ''' creation des fichiers gephi pour visualiser le graphe complet avec sim et rmsd en ponderation des aretes
        et type d'ARN, homologues et clustering perez comme attributs des noeuds '''
        graphe_complet = creation_graphe_complet(liste_num_ARN)

        ### On supprime les aretes en-dessous du seuil ###
        
        with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
            mon_depickler = pickle.Unpickler(fichier_rmsd)
            rmsd = mon_depickler.load()
        
            a_enlever = []
            for u,v,data in graphe_complet.edges(data=True) :
                u_nom = "fichier_%s_%d_taille_4.pdb"%(u[0], u[1])
                v_nom = "fichier_%s_%d_taille_4.pdb"%(v[0], v[1])
                if (u_nom, v_nom) in rmsd.keys() :
                    if rmsd[(u_nom, v_nom)] == None or rmsd[(u_nom, v_nom)] > seuil_rmsd :
                        a_enlever.append((u,v))
                else :
                    if rmsd[(v_nom, u_nom)] == None or rmsd[(v_nom, u_nom)] > seuil_rmsd :
                        a_enlever.append((u,v))
    
            for elt in a_enlever :
                graphe_complet.remove_edge(elt[0], elt[1])
            print(graphe_complet.number_of_nodes())
            print(graphe_complet.number_of_edges())
            
        
                
#                 with open("groupes_homologues_version_4janv_seuil_score_70_rmsd_2.5.pickle", 'rb') as fichier_hom :
#                     mon_depickler = pickle.Unpickler(fichier_hom)
#                     homologues = mon_depickler.load()


                #with open("groupes_homologues_version_5fev_seuil_score_%s_%s_rmsd_%s_%s.pickle"%(60,80,2.5,4.62), 'rb') as fichier_hom :
                #with open("groupes_homologues_version_5fev_seuil_%s_%s_%s.pickle"%(55.5,0.034,0.65), 'rb') as fichier_hom :
            with open("groupes_homologues_60_2.pickle", 'rb') as fichier_pickle :
                    mon_depickler = pickle.Unpickler(fichier_pickle)
                    homologues = mon_depickler.load()
                    print(homologues)

                ### Creation des fichiers Gephi pour visualiser le graphe ###
            with open("/media/coline/Maxtor/fichier_csv_new_data_inter_groupe_seuil_rmsd_%s.csv"%seuil_rmsd,'w') as fichier_csv:
                        csvwriter = csv.writer(fichier_csv)
                        csvwriter.writerow(["source", "target", "weight", "sim"])
                        
                        
                        print(graphe_complet.nodes.data())
                        for u,v,data in graphe_complet.edges(data=True) :
                            u_nom = "fichier_%s_%d_taille_4.pdb"%(u[0], u[1])
                            v_nom = "fichier_%s_%d_taille_4.pdb"%(v[0], v[1])
                            if (u_nom, v_nom) in rmsd.keys() :
                                val_rmsd = rmsd[(u_nom, v_nom)] 
                            else :
                                val_rmsd = rmsd[(v_nom, u_nom)] 
                            
                            if val_rmsd != None :
                                csvwriter.writerow([u,v, round(val_rmsd,2), round(data["sim"],2)])
                            else :
                                csvwriter.writerow([u,v, None, round(data["sim"],2)]) 
                            
                            with open("/media/coline/Maxtor/fichier_csv_new_data_inter_groupe_seuil_rmsd_%s_noeud.csv"%(seuil_rmsd),'w') as fichier_csv:
                                csvwriter2 = csv.writer(fichier_csv)
                                csvwriter2.writerow(["id", "label", "type", "homologues"])
                                     
                                for noeud,data in graphe_complet.nodes(data=True) :
                                    #print(noeud)
                                    num_hom = -1
                                            
                                    compteur = 0
                                    for groupe in homologues :
                                        if noeud in groupe :
                                            num_hom = compteur
                                        compteur += 1 
                                    
                                    
                                    csvwriter2.writerow([noeud, data["nom"], graphe_complet.nodes[noeud]["type"], num_hom])



def creation_fichier_gephi(liste_num_ARN, seuil):
        ''' creation des fichiers gephi pour visualiser le graphe complet avec sim et rmsd en ponderation des aretes
        et type d'ARN, homologues et clustering perez comme attributs des noeuds '''
        graphe_complet = creation_graphe_complet(liste_num_ARN)

        ### On supprime les aretes en-dessous du seuil ###
        
        a_enlever = []
        for u,v,data in graphe_complet.edges(data=True) :
            if data["sim"] < seuil:
                a_enlever.append((u,v))

        for elt in a_enlever :
            graphe_complet.remove_edge(elt[0], elt[1])
        print(graphe_complet.number_of_nodes())
        print(graphe_complet.number_of_edges())
            
            
        with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
            mon_depickler = pickle.Unpickler(fichier_rmsd)
            rmsd = mon_depickler.load() 
            
        with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%3, 'rb') as fichier_rmsd :
            mon_depickler = pickle.Unpickler(fichier_rmsd)
            rmsd.update(mon_depickler.load())
            
        with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%2, 'rb') as fichier_rmsd :
            mon_depickler = pickle.Unpickler(fichier_rmsd)
            rmsd.update(mon_depickler.load()) 
        with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%0, 'rb') as fichier_rmsd :
            mon_depickler = pickle.Unpickler(fichier_rmsd)
            rmsd.update(mon_depickler.load())
                        
            graphe_complet_copie = copy.deepcopy(graphe_complet)
            clusters, dico_relevance = recup_data.clustering_perez.algo_principal(graphe_complet_copie)
            
            
#         with open("/media/coline/Maxtor/clustering_perez_%s_new_data_res_inf_3_avec_modif_030220.pickle"%liste_num_ARN, 'rb') as fichier_clusters :
#             mon_depickler = pickle.Unpickler(fichier_clusters)
#             clusters = mon_depickler.load()
#             
#             with open("/media/coline/Maxtor/clustering_perez_%s_new_data_res_inf_3_avec_modif_relevance_030220.pickle"%liste_num_ARN, 'rb') as fichier_relevance :
#                 mon_depickler = pickle.Unpickler(fichier_relevance)
#                 dico_relevance = mon_depickler.load()
                
#                 with open("groupes_homologues_version_4janv_seuil_score_70_rmsd_2.5.pickle", 'rb') as fichier_hom :
#                     mon_depickler = pickle.Unpickler(fichier_hom)
#                     homologues = mon_depickler.load()


                #with open("groupes_homologues_version_5fev_seuil_score_%s_%s_rmsd_%s_%s.pickle"%(60,80,2.5,4.62), 'rb') as fichier_hom :
                #with open("groupes_homologues_version_5fev_seuil_%s_%s_%s.pickle"%(55.5,0.034,0.65), 'rb') as fichier_hom :
            with open("groupes_homologues_60_2.pickle", 'rb') as fichier_pickle :
                    mon_depickler = pickle.Unpickler(fichier_pickle)
                    homologues = mon_depickler.load()
                    print(homologues)
                    
                    pas_bons = []
#                 with open("groupes_homologues_pas_bon_version_5fev_seuil_%s_%s_%s.pickle"%(55.5,0.034,0.65), 'rb') as fichier_hom_2 :
#                     mon_depickler = pickle.Unpickler(fichier_hom_2)
#                     pas_bons = mon_depickler.load()  
#                     print(pas_bons)
                ### Creation des fichiers Gephi pour visualiser le graphe ###
            with open("/media/coline/Maxtor/fichier_csv_new_data_seuil_%s_inter_groupe_res_3A_030220_60_2.csv"%seuil,'w') as fichier_csv:
                        csvwriter = csv.writer(fichier_csv)
                        csvwriter.writerow(["source", "target", "sim", "rmsd_4", "rmsd_3", "rmsd_2", "rmsd_1"])
                        
                        
                        print(graphe_complet.nodes.data())
                        
                        for u,v,data in graphe_complet.edges(data=True) :
                            val_rmsd = []
                            for i in [4,3,2,0] :
                                u_nom = "fichier_%s_%d_taille_%d.pdb"%(u[0], u[1], i)
                                v_nom = "fichier_%s_%d_taille_%d.pdb"%(v[0], v[1], i)
                                
                                if (u_nom, v_nom) in rmsd.keys() :
                                    val_rmsd.append(rmsd[(u_nom, v_nom)]) 
                                else :
                                    val_rmsd.append(rmsd[(v_nom, u_nom)])
                            
                            csvwriter.writerow([u,v,round(data["sim"],2), round(val_rmsd[0],2) if val_rmsd[0] != None else None, round(val_rmsd[1],2) if val_rmsd[1] != None else None, round(val_rmsd[2],2) if val_rmsd[2] != None else None, round(val_rmsd[3],2) if val_rmsd[3] != None else None]) 
                            
                        with open("/media/coline/Maxtor/fichier_csv_new_data_seuil_%s_inter_groupe_noeud_res_3A_030220_60_2.csv"%(seuil),'w') as fichier_csv:
                                csvwriter2 = csv.writer(fichier_csv)
                                csvwriter2.writerow(["id", "label", "type", "clustering_perez", "relevance", "homologues", "pb_hom"])
                                     
                                for noeud,data in graphe_complet.nodes(data=True) :
                                    #print(noeud)
                                    compteur = 0
                                    num_hom = -1
                                    clustering = []
                                    for cluster in clusters :
                                        #print(clusters)
                                        if noeud in cluster :
                                            clustering.append(compteur)
                                            
                                    
                                            
                                        compteur += 1
                                    compteur = 0
                                    pb_hom = 0
                                    for groupe in homologues :
                                        if noeud in groupe :
                                            num_hom = compteur
                                            
                                            if compteur in pas_bons :
                                                print(noeud)
                                                pb_hom = 1
                                        compteur += 1 
                                    
                                    
                                    csvwriter2.writerow([noeud, data["nom"], graphe_complet.nodes[noeud]["type"], list(clustering), dico_relevance[noeud], num_hom, pb_hom])
                                    del(clustering[:])
''' 16/03/20 '''
def creation_fichier_gephi_sim_rmsd(liste_num_ARN, seuil_sim, seuil_rmsd):
    
    graphe_sim = creation_graphe_complet(liste_num_ARN)
    graphe_rmsd = creation_graphe_complet_rmsd(liste_num_ARN)
    
    graphe_sim_copy = graphe_sim.copy()
    graphe_rmsd_copy = graphe_rmsd.copy()
    genere_graphe_seuil_sim(graphe_sim, seuil_sim)
    genere_graphe_seuil_rmsd(graphe_rmsd, seuil_rmsd)
    
    clustering_perez_sim, dico_relevance_sim = recup_data.clustering_perez.algo_principal(graphe_sim)
    
#     with open("/media/coline/Maxtor/clustering_perez_tot_new_data_100320_sim_par_branche_0.65.pickle", 'rb') as fichier_sortie :
#         mon_depickler = pickle.Unpickler(fichier_sortie)
#         clustering_perez_sim = mon_depickler.load()
#         
#     with open("/media/coline/Maxtor/clustering_perez_dico_relevance_tot_new_data_100320_sim_par_branche_0.65.pickle", 'rb') as fichier_relevance :
#         mon_depickler = pickle.Unpickler(fichier_relevance)
#         dico_relevance_sim = mon_depickler.load()
    
    clustering_perez_rmsd, dico_relevance_rmsd = recup_data.clustering_perez.algo_principal(graphe_rmsd)
    
    with open("groupes_homologues_60_2.pickle", 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        groupes_homologues = mon_depickler.load()
    
    with open("/media/coline/Maxtor/diff_sim_rmsd/fichier_gephi_sim_noeud_%s_%s.csv"%(seuil_sim, seuil_rmsd), 'w') as fichier_sim :
        csvwriter_sim = csv.writer(fichier_sim)
        csvwriter_sim.writerow(["id", "label", "type", "clustering_perez", "homologues", "relevance"])
    
        for noeud, data in graphe_sim.nodes(data = True) :
            num_hom = 0
            compteur_hom = 0
            for groupe in groupes_homologues :
                if noeud in groupe :
                    num_hom = compteur_hom
                compteur_hom += 1
                
            num_cluster = []
            compteur_cluster = 1
            for cluster in clustering_perez_sim :
                if noeud in cluster :
                    num_cluster.append(compteur_cluster)
                compteur_cluster += 1
        
            csvwriter_sim.writerow([noeud, "", graphe_sim_copy.nodes[noeud]["type"], num_cluster, num_hom, round(dico_relevance_sim[noeud],2)])
            
    
    with open("/media/coline/Maxtor/diff_sim_rmsd/fichier_gephi_rmsd_noeud_%s_%s.csv"%(seuil_sim, seuil_rmsd), 'w') as fichier_rmsd :
        csvwriter_rmsd = csv.writer(fichier_rmsd)
        csvwriter_rmsd.writerow(["id", "label", "type", "clustering_perez", "homologues", "relevance"])
    
        for noeud, data in graphe_rmsd.nodes(data = True) :
            num_hom = 0
            compteur_hom = 0
            for groupe in groupes_homologues :
                if noeud in groupe :
                    num_hom = compteur_hom
                compteur_hom += 1
                
            num_cluster = []
            compteur_cluster = 1
            for cluster in clustering_perez_rmsd :
                if noeud in cluster :
                    num_cluster.append(compteur_cluster)
                compteur_cluster += 1
        
            csvwriter_rmsd.writerow([noeud, "", graphe_rmsd_copy.nodes[noeud]["type"], num_cluster, num_hom, round(dico_relevance_rmsd[noeud],2)])
    
    with open("/media/coline/Maxtor/diff_sim_rmsd/fichier_gephi_sim_%s_%s.csv"%(seuil_sim, seuil_rmsd), 'w') as fichier_sim :
        csvwriter_sim = csv.writer(fichier_sim)
        csvwriter_sim.writerow(["source", "target", "sim", "rmsd", "dans_rmsd"])
        
        for u,v,data in graphe_sim.edges(data=True) :
            idem = False
            for cluster in clustering_perez_rmsd :
                if u in cluster and v in cluster :
                    idem = True

            
            if not idem :
                csvwriter_sim.writerow([u,v,round(data["sim"],2), round(graphe_rmsd_copy.edges[u,v]["rmsd"], 2) if graphe_rmsd_copy.edges[u,v]["rmsd"] != None else None, 0])
            else :
                csvwriter_sim.writerow([u,v,round(data["sim"],2), round(graphe_rmsd_copy.edges[u,v]["rmsd"], 2) if graphe_rmsd_copy.edges[u,v]["rmsd"] != None else None, 1])
        
    with open("/media/coline/Maxtor/diff_sim_rmsd/fichier_gephi_rmsd_%s_%s.csv"%(seuil_sim, seuil_rmsd), 'w') as fichier_sim :
        csvwriter_rmsd = csv.writer(fichier_sim)
        csvwriter_rmsd.writerow(["source", "target", "rmsd", "sim", "dans_sim"])
        
        for u,v,data in graphe_rmsd.edges(data=True) :
            idem = False
            for cluster in clustering_perez_sim :
                if u in cluster and v in cluster :
                    idem = True
            
            if not idem :
                csvwriter_rmsd.writerow([u,v,round(data["rmsd"],2) if data["rmsd"] != None else None, round(graphe_sim_copy.edges[u,v]["sim"],2), 0])
            else :
                csvwriter_rmsd.writerow([u,v,round(data["rmsd"],2) if data["rmsd"] != None else None, round(graphe_sim_copy.edges[u,v]["sim"],2), 1])
        

''' execution de toutes les comparaisons des occurrences selectionnees dans liste_tout
et stockage dans un fichier pickle (pour sim par branche > 0.65)'''
def recherche_comparaison_groupe(liste_tout):
    print(len(liste_tout))
    compteur = 0
    dico_new = {}
    
    with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_graphe :    
        mon_depickler = pickle.Unpickler(fichier_graphe)
        dico_graphes = mon_depickler.load()  
    
        compte_diff = 0
        for i in range(len(liste_tout)) :
            elt = liste_tout[i]
            with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt[0]+"_"+str(elt[1])+"_2.pickle", 'rb') as fichier1 :
                mon_depickler_1 = pickle.Unpickler(fichier1)
                graphe1 = mon_depickler_1.load()
                for j in range(i+1, len(liste_tout)) :
                    
                    elt2 = liste_tout[j]
                        
                    print(compteur)
                    
                    if (elt, elt2) in dico_graphes.keys() :
                        sim = dico_graphes[(elt, elt2)]["sim"]
                        graphe = dico_graphes[(elt, elt2)]["graphe"]
                    else :
                        sim = dico_graphes[(elt2, elt)]["sim"]
                        graphe = dico_graphes[(elt2, elt)]["graphe"]
                        
                    if sim >= 0.65 :
                    
                        with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt2[0]+"_"+str(elt2[1])+"_2.pickle", 'rb') as fichier2 :
                                        mon_depickler = pickle.Unpickler(fichier2)
                                        graphe2 = mon_depickler.load()
                                            
                                        graphe_commun_max, sim_max, diff = comparaison(graphe1, graphe2, "petit rat") 
                                        if diff : 
                                            compte_diff += 1
                                        dico_new.update({(elt, elt2) : {"graphe" : graphe_commun_max, "sim" : sim_max} })
                    else :
                        dico_new.update({(elt, elt2) : {"graphe" : graphe, "sim" : sim} })
  
                    compteur += 1
    with open("/media/coline/Maxtor/dico_new_100320_sim_par_branche_0.65.pickle", "wb") as fichier_pickle :
        mon_pickler = pickle.Pickler(fichier_pickle)
        mon_pickler.dump(dico_new)  
    print(compte_diff)   
    
    
''' execution de toutes les comparaisons des occurrences selectionnees dans liste_tout
et stockage dans un fichier pickle (test CWW toutes identiques)'''
def recherche_comparaison_groupe_200320(liste_tout):
    print(len(liste_tout))
    compteur = 0
    dico_new = {}

    compte_diff = 0
    for i in range(len(liste_tout)) :
        elt = liste_tout[i]
        with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt[0]+"_"+str(elt[1])+"_2.pickle", 'rb') as fichier1 :
            mon_depickler_1 = pickle.Unpickler(fichier1)
            graphe1 = mon_depickler_1.load()
            for j in range(i+1, len(liste_tout)) :
                  
                elt2 = liste_tout[j]
                      
                print(compteur)
                  
                with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt2[0]+"_"+str(elt2[1])+"_2.pickle", 'rb') as fichier2 :
                    mon_depickler = pickle.Unpickler(fichier2)
                    graphe2 = mon_depickler.load()
                      
                    graphe_commun_max, sim_max, diff = comparaison(graphe1, graphe2, "petit rat") 
                    if diff : 
                        compte_diff += 1
                    dico_new.update({(elt, elt2) : {"graphe" : graphe_commun_max, "sim" : sim_max} })

    
                    compteur += 1
    with open("/media/coline/Maxtor/dico_new_220320_sim_par_branche_0.65_sans_liaisons_near.pickle", "wb") as fichier_pickle :
        mon_pickler = pickle.Pickler(fichier_pickle)
        mon_pickler.dump(dico_new)  
    print(compte_diff)   

''' recherche si le graphe_commun passe en parametre se trouve ou non dans les graphes selectionnes dans liste_tout
pour cela : on fait la comparaison entre chaque graphe d'extension et le graphe_commun
et si le sous-graphe commun est bien egal au graphe_commun, on ajoute le graphe d'extension a la liste des graphes qui contiennent le graphe_commun
et on retourne la liste a la fin '''     
def recherche_comparaison_groupe_moyen(graphe_commun, liste_tout, pourcentage, num):   
    print(len(liste_tout))
    compteur = 0
    dico_new = {}
    compte = 0
    liste_contient = []

    for i in range(len(liste_tout)) :
        elt = liste_tout[i]
        #if elt == ('1yhq', 17) :
        with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt[0]+"_"+str(elt[1])+"_4.pickle", 'rb') as fichier1 :
                mon_depickler_1 = pickle.Unpickler(fichier1)
                graphe1 = mon_depickler_1.load()
                            
                graphe_commun_max, sim_max, diff = comparaison(graphe1, graphe_commun, "petit rat")
                sim = calcul_sim_aretes_avec_coeff_graphe_moyen(graphe1, graphe_commun, graphe_commun_max, "petit rat",1,1,1)
                
                dico_new.update({(elt) : {"graphe" : graphe_commun_max, "sim" : sim} })
                if sim >= pourcentage :
                    compteur += 1
                    liste_contient.append(elt)
                compte += 1
                print(compte)
            
    print(compteur)
    print(liste_contient)
    for elt in dico_new.keys() :
        if elt in liste_contient :
            print(elt)
            print(dico_new[elt])
    with open(NEW_EXTENSION_PATH_TAILLE+"Graphes_communs_avril_2020/dico_graphes_cluster_%s_0.9.pickle"%num, "wb") as fichier_pickle :
        mon_pickler = pickle.Pickler(fichier_pickle)
        mon_pickler.dump(dico_new)  
         
    with open(NEW_EXTENSION_PATH_TAILLE+"Graphes_communs_avril_2020/liste_contient_%s_0.9.pickle"%num, "wb") as fichier_pickle_2 :
        mon_pickler = pickle.Pickler(fichier_pickle_2)
        mon_pickler.dump(liste_contient)  
    return len(liste_contient)

'''  probleme de doublons a cause de 3 occurrences bizarres '''
def test(liste_tout):
    with open("/media/coline/Maxtor/dico_new_avec_derniere_modif.pickle", "rb") as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        dico_new = mon_depickler.load()
        
        compteur1 = 0
        compteur2 = 0
        for elt in liste_tout :
            if elt != ('6hrm', 27) :
                if (('6hrm', 27), elt) in dico_new.keys() and (elt, ('6hrm', 27)) in dico_new.keys() :
                    del(dico_new[(elt, ('6hrm', 27))])
                    
            if elt != ('6hrm', 11) :
                if (('6hrm', 11), elt) in dico_new.keys() and (elt, ('6hrm', 11)) in dico_new.keys() :
                    del(dico_new[(elt, ('6hrm', 11))])
            
            if elt != ('6hrm', 22) :
                if (('6hrm', 22), elt) in dico_new.keys() and (elt, ('6hrm', 22)) in dico_new.keys() :
                    del(dico_new[(elt, ('6hrm', 22))])
        print(len(dico_new))          
        
        liste = []
        compteur = 0
        for elt in dico_new.keys() :
            if (elt[1], elt[0]) in dico_new.keys() :
                compteur += 1
                print(elt)
                print(dico_new[elt]["sim"])
                liste.append(elt)
                
        for elt in liste :
            del(dico_new[elt])
                
#         for cle in dico_new.keys() :
#             if cle[0] == ('6hrm', 27) :
#                 compteur1 += 1
#             if cle[1] == ('6hrm', 27) :
#                 compteur2 += 1
                
        print(compteur1)
        print(compteur2)
        print(liste)
        print(compteur)
        print(len(dico_new))  
        
        with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_sans_doublons.pickle", "wb") as fichier_pickle :
            mon_pickler = pickle.Pickler(fichier_pickle)
            mon_pickler.dump(dico_new)
if __name__ == '__main__':
#     print("rapalaaaaa")
#     exit()
#     with open("Nouvelles_donnees/fichier_4ybb_6_4.pickle", 'rb') as fichier_graphe1 :
#         mon_depickler = pickle.Unpickler(fichier_graphe1)
#         graphe1 = mon_depickler.load()
#     with open("Nouvelles_donnees/fichier_4y4o_29_2.pickle", 'rb') as fichier_graphe2 :
#         mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
#         graphe2 = mon_depickler_2.load()
#                 
#         sous_graphe_commun, sim, diff = comparaison(graphe1, graphe2, "petit rat")
#         print(diff)
#         print(sim)
#         #exit()
#         dico_test = {(('4ybb', 6), ('4y4o', 29)) : {"graphe" : sous_graphe_commun, "sim" : sim}}
#             
#         with open("test_isomorphisme.pickle", 'wb') as fichier_iso :
#             mon_pickler = pickle.Pickler(fichier_iso)
#             mon_pickler.dump(dico_test)
#         #exit()
#         print("gros tas")
#         voisinage_plus_grand = []
#         dico_rappport_spec = {}
#         liste_rapport_spec = []
#         with open("Nouvelles_donnees/Graphes_communs_avril_2020/fichier_csv_rapport.csv", 'w', newline="") as fichier_csv :
#             csvwriter = csv.writer(fichier_csv)
#             csvwriter.writerow(["Num cluster", "Rapport"])
#             with open("/media/coline/Maxtor/clustering_perez_tot_new_data_100320_sim_par_branche_0.65.pickle", 'rb') as fichier_sortie :
#                 mon_depickler = pickle.Unpickler(fichier_sortie)
#                 clusters = mon_depickler.load()
#                 for fic in os.listdir("Nouvelles_donnees/Graphes_communs_avril_2020") :
#                     if "liste_contient" in fic  and "0.9" in fic :
#                         with open("Nouvelles_donnees/Graphes_communs_avril_2020/"+fic, 'rb') as fichier_liste :
#                             mon_depickler = pickle.Unpickler(fichier_liste)
#                             liste_contient = mon_depickler.load()
#                             
#                             
# #                             if "connexe" not in fic :
# #                                 num_cluster = int(fic.split("_")[2][:len(fic.split("_")[2])-7])
# #                             else :
#                             num_cluster = int(fic.split("_")[2])
#                             print(num_cluster)
#                             if len(liste_contient) > len(clusters[num_cluster-1]) :
#                                 c6hrm = False
#                                 for elt in liste_contient :
#                                     if elt[0] == '6hrm' :
#                                         c6hrm = True
#                                         num_6hrm = elt
#                                 
#                                 if len(liste_contient) - len(clusters[num_cluster-1]) != 1 or not c6hrm : 
#                                     
#                                     print(num_cluster)
#                                     print(clusters[num_cluster-1])
#                                     print(liste_contient)
#                                     print(len(clusters[num_cluster-1]))
#                                     print(len(liste_contient))
#                                     print([x for x in liste_contient if x not in clusters[num_cluster-1]])
#                                     voisinage_plus_grand.append(num_cluster)
#     #                                 print("\n")
#                                 if c6hrm : 
#                                     liste_contient.remove(num_6hrm)
#                                     
#                             for elt in clusters[num_cluster-1] :
#                                 if elt not in liste_contient :
#                                     print(elt)
#                                     print("pas bon")
#                                     print(num_cluster)         
#                                     print(liste_contient)
#                                     liste_contient.append(elt)
#                                     
#                             
#                             print(len(clusters[num_cluster-1])/len(liste_contient))        
#                             dico_rappport_spec.update({num_cluster : len(clusters[num_cluster-1])/len(liste_contient)})
#                             liste_rapport_spec.append(len(clusters[num_cluster-1])/len(liste_contient))
#                             csvwriter.writerow([num_cluster, len(clusters[num_cluster-1])/len(liste_contient)])
#                             
#                                 
#     
#                             
#                 
#     #             print(clusters[13])
#     #             print(clusters[14])                    
#                 print(voisinage_plus_grand)
#                 print(dico_rappport_spec.values())
#                 
#                 sns.distplot(liste_rapport_spec, kde=False)
#                 plt.show()
#                 with open("Nouvelles_donnees/Graphes_communs_avril_2020/dico_graphes_cluster_76.pickle", 'rb') as fichier_dico :
#                     mon_depickler = pickle.Unpickler(fichier_dico)
#                     dico_graphe = mon_depickler.load()
#                      
#                     print(dico_graphe[('4u4r', 18)])
#                     
#                 with open("Nouvelles_donnees/Graphes_communs_avril_2020/liste_contient_76.pickle", 'rb') as fichier_liste :
#                     mon_depickler = pickle.Unpickler(fichier_liste)
#                     liste_contient = mon_depickler.load()
#                     
#                     print(liste_contient)
                
            
                

        #exit()
                    
         
#         dico_test = {(('4ybb', 4), ('6hrm', 3)) : {"graphe" : sous_graphe_commun, "sim" : sim}}
#          
#         with open("/media/coline/Maxtor/dico_test.pickle", 'wb') as fichier_sortie :
#             mon_pickler = pickle.Pickler(fichier_sortie)
#             mon_pickler.dump(dico_test)
#     tps1 = time.time() 
    tps1 = time.time()
    types_arn = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16S_arnm", "arnt_16S"]
    
    obtention_differents_graphes_avec_liaison_near_v2(types_arn)
#     obtention_comparaison_avec_liaison_near_v2(types_arn)
    
    tps2 = time.time()
    print(tps2 - tps1)
    #obtention_graphe_extension_sans_liaison_near(types_arn)
    exit()
#     obtention_comparaison_sans_liaison_near(types_arn)
    #exit()
#     graphe_sim = creation_graphe_complet(types_arn)
#     graphe_sim_copy = graphe_sim.copy()
#     genere_graphe_seuil_sim(graphe_sim, 0.7)
#     clusters, dico_relevance = recup_data.clustering_perez.algo_principal(graphe_sim)
#     distrib_rmsd_clusters(types_arn, clusters, graphe_sim_copy)
    #exit()
#     for seuil_sim in [0.65, 0.7, 0.75, 0.8, 0.85] :
#         for seuil_rmsd in [1.0, 1.5, 2.0, 2.5] :
    #creation_fichier_gephi_sim_rmsd(types_arn, 0.75, 2.5)

    #exit()
    #exit()
#     exit()
    #obtention_comparaison_avec_liaison_near_v2(types_arn)
    
    #exit()
#     obtention_comparaison_avec_liaison_near_v2(types_arn)
#     exit()
# # # #     #creation_graphe_complet(types_arn)
#     creation_fichier_gephi(types_arn, 0.75)
    #creation_fichier_gephi_seuil_rmsd(types_arn, 4)
    
    #distrib_rmsd_clusters(types_arn)

# # #     
# #     
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
#     recherche_comparaison_groupe_200320(liste_tout)
#     
#     tps2 = time.time()
#     print(tps2-tps1)
    #exit()
    tps1 = time.time()
    #recherche_comparaison_groupe_200320(liste_tout)
    tps2 = time.time()
    print("Temps d'execution :")
    print(tps2-tps1)
    #exit()
      
    for fic in os.listdir("Nouvelles_donnees/Graphes_communs_avril_2020/") :
        if "53_4_avec_coord_2" in fic :
            with open("Nouvelles_donnees/Graphes_communs_avril_2020/"+fic, 'rb') as fichier_ecriture :
            #with open(EXTENSION_PATH%"taille_max"+"result_k_max_4_10_toutes_aretes/digraphe_commun_0.7_clustering_perez_groupe_12_taille_4.pickle", 'rb') as fichier_ecriture :
                mon_depickler = pickle.Unpickler(fichier_ecriture)
                graphe_commun = mon_depickler.load()
                
                
                  
                print(graphe_commun.nodes.data())
                
                
                for noeud, data in graphe_commun.nodes(data=True) :
                    print(noeud, data)
                
                #graphe_commun = nx.convert_node_labels_to_integers(graphe_commun, first_label=1)
                
                for noeud, data in graphe_commun.nodes(data=True) :
                    print(noeud, data)
                      
                compter_chaines = [0,0,0,0]
                for noeud, data in graphe_commun.nodes(data=True) :
                    for i in range(1,5) :
                        if i in data["chaine"] and data["type"] != None and data["type"] != -1 :
                            compter_chaines[i-1] += 1
                  
                print(compter_chaines)
                pas_bon = True
                while pas_bon :
                    pas_bon = False
                    chaines = recup_chaines(graphe_commun) 
                      
                    compteur = 1
                    for ch in chaines :
                        if len(ch) < compter_chaines[compteur-1] :
                            pas_bon = True
                            print(ch[len(ch)-1])
                            print(graphe_commun.nodes[2])
                            pos = graphe_commun.nodes[ch[len(ch)-1]]["position"][0]
                            noeud_plus_pres = -1
                            plus_pres = 150
                            for noeud, data in graphe_commun.nodes(data=True) :
                                if noeud not in [1,2,3,4] and noeud != ch[len(ch)-1] and compteur in data["chaine"] and data["type"] != None and data["type"] != -1 and abs(data["position"][0]-pos) < plus_pres :
                                    noeud_plus_pres = noeud
                                    plus_pres = abs(data["position"][0]-pos)
                            if ch[len(ch)-1] == 14 : 
                                print(noeud_plus_pres)
                            #print(noeud_plus_pres)
                            if compteur == 1 or compteur == 4 :
                                graphe_commun.add_edge(ch[len(ch)-1], noeud_plus_pres, label="B53")
                            else :
                                graphe_commun.add_edge(noeud_plus_pres, ch[len(ch)-1], label="B53")
                                  
                        compteur += 1
                        #print(chaines)
                print(chaines)
                  
                graphe_commun.nodes[5]["position"] = [0]
                
                for noeud, data in graphe_commun.nodes(data=True) :
                    print(noeud, data)
                    
                graphe_commun.remove_node(7)
                graphe_commun.remove_node(16)
                graphe_commun.remove_node(9)
                graphe_commun.remove_node(17)
                
                #graphe_commun_clusters.draw_new_data(graphe_commun, "Graphes_communs_avril_2020", 54, 4, True)
                #exit()
                
                
                #exit()
                
                num = fic.split("_")[3]
                print(num)
                #exit()
                nb_idem = []
                for i in np.arange(1.0, -0.1, -0.1) :
                    nb_idem.append(recherche_comparaison_groupe_moyen(graphe_commun, liste_tout, i, num))
                print(nb_idem)
                #print(recherche_comparaison_groupe_moyen(graphe_commun, liste_tout, 0.9, num))
            
            #nb_idem = [12, 12, 20, 79, 175, 260, 321, 373, 413, 418, 418]
                plt.clf()
                ax = plt.gca()
                ax.set_yticks(np.arange(0,421,20))
                ax.set_xticks(np.arange(1.0,-0.02,-0.1))
                ax.set_xlabel("Proportion d'arêtes en commun sur le nombre d'arêtes total du graphe commun moyen")
                ax.set_ylabel("Nombre d'occurrences")
                ax.set_title("Nombre d'occurrences possédant exactement le graphe commun moyen du groupe 49 \n ou une certaine proportion d'arêtes de ce graphe commun ")
                plt.plot(np.arange(1.0, -0.1, -0.1), nb_idem)
                plt.show()  
        
#         
#     
#      print(len(liste_tout))
#      test(liste_tout)
    #recherche_comparaison_groupe(liste_tout)
#     tps2 = time.time()
#     print("temps d'execution : %d secondes"%(tps2-tps1))
#     with open("/media/coline/Maxtor/dico_graphe_sim_rassembles_4.pickle", 'rb') as fichier_sim_sans_doublons_avec_manque :
#             mon_depickler_2 = pickle.Unpickler(fichier_sim_sans_doublons_avec_manque)
#             dico_manque = mon_depickler_2.load()
#             print(len(dico_manque))
#             
#             dico_complet = {}
#             
#             compteur = 0
#             doublons = 0
#             doublons_not_idem = 0
#             doublons_que_sim = 0
#             doublons_que_sim_not_idem = 0
#             for cle in dico_manque.keys() :
#                 if isinstance(cle[0], tuple) :
#                     if (cle[0][0]+"_"+str(cle[0][1]),cle[1][0]+"_"+str(cle[1][1])) in dico_manque.keys() :
#                         doublons += 1
#                         if isinstance(dico_manque[(cle[0][0]+"_"+str(cle[0][1]),cle[1][0]+"_"+str(cle[1][1]))], dict) : 
#                             if dico_manque[cle]["sim"] != dico_manque[(cle[0][0]+"_"+str(cle[0][1]),cle[1][0]+"_"+str(cle[1][1]))]["sim"] :
#                                 doublons_not_idem +=1
#                             else :
#                                 dico_complet.update({cle : {"sim" : dico_manque[cle]["sim"], "graphe" : dico_manque[cle]["graphe"]}})
#                         else :
#                             doublons_que_sim += 1
#                             if dico_manque[cle]["sim"] != dico_manque[(cle[0][0]+"_"+str(cle[0][1]),cle[1][0]+"_"+str(cle[1][1]))] :
#                                 doublons_que_sim_not_idem += 1
#                     else :
#                         dico_complet.update({cle : {"sim" : dico_manque[cle]["sim"], "graphe" : dico_manque[cle]["graphe"]}}) 
#                         
#                 else :
#                     if ((cle[0].split("_")[0], int(cle[0].split("_")[1])), (cle[1].split("_")[0], int(cle[1].split("_")[1]))) not in dico_manque.keys() :
#                         if isinstance(dico_manque[cle], dict) :
#                             compteur += 1
#                             dico_complet.update({((cle[0].split("_")[0], int(cle[0].split("_")[1])), (cle[1].split("_")[0], int(cle[1].split("_")[1]))) : {"sim" : dico_manque[cle]["sim"], "graphe" : dico_manque[cle]["graphe"]}}) 
#             print(compteur)
#             print(doublons)
#             print(doublons_not_idem)
#             print(doublons_que_sim)
#             print(doublons_que_sim_not_idem)
#             
#             compteur_doublons_inverse = 0
#             a_enlever = []
#             for cle in dico_complet.keys() :
#                 if (cle[1], cle[0]) in dico_complet.keys() :
#                     compteur_doublons_inverse += 1
#                     if dico_complet[cle]["sim"] != dico_complet[(cle[1], cle[0])]["sim"] : 
#                         print("rapoulou")
#                     if cle not in a_enlever :
#                         a_enlever.append((cle[1], cle[0]))
#                         
#             
#             print(len(a_enlever))
#             print(a_enlever)
#             
#             for elt in a_enlever :
#                 del(dico_complet[elt])
#             
#             compter_nb_comp = 0
#             
#             compter_nb_comp_pas_la = 0
#             print(len(dico_complet))
#             for i in range(len(liste_tout)) :
#                 for j in range(i+1, len(liste_tout)) :
#                     if (liste_tout[i], liste_tout[j]) not in dico_complet.keys() :
#                         #print("pas la")
#                         compter_nb_comp_pas_la += 1
#                     compter_nb_comp += 1
# #             for cle in dico_complet.keys() :
# #                 print(cle)
#             print(compter_nb_comp)
#             print(compteur_doublons_inverse)
#             print(compter_nb_comp_pas_la)
#             for cle in dico_manque.keys() :
#                 print(cle)
            
#             with open("/media/coline/Maxtor/dico_graphe_sim_rassembles_4.pickle", 'rb') as fichier_sim_sans_doublons_avec_manque :
#                 mon_depickler_2 = pickle.Unpickler(fichier_sim_sans_doublons_avec_manque)
#                 dico_manque = mon_depickler_2.load()
   
   
    #graphe_des_sim_inter_groupes(["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16s_arnm", "arnt_16s"])
#     types_arn = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16s_arnm", "arnt_16s"]
#     graphe_des_sim_inter_groupes_res_3A(types_arn)
    #for type in types_arn :
#     with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_avec_manque_%s_res_3A.pickle"%types_arn, 'rb') as fichier_graphe_sim_sans_doublons_avec_manque :
#             mon_depickler = pickle.Unpickler(fichier_graphe_sim_sans_doublons_avec_manque)
#             graphe_complet = mon_depickler.load() 
#              
#             graphe_complet_new_noms = nx.Graph()
#              
#             for noeud, data in graphe_complet.nodes(data=True) :
#                 graphe_complet_new_noms.add_node(data["nom"])
#                  
#             for u,v,data in graphe_complet.edges(data=True) :
#                 if u != v :
#                     graphe_complet_new_noms.add_edge(graphe_complet.nodes[u]["nom"], graphe_complet.nodes[v]["nom"], **data)
#              
#             print(graphe_complet_new_noms.nodes.data())
#             #print(graphe_complet_new_noms.edges.data())
#              
#             with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_avec_manque_%s_new_noms_res_3A.pickle"%types_arn, 'wb') as fichier_graphe_sim :
#                 mon_pickler = pickle.Pickler(fichier_graphe_sim)
#                 mon_pickler.dump(graphe_complet_new_noms)    
                
#     for typ in types_arn :
#         with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_avec_manque_%s_new_noms.pickle"%typ, 'rb') as fichier_graphe_sim_sans_doublons_avec_manque :
#                 mon_depickler = pickle.Unpickler(fichier_graphe_sim_sans_doublons_avec_manque)
#                 graphe_complet = mon_depickler.load()      
#                  
#                 with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%typ, 'rb') as fichier_repr :
#                     mon_depickler = pickle.Unpickler(fichier_repr)
#                     liste_representant_3a = mon_depickler.load()       
#                     
#                     a_enlever = []
#                     for noeud in graphe_complet.nodes() :
#                         if noeud not in liste_representant_3a :
#                             a_enlever.append(noeud)
#                             
#                     for elt in a_enlever :
#                         graphe_complet.remove_node(elt)
#                     
#                     
#                     print(graphe_complet.number_of_nodes())
#                     
#                     with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_avec_manque_%s_new_noms_res_3a.pickle"%typ, 'wb') as fichier_graphe_sim_sans_doublons_avec_manque_res_3a :
#                         mon_pickler = pickle.Pickler(fichier_graphe_sim_sans_doublons_avec_manque_res_3a)
#                         mon_pickler.dump(graphe_complet)
#                         
#                     ecriture_csv_gephi_graphe(typ)
    
        
#     with open("groupes_%s_homologues_sequences.pickle"%"18S", 'rb') as fichier_homologues :
#             mon_depickler_1 = pickle.Unpickler(fichier_homologues)
#             groupes_homologues = mon_depickler_1.load()
#             nb = 0
#             liste_vu = []
#             print(len(groupes_homologues))
#             compteur = 0
#             for i in range(len(groupes_homologues)) :
#                 for j in range(i+1, len(groupes_homologues)) :
#                     compter = 0
#                     for elt1 in groupes_homologues[i] :
#                         for elt2 in groupes_homologues[j] :
#                             if elt1 == elt2 :
#                                 compter += 1
#                                 print((i,j))
#                                 print(len(groupes_homologues[i]))
#                                 print(len(groupes_homologues[j]))
#                     if compter == len(groupes_homologues[i]) and compter == len(groupes_homologues[j]) :
#                         print("petit rat")
#                     elif compter > 0 :
#                         print("bizarre")
# 
#             print(groupes_homologues)
    
    #obs_resultats("18S")
    
#     big_dico = {}
# #     for elt in os.listdir("/media/coline/Maxtor/Resultats") :
# #         if not os.path.isfile(elt) : 
# #             print(elt)
#     compteur = 0
#     for fic in os.listdir("/media/coline/Maxtor/Resultats/23S") :
#                 with open("/media/coline/Maxtor/Resultats/%s/%s"%("23S",fic), "rb") as fichier_graphe_sim :
#                     mon_depickler = pickle.Unpickler(fichier_graphe_sim)
#                     graphe_sim = mon_depickler.load()
#                      
#                     for cle in graphe_sim.keys() :
#                         if isinstance(graphe_sim[cle], dict): 
#                             if cle not in big_dico.keys() and (cle[1], cle[0]) not in big_dico.keys() :
#                                 big_dico.update({cle : graphe_sim[cle]["sim"]})
#                         else :
#                             if cle not in big_dico.keys() and (cle[1], cle[0]) not in big_dico.keys() :
#                                 big_dico.update({cle : graphe_sim[cle]})
#                 print(compteur)   
#                 compteur += 1
#                                  
#     with open("/media/coline/Maxtor/big_dico_sim_23s.pickle", 'wb') as fichier_dico_sim :
#         mon_pickler = pickle.Pickler(fichier_dico_sim)
#         mon_pickler.dump(big_dico)
        
        
#     with open("/media/coline/Maxtor/Resultats/%s/%s"%("16S","dico_sim_new_algo_16S_homologues_groupe_1_groupe_2_part_1.pickle"), "rb") as fichier_graphe_sim :
#                     mon_depickler = pickle.Unpickler(fichier_graphe_sim)
#                     graphe_sim = mon_depickler.load()
#                     
#                     print(len(graphe_sim))
    
    ''' Comparaison '''
#     with open("groupes_25S_homologues.pickle", "rb") as fichier_homologues :
#         mon_depickler_1 = pickle.Unpickler(fichier_homologues)
#         groupes1_homologues = mon_depickler_1.load()
#         
#     with open("groupes_18S_homologues.pickle", "rb") as fichier2_homologues :
#         mon_depickler_2 = pickle.Unpickler(fichier2_homologues)
#         groupes2_homologues = mon_depickler_2.load()
#         print("petit rat")
#         
#         for groupe in groupes1_homologues :
#             print(len(groupe))
#             
#         for groupe in groupes2_homologues :
#             print(len(groupe))
#             
#         liste_a_faire = []
#         for i in range(len(groupes1_homologues)) :
#         
#             for j in range(len(groupes2_homologues)) :
#                     #if not (i == 1 and j in [2,3,4]) and not (i == 2 and j in [97,98,99,100,101,102,103,104,105,106]) and not (i==3 and j in [4,5]) and not (i==4 and (j in range(41,52) or j in [91,92])) and not (i==6 and j in [89,90,91]) and not (i==7 and j in [89,90,91]) and not (i==8 and j==91) and not (i==10 and j in range(97,107)) and not (i==12 and j == 13) and not (i==14 and j in range(15, 53)) :
#                     liste_a_faire.append((i,j))
#                     #comparaison_homologues((18,19))
#         print(liste_a_faire)
#         for elt in liste_a_faire :
#             comparaison_homologues("28S", "18S", elt[0], elt[1])







#         p = multiprocessing.Pool(8)
#         result = p.map(comparaison_homologues, liste_a_faire)
#         p.close()
#         p.join()
    
    #tps1 = time.time()
#     liste_a_faire = ["23S"]
#     p = multiprocessing.Pool(1)
#     result = p.map(comparaison_homologues, liste_a_faire)
#     p.close()
#     p.join()
#     comparaison_homologues("23S")
#     print(time.time() - tps1)
    
    
    #groupe_entre_deux_chaines("SRP")
    
#     with open("groupes_ARNm_identiques.pickle", 'rb') as fichier_identiques :
#         mon_depickler = pickle.Unpickler(fichier_identiques)
#         groupes_identiques = mon_depickler.load()
#         print(groupes_identiques)
    #distribution_homologues('23S', 11)
#'23S', '28S', '16S', '18S', '25S', '4.5S', '4.8S', "SRP", "ARNt", "Riboswitch", "Ribozyme", "Intron",  "ARNm"
#     liste_types = ["25S"]
#     for elt in liste_types :
#         #comparaison_homologues(elt)
#         obs_resultats(elt)
#         ecarter_doublons(elt)
#         print(elt)
#         time.sleep(5)
#     obs_resultats("Ribozyme")
#     
#     with open(NEW_EXTENSION_PATH_TAILLE+"Resultats/groupes_vraiment_identiques_Ribozyme.pickle", 'rb') as fichier_vraiment_id :
#         mon_depickler = pickle.Unpickler(fichier_vraiment_id)
#         groupes_a_rassembler = mon_depickler.load()
#         
#         print(groupes_a_rassembler)


#     with open("Nouvelles_donnees/Resultats/dico_graphe_sim_new_algo.pickle", 'rb') as fichier_lecture :
#         mon_depickler = pickle.Unpickler(fichier_lecture)
#         dico_graphe = mon_depickler.load()
#                
#         print(len(dico_graphe))
    
    #toutes_comparaison()
    
#     with open("dico_graphe_sim_new_algo_test.pickle", 'wb') as fichier_graphe :
#         mon_pickler = pickle.Pickler(fichier_graphe)
#         dico_graphe_new = []
#              
#         compteur = 0
#         for fic in os.listdir(NEW_EXTENSION_PATH_TAILLE+"/Resultats") :
#             if "graphe" in fic and "dico" not in fic :
#                 print(fic)
#                 with open(NEW_EXTENSION_PATH_TAILLE+"/Resultats/"+fic, 'rb') as fichier :
#                     mon_depickler = pickle.Unpickler(fichier)
#                     dico_graphe_sim = mon_depickler.load()
#                      
# #                     with open("test_%s.pickle"%compteur, 'wb') as fichier_ecriture :
# #                         mon_pickler = pickle.Pickler(fichier_ecriture)
# #                         mon_pickler.dump(dico_graphe_sim["sim"])
#                          
# #                     with open(NEW_EXTENSION_PATH_TAILLE+"/Resultats/test_%s.py"%compteur, 'w') as f:
# #                         f.write(repr(list(dico_graphe_sim["graphe"].edges())))
#                      
#                     dico_graphe_new.append(((fic.split("_")[1]+ "_"+fic.split("_")[2], fic.split("_")[3]+ "_" + fic.split("_")[4][:len(fic.split("_")[4])-7]),dico_graphe_sim["sim"]))
#                            
#                     compteur += 1
#                     #break
#                            
#         print(compteur)
#         print(len(dico_graphe_new)) 
#                 
#         #print(dico_graphe_new)       
#                         
#         mon_pickler.dump(dico_graphe_new)    
    
    