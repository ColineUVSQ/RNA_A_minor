'''
Created on 29 août 2019

@author: coline

Nouvel algo pour comparer les graphes d'extension deux à deux (Methode Branch and Cut)
(version toutes donnees PDB)
'''
import recup_data.clustering_perez
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
    
    return somme_aretes

'''idem que dans calcul_sim '''
def calcul_aretes_communes_avec_coeff(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
    for u,v,data in graphe_commun.edges(data=True) :
        if data["label"] != 'B53' :
            if data["label"] == '0' :
                if coeffa == 1 :
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
                    
    somme_aretes = somme_aretes/2 - 4
    
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
    dico_chaines = [{}]
    for i in range(len(chaines)) :
        if i != 0 :
            dico_chaines.append({})
        for elt in chaines[i] :
            for elt_comp in chaines_comp[i] :
                if graphe.nodes[elt]["type"] == graphe_comp.nodes[elt_comp]["type"] :
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

def meme_type_meme_chaine(noeud1, noeud2, graphe1, graphe2):
    ''' renvoie vrai si les 2 noeuds (issus respectivement de graphe1 et graphe2) sont de meme type et appartiennent a la meme chaine '''
    meme_chaine = False
    for elt in graphe1.nodes[noeud1]["chaine"] :
        if elt in graphe2.nodes[noeud2]["chaine"] :
            meme_chaine = True
    if (graphe1.nodes[noeud1]["type"] == graphe2.nodes[noeud2]["type"] and meme_chaine) or (noeud1 in [1,2,3,4,5] and noeud2 in [1,2,3,4,5]):
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
                voisins = {}
                for voisin1 in graphe1[x] :
                    
                    for edge1 in graphe1[x][voisin1] :
                        if graphe1[x][voisin1][edge1]["label"] != 'B53' :
                            label1 = graphe1[x][voisin1][edge1]["label"]
                            #long_range1 = graphe1[x][voisin1][edge1]["long_range"]
                            for voisin2 in graphe2[y] :
                                for edge2 in graphe2[y][voisin2] :
                                    if graphe2[y][voisin2][edge2]["label"] != 'B53' :
                                        
                                        if label1 == graphe2[y][voisin2][edge2]["label"] : #and long_range1 == graphe2[y][voisin2][edge2]["long_range"] : ## meme label des aretes
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
            
    new_sim_max = verif_sim_branche(graphe_commun_max, sim_max, chaines_1, chaines_2, graphe1, graphe2, 1,1,1, 0.4)
    print(sim_max)
    if new_sim_max != sim_max :
        return graphe_commun_max, new_sim_max, True
    else :
        return graphe_commun_max, new_sim_max, False

def verif_sim_branche(graphe_commun_max, sim_max, chaines_1, chaines_2, graphe1, graphe2, coeffc, coeffa, coeffn, seuil):
    
    
    compteur = 0
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
        if sim <= seuil :
            sim_max = 0.75*sim_max

        compteur += 1
    return sim_max


                

''' Obtention de resultats '''

def nb_liaison_near(graphe):
    nb_near = 0
    for u,v, data in graphe.edges(data=True) :   
        if data["near"] == True and (u,v) not in [(1,2), (2,1), (3,4), (4,3), (1,5), (5,1)] :
            nb_near += 1
            
    return int(nb_near/2)

def enum_graphe_liaison_near(graphe):
    liste_liaison_near = []
    liste_aretes = []
    for u,v, key, data in graphe.edges(keys=True, data=True) :
        if data["near"] and (u,v) not in liste_aretes and (u,v) not in [(1,2), (2,1), (3,4), (4,3), (1,5), (5,1)]:
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
                        if resolutions[element[0]] <= 3 and element not in liste_pbs :
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
 
    
    with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_graphe_sim : ## Fichier contenant les graphes communs max et les sim pour chaque paire
            mon_depickler_2 = pickle.Unpickler(fichier_graphe_sim)
            dico_graphe_sim = mon_depickler_2.load()
            
            for cle in dico_graphe_sim.keys() :
                if cle[0] not in liste_pbs and cle[1] not in liste_pbs :
                    graphe_complet.add_edge(cle[0], cle[1], sim=dico_graphe_sim[cle]["sim"])
                    
            print(list(graphe_complet.edges.data())[0])
            #print(graphe_complet.nodes.data())
            print(graphe_complet.number_of_nodes())
            print(graphe_complet.number_of_edges())
            
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
                        if resolutions[element[0]] <= 3 and element not in liste_pbs :
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
            if cle_0 not in liste_pbs and cle_1 not in liste_pbs and rmsd[cle] != None :
                graphe_complet.add_edge(cle_0, cle_1, rmsd=rmsd[cle])
        print(len(rmsd.keys()))        
        print(list(graphe_complet.edges.data())[0])
        #print(graphe_complet.nodes.data())
        print(graphe_complet.number_of_nodes())
        print(graphe_complet.number_of_edges())
            
    return graphe_complet

    

def distrib_rmsd_clusters(liste_num_ARN):
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load() 
            
    with open("/media/coline/Maxtor/clustering_perez_%s_new_data_res_inf_3_avec_modif_030220.pickle"%liste_num_ARN, 'rb') as fichier_clusters :
        mon_depickler = pickle.Unpickler(fichier_clusters)
        clusters = mon_depickler.load()
    
    with open("groupes_homologues_version_4janv_seuil_score_70_rmsd_2.5.pickle", 'rb') as fichier_hom :
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
            trouve = False
            pas_homologues = True
            for groupe in homologues : 
                if nom1 in groupe and nom2 in groupe :
                    pas_homologues = False
            for c in clusters :
                if nom1 in c and nom2 in c and not trouve :
                    if pas_homologues:
                        distrib_rmsd_intra.append(rmsd[cle])
                    else :
                        print("rapoulou")
                        distrib_rmsd_homologues_intra.append(rmsd[cle])
                    trouve = True
            if not trouve :
                if pas_homologues:
                    distrib_rmsd_inter.append(rmsd[cle])
                else :
                    print("rapili")
                    distrib_rmsd_homologues_inter.append(rmsd[cle])
    
    #print(distrib_rmsd_inter)
    for elt in distrib_rmsd_inter : 
        if elt > 20 :
            print(elt)
    sns.distplot(distrib_rmsd_inter, kde=False, norm_hist=True)
    sns.distplot(distrib_rmsd_intra, kde=False, norm_hist = True)
    sns.distplot(distrib_rmsd_homologues_inter, kde=False, norm_hist= True)
    sns.distplot(distrib_rmsd_homologues_intra, kde=False, norm_hist=True)
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

''' execution de toutes les comparaisons des occurrences selectionnees dans liste_tout
et stockage dans un fichier pickle '''
def recherche_comparaison_groupe(liste_tout):
    print(len(liste_tout))
    compteur = 0
    dico_new = {}

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
                                    
                                graphe_commun_max, sim_max = comparaison(graphe1, graphe2, "petit rat") 
                                dico_new.update({(elt, elt2) : {"graphe" : graphe_commun_max, "sim" : sim_max} })
                compteur += 1
    with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", "wb") as fichier_pickle :
        mon_pickler = pickle.Pickler(fichier_pickle)
        mon_pickler.dump(dico_new)       

''' recherche si le graphe_commun passe en parametre se trouve ou non dans les graphes selectionnes dans liste_tout
pour cela : on fait la comparaison entre chaque graphe d'extension et le graphe_commun
et si le sous-graphe commun est bien egal au graphe_commun, on ajoute le graphe d'extension a la liste des graphes qui contiennent le graphe_commun
et on retourne la liste a la fin '''     
def recherche_comparaison_groupe_moyen(graphe_commun, liste_tout):   
    print(len(liste_tout))
    compteur = 0
    dico_new = {}
    compte = 0
    liste_contient = []

    for i in range(len(liste_tout)) :
        elt = liste_tout[i]
        with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt[0]+"_"+str(elt[1])+"_2.pickle", 'rb') as fichier1 :
            mon_depickler_1 = pickle.Unpickler(fichier1)
            graphe1 = mon_depickler_1.load()
                        
            graphe_commun_max, sim_max = comparaison(graphe1, graphe_commun, "petit rat")
            sim = calcul_sim_aretes_avec_coeff_graphe_moyen(graphe1, graphe_commun, graphe_commun_max, "petit rat",1,1,1)
            
            dico_new.update({(elt) : {"graphe" : graphe_commun_max, "sim" : sim} })
            if sim == 1.0 :
                compteur += 1
                liste_contient.append(elt)
            compte += 1
            print(compte)
    print(compteur)
    print(liste_contient)
    for elt in dico_new.keys() :
        print(elt)
        print(dico_new[elt])
#     with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_612.pickle", "wb") as fichier_pickle :
#         mon_pickler = pickle.Pickler(fichier_pickle)
#         mon_pickler.dump(dico_new)   

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
    
    
    with open("Nouvelles_donnees/fichier_5wfs_11_2.pickle", 'rb') as fichier_graphe1 :
        mon_depickler = pickle.Unpickler(fichier_graphe1)
        graphe1 = mon_depickler.load()
    with open("Nouvelles_donnees/fichier_4u3u_7_2.pickle", 'rb') as fichier_graphe2 :
        mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
        graphe2 = mon_depickler_2.load()
         
        sous_graphe_commun, sim, diff = comparaison(graphe1, graphe2, "petit rat")
        print(diff)
        print(sim)
        exit()
        
#         
#         dico_test = {(('4ybb', 4), ('6hrm', 3)) : {"graphe" : sous_graphe_commun, "sim" : sim}}
#          
#         with open("/media/coline/Maxtor/dico_test.pickle", 'wb') as fichier_sortie :
#             mon_pickler = pickle.Pickler(fichier_sortie)
#             mon_pickler.dump(dico_test)
#     tps1 = time.time() 
    types_arn = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16S_arnm", "arnt_16S"]
    obtention_comparaison_avec_liaison_near(types_arn)
    exit()
# # #     #creation_graphe_complet(types_arn)
    creation_fichier_gephi(types_arn, 0.75)
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
    
    recherche_comparaison_groupe(liste_tout)
    exit()
      
    #with open("Nouvelles_donnees/Graphes_communs/graphe_commun_cluster_38_4_avec_coord.pickle", 'rb') as fichier_ecriture :
    with open(EXTENSION_PATH%"taille_max"+"result_k_max_4_10_toutes_aretes/digraphe_commun_0.7_clustering_perez_groupe_12_taille_4.pickle", 'rb') as fichier_ecriture :
        mon_depickler = pickle.Unpickler(fichier_ecriture)
        graphe_commun = mon_depickler.load()
          
        print(graphe_commun.nodes.data())
        
        
        for noeud, data in graphe_commun.nodes(data=True) :
            print(noeud, data)
        
        graphe_commun = nx.convert_node_labels_to_integers(graphe_commun, first_label=1)
        
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
        
        #exit()
        recherche_comparaison_groupe_moyen(graphe_commun, liste_tout)  
        
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
    
    