'''
Created on 10 déc. 2018

@author: coline

Calculs de toutes les metriques testees
et generation de fichiers pour les stocker (version CaRNAval)

'''
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import csv
import os
import numpy as np
from recup_data.constantes import EXTENSION_PATH, EXTENSION_PATH_TAILLE


def calcul_sim_raymond_avec_motif(graphe1, graphe2, graphe_commun):
    
    sim = ((graphe_commun.number_of_nodes() + graphe_commun.number_of_edges())*(graphe_commun.number_of_nodes() + graphe_commun.number_of_edges()))/((graphe1.number_of_nodes()+graphe1.number_of_edges())*(graphe2.number_of_nodes()+graphe2.number_of_edges()))
    
    return sim

def calcul_sim_raymond_sans_motif(graphe1, graphe2, graphe_commun):
    
    sim = ((graphe_commun.number_of_nodes() - 5 + graphe_commun.number_of_edges()-4)*(graphe_commun.number_of_nodes()-5 + graphe_commun.number_of_edges()-4))/((graphe1.number_of_nodes()-5+graphe1.number_of_edges()-4)*(graphe2.number_of_nodes()-5+graphe2.number_of_edges()-4))
    
    return sim

def calcul_sim_raymond_non_cov_avec_motif(graphe1, graphe2, graphe_commun) :
    compteur_arc = 0
    for (u, v, keys, t) in graphe1.edges(data="label", keys = True) :
        if t == "B53" :
            compteur_arc += 1
    compteur_arete_1 = (graphe1.number_of_edges() - compteur_arc)/2
    
    compteur_arc = 0
    for (u, v, keys, t) in graphe2.edges(data="label", keys = True) :
        if t == "B53" :
            compteur_arc += 1
    compteur_arete_2 = (graphe2.number_of_edges() - compteur_arc)/2
    
    
    compteur_arc = 0
    for (u, v, keys, t) in graphe_commun.edges(data="type", keys = True) :
        if t == "COV" :
            compteur_arc += 1
    compteur_arete_commun = graphe_commun.number_of_edges() - compteur_arc
    
    sim = ((graphe_commun.number_of_nodes() + compteur_arete_commun)*(graphe_commun.number_of_nodes() + compteur_arete_commun))/((graphe1.number_of_nodes()+compteur_arete_1)*(graphe2.number_of_nodes()+compteur_arete_2))
    
    return sim

def calcul_sim_raymond_non_cov_sans_motif(graphe1, graphe2, graphe_commun) :
    compteur_arc = 0
    for (u, v, keys, t) in graphe1.edges(data="label", keys = True) :
        if t == "B53" :
            compteur_arc += 1
    compteur_arete_1 = (graphe1.number_of_edges() - compteur_arc)/2 - 4
    
    compteur_arc = 0
    for (u, v, keys, t) in graphe2.edges(data="label", keys = True) :
        if t == "B53" :
            compteur_arc += 1
    compteur_arete_2 = (graphe2.number_of_edges() - compteur_arc)/2 - 4
    
    
    compteur_arc = 0
    for (u, v, keys, t) in graphe_commun.edges(data="type", keys = True) :
        if t == "COV" :
            compteur_arc += 1
    compteur_arete_commun = graphe_commun.number_of_edges() - compteur_arc -4
    
    sim = ((graphe_commun.number_of_nodes() - 5 + compteur_arete_commun)*(graphe_commun.number_of_nodes()-5 + compteur_arete_commun))/((graphe1.number_of_nodes()-5+compteur_arete_1)*(graphe2.number_of_nodes()-5+compteur_arete_2))
    
    return sim

def calcul_sim_non_cov_avec_motif(graphe1, graphe2, graphe_commun):
    
    compteur_arc = 0
    for (u, v, keys, t) in graphe1.edges(data="label", keys = True) :
        if t == "B53" :
            compteur_arc += 1

    compteur_arete_1 = (graphe1.number_of_edges() - compteur_arc)/2
#     if element1 == "fichier_5DM6_X_197_1" and element2 == "fichier_5DM6_X_48_9" :
#         print(compteur_arc)
    
    
    compteur_arc = 0
    for (u, v, keys, t) in graphe2.edges(data="label", keys = True) :
        if t == "B53" :
            compteur_arc += 1
    compteur_arete_2 = (graphe2.number_of_edges() - compteur_arc)/2
#     if element1 == "fichier_5DM6_X_197_1" and element2 == "fichier_5DM6_X_48_9" :
#         print(compteur_arc)

    compteur_arc = 0
    for (u, v, keys, t) in graphe_commun.edges(data="type", keys = True) :
        if t == "COV" :
            compteur_arc += 1
    compteur_arete_commun = graphe_commun.number_of_edges() - compteur_arc
    

    
    
#     if element1 == "fichier_5DM6_X_197_1" and element2 == "fichier_5DM6_X_48_9" :
#         print(compteur_arete_commun)
#         print(compteur_arete_1)
#         print(compteur_arete_2)
    
#     print(compteur_arete_1)
#     print(compteur_arete_2)
#     print(compteur_arete_commun)
    
    sim = compteur_arete_commun/max(compteur_arete_1, compteur_arete_2)
    return sim, compteur_arete_commun


def calcul_sim_non_cov_sans_motif(graphe1, graphe2, graphe_commun):
    compteur_arc = 0
    for (u, v, keys, t) in graphe1.edges(data="label", keys = True) :
        if t == "B53" :
            compteur_arc += 1

    compteur_arete_1 = (graphe1.number_of_edges() - compteur_arc)/2 - 4

    
    compteur_arc = 0
    for (u, v, keys, t) in graphe2.edges(data="label", keys = True) :
        if t == "B53" :
            compteur_arc += 1
    compteur_arete_2 = (graphe2.number_of_edges() - compteur_arc)/2 - 4


    compteur_arc = 0
    for (u, v, keys, t) in graphe_commun.edges(data="type", keys = True) :
        if t == "COV" :
            compteur_arc += 1
    compteur_arete_commun = graphe_commun.number_of_edges() - compteur_arc - 4

    if max(compteur_arete_1, compteur_arete_2) > 0 :
        sim = compteur_arete_commun/max(compteur_arete_1, compteur_arete_2)
    else :
        sim = 0.0
    return sim, compteur_arete_commun, compteur_arete_1, compteur_arete_2

def calcul_sim_non_cov_sans_motif_par_chaine(graphe1, graphe2, graphe_commun, chaines_1, chaines_2, chaines_commun, i):
    
    compteur_arete_1 = 0
    for (u, v, t) in graphe1.edges(data="label") :
        if (u in chaines_1[i] and v in chaines_1[i]) :
            if t != 'B53' :
                print((u,v))
                compteur_arete_1 += 1
                
#     if i == 2 or i == 3 :
#         compteur_arete_1 = compteur_arete_1 - 2
#     else :
#         compteur_arete_1 = compteur_arete_1 - 4
    compteur_arete_1 = compteur_arete_1/2
#     if element1 == "fichier_5DM6_X_197_1" and element2 == "fichier_5DM6_X_48_9" :

    compteur_arete_2 = 0
    for (u, v, t) in graphe2.edges(data="label") :
        if (u in chaines_2[i] and v in chaines_2[i]) :
            if t != 'B53' :
                compteur_arete_2 += 1
#     if i == 2 or i == 3 :
#         compteur_arete_2 = compteur_arete_2 - 2
#     else :
#         compteur_arete_2 = compteur_arete_2 - 4
    compteur_arete_2 = compteur_arete_2/2
#     if element1 == "fichier_5DM6_X_197_1" and element2 == "fichier_5DM6_X_48_9" :
    
    compteur_arete_commun = 0
    for (u, v, t) in graphe_commun.edges(data="type") :
        if (u in chaines_commun[i] and v in chaines_commun[i]) :
            if t != 'COV' :
                compteur_arete_commun += 1
#     if i == 2 or i == 3 :
#         compteur_arete_commun = compteur_arete_commun - 1
#     else :
#         compteur_arete_commun = compteur_arete_commun - 2
        
    #compteur_arete_commun = compteur_arete_commun/2
    
#     if element1 == "fichier_5DM6_X_197_1" and element2 == "fichier_5DM6_X_48_9" :
#         print(compteur_arete_commun)
#         print(compteur_arete_1)
#         print(compteur_arete_2)
    
#     print(compteur_arete_1)
#     print(compteur_arete_2)
#     print(compteur_arete_commun)
    if max(compteur_arete_1, compteur_arete_2) > 0 :
        sim = compteur_arete_commun/max(compteur_arete_1, compteur_arete_2)
    else :
        sim = 0.0
    return sim, compteur_arete_1, compteur_arete_2, compteur_arete_commun

def calcul_sommets_aretes_grands_graphes(graphe):
    
    somme_aretes = 0
    for u,v,data in graphe.edges(data=True) :
        if data["label"] != "B53" :
            somme_aretes += 1
    somme_aretes = somme_aretes/2 - 4
    
    return somme_aretes + graphe.number_of_nodes() - 5

def calcul_sommets_aretes_grands_graphes_commun(graphe):
    somme_aretes = 0
    for u,v,data in graphe.edges(data=True) :
        if data["type"] != "COV" :
            somme_aretes += 1
    somme_aretes = somme_aretes - 4
    
    return somme_aretes + graphe.number_of_nodes() - 5
    

'''Calcul metrique sommets_aretes_ponderees V2
pour les sommets de type 1, si le sommet n'a pas de voisin non artificiel, on met un coeff qui vaut 3, pour compter les deux nucleotides 
et l' arete, et si le sommet a un voisin non artificiel dans l'extension, on met un coeff qui vaut 1 pour compter le sommet uniquement
(l'arete et l'autre sommet seront comptes separement, chacun de leur cote) '''
def calcul_sommets_aretes_graphe(graphe, cle):
    somme_sommets = 0
    for noeud, data in graphe.nodes(data=True) :
        if data["type"] == 1 : ## si de type 1 on multiplie par le coeff 3
            if len(graphe[noeud]) <= 1 :
                somme_sommets += 3*data["poids"] ## cas ou le sommet de type 1 est lie a un nt qui nest pas present dans le graphe
            else  :
                voisin = -1
                for u,v, label in graphe.edges(data="label") :
                    if noeud == u or noeud == v :
                        if label != 'B53' :
                            if u == noeud :
                                voisin = v
                            if v == noeud :
                                voisin = u
                if voisin != -1 :
                    if graphe.nodes[voisin]["type"] == 1 :            
                        somme_sommets += 3*data["poids"]/2 ## cas ou le sommet de type 1 est lie a un sommet qui est present dans le graphe et qui est aussi de type 1
                    else :
                        somme_sommets += data["poids"] 
                        
                else :
                    somme_sommets += data["poids"] 
#                 except KeyError :
#                     with open("fichier_pas_bon.txt", 'a') as f :
#                         f.write(str(cle) + "\n")
        else :
            somme_sommets += data["poids"]
            
    somme_sommets -= 5

    compteur_arete = 0
    for u,v, data in graphe.edges(data=True) :
        if data["label"] != 'B53' and (graphe.nodes[u]["type"] != 1 or graphe.nodes[v]["type"] != 1):
            compteur_arete += 1
    compteur_arete = compteur_arete/2 - 4
    

#     print(somme_sommets)
#     print(compteur_arete)
    
    return somme_sommets + compteur_arete

'''Calcul metrique sommets_aretes_ponderees V2 pour graphe commun
meme principe que pour les extensions '''
def calcul_sommets_graphe_commun(graphe1, graphe2, graphe_commun, cle):
    somme_sommets = 0
    for noeud in graphe_commun.nodes() :
        
        if graphe1.nodes[noeud[0]]["type"] == 1 and graphe2.nodes[noeud[1]]["type"] == 1 :
            if len(graphe1[noeud[0]]) <= 1 and len(graphe2[noeud[1]]) <= 1  : ## cas ou les deux sont des helices isolees du reste du graphe
                somme_sommets += 3*min(graphe1.nodes[noeud[0]]["poids"], graphe2.nodes[noeud[1]]["poids"])
#                 if cle[0] == "fichier_5J7L_DA_50_21" and cle[1] == "fichier_5J7L_DA_48_20" :
#                     print(noeud)
            else :
                if len(graphe1[noeud[0]]) == 2 and len(graphe2[noeud[1]]) == 2 : 
                    voisin = None
                    for u,v, data in graphe_commun.edges(data=True) :
                        if data["type"] == "CAN" :
                            if u == noeud :
                                voisin = v
                            if v == noeud :
                                voisin = u        
                    if voisin != None : ## cas ou les deux sont des helices au sein du graphe, entre les memes chaines
                        somme_sommets += 3*min(graphe1.nodes[noeud[0]]["poids"], graphe2.nodes[noeud[1]]["poids"])/2
                    else : ## cas ou les deux sont des helices au sein du graphe, pas entre les memes chaines
                        somme_sommets += min(graphe1.nodes[noeud[0]]["poids"], graphe2.nodes[noeud[1]]["poids"])
#                         if cle[0] == "fichier_4V9F_0_30_4" and cle[1] == "fichier_4V9F_0_48_21" :
                        #print("probleme")
#                             print(noeud)
#                             print(len(graphe1[noeud[0]]))
#                             print(len(graphe2[noeud[1]]))
#                         with open("fichier_pas_bon.txt", 'a') as f :
#                             f.write(str(cle) + "\n")
                        
                else : ## cas ou l une est une helice isolee l autre une helice au sein du graphe
                    somme_sommets += min(graphe1.nodes[noeud[0]]["poids"], graphe2.nodes[noeud[1]]["poids"])
            
        elif graphe1.nodes[noeud[0]]["type"] == 0 and graphe2.nodes[noeud[1]]["type"] == 0 :
            somme_sommets += min(graphe1.nodes[noeud[0]]["poids"], graphe2.nodes[noeud[1]]["poids"])
            
        else :
            somme_sommets += 1
    somme_sommets -= 5
    
    
    compteur_arete = 0    
    for noeud1,noeud2,typ in graphe_commun.edges(data="type") :
        if typ != 'COV' and (graphe1.nodes[noeud1[0]]["type"] != 1 and graphe1.nodes[noeud2[0]]["type"] != 1 and graphe2.nodes[noeud1[1]]["type"] != 1 and graphe2.nodes[noeud2[1]]["type"] != 1):
            compteur_arete += 1
#             if cle[0] == "fichier_5J7L_DA_50_21" and cle[1] == "fichier_5J7L_DA_48_20" :
#                 print(noeud1)
#                 print(noeud2)
    compteur_arete -= 4    
    
#     print(somme_sommets)
#     print(compteur_arete)
    
    return somme_sommets + compteur_arete

'''Calcul sim pour metrique sommets_aretes_ponderees V2'''
def calcul_sim_avec_poids(graphe1, graphe2, graphe_commun, cle):
    poids_sommets_aretes_1 = calcul_sommets_aretes_graphe(graphe1, cle)
    poids_sommets_aretes_2 = calcul_sommets_aretes_graphe(graphe2, cle)
    poids_sommets_aretes_commun = calcul_sommets_graphe_commun(graphe1, graphe2, graphe_commun, cle)

    print(cle)
    print(poids_sommets_aretes_1)
    print(poids_sommets_aretes_2)
    print(poids_sommets_aretes_commun)
    
    sim = (poids_sommets_aretes_commun)/max(poids_sommets_aretes_1, poids_sommets_aretes_2)  
    return sim

'''Calcul metrique aretes ponderees uniquement (toutes aretes)'''
def calcul_aretes_avec_coeff(graphe, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0

    for u,v,data in graphe.edges(data=True) :
        #if (u,v) not in [(1,2), (1,5), (2,1), (3,4), (4,3), (5,1), (2,5), (5,2)] :
            if data["label"] != 'B53' :
                if data["label"] == '0' :
                    if coeffa == 1 :
                        somme_aretes += graphe.nodes[u]["poids"]
                elif (graphe.nodes[u]["type"] == 1 and graphe.nodes[v]["type"] == 1) or (data["label"] == 'CWW' and data["long_range"] == False) :
                    if coeffc == 1 :
                        somme_aretes += graphe.nodes[u]["poids"]
                else :
                    if coeffn == 1 :
                        somme_aretes += graphe.nodes[u]["poids"]
    somme_aretes = somme_aretes/2 - 4
    
    return somme_aretes

'''Calcul metrique aretes ponderees uniquement (toutes aretes)'''
def calcul_aretes_communes_avec_coeff(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
    liste_aretes = []
    for u,v,data in graphe_commun.edges(data=True) :
        #if (u,v) not in [((1, 1), (2, 2)), ((1, 1), (5, 5)), ((3, 3), (4, 4)), ((2,2), (5,5))] :
            if data["type"] != 'B53' and (u,v) not in liste_aretes:
                print(u[0])
                print(u[1])
                print(graphe1.nodes.data())
                print(graphe2.nodes.data())
                print(cle)
                if data["type"] == '0' :
                    if coeffa == 1 :
                        somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
                elif (graphe1.nodes[u[0]]["type"] == 1 and graphe1.nodes[v[0]]["type"] == 1) or (data["type"] == 'CWW' and data["long_range"] == False):
                    if coeffc == 1 :
                        somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
                else :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
                liste_aretes.append((u,v))
    
    return somme_aretes/2 - 4
       
'''Calcul sim avec metrique aretes ponderees uniquement (toutes aretes)'''
def calcul_sim_aretes_avec_coeff(graphe1, graphe2, graphe_commun, cle, coeffc, coeffa, coeffn):
    aretes_1 = calcul_aretes_avec_coeff(graphe1, cle, coeffc, coeffa, coeffn)
    aretes_2 = calcul_aretes_avec_coeff(graphe2, cle, coeffc, coeffa, coeffn)
    aretes_commun = calcul_aretes_communes_avec_coeff(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn)
    
    #if cle == ('fichier_1FJG_A_58_23', 'fichier_5J5B_BA_58_3') or cle == ('fichier_5J5B_BA_58_3', 'fichier_1FJG_A_58_23') :
    print(aretes_1)
    print(aretes_2)
    print(aretes_commun)
    
    print(graphe_commun.edges.data())
    
    return aretes_commun/max(aretes_1, aretes_2)

'''Calcul metrique aretes ponderees uniquement (toutes aretes) en comptant le motif'''
def calcul_aretes_avec_coeff_avec_motif(graphe, cle, coeffc, coeffa, coeffn):
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
    somme_aretes = somme_aretes/2
    
    return somme_aretes

'''Calcul metrique aretes ponderees uniquement (toutes aretes) en comptant le motif'''
def calcul_aretes_communes_avec_coeff_avec_motif(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
    for u,v,data in graphe_commun.edges(data=True) :
        if data["type"] != 'COV' :
            if data["type"] == '0' :
                if coeffa == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
            elif graphe1.nodes[u[0]]["type"] == 1 and graphe1.nodes[v[0]]["type"] == 1 :
                if coeffc == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
            else :
                if coeffn == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
                    
    somme_aretes = somme_aretes
    
    return somme_aretes
       
'''Calcul sim avec metrique aretes ponderees uniquement (toutes aretes) en comptant le motif'''
def calcul_sim_aretes_avec_coeff_avec_motif(graphe1, graphe2, graphe_commun, cle, coeffc, coeffa, coeffn):
    aretes_1 = calcul_aretes_avec_coeff_avec_motif(graphe1, cle, coeffc, coeffa, coeffn)
    aretes_2 = calcul_aretes_avec_coeff_avec_motif(graphe2, cle, coeffc, coeffa, coeffn)
    aretes_commun = calcul_aretes_communes_avec_coeff_avec_motif(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn)
    
    return aretes_commun/max(aretes_1, aretes_2)


# def calcul_aretes_avec_coeff_par_chaine(graphe, chaines, cle, coeffc, coeffa, coeffn, chaines_a_regarder):
#     somme_aretes = 0
#     
#     print(chaines)
#     
#     for u,v,data in graphe.edges(data=True) :
#         dans_chaines_a_regarder = False
#         for elt in chaines_a_regarder :
#             if u in chaines[elt] or v in chaines[elt] :
#                 dans_chaines_a_regarder = True
#         
#         if (u,v) not in [(1,2), (1,5), (2,1), (3,4), (4,3), (5,1), (2,5), (5,2)] and dans_chaines_a_regarder :
#             print((u,v))
#             if data["label"] != 'B53' :
#                 if data["label"] == '0' :
#                     if coeffa == 1 :
#                         somme_aretes += graphe.nodes[u]["poids"]
#                 elif (graphe.nodes[u]["type"] == 1 and graphe.nodes[v]["type"] == 1) or (data["label"] == 'CWW' and data["long_range"] == False) :
#                     if coeffc == 1 :
#                         somme_aretes += graphe.nodes[u]["poids"]
#                 else :
#                     if coeffn == 1 :
#                         somme_aretes += graphe.nodes[u]["poids"]
#     somme_aretes = somme_aretes/2
#     
#     return somme_aretes

'''Calcul metrique aretes ponderees uniquement (toutes aretes) par chaine'''
def calcul_aretes_avec_coeff_par_chaine(graphe, chaines, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
    
    #print(chaines)
    
    for u,v,data in graphe.edges(data=True) :
        a_prendre_u = False
        for num_chaine in graphe.nodes[u]["chaine"] :
            if num_chaine in chaines :
                a_prendre_u = True
        
        a_prendre_v = False
        for num_chaine in graphe.nodes[v]["chaine"] :
            if num_chaine in chaines :
                    a_prendre_v = True
                    
        if (a_prendre_u or a_prendre_v) and not (u in [1,2,3,4,5] and v in [1,2,3,4,5]) :
            #print((u,v))
            if data["label"] != 'B53' :
                if data["label"] == '0' :
                    if coeffa == 1 :
                        somme_aretes += graphe.nodes[u]["poids"]
                elif (graphe.nodes[u]["type"] == 1 and graphe.nodes[v]["type"] == 1) or (data["label"] == 'CWW' and data["long_range"] == False) :
                    if coeffc == 1 :
                        somme_aretes += graphe.nodes[u]["poids"]
                else :
                    if coeffn == 1 :
                        somme_aretes += graphe.nodes[u]["poids"]
    somme_aretes = somme_aretes/2
    
    return somme_aretes

# def calcul_aretes_communes_avec_coeff_par_chaine(graphe_commun, chaines_1, chaines_2, chaines_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn, chaines_a_regarder):
#     somme_aretes = 0
#     
#     for u,v,data in graphe_commun.edges(data=True) :
#         
#         dans_chaines_a_regarder = False
#         for elt in chaines_a_regarder :
#             if u in chaines_commun[elt] or v in chaines_commun[elt] :
#                 dans_chaines_a_regarder = True
#             
#         if (u,v) not in [((1, 1), (2, 2)), ((1, 1), (5, 5)), ((3, 3), (4, 4)), ((2,2), (5,5))] and dans_chaines_a_regarder:
#             if data["type"] != 'COV' :
#                 if data["type"] == 'ART' :
#                     if coeffa == 1 :
#                         somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
#                 elif (graphe1.nodes[u[0]]["type"] == 1 and graphe1.nodes[v[0]]["type"] == 1) or (data["type"] == 'CAN' and data["long_range"] == False):
#                     if coeffc == 1 :
#                         somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
#                 else :
#                     if coeffn == 1 :
#                         somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
#     
#     return somme_aretes

'''Calcul metrique aretes ponderees uniquement (toutes aretes) par chaine'''
def calcul_aretes_communes_avec_coeff_par_chaine(graphe_commun, graphe1, graphe2, chaines, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
     
    for u,v,data in graphe_commun.edges(data=True) :
        a_prendre_u = False
        for num_chaine in graphe1.nodes[u[0]]["chaine"] :
            if num_chaine in chaines :
                a_prendre_u = True
        
        a_prendre_v = False
        for num_chaine in graphe1.nodes[v[0]]["chaine"] :
            if num_chaine in chaines :
                    a_prendre_v = True
                    
        if (a_prendre_u or a_prendre_v) and not (u in [(1,1), (2,2), (3,3), (4,4), (5,5)] and v in [(1,1), (2,2), (3,3), (4,4), (5,5)])  :
            if data["type"] != 'COV' :
                if data["type"] == 'ART' :
                    if coeffa == 1 :
                        somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
                elif (graphe1.nodes[u[0]]["type"] == 1 and graphe1.nodes[v[0]]["type"] == 1) or (data["type"] == 'CAN' and data["long_range"] == False):
                    if coeffc == 1 :
                        somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
                else :
                    if coeffn == 1 :
                        somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
     
    return somme_aretes

# def calcul_sim_aretes_avec_coeff_par_chaine(graphe1, graphe2, graphe_commun, chaines_1, chaines_2, chaines_commun, cle, coeffc, coeffa, coeffn, i):
#     aretes_1 = calcul_aretes_avec_coeff_par_chaine(graphe1, chaines_1, cle, coeffc, coeffa, coeffn, i)
#     aretes_2 = calcul_aretes_avec_coeff_par_chaine(graphe2, chaines_2, cle, coeffc, coeffa, coeffn, i)
#     aretes_commun = calcul_aretes_communes_avec_coeff_par_chaine(graphe_commun, chaines_1, chaines_2, chaines_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn, i)
#     
#     print(aretes_1)
#     print(aretes_2)
#     print(aretes_commun)
#     
#     return aretes_commun/max(aretes_1, aretes_2)
    

'''Calcul sim avec metrique aretes ponderees uniquement (toutes aretes) par chaine'''
def calcul_sim_aretes_avec_coeff_par_chaine(graphe1, graphe2, graphe_commun, chaines, cle, coeffc, coeffa, coeffn):
    aretes_1 = calcul_aretes_avec_coeff_par_chaine(graphe1, chaines, cle, coeffc, coeffa, coeffn)
    aretes_2 = calcul_aretes_avec_coeff_par_chaine(graphe2, chaines, cle, coeffc, coeffa, coeffn)
    aretes_commun = calcul_aretes_communes_avec_coeff_par_chaine(graphe_commun, graphe1, graphe2, chaines, cle, coeffc, coeffa, coeffn)
    
    if cle == ('fichier_5DM6_X_134_2', 'fichier_4V88_A5_48_3') :
        print(aretes_1)
        print(aretes_2)
        print(aretes_commun)
    
    return aretes_commun/max(aretes_1, aretes_2)
    
    
  
'''Calcul sim pour metrique sommets_aretes_ponderees V2 par chaine'''
def calcul_sim_avec_poids_par_chaine(graphe1, graphe2, graphe_commun, chaines_1, chaines_2, chaines_commun, i):
    
    poids_sommets_1 = 0
    for u, poids in graphe1.nodes(data="poids") :
        if u in chaines_1[i] :
            poids_sommets_1 += poids
    poids_sommets_1 -= 1
    
    poids_aretes_1 = 0
    for u,v,typ in graphe1.edges(data="label") :
        if typ != 'B53' and (u in chaines_1[i] and v in chaines_1[i]):
            poids_aretes_1 += 1
    poids_aretes_1 = poids_aretes_1/2
#     if i == 0 or i == 1 :
#         poids_aretes_1 = poids_aretes_1/2 - 2
#     else :
#         poids_aretes_1 = poids_aretes_1/2 - 1
    
    poids_sommets_2 = 0
    for u, poids in graphe2.nodes(data="poids") :
        if u in chaines_2[i] :
            poids_sommets_2 += poids
    poids_sommets_2 -= 1
    
    poids_aretes_2 = 0
    for u,v,typ in graphe2.edges(data="label") :
        if typ != 'B53' and (u in chaines_2[i] and v in chaines_2[i]) :
            poids_aretes_2 += 1   
    poids_aretes_2 = poids_aretes_2/2   
#     if i == 0 or i == 1 :
#         poids_aretes_2 = poids_aretes_2/2 - 2
#     else :
#         poids_aretes_2 = poids_aretes_2/2 - 1
            
    poids_sommets_communs = 0
    for noeud in graphe_commun.nodes() :
        if noeud in chaines_commun[i] :
            poids_sommets_communs += min(graphe1.nodes[noeud[0]]["poids"], graphe2.nodes[noeud[1]]["poids"])
    poids_sommets_communs -= 1
    
    poids_aretes_communes = 0    
    for u,v,typ in graphe_commun.edges(data="type") :
        if typ != 'COV' and (u in chaines_commun[i] and v in chaines_commun[i]):
            poids_aretes_communes += 1
#     if i == 0 or i == 1 :
#         poids_aretes_communes -= 2
#     else :
#         poids_aretes_communes -= 1
    
#     print(poids_sommets_1)
#     print(poids_sommets_2)
#     print(poids_aretes_1)
#     print(poids_aretes_2)
#     print(poids_sommets_communs)
#     print(poids_aretes_communes)
    sim = (poids_sommets_communs + poids_aretes_communes)/max(poids_sommets_1 + poids_aretes_1, poids_sommets_2, poids_aretes_2)  
    return sim, poids_sommets_1, poids_aretes_1, poids_sommets_2, poids_aretes_2, poids_sommets_communs, poids_aretes_communes

'''Calcul sim pour metrique nb_aretes / k'''
def calcul_sim_nb_aretes_par_k(graphe_commun, graphe1, graphe2, cle, taille_ext, coeffc, coeffa, coeffn):
    somme_aretes = calcul_aretes_communes_avec_coeff(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn)
    return somme_aretes/taille_ext

liste = ['5J7L_DA_191_3', '5J7L_DA_191_4', '5FDU_1A_301_1', '5J7L_DA_301_2', '5DM6_X_334_1', '5FDU_1A_334_2', '4V9F_0_335_1', '5J7L_DA_335_2', '3JCS_1_137_4', '4V88_A5_290_1', '4V88_A6_314_2', '5J7L_DA_218_3', '4V9F_0_251_2', '1FJG_A_62_8', '5J7L_DA_137_1', '4V9F_0_118_1', '4V9F_0_62_2', '5J7L_DA_271_2', '4V9F_0_224_1', '5DM6_X_197_1', '3GX5_A_138_6', '1FJG_A_317_2', '5J5B_BA_317_1', '1FJG_A_326_1', '5DM6_X_137_3', '5J5B_BA_314_1', '4V9F_0_134_6', '4V9F_0_328_1', '4V9F_0_197_2', '4V9F_0_62_16', '5J7L_DA_282_2', '4V88_A5_137_2', '5FDU_1A_224_3', '5J7L_DA_326_2']

'''Generation d'un fichier csv pour obtenir la distribution des valeurs de sim '''
def generation_fichier_csv_sim(version, type_sim, fichier, type_fic) :
    with open(fichier, 'rb') as fichier_graphe :
            mon_depickler = pickle.Unpickler(fichier_graphe)
            dico_graphe = mon_depickler.load()
            print(dico_graphe.keys())
            
            tot_sim = []
#             with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/sim_extensions_toutes_aretes_coeffn1_a1_c1_taille_4_vrai.pickle", 'rb') as fichier_sim:
#                 mon_depickler_sim = pickle.Unpickler(fichier_sim)
#                 dico_sim = mon_depickler_sim.load()
#                 for cle in dico_sim.keys() :
#                     tot_sim.append(dico_sim[cle])
#                 for fic in os.listdir("graphes_extension/fichiers_couples_epures") :
#                     if "pickle" in fic :
#                         elt1 = fic.split('_')[2] + '_' + fic.split('_')[3] + '_' + fic.split('_')[4] + '_' + fic.split('_')[5] + '_' + fic.split('_')[6]
#                         elt2 = fic.split('_')[7] + '_' + fic.split('_')[8] + '_' + fic.split('_')[9] + '_' + fic.split('_')[10] + '_' + fic.split('_')[11][:len(fic.split('_')[11])-7]
#                         
#                         if (elt1, elt2) not in dico_graphe.keys() :
#                             print(elt1)
#                             print(elt2)  
                
                
                
                
#                 print(len(dico_sim))
                

            for cle in dico_graphe.keys() :
                    element1 = cle[0]
                    element2 = cle[1]
                    print(element1)
                    print(element2)
                    if type_fic == 'structure' :
                        element1_1 = str(element1).split(",")[0][2:len(str(element1).split(",")[0])-1] + "_" + str(element1).split(",")[1][2:len(str(element1).split(",")[1])-1] + "_" + str(element1).split(",")[2][1:] + "_" + str(element1).split(",")[3][1:len(str(element1).split(",")[3])-1] 
                        element2_1 = str(element2).split(",")[0][2:len(str(element2).split(",")[0])-1] + "_" + str(element2).split(",")[1][2:len(str(element2).split(",")[1])-1] + "_" + str(element2).split(",")[2][1:] + "_" + str(element2).split(",")[3][1:len(str(element2).split(",")[3])-1]                 
                        element1 = "fichier_" + element1_1
                        element2 = "fichier_" + element2_1
                    print(element1)
                    print(element2)
                    enlever = False 
                    for elt in liste : 
                        if elt in element1 or elt in element2 :
                            enlever = True
                    #print(enlever)
                     
                    if enlever == False :
  
                        if "crible_taille_extensions_%s_%s.pickle"%(element1[8:], element2[8:]) in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles") :
                            with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/crible_taille_extensions_%s_%s.pickle"%(element1[8:], element2[8:]), 'rb') as fichier_crible:
                                mon_depickler_1 = pickle.Unpickler(fichier_crible)
                                crible = mon_depickler_1.load()
                                for cle in crible.keys() :
                                    tot_sim.append(crible[cle][4])
                             
                        with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_5/"+element1+".pickle", 'rb') as fichier_graphe_1 :
                            mon_depickler_1 = pickle.Unpickler(fichier_graphe_1)
                            graphe1 = mon_depickler_1.load()
                              
                            with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_5/"+element2+".pickle", 'rb') as fichier_graphe_2 :
                                mon_depickler_2 = pickle.Unpickler(fichier_graphe_2)
                                graphe2 = mon_depickler_2.load()
                                if version == 'sans_motif' :
                                    if type_sim == 'raymond' :
                                        sim = calcul_sim_raymond_sans_motif(graphe1, graphe2, dico_graphe[cle])
                                        tot_sim.append(sim)
                                    if type_sim == 'raymond_non_cov' :
                                        sim = calcul_sim_raymond_non_cov_sans_motif(graphe1, graphe2, dico_graphe[cle])
                                        tot_sim.append(sim)
                                    if type_sim == 'non_cov' :
                                        sim = calcul_sim_non_cov_sans_motif(graphe1, graphe2, dico_graphe[cle])
                                        print(sim)
                                        tot_sim.append(sim[0])
                                    if type_sim == 'sommets_aretes' :
                                        sim = calcul_sim_avec_poids(graphe1, graphe2, dico_graphe[cle], cle)
                                        tot_sim.append(sim)
                                    if type_sim == 'toutes_aretes' :
                                        sim = calcul_sim_aretes_avec_coeff(graphe1, graphe2, dico_graphe[cle], cle, 0,0,1)
                                        tot_sim.append(sim)
                                else :
                                    if type_sim == 'raymond' :
                                        sim = calcul_sim_raymond_avec_motif(graphe1, graphe2, dico_graphe[cle])
                                        tot_sim.append(sim)
                                    if type_sim == 'raymond_non_cov' :
                                        sim = calcul_sim_raymond_non_cov_avec_motif(graphe1, graphe2, dico_graphe[cle])
                                        tot_sim.append(sim)
                                    if type_sim == 'non_cov' :
                                        sim = calcul_sim_non_cov_avec_motif(graphe1, graphe2, dico_graphe[cle])
                                        tot_sim.append(sim)
                               #                             if sim < float(dico_sim[cle]) :
    #                                 print("plus petit")
    #                                 print(dico_sim[cle])
    #                                 print(sim)
    #                             else :
    #                                 print("plus grand")
    #                                 print(dico_sim[cle])
    #                                 print(sim)
                                
    #                 else :
    #                     print(element1)
    #                     print(element2)
                            
            #print(len(tot_sim_complet))   
            print(len(tot_sim))
            with open(fichier[:len(fichier)-7]+"_sim_csv_"+type_sim+"_"+version+"_taille5.csv", 'w', newline='') as fichier_csv :
                    csvwriter = csv.writer(fichier_csv, delimiter=',')  
                    i = 0
                    while i+1024 < len(tot_sim) :
                        csvwriter.writerow(tot_sim[i:i+1024])
                        i = i+1024
                    csvwriter.writerow(tot_sim[i:])
#                 with open(fichier[:len(fichier)-7]+"_sim_csv_tot_new_"+version+".csv", 'w', newline='') as fichier_csv_2 :
#                     csvwriter_2 = csv.writer(fichier_csv_2, delimiter=',')    
#                     i = 0
#                     while i+1024 < len(tot_sim) :
#                         csvwriter_2.writerow(tot_sim[i:i+1024])
#                         i = i+1024
#                     csvwriter_2.writerow(tot_sim[i:])
#                 
#                 with open(fichier[:len(fichier)-7]+"_sim_csv_tot_non_raymond_"+version+".csv", 'w', newline='') as fichier_csv_3 :
#                     csvwriter_3 = csv.writer(fichier_csv_3, delimiter=',')
#                     i = 0
#                     while i+1024 < len(tot_sim_new) :
#                         csvwriter_3.writerow(tot_sim_new[i:i+1024])
#                         i = i+1024
#                     csvwriter_3.writerow(tot_sim_new[i:])
                
                
                
#             sns.distplot(tot_sim_new, bins=20)
#             plt.show()

''' Generation d'un fichier pickle comprenant un dictionnaire où une cle
 de comparaison est associee à une valeur de sim'''
def generation_fichier_pickle(version, type):
    with open("dico_graphe_epure_en_tout.pickle", 'rb') as fichier_graphe :
        mon_depickler = pickle.Unpickler(fichier_graphe)
        dico_graphe = mon_depickler.load()
        
        tot_sim = {}
        for cle in dico_graphe.keys() :
            element1 = cle[0]
            element2 = cle[1]
            #print(element1)
            #print(element2)
            
            enlever = False 
            for elt in liste : 
                if elt in element1 or elt in element2 :
                    enlever = True
            #print(enlever)
            
            if enlever == False :
                        
                with open("graphes_extension/"+element1+".pickle", 'rb') as fichier_graphe_1 :
                    mon_depickler_1 = pickle.Unpickler(fichier_graphe_1)
                    graphe1 = mon_depickler_1.load()
                    
                    with open("graphes_extension/"+element2+".pickle", 'rb') as fichier_graphe_2 :
                        mon_depickler_2 = pickle.Unpickler(fichier_graphe_2)
                        graphe2 = mon_depickler_2.load()
                        if version == 'sans_motif' :
                            if type == 'raymond' : 
                                sim = calcul_sim_raymond_sans_motif(graphe1, graphe2, dico_graphe[cle])
                            elif type == 'raymond_non_cov' :
                                sim = calcul_sim_raymond_non_cov_sans_motif(graphe1, graphe2, dico_graphe[cle])
                            else : 
                                sim,compteur_arete_commun = calcul_sim_non_cov_sans_motif(graphe1, graphe2, dico_graphe[cle])
                        else :
                            if type == 'raymond' : 
                                sim = calcul_sim_raymond_avec_motif(graphe1, graphe2, dico_graphe[cle])
                            elif type == 'raymond_non_cov' :
                                sim = calcul_sim_raymond_non_cov_avec_motif(graphe1, graphe2, dico_graphe[cle])
                            else :
                                sim,compteur_arete_commun = calcul_sim_non_cov_avec_motif(graphe1, graphe2, dico_graphe[cle])
                        
                        tot_sim.update({(element1, element2) : sim})
                        
        with open("dico_sim_"+version+"_"+type+"en_tout.pickle", 'wb') as fichier_pickle :
            mon_pickler = pickle.Pickler(fichier_pickle)
            mon_pickler.dump(tot_sim)
            
            
def tri_par_insertion(deb, tab):
    pos_mini = deb
    mini = tab[deb][1]
    for i in range(deb+1, len(tab)) :
        if tab[i][1] < mini :
            mini = tab[i][1]
            pos_mini = i
    
    return pos_mini

def echange(tab, pos_mini, deb):
    temp = tab[deb]
    tab[deb] = tab[pos_mini]
    tab[pos_mini] = temp



'''Obtention d'un fichier csv comportant toutes les valeurs de sim entre un element en entree et tous les autres
pour distribution'''
def calcul_sim_par_element(element, version):
    
    with open("dico_graphe.pickle", 'rb') as fichier_graphe :
            mon_depickler = pickle.Unpickler(fichier_graphe)
            dico_graphe = mon_depickler.load()
            
            tot_sim_new = []
            tot_sim = []
            tab_cle = []
            for cle in dico_graphe.keys() :
                if cle[0] == element or cle[1] == element :
                    enlever = False 
                    for elt in liste : 
                        if elt in cle[0] or elt in cle[1] :
                            enlever = True
                    if enlever == False :
                        with open("graphes_extension/"+cle[0]+".pickle", 'rb') as fichier_graphe_1 :
                            mon_depickler_1 = pickle.Unpickler(fichier_graphe_1)
                            graphe1 = mon_depickler_1.load()
                            
                            with open("graphes_extension/"+cle[1]+".pickle", 'rb') as fichier_graphe_2 :
                                mon_depickler_2 = pickle.Unpickler(fichier_graphe_2)
                                graphe2 = mon_depickler_2.load()
                                
                                if version == "avec_motif" :
                                    sim, compteur_arete_commun = calcul_sim_non_cov_avec_motif(graphe1, graphe2, dico_graphe[cle])
                                else :
                                    sim, compteur_arete_commun = calcul_sim_non_cov_sans_motif(graphe1, graphe2, dico_graphe[cle])
                                
                                if cle[0] == element and cle[1] not in tab_cle :
                                    tot_sim_new.append((cle[1], sim, compteur_arete_commun))
                                    tab_cle.append(cle[1])
                                    tot_sim.append(sim)
                                elif cle[1] == element and cle[0] not in tab_cle :
                                    tot_sim_new.append((cle[0], sim, compteur_arete_commun))
                                    tab_cle.append(cle[0])
                                    tot_sim.append(sim)
            print(len(tot_sim_new))
            print(tot_sim_new)
            
            
            pos_mini = 0
            for i in range(0, len(tot_sim_new)) :
                pos_mini = tri_par_insertion(i, tot_sim_new)
                print(pos_mini)
                echange(tot_sim_new, pos_mini, i)
            
            print(tot_sim_new)
            
            with open("ordre_"+element+"_"+version+".txt", 'w') as fichier :
                for elt in tot_sim_new :
                    fichier.write(elt[0] + " " + str(elt[1]) + " "+ str(elt[2]) +"\n")
            
            with open("sim_"+element+"_"+version+".csv", 'w', newline="") as fichier_csv :
                csvwriter = csv.writer(fichier_csv, delimiter=',') 
                csvwriter.writerow(sorted(tot_sim))
                
    

'''Recherche de la similarite maximum '''
def calcul_max_sim() :            
    with open("dico_graphe.pickle", 'rb') as fichier_graphe :
            mon_depickler = pickle.Unpickler(fichier_graphe)
            dico_graphe = mon_depickler.load()
            
            with open("dico_sim.pickle", 'rb') as fichier_sim :
                mon_depickler_sim = pickle.Unpickler(fichier_sim)
                dico_sim = mon_depickler_sim.load()
 
                #print(len(dico_sim))
                
                for cle in dico_graphe.keys() :
                    element1 = cle[0]
                    element2 = cle[1]
#                     print(element1)
#                     print(element2)
                    
                    enlever = False 
                    for elt in liste : 
                        if elt in element1 or elt in element2 :
                            enlever = True
                    #print(enlever)
                    
                    if enlever == False :        
                        with open("graphes_extension/"+element1+".pickle", 'rb') as fichier_graphe_1 :
                            mon_depickler_1 = pickle.Unpickler(fichier_graphe_1)
                            graphe1 = mon_depickler_1.load()
                            
                            with open("graphes_extension/"+element2+".pickle", 'rb') as fichier_graphe_2 :
                                mon_depickler_2 = pickle.Unpickler(fichier_graphe_2)
                                graphe2 = mon_depickler_2.load()
                                
                                sim, compteur_arete_commun = calcul_sim_non_cov_sans_motif(graphe1, graphe2, dico_graphe[cle])
                                
                            if float(sim) >= 0.8:
                                print(element1)
                                print(element2)
                                print(sim)

'''Creation d'un dictionnaire avec la similarite de chaque paire selon les criteres en parametres
(comme generation fichier pickle mais plus recent pour la version extension)'''
def stockage_sim(typ_graphe, typ_sim, depart, taille_ext):
    if depart == "extensions" :
        with open(EXTENSION_PATH%taille_ext+"dico_comp_complet_metrique_%s_taille_%s.pickle"%(typ_graphe,taille_ext), 'rb') as fichier_graphe :
            mon_depickler = pickle.Unpickler(fichier_graphe)
            dico_graphe = mon_depickler.load()

#         with open("result_graphes_comp/graphe_comp_couples_possibles_fichier_1FJG_A_48_8_fichier_1FJG_A_138_3.pickle", 'rb') as fichier_graphe :
#             mon_depickler = pickle.Unpickler(fichier_graphe)
#             dico_graphe = mon_depickler.load()
             
            dico_sim = {}
            for cle in dico_graphe.keys() : 
                print(cle)
                with open(EXTENSION_PATH_TAILLE%taille_ext+cle[0]+".pickle", 'rb') as fichier_graphe_1 :
                    mon_depickler_1 = pickle.Unpickler(fichier_graphe_1)
                    graphe1 = mon_depickler_1.load()
                
                    with open(EXTENSION_PATH_TAILLE%taille_ext+cle[1]+".pickle", 'rb') as fichier_graphe_2 :
                        mon_depickler_2 = pickle.Unpickler(fichier_graphe_2)
                        graphe2 = mon_depickler_2.load()
                
                        if typ_sim == "raymond" :
                            sim = calcul_sim_raymond_sans_motif(graphe1, graphe2, dico_graphe[cle])
                            
                        elif typ_sim == "longue_distance" :
                            sim = calcul_sim_non_cov_sans_motif(graphe1, graphe2, dico_graphe[cle])
                        
                        elif typ_sim == "non_cov"  :
                            sim = calcul_sim_avec_poids(graphe1, graphe2, dico_graphe[cle], cle)   
#                             if cle[0] == "fichier_4V9F_0_30_4" and cle[1] == "fichier_4V9F_0_48_21" :
#                             print(cle[0])
#                             print(cle[1])
#                             print(sim)
                            print(dico_graphe[cle].edges.data())
                        elif typ_sim == "toutes_aretes_coeff_all1" :
#                             print(cle)
                            sim = calcul_sim_aretes_avec_coeff(graphe1, graphe2, dico_graphe[cle], cle, 1, 1, 1)
#                             print(sim)
                        elif typ_sim == "toutes_aretes_coeff_all1_chaines_1_3" :
                            sim = calcul_sim_aretes_avec_coeff_par_chaine(graphe1, graphe2, dico_graphe[cle], [1,3], cle, 1,1,1)
                        elif typ_sim == "toutes_aretes_coeff_all1_chaines_2_4" :
                            sim = calcul_sim_aretes_avec_coeff_par_chaine(graphe1, graphe2, dico_graphe[cle], [2,4], cle, 1,1,1)   
                        elif typ_sim == "toutes_aretes_coeff_all1_par_k"  :  
                            sim = calcul_sim_nb_aretes_par_k(dico_graphe[cle], graphe1, graphe2, cle, taille_ext, 1, 1, 1)              
                        dico_sim.update({(cle[0][8:], cle[1][8:]) : sim})

    else :
        with open("fichier_comp_grands_graphes_V2.pickle", 'rb') as fichier_graphe :
            mon_depickler = pickle.Unpickler(fichier_graphe)
            dico_graphe = mon_depickler.load()
             
            dico_sim = {}
            for cle in dico_graphe.keys() : 
                element1 = str(cle[0]).split(",")
                element2 = str(cle[1]).split(",")
                elt1 = element1[0][2:len(element1[0])-1] + "_" + element1[1][2:len(element1[1])-1] + "_" + element1[2][1:] + "_" + element1[3][1:len(element1[3])-1]
                elt2 = element2[0][2:len(element2[0])-1] + "_" + element2[1][2:len(element2[1])-1] + "_" + element2[2][1:] + "_" + element2[3][1:len(element2[3])-1]
#                 print(elt1)
#                 print(elt2)
                with open("grands_graphes.pickle", 'rb') as fichier_grands_graphes :
                    mon_depickler_grands_graphes = pickle.Unpickler(fichier_grands_graphes)
                    dico_grands_graphes = mon_depickler_grands_graphes.load()
                    
                    
                    if typ_sim == "raymond" : 
                        sim = calcul_sim_raymond_sans_motif(dico_grands_graphes[cle[0]], dico_grands_graphes[cle[1]], dico_graphe[cle])
                        
                    elif typ_sim == "non_cov" :
                        sim = calcul_sim_non_cov_sans_motif(dico_grands_graphes[cle[0]], dico_grands_graphes[cle[1]], dico_graphe[cle])
                        if sim[0] > 1.0 :
                            print(elt1)
                            print(elt2)
                            print(sim)
                            
                        if elt1 == "5J7L_DA_50_21" and elt2 == "5J7L_DA_48_20" :
                            print("ramou")
                            print(sim)
               
                        
                        
                    else : ## type=="non_cov_sommets"
                        sim = calcul_sommets_aretes_grands_graphes_commun(dico_graphe[cle])/max(calcul_sommets_aretes_grands_graphes(dico_grands_graphes[cle[0]]), calcul_sommets_aretes_grands_graphes(dico_grands_graphes[cle[1]]))
                    dico_sim.update({(elt1, elt2) : sim})    
    
    #print(dico_sim)       
                 
    with open(EXTENSION_PATH%taille_ext+"sim_"+depart+"_"+typ_sim+"_taille_%s.pickle"%taille_ext, 'wb') as fichier_sim :
        mon_pickler = pickle.Pickler(fichier_sim)
        mon_pickler.dump(dico_sim)
    return dico_sim

'''Observation des valeurs de sim entre homologues et entre non homologues'''
def repartition_homologues():
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
                  ['4YAZ_R_36_25', '3UCZ_R_62_15']]# '5FJC_A_138_1', '4L81_A_25_77']]
    
    dico_homologues = {}
    dico_non_homologues = {}
    with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/sim_extensions_toutes_aretes_coeffn1_a1_c1_taille_4_vrai.pickle", 'rb') as fichier_sim:
        mon_depickler_sim = pickle.Unpickler(fichier_sim)
        dico_sim = mon_depickler_sim.load()
        
        for cle in dico_sim.keys() :
            est_homologue = False
            for elt in homologues :
                if cle[0] in elt and cle[1] in elt :
                    est_homologue = True
 
            if est_homologue :
                dico_homologues.update({cle : dico_sim[cle]})
            else :       
                dico_non_homologues.update({cle : dico_sim[cle]})    
    print(len(dico_homologues))
    #print(dico_homologues)
    print(len(dico_non_homologues))
    #print(dico_non_homologues)
    
    somme_sim = 0
    maxi_sim = -0.1
    mini_sim = 1.1
    for cle in dico_homologues.keys() :
        somme_sim += dico_homologues[cle]
        if dico_homologues[cle] < mini_sim :
            mini_sim = dico_homologues[cle]    
        if dico_homologues[cle] > maxi_sim :
            maxi_sim = dico_homologues[cle]
            print(cle)
            
    print(mini_sim)
    print(maxi_sim)
    print(somme_sim/len(dico_homologues))


'''Generation d'un fichier csv avec directement le nombre de valeurs de sim entre intervalles pour distribution
(en fait on pourrait le faire completement avec python) '''
def distribution_similarite(depart_sim, typ_sim, taille_ext):
    with open(EXTENSION_PATH%taille_ext+"sim_%s_%s_taille_%s.pickle"%(depart_sim, typ_sim, taille_ext), 'rb') as fichier_sim :
        mon_depickler = pickle.Unpickler(fichier_sim)
        dico_sim = mon_depickler.load()
        print(dico_sim)
        
        nombre_intervalles = [0]*20
        tab_intervalles = [0]*20
        
        print(max(dico_sim.values()))
        print(min(dico_sim.values()))
        if max(dico_sim.values()) > 1.0 :
            intervalle = (max(dico_sim.values())-min(dico_sim.values()))/20
        else :
            intervalle = 0.05
        print(intervalle)
        
        if max(dico_sim.values()) > 1.0 :
            for i in range(20) :
                tab_intervalles[i] = i*intervalle+min(dico_sim.values())
        else :
            for i in range(20) :
                tab_intervalles[i] = i*intervalle
        
        for cle in dico_sim.keys() :
            if dico_sim[cle] == max(dico_sim.values()) :
                nombre_intervalles[19] += 1
            else :
                print(dico_sim[cle])
                print(int(dico_sim[cle]/intervalle - min(dico_sim.values())/intervalle))
                if max(dico_sim.values()) > 1.0 :
                    nombre_intervalles[int(dico_sim[cle]/intervalle - min(dico_sim.values())/intervalle)] += 1
                else :
                    nombre_intervalles[int(dico_sim[cle]/intervalle)] += 1
        
        proportion_intervalles = [0]*20
        for i in range(20) :
            proportion_intervalles[i] = nombre_intervalles[i]/4005
        
        print(nombre_intervalles) 
        print(proportion_intervalles)  
        
        
        
        with open(EXTENSION_PATH%taille_ext+"csv_distrib_sim_%s.csv"%typ_sim, 'w', newline='') as fichier_csv :
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(tab_intervalles)
            csvwriter.writerow(proportion_intervalles) 
            
            
''' recherche difference entre les trois metriques (non cov dans graphe global, non cov dans extension, non canoniques + canoniques longue distance dans extension)'''
def diff_metrique_non_cov() :    
    with open("sim_extensions_longue_distance.pickle", 'rb') as fichier_sim_1 :
        mon_depickler_1 = pickle.Unpickler(fichier_sim_1)
        dico_sim_extensions_longue_distance = mon_depickler_1.load()
          
        with open("sim_extensions_non_cov_nouvelle_metrique.pickle", 'rb') as fichier_sim_2 :
            mon_depickler_2 = pickle.Unpickler(fichier_sim_2)
            dico_sim_extensions_non_cov = mon_depickler_2.load()
              
            with open("sim_grands_graphes_non_cov_sommets_nouvelle_metrique.pickle", 'rb') as fichier_sim_3 :
                mon_depickler_3 = pickle.Unpickler(fichier_sim_3)
                dico_sim_graphe_global_non_cov = mon_depickler_3.load()
                 
                print(dico_sim_graphe_global_non_cov)
                  
                ordres = ["el,ec,gc", "el,gc,ec", "ec,el,gc", "ec,gc,el", "gc,el,ec", "gc,ec,el"]
                nombre_ordres = dict((el,0) for el in ordres)
                moyenne_ordres = dict((el,[0,0]) for el in ordres)
                print(nombre_ordres.keys())
                print(moyenne_ordres)
#                 print(dico_sim_extensions_longue_distance.keys())
#                 print(dico_sim_extensions_non_cov.keys())
#                 print(dico_sim_graphe_global_non_cov.keys())
                for cle in dico_sim_extensions_non_cov.keys() :
                    if cle not in dico_sim_extensions_longue_distance.keys() :
                        cle_ex_longue_distance = (cle[1], cle[0])
                    else :
                        cle_ex_longue_distance = cle
                    if cle not in dico_sim_graphe_global_non_cov.keys() :
                        cle_gg_non_cov = (cle[1], cle[0])
                    else :
                        cle_gg_non_cov = cle
                              
                    if dico_sim_extensions_longue_distance[cle_ex_longue_distance][0] <= dico_sim_extensions_non_cov[cle] and dico_sim_extensions_longue_distance[cle_ex_longue_distance][0] <= dico_sim_graphe_global_non_cov[cle_gg_non_cov] :
                        if dico_sim_extensions_non_cov[cle] < dico_sim_graphe_global_non_cov[cle_gg_non_cov] :
                            nombre_ordres["el,ec,gc"] += 1
                            moyenne_ordres["el,ec,gc"][0] += dico_sim_extensions_non_cov[cle] - dico_sim_extensions_longue_distance[cle_ex_longue_distance][0]
                            moyenne_ordres["el,ec,gc"][1] += dico_sim_graphe_global_non_cov[cle_gg_non_cov] - dico_sim_extensions_non_cov[cle]
                              
                        else :
                            nombre_ordres["el,gc,ec"] += 1
                            moyenne_ordres["el,gc,ec"][0] += dico_sim_graphe_global_non_cov[cle_gg_non_cov] - dico_sim_extensions_longue_distance[cle_ex_longue_distance][0]
                            moyenne_ordres["el,gc,ec"][1] += dico_sim_extensions_non_cov[cle] - dico_sim_graphe_global_non_cov[cle_gg_non_cov]
                              
                              
                            print(cle) 
#                             print(cle_ex_longue_distance)
#                             print(dico_sim_extensions_longue_distance[cle_ex_longue_distance][0])
#                             print(dico_sim_extensions_non_cov[cle])
#                             print(dico_sim_graphe_global_non_cov[cle_gg_non_cov])
                              
                    elif dico_sim_extensions_non_cov[cle] <= dico_sim_extensions_longue_distance[cle_ex_longue_distance][0] and dico_sim_extensions_non_cov[cle] <= dico_sim_graphe_global_non_cov[cle_gg_non_cov] :
                        if dico_sim_extensions_longue_distance[cle_ex_longue_distance][0] < dico_sim_graphe_global_non_cov[cle_gg_non_cov] :
                            nombre_ordres["ec,el,gc"] += 1
                            moyenne_ordres["ec,el,gc"][0] += dico_sim_extensions_longue_distance[cle_ex_longue_distance][0] - dico_sim_extensions_non_cov[cle]
                            moyenne_ordres["ec,el,gc"][1] += dico_sim_graphe_global_non_cov[cle_gg_non_cov] - dico_sim_extensions_longue_distance[cle_ex_longue_distance][0]
                              
                              
                        else :
                            nombre_ordres["ec,gc,el"] += 1
                            moyenne_ordres["ec,gc,el"][0] += dico_sim_graphe_global_non_cov[cle_gg_non_cov] - dico_sim_extensions_non_cov[cle]
                            moyenne_ordres["ec,gc,el"][1] += dico_sim_extensions_longue_distance[cle_ex_longue_distance][0] - dico_sim_graphe_global_non_cov[cle_gg_non_cov]
                              
                              
                    elif dico_sim_graphe_global_non_cov[cle_gg_non_cov] <= dico_sim_extensions_longue_distance[cle_ex_longue_distance][0] and dico_sim_graphe_global_non_cov[cle_gg_non_cov] <= dico_sim_extensions_non_cov[cle] :
                        if dico_sim_extensions_longue_distance[cle_ex_longue_distance][0] < dico_sim_extensions_non_cov[cle] :
                            nombre_ordres["gc,el,ec"] += 1
  
                            moyenne_ordres["gc,el,ec"][0] += dico_sim_extensions_longue_distance[cle_ex_longue_distance][0] - dico_sim_graphe_global_non_cov[cle_gg_non_cov]
                            moyenne_ordres["gc,el,ec"][1] += dico_sim_extensions_non_cov[cle] - dico_sim_extensions_longue_distance[cle_ex_longue_distance][0]
                              
                        else :
                            nombre_ordres["gc,ec,el"] += 1
                            moyenne_ordres["gc,ec,el"][0] += dico_sim_extensions_non_cov[cle] - dico_sim_graphe_global_non_cov[cle_gg_non_cov]
                            moyenne_ordres["gc,ec,el"][1] += dico_sim_extensions_longue_distance[cle_ex_longue_distance][0] - dico_sim_extensions_non_cov[cle]
                              
                              
                    else :
                        print("autre cas?")
                        print(dico_sim_graphe_global_non_cov[cle_gg_non_cov])
                        print(dico_sim_extensions_longue_distance[cle_ex_longue_distance][0])
                        print(dico_sim_extensions_non_cov[cle])
                          
                print(nombre_ordres)
                  
                moyenne_ordres["el,ec,gc"][0] = moyenne_ordres["el,ec,gc"][0]/nombre_ordres["el,ec,gc"]
                moyenne_ordres["el,ec,gc"][1] = moyenne_ordres["el,ec,gc"][1]/nombre_ordres["el,ec,gc"]
                  
                moyenne_ordres["ec,el,gc"][0] = moyenne_ordres["ec,el,gc"][0]/nombre_ordres["ec,el,gc"]
                moyenne_ordres["ec,el,gc"][1] = moyenne_ordres["ec,el,gc"][1]/nombre_ordres["ec,el,gc"]
                  
                moyenne_ordres["ec,gc,el"][0] = moyenne_ordres["ec,gc,el"][0]/nombre_ordres["ec,gc,el"]
                moyenne_ordres["ec,gc,el"][1] = moyenne_ordres["ec,gc,el"][1]/nombre_ordres["ec,gc,el"]
                  
                moyenne_ordres["el,gc,ec"][0] = moyenne_ordres["el,gc,ec"][0]/nombre_ordres["el,gc,ec"]
                moyenne_ordres["el,gc,ec"][1] = moyenne_ordres["el,gc,ec"][1]/nombre_ordres["el,gc,ec"]
                  
                moyenne_ordres["gc,el,ec"][0] = moyenne_ordres["gc,el,ec"][0]/nombre_ordres["gc,el,ec"]
                moyenne_ordres["gc,el,ec"][1] = moyenne_ordres["gc,el,ec"][1]/nombre_ordres["gc,el,ec"]
                  
                moyenne_ordres["gc,ec,el"][0] = moyenne_ordres["gc,ec,el"][0]/nombre_ordres["gc,ec,el"]
                moyenne_ordres["gc,ec,el"][1] = moyenne_ordres["gc,ec,el"][1]/nombre_ordres["gc,ec,el"]
                  
                print(moyenne_ordres)
                 
                print(dico_sim_graphe_global_non_cov[('5DM6_X_48_9', '1FJG_A_109_6')])
                
                

''' Comparaison des valeurs de sim  aretes ponderees uniquement version n1a1c1 et version a0c0n1 '''

def diff_valeurs_sim_111_001():
    compteur_all_a0_c0 = 0
    compteur_a0_c0_all = 0
    compteur_egalite = 0
    somme_all_a0_c0 = 0
    somme_a0_c0_all = 0
      
    mini_all_a0_c0 = 1.1
    maxi_all_a0_c0 = -1
      
    mini_a0_c0_all = 1.1
    maxi_a0_c0_all = -1
      
    maxi_val_all_a0_c0 = -1
    mini_val_all_a0_c0 = 1.1
      
    compteur_diff_sup_0_1 = 0
    compteur = 0
    somme = 0
      
    intervalles = [{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}]
      
    with open("Extensions/Metrique_toutes_aretes/sim_extensions_toutes_aretes_coeffn1_a1_c1.pickle", 'rb') as fichier_sim_1 :
        mon_depickler_1 = pickle.Unpickler(fichier_sim_1)
        dico_sim_extensions_toutes_aretes_coeffall1 = mon_depickler_1.load()
            
        with open("Extensions/Metrique_toutes_aretes/sim_extensions_toutes_aretescoeff_n1_a0_c0.pickle", 'rb') as fichier_sim_2 :
            mon_depickler_2 = pickle.Unpickler(fichier_sim_2)
            dico_sim_extensions_toutes_aretes_coeffa0_c0 = mon_depickler_2.load() 
              
            for cle in dico_sim_extensions_toutes_aretes_coeffall1.keys() :
                if cle not in dico_sim_extensions_toutes_aretes_coeffa0_c0.keys() :
                    cle_ex_toutes_aretes_coeffa0_c0 = (cle[1], cle[0])
                else :
                    cle_ex_toutes_aretes_coeffa0_c0 = cle   
                      
                      
                #print(dico_sim_extensions_toutes_aretes_coeffall1[cle])
                #print(dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0])
                  
                  
                if dico_sim_extensions_toutes_aretes_coeffall1[cle] < dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] :
                    compteur_all_a0_c0 += 1
                      
                    somme_all_a0_c0 += dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] - dico_sim_extensions_toutes_aretes_coeffall1[cle]
                    if dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] - dico_sim_extensions_toutes_aretes_coeffall1[cle] < mini_all_a0_c0 :
                        mini_all_a0_c0 = dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] - dico_sim_extensions_toutes_aretes_coeffall1[cle]
                    if dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] - dico_sim_extensions_toutes_aretes_coeffall1[cle] > maxi_all_a0_c0 :
                        maxi_all_a0_c0 = dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] - dico_sim_extensions_toutes_aretes_coeffall1[cle]
                        #print(cle)
                    if dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] > maxi_val_all_a0_c0 : 
                        maxi_val_all_a0_c0 = dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0]
                        #print(cle)
                    if dico_sim_extensions_toutes_aretes_coeffall1[cle] < mini_val_all_a0_c0 :
                        mini_val_all_a0_c0 = dico_sim_extensions_toutes_aretes_coeffall1[cle]
                      
                    if dico_sim_extensions_toutes_aretes_coeffall1[cle] > 0.5 and dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] > 0.5 :
                        print(cle)
                        print(dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] - dico_sim_extensions_toutes_aretes_coeffall1[cle])
                        compteur += 1
                        somme += dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] - dico_sim_extensions_toutes_aretes_coeffall1[cle]
                      
                    num_classe = (dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] - dico_sim_extensions_toutes_aretes_coeffall1[cle]) / 0.05
                    intervalles[int(num_classe)].update({(dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0], dico_sim_extensions_toutes_aretes_coeffall1[cle]): dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] - dico_sim_extensions_toutes_aretes_coeffall1[cle]})
                  
                elif dico_sim_extensions_toutes_aretes_coeffall1[cle] > dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] :
                    compteur_a0_c0_all += 1
#                     print(cle)
#                     print(dico_sim_extensions_toutes_aretes_coeffall1[cle])
#                     print(dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0])
                    somme_a0_c0_all += dico_sim_extensions_toutes_aretes_coeffall1[cle] - dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0]
                    if dico_sim_extensions_toutes_aretes_coeffall1[cle] - dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] < mini_a0_c0_all :
                        mini_a0_c0_all = dico_sim_extensions_toutes_aretes_coeffall1[cle] - dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0]
                    if dico_sim_extensions_toutes_aretes_coeffall1[cle] - dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] > maxi_a0_c0_all :
                        maxi_a0_c0_all = dico_sim_extensions_toutes_aretes_coeffall1[cle] - dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0]
                        #print(cle)
                    if dico_sim_extensions_toutes_aretes_coeffall1[cle] - dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0] > 0.3 :
                        compteur_diff_sup_0_1 += 1
                      
                    num_classe = (dico_sim_extensions_toutes_aretes_coeffall1[cle] - dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0]) / 0.05
                    intervalles[int(num_classe)].update({(dico_sim_extensions_toutes_aretes_coeffall1[cle], dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0]) : dico_sim_extensions_toutes_aretes_coeffall1[cle] - dico_sim_extensions_toutes_aretes_coeffa0_c0[cle_ex_toutes_aretes_coeffa0_c0]})
     
                      
                else :
                    compteur_egalite += 1
                      
    print(compteur_all_a0_c0)
    print(compteur_a0_c0_all)
    print(compteur_egalite)
      
    print(compteur_diff_sup_0_1)
      
    #print(somme_all_a0_c0)
    if compteur_all_a0_c0 > 0 :
        print(somme_all_a0_c0/compteur_all_a0_c0)
    if compteur_a0_c0_all > 0 :
        print(somme_a0_c0_all/compteur_a0_c0_all)
      
    print(mini_all_a0_c0)       
    print(maxi_all_a0_c0)
      
    print(mini_a0_c0_all)
    print(maxi_a0_c0_all)
      
    print(mini_val_all_a0_c0)
    print(maxi_val_all_a0_c0)
      
    print(somme/compteur)
     
    liste_max = {}
    liste_moyenne = {}
    for i in range(20) :
        print(len(intervalles[i]))
        print(intervalles[i])
         
        if len(intervalles[i]) > 0 :
            maxi = -1
            somme = 0
            for elt in intervalles[i].keys() :
                print(elt[1])
                if elt[1] > maxi :
                    maxi = elt[1]
                somme += elt[1]
            liste_max.update({i*0.05 : maxi})
            liste_moyenne.update({i*0.05 :somme/len(intervalles[i])})
        else :
            liste_max.update({i*0.05 : 0})
            liste_moyenne.update({i*0.05 : 0})
             
     
    print(liste_max)
    print(liste_moyenne)
     
    x = []
    y = []
    for elt in liste_max.keys() :
        x.append(elt)
        y.append(liste_max[elt])
     
    print(x)
    print(y)
     
    figure = plt.figure(figsize=(10,7))
    plt.title('Valeur du maximum des deux valeurs en fonction de la valeur des différences entre les deux valeurs')
    axes = plt.gca()
    axes.set_xlabel('Catégorie de différences de valeurs')
    axes.set_ylabel('Valeur du max ')
    axes.xaxis.set_ticks(np.arange(0,1,0.05))
    print(np.arange(0,1,0.05))
     
    plt.scatter(x,y)
    plt.savefig("Extensions/Metrique_toutes_aretes/comparaison_diff_valeurs_all1_n1a0c0_max.png")
    plt.show()  
     
    plt.close()
     
     
    x = []
    y = []
    for elt in liste_max.keys() :
        x.append(elt)
        y.append(liste_moyenne[elt])
     
    print(x)
    print(y)
     
    figure = plt.figure(figsize=(10,7))
    plt.title('Valeur de la moyenne des maximums des deux valeurs en fonction de la valeur des différences entre les deux valeurs')
    axes = plt.gca()
    axes.set_xlabel('Catégorie de différences de valeurs')
    axes.set_ylabel('Valeur de la moyenne des maximums ')
    axes.xaxis.set_ticks(np.arange(0,1,0.05))
    print(np.arange(0,1,0.05))
     
    plt.scatter(x,y)
    plt.savefig("Extensions/Metrique_toutes_aretes/comparaison_diff_valeurs_all1_n1a0c0_moy.png")
    plt.show()  
     
    plt.close()

#if __name__ == '__main__':
    
    
    
                
    
    

                