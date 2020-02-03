'''
Created on 4 nov. 2019

@author: coline

Methode heuristique pour faire les comparaisons de graphes d'extension deux a deux
(version toutes donnees PDB)
'''

import networkx as nx
import pickle
import time
import copy
from recup_data.new_algo_comparaison import calcul_sim_aretes_avec_coeff,\
    recup_chaines, chaines_reliees
from recup_data.constantes import NEW_EXTENSION_PATH_TAILLE
from recup_data.recup_new_data import recherche_composante_connexe

''' calcul du nombre d'aretes dans le graphe passe en parametre en comptant le motif '''
def calcul_aretes_avec_coeff_heuristique(graphe, cle, coeffc, coeffa, coeffn):
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

''' calcul du nombre d'aretes dans le graphe commun passe en parametre en comptant le motif '''
def calcul_aretes_communes_avec_coeff_heuristique(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
    for u,v,data in graphe_commun.edges(data=True) :
        if data["label"] != 'B53' :
            if data["label"] == '0' :
                if coeffa == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
                    #print((u,v))
            elif graphe1.nodes[u[0]]["type"] == 1 and graphe1.nodes[v[0]]["type"] == 1 :
                if coeffc == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
                    #print((u,v))
            else :
                if coeffn == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
                    #print((u,v))
                    
    somme_aretes = somme_aretes/2
    
    return somme_aretes
       
''' calcul de la similarite du graphe commun a graphe1 et graphe2 passes en parametre en comptant le motif (version toutes aretes all1) '''
def calcul_sim_aretes_avec_coeff_heuristique(graphe1, graphe2, graphe_commun, cle, coeffc, coeffa, coeffn):
    aretes_1 = calcul_aretes_avec_coeff_heuristique(graphe1, cle, coeffc, coeffa, coeffn)
    aretes_2 = calcul_aretes_avec_coeff_heuristique(graphe2, cle, coeffc, coeffa, coeffn)
    aretes_commun = calcul_aretes_communes_avec_coeff_heuristique(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn)
    
#     print(aretes_1)
#     print(aretes_2)
#     print(aretes_commun)
    
    return aretes_commun/max(aretes_1, aretes_2)

''' renvoie vrai si le noeud1 du graphe1 et le noeud2 du graphe2 sont de meme type et appartiennent a la meme chaine 
et faux sinon'''
def meme_type_meme_chaine(noeud1, noeud2, graphe1, graphe2):
    meme_chaine = False
    for elt in graphe1.nodes[noeud1]["chaine"] :
        if elt in graphe2.nodes[noeud2]["chaine"] :
            meme_chaine = True
    if (graphe1.nodes[noeud1]["type"] == graphe2.nodes[noeud2]["type"] and meme_chaine) :
        return True
    else :
        return False

''' renvoie vrai si l'un des noeuds du couples_a_chercher se trouve deja en paire avec un autre noeud du graphe
et faux sinon '''
def dans_graphe(graphe, couple_a_chercher):
    for noeud in graphe.nodes() :
        #print(noeud)
        if (noeud[0] == couple_a_chercher[0] and noeud[1] != couple_a_chercher[1]) or (noeud[1] == couple_a_chercher[1] and noeud[0] != couple_a_chercher[0])  :
            return True
    return False

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
            
                   
#             if noeud == (1033,1034) :
#                 print(graphe1.nodes[noeud[0]]["chaine"])
            if noeud != noeud2 and meme_chaine and graphe1.nodes[noeud[0]]["type"] != None and graphe1.nodes[noeud2[0]]["type"] != None and noeud[0] not in [1,2,3,4,5] and noeud2[0] not in [1,2,3,4,5]:
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
                    if (abs(noeud[0]-modele) < abs(noeud2[0]-modele) and abs(noeud[1] - modele) > abs(noeud2[1]-modele)) or (abs(noeud[0]-modele) > abs(noeud2[0]-modele) and abs(noeud[1] - modele) < abs(noeud2[1]-modele)) :
                        return False
    
    return True

''' si l'un des noeuds du couples_a_chercher se trouve deja en paire avec un autre noeud du graphe, renvoie ce noeud
et faux sinon '''
def dans_graphe_renvoie_noeud(graphe, couple_a_chercher):
    for noeud in graphe.nodes() :
        #print(noeud)
        if (noeud[0] == couple_a_chercher[0] and noeud[1] != couple_a_chercher[1]) or (noeud[1] == couple_a_chercher[1] and noeud[0] != couple_a_chercher[0])  :
            return noeud
    return False

def test_compatibilite_renvoie_noeud(graphe_commun, noeud, graphe1, graphe2):
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
            
                   
#             if noeud == (1033,1034) :
#                 print(graphe1.nodes[noeud[0]]["chaine"])
            if noeud != noeud2 and meme_chaine and graphe1.nodes[noeud[0]]["type"] != None and graphe1.nodes[noeud2[0]]["type"] != None and noeud[0] not in [1,2,3,4,5] and noeud2[0] not in [1,2,3,4,5]:
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
                    if (abs(noeud[0]-modele) < abs(noeud2[0]-modele) and abs(noeud[1] - modele) > abs(noeud2[1]-modele)) or (abs(noeud[0]-modele) > abs(noeud2[0]-modele) and abs(noeud[1] - modele) < abs(noeud2[1]-modele)) :
                        return noeud2
    
    return True

def liste_meme_type(chaines, graphe, chaines_comp, graphe_comp):
    '''renvoie une liste de dictionnaires dans lequel sont indiques pour chaque sommet du graphe (cle du dictionnaire)
    les sommets (de graphe_comp) de la meme chaine qui sont de meme type '''
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
                dico_chaines[i].update({elt : []})
            
        
    return dico_chaines 


''' on veut essayer d'etendre la superposition de sommet1 et sommet2 aux sommets voisins
sous_graphe_commun :  le sous-graphe commun deja construit avant cette iteration
sous_graphe_commun_1 : la partie du graphe1 du sous-commun deja construit
sous_graphe_commun_2 : la partie du graphe2 du sous-commun deja construit
'''
def extension_sous_graphe_dist_1(graphe1, graphe2, sommet1, sommet2, sous_graphe_commun, sous_graphe_commun_1, sous_graphe_commun_2):
    sous_graphe = nx.MultiDiGraph()
    sous_graphe1 = nx.MultiDiGraph()
    sous_graphe2 = nx.MultiDiGraph()
#     print("gros rat")
#     print(sous_graphe_commun.nodes.data())
    if test_compatibilite(sous_graphe_commun, (sommet1, sommet2), graphe1, graphe2) : ## les voisins ne pourront pas etre compatibles si les sommets de depart ne le sont pas
        
        
        
        sous_graphe.add_node((sommet1, sommet2))
        sous_graphe1.add_node(sommet1)
        sous_graphe2.add_node(sommet2)
    
        c = 0
        liste_noeuds_a_voir = [(sommet1, sommet2)]
        while  c < len(liste_noeuds_a_voir) :
            #print(liste_noeuds_chaines)
            noeud = liste_noeuds_a_voir[c]
            #print(noeud)
            #voisins directs : on ajoute au sous-graphe
            for voisin_1 in graphe1[noeud[0]] :                                                                                                                                                                                    
                if voisin_1 not in sous_graphe_commun_1.nodes() :
                    paire_avec_voisin_1 = 0
                    for voisin_2 in graphe2[noeud[1]] :
                        if voisin_2 not in sous_graphe_commun_2.nodes() :
                            if meme_type_meme_chaine(voisin_1, voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe, (voisin_1, voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe, (voisin_1, voisin_2)) and test_compatibilite(sous_graphe_commun, (voisin_1, voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_1, voisin_2)) :
                                for edge_1 in graphe1[noeud[0]][voisin_1] :
                                    for edge_2 in graphe2[noeud[1]][voisin_2] :
                                        if graphe1[noeud[0]][voisin_1][edge_1]["label"] == graphe2[noeud[1]][voisin_2][edge_2]["label"] :
                                            

                                            
                                            paire_avec_voisin_1 += 1
                                            
                                            if (voisin_1, voisin_2) not in sous_graphe.nodes() :
                                                        sous_graphe.add_node((voisin_1, voisin_2))
                                                                                            
                                                        sous_graphe1.add_node(voisin_1)
                                                        sous_graphe2.add_node(voisin_2)
                                                    
                                            existe = False
                                            if ((noeud[0], noeud[1]), (voisin_1, voisin_2)) in sous_graphe.edges() :
                                                    for edge in sous_graphe[(noeud[0], noeud[1])][(voisin_1, voisin_2)] :
                                                        if sous_graphe[(noeud[0], noeud[1])][(voisin_1, voisin_2)][edge]["label"] == graphe1[noeud[0]][voisin_1][edge_1]["label"] :
                                                            existe = True
                                                
                                            if ((noeud[0], noeud[1]), (voisin_1, voisin_2)) not in sous_graphe.edges() or not existe :

                                            #if ((noeud[0], noeud[1]), (voisin_1, voisin_2), {'label' : graphe1[noeud[0]][voisin_1][edge_1]["label"]} ) not in sous_graphe.nodes.data() :
                                            
                                                sous_graphe.add_edge((noeud[0], noeud[1]), (voisin_1, voisin_2), label=graphe1[noeud[0]][voisin_1][edge_1]["label"])
        
                                                sous_graphe1.add_edge(noeud[0], voisin_1, label=graphe1[noeud[0]][voisin_1][edge_1]["label"])
                                                sous_graphe2.add_edge(noeud[1], voisin_2, label=graphe1[noeud[0]][voisin_1][edge_1]["label"])
                                            
                                            if graphe1[noeud[0]][voisin_1][edge_1]["label"] != "B53" :
                                                if len(graphe1[noeud[0]][voisin_1][edge_1]["label"]) == 3 :
                                                    label_inv = graphe1[noeud[0]][voisin_1][edge_1]["label"][0] + graphe1[noeud[0]][voisin_1][edge_1]["label"][2] + graphe1[noeud[0]][voisin_1][edge_1]["label"][1]
                                                else :
                                                    label_inv = graphe1[noeud[0]][voisin_1][edge_1]["label"]

                                                existe = False
                                                if ((voisin_1, voisin_2), (noeud[0], noeud[1])) in sous_graphe.edges() :
                                                        for edge in sous_graphe[(voisin_1, voisin_2)][(noeud[0], noeud[1])] :
                                                            if sous_graphe[(voisin_1, voisin_2)][(noeud[0], noeud[1])][edge]["label"] == label_inv :
                                                                existe = True
                                                    
                                                if ((voisin_1, voisin_2), (noeud[0], noeud[1])) not in sous_graphe.edges() or not existe :
           
                                                 #:and ((voisin_1, voisin_2), (noeud[0], noeud[1]), {'label' : graphe1[voisin_1][noeud[0]][edge_1]["label"]} ) not in sous_graphe.nodes.data() :
                                                    sous_graphe.add_edge((voisin_1, voisin_2), (noeud[0], noeud[1]), label=label_inv)
                                                
                                                    sous_graphe1.add_edge(voisin_1, noeud[0], label=label_inv)
                                                    sous_graphe2.add_edge(voisin_2, noeud[1], label=label_inv)

                                            if (voisin_1, voisin_2) not in liste_noeuds_a_voir :
                                            #and (voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :
                                                liste_noeuds_a_voir.append((voisin_1, voisin_2))
                                            #dico_ajout_aretes.update({(voisin_1, voisin_2) : ((noeud[0], noeud[1]), (voisin_1, voisin_2), graphe1[noeud[0]][voisin_1][edge_1]["label"] ) })
#                                                     
                                if paire_avec_voisin_1 > len(graphe1[noeud[0]][voisin_1]) :
                                    print("bizarre")
                                    print(paire_avec_voisin_1)
                                    print(sous_graphe.nodes.data())
                                            #exit(0)
                for voisin_1 in graphe1.predecessors(noeud[0]) :
                            if voisin_1 not in sous_graphe_commun_1.nodes() :
                                paire_avec_voisin_1 = 0
                                for voisin_2 in graphe2.predecessors(noeud[1]) :
                                    if voisin_2 not in sous_graphe_commun_2.nodes() :
                                        if meme_type_meme_chaine(voisin_1, voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe, (voisin_1, voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe, (voisin_1, voisin_2)) and test_compatibilite(sous_graphe_commun, (voisin_1, voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_1, voisin_2)) :
                                            for edge_1 in graphe1[voisin_1][noeud[0]] :
                                                for edge_2 in graphe2[voisin_2][noeud[1]] :
                                                    if graphe1[voisin_1][noeud[0]][edge_1]["label"] == graphe2[voisin_2][noeud[1]][edge_2]["label"] and graphe1[voisin_1][noeud[0]][edge_1]["label"] == 'B53':
                                                        #if (voisin_1, voisin_2) not in liste_noeuds_chaines :
                                                        paire_avec_voisin_1 += 1
                                                        
                                                        
                                                        if (voisin_1, voisin_2) not in sous_graphe.nodes() :
                                                                sous_graphe.add_node((voisin_1, voisin_2))
                                                                                                    
                                                                sous_graphe1.add_node(voisin_1)
                                                                sous_graphe2.add_node(voisin_2)
                                                            
                                                        #if ((voisin_1, voisin_2), (noeud[0], noeud[1]), {'label' : graphe1[voisin_1][noeud[0]][edge_1]["label"]} ) not in sous_graphe.nodes.data() :
                                                        existe = False
                                                        if ((voisin_1, voisin_2), (noeud[0], noeud[1])) in sous_graphe.edges() :
                                                                for edge in sous_graphe[(voisin_1, voisin_2)][(noeud[0], noeud[1])] :
                                                                    if sous_graphe[(voisin_1, voisin_2)][(noeud[0], noeud[1])][edge]["label"] == graphe1[voisin_1][noeud[0]][edge_1]["label"] :
                                                                        existe = True
                                                            
                                                        if ((voisin_1, voisin_2), (noeud[0], noeud[1])) not in sous_graphe.edges() or not existe :
                   
                                                            sous_graphe.add_edge((voisin_1, voisin_2), (noeud[0], noeud[1]), label=graphe1[voisin_1][noeud[0]][edge_1]["label"])
                    
                                                            sous_graphe1.add_edge(voisin_1, noeud[0], label=graphe1[voisin_1][noeud[0]][edge_1]["label"])
                                                            sous_graphe2.add_edge(voisin_2, noeud[1], label=graphe1[voisin_1][noeud[0]][edge_1]["label"])
                                                        
                                                        if (voisin_1, voisin_2) not in liste_noeuds_a_voir :
                                                            #if (voisin_1, voisin_2) not in liste_noeuds_a_voir  :#and (voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :                                                    
                                                            liste_noeuds_a_voir.append((voisin_1, voisin_2))
                                                        #dico_ajout_aretes.update({(voisin_1, voisin_2) : ((voisin_1, voisin_2), (noeud[0], noeud[1]), graphe1[voisin_1][noeud[0]][edge_1]["label"] ) })
    
                                            if paire_avec_voisin_1 > 1 :
                                                print("bizarre")
                                                print(paire_avec_voisin_1)
                                                print(sous_graphe.nodes.data())
                                                #exit(0)
            c += 1

        a_enlever = []
        for noeud in sous_graphe.nodes() :
            if len(sous_graphe[noeud]) == 0 :
                diff_b53 = False
            elif len(sous_graphe[noeud]) == 1 :
                diff_b53 = False
                for voisin in sous_graphe[noeud] :
                    for edge in sous_graphe[noeud][voisin] :
                        if sous_graphe[noeud][voisin][edge]["label"] != 'B53' :
                            diff_b53 = True
            else :
                diff_b53 = True
            if not diff_b53 :
                a_enlever.append(noeud)
                
    #     print(a_enlever)
    #     print(sous_graphe.nodes())
    #     print(sous_graphe.edges.data())
    
        for elt in a_enlever : 
            sous_graphe.remove_node(elt)    
            sous_graphe1.remove_node(elt[0])
            sous_graphe2.remove_node(elt[1])
            
        #exit(0)
#     if sommet1 == 9 and sommet2 == 17 :
#             print("gros rat") 
#             print(sous_graphe.nodes.data())    
        #exit(0)               
    return sous_graphe, sous_graphe1, sous_graphe2  


''' on veut essayer d'etendre la superposition de sommet1 et sommet2 aux sommets voisins (distance 1 et 2)
sous_graphe_commun :  le sous-graphe commun deja construit avant cette iteration
sous_graphe_commun_1 : la partie du graphe1 du sous-commun deja construit
sous_graphe_commun_2 : la partie du graphe2 du sous-commun deja construit '''
def extension_sous_graphe(graphe1, graphe2, sommet1, sommet2, sous_graphe_commun, sous_graphe_commun_1, sous_graphe_commun_2):
    sous_graphe = nx.MultiDiGraph()
    sous_graphe1 = nx.MultiDiGraph()
    sous_graphe2 = nx.MultiDiGraph()
    liste_noeuds = []
#     print("gros rat")
#     print(sous_graphe_commun.nodes.data())
    if test_compatibilite(sous_graphe_commun, (sommet1, sommet2), graphe1, graphe2) : ## les voisins ne pourront pas etre compatibles si les sommets de depart ne le sont pas
        
        
        
        sous_graphe.add_node((sommet1, sommet2))
        sous_graphe1.add_node(sommet1)
        sous_graphe2.add_node(sommet2)
    
    
        liste_noeuds_a_voir = [(sommet1, sommet2)]
        liste_noeuds_chaines = []
        liste_sous_graphe = []
        liste_sous_graphe1 = []
        liste_sous_graphe2 = []
#         if (sommet1, sommet2) in [(1,1), (2,2), (3,3), (4,4)] :
#             noeud = liste_noeuds_a_voir[0]
#             del(liste_noeuds_a_voir[0])
#             compteur_noeud = 0
#             compteur_arete = 0
#             
#             liste_noeuds_a_voir.append(noeud)
#             # voisins directs 
#             for voisin_1 in graphe1[noeud[0]] :                                                                                                                                                                                    
#                         if voisin_1 not in sous_graphe_commun_1.nodes() :
#                             paire_avec_voisin_1 = 0
#                             for voisin_2 in graphe2[noeud[1]] :
#                                 if voisin_2 not in sous_graphe_commun_2.nodes() :
#                                     if meme_type_meme_chaine(voisin_1, voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_1, voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe, (voisin_1, voisin_2)):
#                                         for edge_1 in graphe1[noeud[0]][voisin_1] :
#                                             for edge_2 in graphe2[noeud[1]][voisin_2] :
#                                                 if graphe1[noeud[0]][voisin_1][edge_1]["label"] == graphe2[noeud[1]][voisin_2][edge_2]["label"] :
#                                                     paire_avec_voisin_1 += 1
#                                                         
#                                                     if (voisin_1, voisin_2) not in liste_noeuds_a_voir :#and (voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :
#                                                         liste_noeuds_a_voir.append((voisin_1, voisin_2))
#                                                     #dico_ajout_aretes.update({(voisin_1, voisin_2) : ((noeud[0], noeud[1]), (voisin_1, voisin_2), graphe1[noeud[0]][voisin_1][edge_1]["label"] ) })
# #                                                     
#                                         if paire_avec_voisin_1 > len(graphe1[noeud[0]][voisin_1]) :
#                                             print("bizarre")
#                                             print(paire_avec_voisin_1)
#                                             print(sous_graphe.nodes.data())
#                                             #exit(0)
#             for voisin_1 in graphe1.predecessors(noeud[0]) :
#                         if voisin_1 not in sous_graphe_commun_1.nodes() :
#                             paire_avec_voisin_1 = 0
#                             for voisin_2 in graphe2.predecessors(noeud[1]) :
#                                 if voisin_2 not in sous_graphe_commun_2.nodes() :
#                                     if meme_type_meme_chaine(voisin_1, voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_1, voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe, (voisin_1, voisin_2)) :
#                                         for edge_1 in graphe1[voisin_1][noeud[0]] :
#                                             for edge_2 in graphe2[voisin_2][noeud[1]] :
#                                                 if graphe1[voisin_1][noeud[0]][edge_1]["label"] == graphe2[voisin_2][noeud[1]][edge_2]["label"] and graphe1[voisin_1][noeud[0]][edge_1]["label"] == 'B53':
#                                                     paire_avec_voisin_1 += 1
#                                                     
#                     
#                                                     compteur_arete += 1
#                                                     if (voisin_1, voisin_2) not in liste_noeuds_a_voir  :#and (voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :                                                    
#                                                         liste_noeuds_a_voir.append((voisin_1, voisin_2))
#                                                     #dico_ajout_aretes.update({(voisin_1, voisin_2) : ((voisin_1, voisin_2), (noeud[0], noeud[1]), graphe1[voisin_1][noeud[0]][edge_1]["label"] ) })
# 
#                                         if paire_avec_voisin_1 > 1 :
#                                             print("bizarre")
#                                             print(paire_avec_voisin_1)
#                                             print(sous_graphe.nodes.data())
#                                             #exit(0)
#                                             
#             ## voisins indirects dans g1
#             for voisin_1 in graphe1[noeud[0]] :
#                 if voisin_1 not in sous_graphe_commun_1.nodes() :
#                     for voisin_2 in graphe2[noeud[1]] :
#                         #if voisin_2 not in [2,3,4] :
#                             for voisin_voisin_2 in graphe2[voisin_2] :
#                                 if voisin_voisin_2 not in sous_graphe_commun_2.nodes() :
#                                         if meme_type_meme_chaine(voisin_1, voisin_voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_1, voisin_voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_1, voisin_voisin_2)) :
#                                             if (voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir  and (voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:
#                                                 liste_noeuds_a_voir.append((voisin_1, voisin_voisin_2))
#                             for voisin_voisin_2 in graphe2.predecessors(voisin_2) :                    
#                                 if voisin_voisin_2 not in sous_graphe_commun_2.nodes() :
#                                             if meme_type_meme_chaine(voisin_1, voisin_voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_1, voisin_voisin_2), graphe1, graphe2) and  not dans_graphe(sous_graphe_commun, (voisin_1, voisin_voisin_2))  :                                                 
#                                                 if (voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir  and (voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:
#                                                     liste_noeuds_a_voir.append((voisin_1, voisin_voisin_2))
#             for voisin_1 in graphe1.predecessors(noeud[0]) :
#                         if voisin_1 not in sous_graphe_commun_1.nodes() :
#                             paire_avec_voisin_1 = 0
#                             for voisin_2 in graphe2.predecessors(noeud[1]) :
#                                 #if voisin_2 not in [2,3,4] :
#                                     for voisin_voisin_2 in graphe2[voisin_2] :
#                                         
#                                         if voisin_voisin_2 not in sous_graphe_commun_2.nodes() :
#                                             if meme_type_meme_chaine(voisin_1, voisin_voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_1, voisin_voisin_2), graphe1, graphe2) and  not dans_graphe(sous_graphe_commun, (voisin_1, voisin_voisin_2)) :                                                 
#                                                 if (voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir  and (voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:
#                                                     liste_noeuds_a_voir.append((voisin_1, voisin_voisin_2))
#                                     for voisin_voisin_2 in graphe2.predecessors(voisin_2) :
#                                         if voisin_voisin_2 not in sous_graphe_commun_2.nodes() :
#                                             if meme_type_meme_chaine(voisin_1, voisin_voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_1, voisin_voisin_2), graphe1, graphe2) and  not dans_graphe(sous_graphe_commun, (voisin_1, voisin_voisin_2)) :                                                 
#                                                 if (voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir  and (voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :
#                                                     liste_noeuds_a_voir.append((voisin_1, voisin_voisin_2))
#                                                 
#             ## voisins indirects dans g2
#             for voisin_2 in graphe2[noeud[1]] :
#                 if voisin_2 not in sous_graphe_commun_2.nodes() :
#                     for voisin_1 in graphe1[noeud[0]] :
#                         #if voisin_1 not in [2,3,4] :
#                             
#                             for voisin_voisin_1 in graphe1[voisin_1] :
#                                 if voisin_voisin_1 not in sous_graphe_commun_1.nodes() :
#                                     
#                                     if meme_type_meme_chaine(voisin_2, voisin_voisin_1, graphe2, graphe1)  and test_compatibilite(sous_graphe_commun, (voisin_voisin_1, voisin_2), graphe1, graphe2) and  not dans_graphe(sous_graphe_commun, (voisin_voisin_1, voisin_2)):
#                                             if (voisin_voisin_1, voisin_2) not in liste_noeuds_a_voir and (voisin_voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:
#                                                 liste_noeuds_a_voir.append((voisin_voisin_1, voisin_2))
#                             for voisin_voisin_1 in graphe1.predecessors(voisin_1) :    
#                                 if voisin_voisin_1 not in sous_graphe_commun_1.nodes() :
#                                             if meme_type_meme_chaine(voisin_voisin_1, voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_voisin_1, voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_voisin_1, voisin_2)) :                                                 
#                                                 if (voisin_voisin_1, voisin_2) not in liste_noeuds_a_voir  and (voisin_voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :
#                                                     liste_noeuds_a_voir.append((voisin_voisin_1, voisin_2))                    
#                             
#             for voisin_2 in graphe2.predecessors(noeud[1]) :
#                         if voisin_2 not in sous_graphe_commun_2.nodes() :
#                             for voisin_1 in graphe1.predecessors(noeud[0]) :
#                                 #if voisin_1 not in [2,3,4] :
#                                     
#                                     for voisin_voisin_1 in graphe1[voisin_1] :
#                                         if voisin_voisin_1 not in sous_graphe_commun_1.nodes() :
#                                             if meme_type_meme_chaine(voisin_voisin_1, voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_voisin_1, voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_voisin_1, voisin_2)) :                                                 
#                                                 if (voisin_voisin_1, voisin_2) not in liste_noeuds_a_voir  and (voisin_voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :
#                                                     liste_noeuds_a_voir.append((voisin_voisin_1, voisin_2))
#                                     for voisin_voisin_1 in graphe1.predecessors(voisin_1) :
#                                         if voisin_voisin_1 not in sous_graphe_commun_1.nodes() :
#                                             if meme_type_meme_chaine(voisin_voisin_1, voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_voisin_1, voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_voisin_1, voisin_2)) :                                                 
#                                                 if (voisin_voisin_1, voisin_2) not in liste_noeuds_a_voir and (voisin_voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:
#                                                     liste_noeuds_a_voir.append((voisin_voisin_1, voisin_2))
    
            
#         print("gros tas")
#         print(liste_noeuds_a_voir)    
        #exit(0)   
        c = 0
        cc = 0
#         tps1 = time.time()
        while  c < len(liste_noeuds_a_voir) :
            #print(liste_noeuds_chaines)
            if cc < len(liste_noeuds_chaines) :
                noeud = liste_noeuds_chaines[cc]
                cc += 1
    #         print(liste_noeuds_a_voir)
            #noeud_ancien = noeud
            else :
                noeud = liste_noeuds_a_voir[c]
                c += 1
                cc = 0
                del(liste_noeuds_chaines[:])
#             print(noeud)
#             print(noeud)
#             if noeud == (2,2) :
#                 print("avant")
#                 print(sous_graphe.nodes())
#                 print(liste_noeuds_a_voir)
#                 print(liste_noeuds_chaines)
#                 print(graphe1[noeud[0]])
#                 print(graphe2[noeud[0]])
#                 print(graphe1[5])
#             print(sous_graphe.nodes()) 
            
#             if noeud == (6,6) :
#                 print("rasmounif")
#                 print(sous_graphe.nodes.data())
             
            if dans_graphe(sous_graphe, noeud) or not test_compatibilite(sous_graphe, noeud, graphe1, graphe2) :
#                 if noeud == (6,7) : 
#                     print("ramousnif")
                    
                liste_sous_graphe.append(copy.deepcopy(sous_graphe))
                sous_graphe.clear()
                sous_graphe.add_node(noeud)
                liste_sous_graphe1.append(copy.deepcopy(sous_graphe1))
                sous_graphe1.clear()
                sous_graphe1.add_node(noeud[0])
                liste_sous_graphe2.append(copy.deepcopy(sous_graphe2))
                sous_graphe2.clear()
                sous_graphe2.add_node(noeud[1])
#                 for n in sous_graphe.nodes() :
#                     if not dans_graphe(new_sg, n) and test_compatibilite(new_sg, n, graphe1, graphe2):
#                         new_sg.add_node(n)
#                         new_sg1.add_node(n[0])
#                         new_sg2.add_node(n[1])
#                 for u,v,data in sous_graphe.edges(data=True) :
#                     if u in new_sg.nodes() and v in new_sg.nodes() :
#                         new_sg.add_edge(u,v,**data)
#                         new_sg1.add_edge(u[0], v[0], **data)
#                         new_sg2.add_edge(u[1], v[1], **data)
            
#                     print(n)
#                     print(dans_graphe(sous_graphe, (noeud, n)))
#                     if dans_graphe(sous_graphe, (noeud, n)) :
#                 print(new_sg.nodes.data())      
#                 liste_sous_graphe.append(copy.deepcopy(sous_graphe))
#                 sous_graphe.clear()
#                 sous_graphe = copy.deepcopy(new_sg)
#                 liste_sous_graphe1.append(copy.deepcopy(sous_graphe1))
#                 sous_graphe1.clear()
#                 sous_graphe1 = copy.deepcopy(new_sg1)
#                 liste_sous_graphe2.append(copy.deepcopy(sous_graphe2))
#                 sous_graphe2.clear()
#                 sous_graphe2 = copy.deepcopy(new_sg2)
#                 print(sous_graphe.nodes.data()) 
#                 print(liste_sous_graphe[len(liste_sous_graphe)-1].nodes.data())
#                         break
            
            #if noeud in [(1,1), (2,2), (3,3), (4,4)] :
            #voisins directs : on ajoute au sous-graphe
            for voisin_1 in graphe1[noeud[0]] :                                                                                                                                                                                    
                        if voisin_1 not in sous_graphe_commun_1.nodes() :
                            paire_avec_voisin_1 = 0
                            for voisin_2 in graphe2[noeud[1]] :
                                if voisin_2 not in sous_graphe_commun_2.nodes() : # ((voisin_1, voisin_2) not in liste_noeuds_a_voir  and (voisin_1, voisin_2) not in liste_noeuds_chaines) :
                                    if meme_type_meme_chaine(voisin_1, voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe, (voisin_1, voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe, (voisin_1, voisin_2)):
                                        for edge_1 in graphe1[noeud[0]][voisin_1] :
                                            for edge_2 in graphe2[noeud[1]][voisin_2] :
                                                if graphe1[noeud[0]][voisin_1][edge_1]["label"] == graphe2[noeud[1]][voisin_2][edge_2]["label"] :
                                                    

                                                    
                                                    paire_avec_voisin_1 += 1
                                                    
                                                    if (voisin_1, voisin_2) not in sous_graphe.nodes() :
                                                                sous_graphe.add_node((voisin_1, voisin_2))
                                                                                                    
                                                                sous_graphe1.add_node(voisin_1)
                                                                sous_graphe2.add_node(voisin_2)
                                                            
                                                    existe = False
                                                    if ((noeud[0], noeud[1]), (voisin_1, voisin_2)) in sous_graphe.edges() :
                                                            for edge in sous_graphe[(noeud[0], noeud[1])][(voisin_1, voisin_2)] :
                                                                if sous_graphe[(noeud[0], noeud[1])][(voisin_1, voisin_2)][edge]["label"] == graphe1[noeud[0]][voisin_1][edge_1]["label"] :
                                                                    existe = True
                                                        
                                                    if ((noeud[0], noeud[1]), (voisin_1, voisin_2)) not in sous_graphe.edges() or not existe :

                                                    #if ((noeud[0], noeud[1]), (voisin_1, voisin_2), {'label' : graphe1[noeud[0]][voisin_1][edge_1]["label"]} ) not in sous_graphe.nodes.data() :
                                                    
                                                        sous_graphe.add_edge((noeud[0], noeud[1]), (voisin_1, voisin_2), label=graphe1[noeud[0]][voisin_1][edge_1]["label"])
                
                                                        sous_graphe1.add_edge(noeud[0], voisin_1, label=graphe1[noeud[0]][voisin_1][edge_1]["label"])
                                                        sous_graphe2.add_edge(noeud[1], voisin_2, label=graphe1[noeud[0]][voisin_1][edge_1]["label"])
                                                    
                                                    if graphe1[noeud[0]][voisin_1][edge_1]["label"] != "B53" :
                                                        if len(graphe1[noeud[0]][voisin_1][edge_1]["label"]) == 3 :
                                                            label_inv = graphe1[noeud[0]][voisin_1][edge_1]["label"][0] + graphe1[noeud[0]][voisin_1][edge_1]["label"][2] + graphe1[noeud[0]][voisin_1][edge_1]["label"][1]
                                                        else :
                                                            label_inv = graphe1[noeud[0]][voisin_1][edge_1]["label"]
        
                                                        existe = False
                                                        if ((voisin_1, voisin_2), (noeud[0], noeud[1])) in sous_graphe.edges() :
                                                                for edge in sous_graphe[(voisin_1, voisin_2)][(noeud[0], noeud[1])] :
                                                                    if sous_graphe[(voisin_1, voisin_2)][(noeud[0], noeud[1])][edge]["label"] == label_inv :
                                                                        existe = True
                                                            
                                                        if ((voisin_1, voisin_2), (noeud[0], noeud[1])) not in sous_graphe.edges() or not existe :
                   
                                                         #:and ((voisin_1, voisin_2), (noeud[0], noeud[1]), {'label' : graphe1[voisin_1][noeud[0]][edge_1]["label"]} ) not in sous_graphe.nodes.data() :
                                                            sous_graphe.add_edge((voisin_1, voisin_2), (noeud[0], noeud[1]), label=label_inv)
                                                        
                                                            sous_graphe1.add_edge(voisin_1, noeud[0], label=label_inv)
                                                            sous_graphe2.add_edge(voisin_2, noeud[1], label=label_inv)
        
                                                    if (voisin_1, voisin_2) not in liste_noeuds_chaines and (voisin_1, voisin_2) not in liste_noeuds_a_voir :
                                                    #and (voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :
                                                        liste_noeuds_chaines.append((voisin_1, voisin_2))
                                                    #dico_ajout_aretes.update({(voisin_1, voisin_2) : ((noeud[0], noeud[1]), (voisin_1, voisin_2), graphe1[noeud[0]][voisin_1][edge_1]["label"] ) })
#                                                     
                                        if paire_avec_voisin_1 > len(graphe1[noeud[0]][voisin_1]) :
                                            print("bizarre")
                                            print(paire_avec_voisin_1)
                                            print(sous_graphe.nodes.data())
                                            #exit(0)
            for voisin_1 in graphe1.predecessors(noeud[0]) :
                            if voisin_1 not in sous_graphe_commun_1.nodes() :
                                paire_avec_voisin_1 = 0
                                for voisin_2 in graphe2.predecessors(noeud[1]) :
                                    if voisin_2 not in sous_graphe_commun_2.nodes() : # ((voisin_1, voisin_2) not in liste_noeuds_a_voir  and (voisin_1, voisin_2) not in liste_noeuds_chaines) :
                                        if meme_type_meme_chaine(voisin_1, voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe, (voisin_1, voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe, (voisin_1, voisin_2)) :
                                            for edge_1 in graphe1[voisin_1][noeud[0]] :
                                                for edge_2 in graphe2[voisin_2][noeud[1]] :
                                                    if graphe1[voisin_1][noeud[0]][edge_1]["label"] == graphe2[voisin_2][noeud[1]][edge_2]["label"] and graphe1[voisin_1][noeud[0]][edge_1]["label"] == 'B53':
                                                        #if (voisin_1, voisin_2) not in liste_noeuds_chaines :
                                                        paire_avec_voisin_1 += 1
                                                        
                                                        
                                                        if (voisin_1, voisin_2) not in sous_graphe.nodes() :
                                                                sous_graphe.add_node((voisin_1, voisin_2))
                                                                                                    
                                                                sous_graphe1.add_node(voisin_1)
                                                                sous_graphe2.add_node(voisin_2)
                                                            
                                                        #if ((voisin_1, voisin_2), (noeud[0], noeud[1]), {'label' : graphe1[voisin_1][noeud[0]][edge_1]["label"]} ) not in sous_graphe.nodes.data() :
                                                        existe = False
                                                        if ((voisin_1, voisin_2), (noeud[0], noeud[1])) in sous_graphe.edges() :
                                                                for edge in sous_graphe[(voisin_1, voisin_2)][(noeud[0], noeud[1])] :
                                                                    if sous_graphe[(voisin_1, voisin_2)][(noeud[0], noeud[1])][edge]["label"] == graphe1[voisin_1][noeud[0]][edge_1]["label"] :
                                                                        existe = True
                                                            
                                                        if ((voisin_1, voisin_2), (noeud[0], noeud[1])) not in sous_graphe.edges() or not existe :
                   
                                                            sous_graphe.add_edge((voisin_1, voisin_2), (noeud[0], noeud[1]), label=graphe1[voisin_1][noeud[0]][edge_1]["label"])
                    
                                                            sous_graphe1.add_edge(voisin_1, noeud[0], label=graphe1[voisin_1][noeud[0]][edge_1]["label"])
                                                            sous_graphe2.add_edge(voisin_2, noeud[1], label=graphe1[voisin_1][noeud[0]][edge_1]["label"])
                                                        
                                                        if (voisin_1, voisin_2) not in liste_noeuds_chaines and (voisin_1, voisin_2) not in liste_noeuds_a_voir :
                                                            #if (voisin_1, voisin_2) not in liste_noeuds_a_voir  :#and (voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :                                                    
                                                            liste_noeuds_chaines.append((voisin_1, voisin_2))
                                                        #dico_ajout_aretes.update({(voisin_1, voisin_2) : ((voisin_1, voisin_2), (noeud[0], noeud[1]), graphe1[voisin_1][noeud[0]][edge_1]["label"] ) })
    
                                            if paire_avec_voisin_1 > 1 :
                                                print("bizarre")
                                                print(paire_avec_voisin_1)
                                                print(sous_graphe.nodes.data())
                                                #exit(0)
                
            ## voisins indirects : on n'ajoute pas encore au sous-graphe                                
            ## voisins indirects dans g1
            for voisin_1 in graphe1[noeud[0]] :
                    if voisin_1 not in sous_graphe_commun_1.nodes() :
                        for voisin_2 in graphe2[noeud[1]] :
                            #if voisin_2 not in [2,3,4] :
                                for voisin_voisin_2 in graphe2[voisin_2] :
                                    if voisin_voisin_2 not in sous_graphe_commun_2.nodes() and (voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir  and (voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:
                                        if meme_type_meme_chaine(voisin_1, voisin_voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_1, voisin_voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_1, voisin_voisin_2)) :
                                                #if (voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir  and (voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:
                                                liste_noeuds_a_voir.append((voisin_1, voisin_voisin_2))
                                for voisin_voisin_2 in graphe2.predecessors(voisin_2) :                    
                                    if voisin_voisin_2 not in sous_graphe_commun_2.nodes() :
                                                if (voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir  and (voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:

                                                    if meme_type_meme_chaine(voisin_1, voisin_voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_1, voisin_voisin_2), graphe1, graphe2) and  not dans_graphe(sous_graphe_commun, (voisin_1, voisin_voisin_2))  :                                                 
                                                        liste_noeuds_a_voir.append((voisin_1, voisin_voisin_2))
            for voisin_1 in graphe1.predecessors(noeud[0]) :
                            if voisin_1 not in sous_graphe_commun_1.nodes() :
                                paire_avec_voisin_1 = 0
                                for voisin_2 in graphe2.predecessors(noeud[1]) :
                                    #if voisin_2 not in [2,3,4] :
                                        for voisin_voisin_2 in graphe2[voisin_2] :
                                            
                                            if voisin_voisin_2 not in sous_graphe_commun_2.nodes() :
                                                if (voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir  and (voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:
                                                    if meme_type_meme_chaine(voisin_1, voisin_voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_1, voisin_voisin_2), graphe1, graphe2) and  not dans_graphe(sous_graphe_commun, (voisin_1, voisin_voisin_2)) :                                                 
                                                        liste_noeuds_a_voir.append((voisin_1, voisin_voisin_2))
                                        for voisin_voisin_2 in graphe2.predecessors(voisin_2) :
                                            if voisin_voisin_2 not in sous_graphe_commun_2.nodes() :
                                                if (voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir  and (voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :
                                                    if meme_type_meme_chaine(voisin_1, voisin_voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_1, voisin_voisin_2), graphe1, graphe2) and  not dans_graphe(sous_graphe_commun, (voisin_1, voisin_voisin_2)) :                                                 
                                                        liste_noeuds_a_voir.append((voisin_1, voisin_voisin_2))
                                                    
            ## voisins indirects dans g2
            for voisin_2 in graphe2[noeud[1]] :
                    if voisin_2 not in sous_graphe_commun_2.nodes() :
                        for voisin_1 in graphe1[noeud[0]] :
                            #if voisin_1 not in [2,3,4] :
                                
                                for voisin_voisin_1 in graphe1[voisin_1] :
                                    if voisin_voisin_1 not in sous_graphe_commun_1.nodes() :
                                        if (voisin_voisin_1, voisin_2) not in liste_noeuds_a_voir and (voisin_voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:
                                            if meme_type_meme_chaine(voisin_2, voisin_voisin_1, graphe2, graphe1)  and test_compatibilite(sous_graphe_commun, (voisin_voisin_1, voisin_2), graphe1, graphe2) and  not dans_graphe(sous_graphe_commun, (voisin_voisin_1, voisin_2)):
                                                liste_noeuds_a_voir.append((voisin_voisin_1, voisin_2))
                                for voisin_voisin_1 in graphe1.predecessors(voisin_1) :    
                                    if voisin_voisin_1 not in sous_graphe_commun_1.nodes() :
                                                if (voisin_voisin_1, voisin_2) not in liste_noeuds_a_voir  and (voisin_voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :
                                                    if meme_type_meme_chaine(voisin_voisin_1, voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_voisin_1, voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_voisin_1, voisin_2)) :                                                 
                                                        liste_noeuds_a_voir.append((voisin_voisin_1, voisin_2))                    
                                
            for voisin_2 in graphe2.predecessors(noeud[1]) :
                            if voisin_2 not in sous_graphe_commun_2.nodes() :
                                for voisin_1 in graphe1.predecessors(noeud[0]) :
                                    #if voisin_1 not in [2,3,4] :
                                        
                                        for voisin_voisin_1 in graphe1[voisin_1] :
                                            if voisin_voisin_1 not in sous_graphe_commun_1.nodes() :
                                                if (voisin_voisin_1, voisin_2) not in liste_noeuds_a_voir  and (voisin_voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :
                                                    if meme_type_meme_chaine(voisin_voisin_1, voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_voisin_1, voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_voisin_1, voisin_2)) :                                                 
                                                        liste_noeuds_a_voir.append((voisin_voisin_1, voisin_2))
                                        for voisin_voisin_1 in graphe1.predecessors(voisin_1) :
                                            if voisin_voisin_1 not in sous_graphe_commun_1.nodes() :
                                                if (voisin_voisin_1, voisin_2) not in liste_noeuds_a_voir and (voisin_voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:
                                                    if meme_type_meme_chaine(voisin_voisin_1, voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_voisin_1, voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_voisin_1, voisin_2)) :                                                 
                                                        liste_noeuds_a_voir.append((voisin_voisin_1, voisin_2))
            
            ## voisins indirects dans g1 et g2  
            for voisin_2 in graphe2[noeud[1]] :
                for voisin_voisin_2 in graphe2[voisin_2] :
                    if voisin_voisin_2 not in sous_graphe_commun_2.nodes() :
                        for voisin_1 in graphe1[noeud[0]] :
                            #if voisin_1 not in [2,3,4] :
                                for voisin_voisin_1 in graphe1[voisin_1] :
                                    if voisin_voisin_1 not in sous_graphe_commun_1.nodes() :
                                        if (voisin_voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir and (voisin_voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:
                                            if meme_type_meme_chaine(voisin_voisin_2, voisin_voisin_1, graphe2, graphe1)  and test_compatibilite(sous_graphe_commun, (voisin_voisin_1, voisin_voisin_2), graphe1, graphe2) and  not dans_graphe(sous_graphe_commun, (voisin_voisin_1, voisin_voisin_2)):
                                                liste_noeuds_a_voir.append((voisin_voisin_1, voisin_voisin_2))
                                                    
                                for voisin_voisin_1 in graphe1.predecessors(voisin_1) :    
                                    if voisin_voisin_1 not in sous_graphe_commun_1.nodes() :
                                        if (voisin_voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir  and (voisin_voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :
                                            if meme_type_meme_chaine(voisin_voisin_1, voisin_voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_voisin_1, voisin_voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_voisin_1, voisin_voisin_2)) :                                                 
                                                liste_noeuds_a_voir.append((voisin_voisin_1, voisin_voisin_2))                    
            ## voisins indirects 1 et 2 NEW                                            
            for voisin_2 in graphe2.predecessors(noeud[1]) :
                for voisin_voisin_2 in graphe2.predecessors(voisin_2) :
                            if voisin_voisin_2 not in sous_graphe_commun_2.nodes() :
                                for voisin_1 in graphe1.predecessors(noeud[0]) :
                                    #if voisin_1 not in [2,3,4] :
                                         
                                        for voisin_voisin_1 in graphe1[voisin_1] :
                                            if voisin_voisin_1 not in sous_graphe_commun_1.nodes() :
                                                if (voisin_voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir  and (voisin_voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :
                                                    if meme_type_meme_chaine(voisin_voisin_1, voisin_voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_voisin_1, voisin_voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_voisin_1, voisin_voisin_2)) :                                                 
                                                        liste_noeuds_a_voir.append((voisin_voisin_1, voisin_voisin_2))
                                        for voisin_voisin_1 in graphe1.predecessors(voisin_1) :
                                            if voisin_voisin_1 not in sous_graphe_commun_1.nodes() :
                                                if (voisin_voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir and (voisin_voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:
                                                    if meme_type_meme_chaine(voisin_voisin_1, voisin_voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_voisin_1, voisin_voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_voisin_1, voisin_voisin_2)) :                                                 
                                                        liste_noeuds_a_voir.append((voisin_voisin_1, voisin_voisin_2))
#                                                         print(voisin_1)
#                                                         print(voisin_2)
    
            
#             if noeud == (2,2) :
#                 print("apres")
#                 print(sous_graphe.nodes())
#                 print(liste_noeuds_a_voir)
#                 print(liste_noeuds_chaines)
                #exit(0)
#                 else :
#                 for voisin_1 in graphe1[noeud[0]] :                                                                                                                                                                                    
#                             if voisin_1 not in sous_graphe_commun_1.nodes() :
#                                 paire_avec_voisin_1 = 0
#                                 for voisin_2 in graphe2[noeud[1]] :
#                                     if voisin_2 not in sous_graphe_commun_2.nodes() :
#                                         if meme_type_meme_chaine(voisin_1, voisin_2, graphe1, graphe2) and not dans_graphe(sous_graphe, (voisin_1, voisin_2)) and test_compatibilite(sous_graphe, (voisin_1, voisin_2), graphe1, graphe2) : #and noeud_ancien != (voisin_1, voisin_2):#and (voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:
#                                             for edge_1 in graphe1[noeud[0]][voisin_1] :
#                                                 for edge_2 in graphe2[noeud[1]][voisin_2] :
#                                                     if graphe1[noeud[0]][voisin_1][edge_1]["label"] == graphe2[noeud[1]][voisin_2][edge_2]["label"] :
#                                                         #if (voisin_1, voisin_2) not in liste_noeuds_chaines :
#                                                         #if graphe1[noeud[0]][voisin_1][edge_1]["label"] != 'B53' or noeud not in [(1,1), (2,2), (3,3), (4,4)] : #or (noeud in [(1,1), (2,2), (3,3), (4,4)] and (voisin_1, voisin_2) in [(1,1), (2,2), (3,3), (4,4)]): 
# #                                                         if noeud == (6,6) :
# #                                                             print("rat")
# #                                                             print(voisin_1)
# #                                                             print(voisin_2)
#                                                         paire_avec_voisin_1 += 1
#                                                         if (voisin_1, voisin_2) not in sous_graphe.nodes() :
#                                                             sous_graphe.add_node((voisin_1, voisin_2))
#                                                                                                 
#                                                             sous_graphe1.add_node(voisin_1)
#                                                             sous_graphe2.add_node(voisin_2)
#                                                         
#                                                         existe = False
#                                                         if ((noeud[0], noeud[1]), (voisin_1, voisin_2)) in sous_graphe.edges() :
#                                                                 for edge in sous_graphe[(noeud[0], noeud[1])][(voisin_1, voisin_2)] :
#                                                                     if sous_graphe[(noeud[0], noeud[1])][(voisin_1, voisin_2)][edge]["label"] == graphe1[noeud[0]][voisin_1][edge_1]["label"] :
#                                                                         existe = True
#                                                             
#                                                         if ((noeud[0], noeud[1]), (voisin_1, voisin_2)) not in sous_graphe.edges() or not existe :
#     
#                                                         #if ((noeud[0], noeud[1]), (voisin_1, voisin_2), {'label' : graphe1[noeud[0]][voisin_1][edge_1]["label"]} ) not in sous_graphe.nodes.data() :
#                                                         
#                                                             sous_graphe.add_edge((noeud[0], noeud[1]), (voisin_1, voisin_2), label=graphe1[noeud[0]][voisin_1][edge_1]["label"])
#                     
#                                                             sous_graphe1.add_edge(noeud[0], voisin_1, label=graphe1[noeud[0]][voisin_1][edge_1]["label"])
#                                                             sous_graphe2.add_edge(noeud[1], voisin_2, label=graphe1[noeud[0]][voisin_1][edge_1]["label"])
#                                                         
#                                                         if graphe1[noeud[0]][voisin_1][edge_1]["label"] != "B53" :
#                                                             if len(graphe1[noeud[0]][voisin_1][edge_1]["label"]) == 3 :
#                                                                     label_inv = graphe1[noeud[0]][voisin_1][edge_1]["label"][0] + graphe1[noeud[0]][voisin_1][edge_1]["label"][2] + graphe1[noeud[0]][voisin_1][edge_1]["label"][1]
#                                                             else :
#                                                                     label_inv = graphe1[noeud[0]][voisin_1][edge_1]["label"]
#                                                             existe = False
#                                                             if ((voisin_1, voisin_2), (noeud[0], noeud[1])) in sous_graphe.edges() :
#                                                                     for edge in sous_graphe[(voisin_1, voisin_2)][(noeud[0], noeud[1])] :
#                                                                         if sous_graphe[(voisin_1, voisin_2)][(noeud[0], noeud[1])][edge]["label"] == label_inv :
#                                                                             existe = True
#                                                                 
#                                                             if ((voisin_1, voisin_2), (noeud[0], noeud[1])) not in sous_graphe.edges() or not existe :
# #                                                                 print(voisin_1, voisin_2)
# #                                                                 print(noeud)
#                                                              #:and ((voisin_1, voisin_2), (noeud[0], noeud[1]), {'label' : graphe1[voisin_1][noeud[0]][edge_1]["label"]} ) not in sous_graphe.nodes.data() :
#                                                                 
#                                                                 sous_graphe.add_edge((voisin_1, voisin_2), (noeud[0], noeud[1]), label=label_inv)
#                                                             
#                                                                 sous_graphe1.add_edge(voisin_1, noeud[0], label=label_inv)
#                                                                 sous_graphe2.add_edge(voisin_2, noeud[1], label=label_inv)
#                                                                 
#                                                         if (voisin_1, voisin_2) not in liste_noeuds_chaines and (voisin_1, voisin_2) not in liste_noeuds_a_voir :
#                                                             #if noeud not in [(1,1), (2,2), (3,3), (4,4)] or (voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] and (voisin_1, voisin_2) not in liste_noeuds_a_voir :
#                                                             liste_noeuds_chaines.append((voisin_1, voisin_2))
#                                             if paire_avec_voisin_1 > len(graphe1[noeud[0]][voisin_1]) :
#                                                 print("bizarre")
#                                                 exit(0)
#                 for voisin_1 in graphe1.predecessors(noeud[0]) :
#                             if voisin_1 not in sous_graphe_commun_1.nodes() :
#                                 paire_avec_voisin_1 = 0
#                                 for voisin_2 in graphe2.predecessors(noeud[1]) :
#                                     if voisin_2 not in sous_graphe_commun_2.nodes() :
#                                         if meme_type_meme_chaine(voisin_1, voisin_2, graphe1, graphe2) and not dans_graphe(sous_graphe, (voisin_1, voisin_2)) and test_compatibilite(sous_graphe, (voisin_1, voisin_2), graphe1, graphe2) : #and noeud_ancien != (voisin_1, voisin_2):# and ((voisin_1, voisin_2), (noeud[0], noeud[1])) not in sous_graphe.edges() :#and (voisin_1, voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :
#                                             for edge_1 in graphe1[voisin_1][noeud[0]] :
#                                                 for edge_2 in graphe2[voisin_2][noeud[1]] :
#                                                     if graphe1[voisin_1][noeud[0]][edge_1]["label"] == graphe2[voisin_2][noeud[1]][edge_2]["label"] and graphe1[voisin_1][noeud[0]][edge_1]["label"] == 'B53':
#                                                         #if noeud not in [(1,1), (2,2), (3,3), (4,4)] or (noeud in [(1,1), (2,2), (3,3), (4,4)] and (voisin_1, voisin_2) in [(1,1), (2,2), (3,3), (4,4)]) : 
#                                                         #if (voisin_1, voisin_2) not in liste_noeuds_chaines :
#                                                         paire_avec_voisin_1 += 1
#                                                         if (voisin_1, voisin_2) not in sous_graphe.nodes() :
#                                                             sous_graphe.add_node((voisin_1, voisin_2))
#                                                                                                 
#                                                             sous_graphe1.add_node(voisin_1)
#                                                             sous_graphe2.add_node(voisin_2)
#                                                         
#                                                         #if ((voisin_1, voisin_2), (noeud[0], noeud[1]), {'label' : graphe1[voisin_1][noeud[0]][edge_1]["label"]} ) not in sous_graphe.nodes.data() :
#                                                         existe = False
#                                                         if ((voisin_1, voisin_2), (noeud[0], noeud[1])) in sous_graphe.edges() :
#                                                                 for edge in sous_graphe[(voisin_1, voisin_2)][(noeud[0], noeud[1])] :
#                                                                     if sous_graphe[(voisin_1, voisin_2)][(noeud[0], noeud[1])][edge]["label"] == graphe1[voisin_1][noeud[0]][edge_1]["label"] :
#                                                                         existe = True
#                                                             
#                                                         if ((voisin_1, voisin_2), (noeud[0], noeud[1])) not in sous_graphe.edges() or not existe :
#                    
#                                                             sous_graphe.add_edge((voisin_1, voisin_2), (noeud[0], noeud[1]), label=graphe1[voisin_1][noeud[0]][edge_1]["label"])
#                     
#                                                             sous_graphe1.add_edge(voisin_1, noeud[0], label=graphe1[voisin_1][noeud[0]][edge_1]["label"])
#                                                             sous_graphe2.add_edge(voisin_2, noeud[1], label=graphe1[voisin_1][noeud[0]][edge_1]["label"])
#                                                         
#                                                         if (voisin_1, voisin_2) not in liste_noeuds_chaines and (voisin_1, voisin_2) not in liste_noeuds_a_voir :
#                                                             #if (voisin_1, voisin_2) not in liste_noeuds_a_voir :
#                                                             #if (voisin_1, voisin_2) not in liste_noeuds_a_voir : #not (noeud in [(1,1), (2,2), (3,3), (4,4)] and (voisin_1, voisin_2) in [(1,1), (2,2), (3,3), (4,4)]) and (voisin_1, voisin_2) not in liste_noeuds_a_voir:
#                                                             liste_noeuds_chaines.append((voisin_1, voisin_2))
#                                             if paire_avec_voisin_1 > 1 :
#                                                 print("bizarre")
#                                                 exit(0)
#             print(sous_graphe.nodes.data())
#             print(sous_graphe.edges.data())
#             print(liste_noeuds_a_voir)
#             if c == 3 :
#                 exit(0)
            
        #exit(0)    
#         print(liste_sous_graphe)
#         print("ramou")
#         
#         print(liste_noeuds_a_voir)
#         print(liste_noeuds_chaines)
        #if len(liste_sous_graphe) == 0 :
#         tps2 = time.time()
#         print(tps2-tps1)
        liste_sous_graphe.append(sous_graphe)
        liste_sous_graphe1.append(sous_graphe1)
        liste_sous_graphe2.append(sous_graphe2)           
        maxi_sim = 0
        sous_graphe_maxi_sim = -1
        compteur = 0

        #print("petit rat")
        sous_graphe_max = nx.MultiDiGraph()
        compteur = 0
#         print("rasmo")
# #         
#         for sous_graphe in liste_sous_graphe :
#             print(sous_graphe.nodes())
# #             print(sous_graphe.edges.data())
#         
#         print("rasmo2")
#         tps1 = time.time()
        for sous_graphe in liste_sous_graphe :
#             print(sous_graphe.nodes())
            if sous_graphe_max.number_of_nodes() == 0 :
                sous_graphe_max = nx.union(sous_graphe_max,sous_graphe)
            else :
                    couples_pbs = []
                
                    liste_noeuds_traites = []
                    
                    #del(couples_pbs[:])
                    for noeud in sous_graphe.nodes() :
                            #if noeud not in liste_noeuds_traites : 
                                res_compa_seq = test_compatibilite_renvoie_noeud(sous_graphe_max, noeud, graphe1, graphe2)
                                res_compa_doub = dans_graphe_renvoie_noeud(sous_graphe_max, noeud)
                                if res_compa_seq != True :
                                    couples_pbs.append((noeud, res_compa_seq))
                                if res_compa_doub != False :
                                    couples_pbs.append((noeud, res_compa_doub))
                    
                    while len(couples_pbs) > 0 :
                        
#                         print("gros tas")
#                         print(sous_graphe.nodes())
#                         print(sous_graphe_max.nodes())
#                         print(couples_pbs)
#                         print(liste_noeuds_traites)
                        #for elt in couples_pbs :
                        elt = couples_pbs[0]
                        liste_a_enlever_1 = [elt[0]]
                        subgraph1 = nx.MultiDiGraph()
                        liste_sommets = [elt[0]]
#                         if elt[0] not in liste_noeuds_traites :
#                             liste_noeuds_traites.append(elt[0])
                        #liste_noeuds_pb.append(elt[0])
                        for voisin in sous_graphe[elt[0]] :
                            if voisin not in liste_sommets :
                                liste_sommets.append(voisin)
                            for edge in sous_graphe[elt[0]][voisin] :
                                if sous_graphe[elt[0]][voisin][edge]["label"] != 'B53' :
                                    if voisin not in liste_a_enlever_1 :
                                        liste_a_enlever_1.append(voisin)
#                             if voisin not in liste_noeuds_traites :
#                                 liste_noeuds_traites.append(voisin)
                                #liste_noeuds_pb.append(voisin)
                        for voisin in sous_graphe.predecessors(elt[0]) :
                            if voisin not in liste_sommets :
                                liste_sommets.append(voisin)
                            for edge in sous_graphe[voisin][elt[0]] :
                                if sous_graphe[voisin][elt[0]][edge]["label"] != 'B53' :
                                    if voisin not in liste_a_enlever_1 :
                                        liste_a_enlever_1.append(voisin)
#                             if voisin not in liste_noeuds_traites :
#                                 liste_noeuds_traites.append(voisin)
                                #liste_noeuds_pb.append(voisin)
                        subgraph1 = sous_graphe.subgraph(liste_sommets)   
                        
                        sim1 = calcul_sim_aretes_avec_coeff_heuristique(graphe1, graphe2, subgraph1, "petit rat", 1, 1, 1)
                        
                        subgraph2 = nx.MultiDiGraph()
                        
                        liste_a_enlever_2 = [elt[1]]
                        liste_sommets = [elt[1]]
#                         print(sous_graphe.nodes())
#                         print(elt[0])
#                         print(sous_graphe_max.nodes())
#                         print(elt[1])
                        
                        for voisin in sous_graphe_max[elt[1]] :
                            if voisin not in liste_sommets :
                                liste_sommets.append(voisin)
                            for edge in sous_graphe_max[elt[1]][voisin] :
                                if sous_graphe_max[elt[1]][voisin][edge]["label"] != 'B53' :
                                    if voisin not in liste_a_enlever_2 :
                                        liste_a_enlever_2.append(voisin)
                                    
                        for voisin in sous_graphe_max.predecessors(elt[1]) :
                            if voisin not in liste_sommets :
                                liste_sommets.append(voisin)
                            for edge in sous_graphe_max[voisin][elt[1]] :
                                if sous_graphe_max[voisin][elt[1]][edge]["label"] != 'B53' :
                                    if voisin not in liste_a_enlever_2 :
                                        liste_a_enlever_2.append(voisin)
                
                        subgraph2 = sous_graphe_max.subgraph(liste_sommets)   
                        
                        sim2 = calcul_sim_aretes_avec_coeff_heuristique(graphe1, graphe2, subgraph2, "petit rat", 1, 1, 1)
#                         if elt == ((14, 12), (10, 12)) :
#                             print(sim1)
#                             print(sim2) 
                        
                        if sim1 > sim2 :
                                
                            for noeud in liste_a_enlever_2 :
                                if noeud in sous_graphe_max.nodes() :
                                    sous_graphe_max.remove_node(noeud)
                            
                            #sous_graphe_max = nx.union(a_garder_sous_graphes[c], sous_graphe_max)
                            
#                                 for noeud in subgraph1.nodes() :
#                                     if noeud not in sous_graphe_max.nodes()  :
#                                         sous_graphe_max.add_node(noeud)
#                                 for u,v,data in subgraph1.edges(data=True) :
#                                     
#                                     if (u,v) not in sous_graphe_max.edges() :
#                                                 sous_graphe_max.add_edge(u,v,**data)
#                                     else :
#                                                 existe = False
#                                                 for edge in sous_graphe_max[u][v] :
#                                                     if sous_graphe_max[u][v][edge]["label"] == data["label"] :
#                                                         existe = True
#                                                 if not existe :
#                                                     sous_graphe_max.add_edge(u,v,**data) 
#                             if elt == ((14, 12), (10, 12)) :
#                                 print(sim1)
#                                 print(liste_a_enlever_2)
#                                 print(sous_graphe_max.nodes())
    
                            
                        else :
                            for noeud in liste_a_enlever_1 :
                                if noeud in sous_graphe.nodes() :
                                    sous_graphe.remove_node(noeud)
#                                 a_garder.append(1)
#                                 a_garder_sous_graphes.append(subgraph2)
#                                 a_enlever_sous_graphes.append(subgraph1)
                        del(couples_pbs[:])
                        for noeud in sous_graphe.nodes() :
                                #if noeud not in liste_noeuds_traites : 
                                    res_compa_seq = test_compatibilite_renvoie_noeud(sous_graphe_max, noeud, graphe1, graphe2)
                                    res_compa_doub = dans_graphe_renvoie_noeud(sous_graphe_max, noeud)
                                    if res_compa_seq != True :
                                        couples_pbs.append((noeud, res_compa_seq))
                                    if res_compa_doub != False :
                                        couples_pbs.append((noeud, res_compa_doub))

#                     c = 0
#                     for elt in couples_pbs :
#                         if a_garder[c] == 0 :
#                             print("gros rat")
#                             print(a_enlever_sous_graphes[c])
#                             print(a_garder_sous_graphes[c].nodes())
                            
#                         else :
#                             for noeud in sous_graphe.nodes() :
#                                 if noeud not in a_garder_sous_graphes[c].nodes() :
#                                     sous_graphe_max.add_node(noeud)
#                                     
#                             for u,v,data in sous_graphe.edges(data=True) :
#                                 if (u,v) not in subgraph1.edges() :
#                                     sous_graphe_max.add_edge(u,v,**data)
#                                 else :
#                                     existe = False
#                                     for edge in subgraph1[u][v] :
#                                         if subgraph1[u][v][edge]["label"] == data["label"] :
#                                             existe = True
#                                     if not existe :
#                                         sous_graphe_max.add_edge(u,v,**data)
                    
                        #c += 1
                    
#                     print(sous_graphe_max.nodes())
#                     print("rat")
#                     print(liste_noeuds_pb)
                    for noeud in sous_graphe.nodes() :
                        if noeud not in sous_graphe_max.nodes():
                            sous_graphe_max.add_node(noeud)
                    
                    for u,v,data in sous_graphe.edges(data=True) :
                        #if u not in liste_noeuds_traites and v not in liste_noeuds_traites :
                            if (u,v) not in sous_graphe_max.edges() :
                                    sous_graphe_max.add_edge(u,v,**data)
                            else :
                                    existe = False
                                    for edge in sous_graphe_max[u][v] :
                                        if sous_graphe_max[u][v][edge]["label"] == data["label"] :
                                            existe = True
                                    if not existe :
                                        sous_graphe_max.add_edge(u,v,**data)      
                    
#                 print(liste_sous_graphe_max[0].nodes())
#                 print(liste_sous_graphe_max[1].nodes())
#                 print(liste_sous_graphe_max[2].nodes())
#             if compteur == 4 :
#                 exit(0)
            compteur += 1
#         print("gros tas garde")
#         print(sous_graphe_max.nodes())
#         print(sous_graphe_max.edges())
        #exit(0)            
#         maxi_sim = 0
#         sous_graphe_maxi_sim = -1
#         tps2 = time.time()
#         print(tps2-tps1)
        compteur = 0
#         for sous_graphe in liste_sous_graphe_max :
#                 
#             print(sous_graphe.nodes.data())
#             sim = calcul_sim_aretes_avec_coeff_heuristique(graphe1, graphe2, sous_graphe, "petit rat", 1, 1, 1)
#             print("tralala")
#             print(sous_graphe.nodes.data())
#             print(sous_graphe.edges.data())
#             print(sim)
#             if sim > maxi_sim :
#                 
#                 maxi_sim = sim
#                 sous_graphe_maxi_sim = compteur
#             compteur += 1
        
        
        
        new_sous_graphe = copy.deepcopy(sous_graphe_max)
        
        new_sous_graphe1 = nx.MultiDiGraph()
        new_sous_graphe2 = nx.MultiDiGraph()
        for noeud in new_sous_graphe.nodes() :
            liste_noeuds.append(noeud)
            new_sous_graphe1.add_node(noeud[0], **graphe1.nodes[noeud[0]])
            new_sous_graphe2.add_node(noeud[1], **graphe2.nodes[noeud[1]])
        
        for u,v,data in new_sous_graphe.edges(data=True) :
            new_sous_graphe1.add_edge(u[0], v[0], **data)
            new_sous_graphe2.add_edge(u[1], v[1], **data)
        
#         new_sous_graphe1 = copy.deepcopy(liste_sous_graphe1[sous_graphe_maxi_sim])
#         new_sous_graphe2 = copy.deepcopy(liste_sous_graphe2[sous_graphe_maxi_sim])
#         print(new_sous_graphe.nodes.data())
#         print(new_sous_graphe.edges.data())
#         print(new_sous_graphe1.nodes.data())
#         print(new_sous_graphe2.nodes.data())
        ## recherche des noeuds en bout de chaine qui ne sont relies que par des liaisons b53
        
        #exit(0)
        a_enlever = []
        for noeud in new_sous_graphe.nodes() :
            if len(new_sous_graphe[noeud]) == 0 :
                diff_b53 = False
            elif len(new_sous_graphe[noeud]) == 1 :
                diff_b53 = False
                for voisin in new_sous_graphe[noeud] :
                    for edge in new_sous_graphe[noeud][voisin] :
                        if new_sous_graphe[noeud][voisin][edge]["label"] != 'B53' :
                            diff_b53 = True
            else :
                diff_b53 = True
            if not diff_b53 :
                a_enlever.append(noeud)
                
    #     print(a_enlever)
    #     print(sous_graphe.nodes())
    #     print(sous_graphe.edges.data())
    
        for elt in a_enlever :
            liste_noeuds.remove(elt)
            new_sous_graphe.remove_node(elt)    
            new_sous_graphe1.remove_node(elt[0])
            new_sous_graphe2.remove_node(elt[1])
        
           
        #exit(0)
#     if sommet1 == 9 and sommet2 == 17 :
#             print("gros rat") 
#             print(sous_graphe.nodes.data())    
        #exit(0)               
        return new_sous_graphe, new_sous_graphe1, new_sous_graphe2, liste_noeuds 
    return sous_graphe, sous_graphe1, sous_graphe2, liste_noeuds    
                    
''' renvoie vrai si les deux sous-graphes sont identiques et faux sinon '''
def sous_graphe_identique(sous_graphe1, sous_graphe2):
    if sous_graphe1.number_of_nodes() != sous_graphe2.number_of_nodes() or sous_graphe1.number_of_edges() != sous_graphe2.number_of_edges() :
        return False
    
    noeuds_vus = []
    for noeud1 in sous_graphe1.nodes() :
        ok = False
        for noeud2 in sous_graphe2.nodes() :
            if noeud1 == noeud2 and noeud2 not in noeuds_vus :
                noeuds_vus.append(noeud2)
                ok = True
        if not ok :
            return False
    if len(noeuds_vus) != sous_graphe2.number_of_nodes() :
        return False
    
    aretes_vus = []
    for arete1, data1 in sous_graphe1.edges(data=True) :
        ok = False
        for arete2, data2 in sous_graphe2.edges(data=True) :
            if arete1 == arete2 and data1["label"] == data2["label"] and arete2 not in aretes_vus :
                aretes_vus.append(arete2)
                ok = True
        if not ok :
            return False
    if len(aretes_vus) != sous_graphe2.number_of_edges() :
        return False
    
    return True

''' algo principal '''
def heuristique(graphe1, graphe2):
    liste_sous_graphe_meilleur = []
    noeuds1_vus = []
    sous_graphe_commun = nx.MultiDiGraph()
    sous_graphe_commun1 = nx.MultiDiGraph()
    sous_graphe_commun2 = nx.MultiDiGraph()
    
    ancien_nb_noeuds = -1
    
    
    chaines_1 = recup_chaines(graphe1)
    chaines_2 = recup_chaines(graphe2)

    print(chaines_1)
    print(chaines_2)
    
    dico_sous_graphe = {}
    liste_noeuds_vus = []
    
    dico_chaines_1 = liste_meme_type(chaines_1, graphe1, chaines_2, graphe2)
#     print(dico_chaines_1)
    while ancien_nb_noeuds < sous_graphe_commun.number_of_nodes() :
#         print("ramou")
        dico_sous_graphe.clear()
        del(liste_noeuds_vus[:])
        compteur = 0
        for ch in chaines_1 :
            for noeud1 in ch :
                if noeud1 not in sous_graphe_commun1.nodes() and noeud1 not in [2,3,4] :
                    if noeud1 not in dico_sous_graphe.keys() :
                        dico_sous_graphe.update({noeud1 : []})
                    for noeud2 in dico_chaines_1[compteur][noeud1] :
#                         print(noeud1)
#                         print(noeud2)
#                         print(sous_graphe_commun1.nodes.data())
#                         print(sous_graphe_commun2.nodes.data())
#                         print(sous_graphe_commun.nodes.data())
                        if noeud2 not in sous_graphe_commun2.nodes()and noeud2 not in [2,3,4] :
                            
                            #if (noeud1, noeud2) not in liste_noeuds_vus and meme_type_meme_chaine(noeud1, noeud2, graphe1, graphe2) : 
                            if meme_type_meme_chaine(noeud1, noeud2, graphe1, graphe2) :
#                                 print(noeud1)
#                                 print(noeud2)
    #                             print(noeud1, noeud2)
#                                 print("tout petit rat")
#                                 print(sous_graphe_commun.nodes.data())
                                sous_graphe, sous_graphe1, sous_graphe2, liste_noeuds = extension_sous_graphe(graphe1, graphe2, noeud1, noeud2, sous_graphe_commun, sous_graphe_commun1, sous_graphe_commun2)
#                                 print(sous_graphe.nodes())
                                liste_noeuds_vus.extend(liste_noeuds)
              
                                #print(sous_graphe.edges.data())
                                
                                dico_sous_graphe[noeud1].append((sous_graphe, sous_graphe1, sous_graphe2))
#                                 if noeud1 == 9 :
#                                     print("petit rat")
#                                     print(sous_graphe.nodes.data())
#                                     print(dico_sous_graphe[noeud1])
#                                     print(dico_sous_graphe[noeud1][0][0].nodes.data())
                                    
                                #sim = calcul_sim_aretes_avec_coeff_heuristique(graphe1, graphe2, sous_graphe, "petit rat", 1, 1, 1)
#                             print(sim)
        #                 trouve = False
        #                 for elt in dico_sous_graphe[noeud1].keys() :
        #                     if sous_graphe_identique(sous_graphe1, elt) :
        #                         dico_sous_graphe[noeud1][elt].append((sous_graphe2, sous_graphe))
        #                         trouve = True
        #                 if not trouve :
        #                     dico_sous_graphe[noeud1].update({sous_graphe1 : (sous_graphe2, sous_graphe)})
        #print(dico_sous_graphe)
            compteur += 1
        sim_max = 0
    #     for noeud in dico_sous_graphe.keys():
    #         for sous_graphe_par_noeud in dico_sous_graphe[noeud] :
        sous_graphe_meilleur = (-1, nx.MultiDiGraph(), nx.MultiDiGraph(), nx.MultiDiGraph())
        for noeud in dico_sous_graphe.keys() :
#             print(noeud)
#             print(dico_sous_graphe[noeud])
            for sous_graphe_par_noeud in dico_sous_graphe[noeud] :
#                 print(sous_graphe_par_noeud)
                print(sous_graphe_par_noeud[0].nodes.data())
                #print(sous_graphe_par_noeud[0].nodes.data())
                sim = calcul_sim_aretes_avec_coeff_heuristique(graphe1, graphe2, sous_graphe_par_noeud[0], "petit rat", 1, 1, 1)
#                 print(sim)
    #             if sim == sim_max :
    #                 liste_sous_graphe_meilleur.append(sous_graphe_par_noeud)
    #             elif sim > sim_max :
    #                 del(liste_sous_graphe_meilleur[:])
    #                 sim_max = sim
    #                 liste_sous_graphe_meilleur.append(sous_graphe_par_noeud)    
                
                
                if sim > sim_max :
#                     compa = True
#                     for noeud in sous_graphe_par_noeud[0].nodes() :
#                         compa = test_compatibilite(sous_graphe_commun, noeud, graphe1, graphe2)
#                         if not compa : break
#                     if compa :
                    sous_graphe_meilleur = (noeud, sous_graphe_par_noeud[0].copy(), sous_graphe_par_noeud[1].copy(), sous_graphe_par_noeud[2].copy())
                    sim_max = sim
                elif sim == sim_max and sous_graphe_par_noeud[0].number_of_edges() > sous_graphe_meilleur[1].number_of_edges() :
                    sous_graphe_meilleur = (noeud, sous_graphe_par_noeud[0].copy(), sous_graphe_par_noeud[1].copy(), sous_graphe_par_noeud[2].copy())
                    sim_max = sim
                #print("petit rat")
#         print("gros rat debile")
#         print(sim_max)
#         print(sous_graphe_meilleur[1].nodes())
#         print(sous_graphe_meilleur[1].nodes.data())
#         print(sous_graphe_meilleur[2].nodes.data())
        ancien_nb_noeuds = sous_graphe_commun.number_of_nodes()
#         print(sous_graphe_commun.nodes.data())
#         print(sous_graphe_meilleur[1].nodes.data())
        
        sous_graphe_commun = nx.union(sous_graphe_commun, sous_graphe_meilleur[1])
        sous_graphe_commun1 = nx.union(sous_graphe_commun1, sous_graphe_meilleur[2])
        sous_graphe_commun2 = nx.union(sous_graphe_commun2, sous_graphe_meilleur[3])
#         print("encore plus petit rat")
#         print(sous_graphe_commun.nodes.data())
#         print(sous_graphe_commun.edges.data())
        #print(sous_graphe_commun.nodes.data())
#         print(sous_graphe_meilleur[1].nodes.data())
#         print(sous_graphe_commun.nodes.data())
#         print(sous_graphe_commun.edges.data())
#         print(sous_graphe_commun1.nodes.data())
#         print(sous_graphe_commun2.nodes.data())
        #exit(0)
        
    return sous_graphe_commun, sous_graphe_commun1, sous_graphe_commun2
        
        #liste_sous_graphe_meilleur.append(sous_graphe_meilleur)
    


''' new heuristique 24/01/20 '''

''' renvoie le nombre de voisins non covalents du noeud passe en parametre '''
def nb_voisins_non_cov(graphe, noeud):
    nombre = 0
    for voisin in graphe[noeud] :
        for edge in graphe[noeud][voisin] :
            if graphe[noeud][voisin][edge]["label"] != 'B53' :
                nombre += graphe.nodes[noeud]["poids"]
                
    return nombre

''' ajoute a tous les noeuds un attribut permettant de stocker le nombre de voisins non covalents (nb_vnc) '''
def ajout_nb_voisins_non_cov(graphe):
    nx.set_node_attributes(graphe, -1, 'nb_vnc')
    for noeud, data in graphe.nodes(data=True):
        if data["type"] != None and data["type"] != -1 :
            graphe.nodes[noeud]["nb_vnc"] = nb_voisins_non_cov(graphe, noeud)

''' renvoie une liste de listes contenant les numeros des chaines qui ont des liens les unes avec les autres
(si une chaine n'a de lien avec aucune autre chaine, la liste correspondante ne contiendra que le numero de la chaine)  '''
def chaines_ensemble(graphe):   
    
    chaine_reliee = chaines_reliees(graphe)
    

    petit_dico = {}
    for i in range(1,5) :
        for elt in chaine_reliee:
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
    
                
    
    return chaines_a_mettre_ensemble

''' ajoute un attribut a tous les noeuds du graphe permettant de stocker le numero de la composante a laquelle appartient le noeud 
(chaine ou groupe de chaines s'il y a des liens entre chaines)
chaines : une liste de listes stockant les numeros de noeuds de chacune des chaines (dans la liste d'indice 0, il y aura les noeuds de la chaine 1, etc...)
chaines_a_mettre_ensemble : une liste de listes contenant les liens entre chaines : ex : [[1,3], [2], [4]] s'il y a un lien entre la ch 1 et la ch 3 et si les autres sont independantes
(issue de chaines_ensemble)  '''
def reetiquetage(graphe, chaines, chaines_a_mettre_ensemble):
    composantes = {}
    
    nx.set_node_attributes(graphe, -1, 'num_composante')
    for ch in chaines_a_mettre_ensemble :
        if len(ch) > 1 :
            liste_chaines = []
            for elt in ch :
                liste_chaines.extend(chaines[elt-1])
            for noeud in liste_chaines :
                graphe.nodes[noeud]["num_composante"] = list(ch)
                for voisin in graphe[noeud] :
                    if graphe.nodes[voisin]["num_composante"] == -1 :
                        graphe.nodes[voisin]["num_composante"] = list(ch)
            composantes.update({tuple(ch) :liste_chaines})
        else :
            for noeud in chaines[ch[0]-1] : 
                graphe.nodes[noeud]["num_composante"] = list(ch)
                for voisin in graphe[noeud] :
                    if graphe.nodes[voisin]["num_composante"] == -1 :
                        graphe.nodes[voisin]["num_composante"] = list(ch)
            composantes.update({tuple(ch) :chaines[ch[0]-1]})      
    return composantes

def maximum_voisin_non_cov(liste_noeuds, graphe):
    maxi = -1
    noeud_maxi = -1
    for noeud in liste_noeuds :
        if graphe.nodes[noeud]["marque"] == 0 :
            if graphe.nodes[noeud]["nb_vnc"] > maxi :
                maxi = graphe.nodes[noeud]["nb_vnc"]
                noeud_maxi = noeud
            elif graphe.nodes[noeud]["nb_vnc"] == maxi :
#                 print("houhou")
#                 print(noeud)
#                 print(graphe.nodes[noeud]["type"])
#                 print(graphe.nodes[noeud_maxi]["type"])
                if (graphe.nodes[noeud]["type"] in [0,1] and graphe.nodes[noeud_maxi]["type"] in [2,3]) or (abs(graphe.nodes[noeud]["position"][0] - graphe.nodes[graphe.nodes[noeud]["chaine"][0]]["position"][0]) < abs(graphe.nodes[noeud_maxi]["position"][0] - graphe.nodes[graphe.nodes[noeud_maxi]["chaine"][0]]["position"][0])) : 
                #if abs(graphe.nodes[noeud]["position"][0] - graphe.nodes[graphe.nodes[noeud]["chaine"][0]]["position"][0]) < abs(graphe.nodes[noeud_maxi]["position"][0] - graphe.nodes[graphe.nodes[noeud_maxi]["chaine"][0]]["position"][0]) : 
                #      print("oh la la")
                    noeud_maxi = noeud 
    return noeud_maxi 


def test_compatibilite_new(noeud1, noeud2, graphe1, graphe2, liste_superposes):
    ''' renvoie vrai si ajouter le noeud au graphe_commun ne provoquera pas l'apparition d'une incompatibilité 
    de séquence dans le graphe_commun (deux superpositions de noeuds qui sont dans un ordre dans le 1er graphe
    et dans l'autre ordre dans le 2e graphe)'''

    if graphe1.nodes[noeud1]["type"] not in [None, -1] : 
        for noeud in liste_superposes : 
                meme_chaine = False
                    
                modele = graphe1.nodes[1]["position"][0]
                    
                for elt in graphe1.nodes[noeud1]["chaine"] :
                        if elt in graphe1.nodes[noeud[0]]["chaine"] :
                            meme_chaine = True
                            
                for elt in graphe2.nodes[noeud2]["chaine"] :
                        if elt in graphe2.nodes[noeud[1]]["chaine"] :
                            meme_chaine = True
        #             if noeud == (24,27) :
        #                 print(noeud)
        #                 print(graphe1.nodes[noeud[0]]["type"])
        #                 print(graphe1.nodes.data())
        #                 print(graphe2.nodes.data())
                           
        #             if noeud == (1033,1034) :
        #                 print(graphe1.nodes[noeud[0]]["chaine"])
                if (noeud1, noeud2) != noeud and meme_chaine and graphe1.nodes[noeud1]["type"] not in [None,-1] and graphe1.nodes[noeud[0]]["type"] not in [None,-1] and noeud1 not in [1,2,3,4,5] and noeud[0] not in [1,2,3,4,5]:
        #                 print("test compa")
        #                 print(graphe1.nodes[noeud[0]]["position"])
        #                 print(graphe1.nodes[noeud2[0]]["position"])
        #                 
        #                 print(graphe1.nodes[noeud[0]]["position"])
        #                 print(graphe1.nodes[noeud2[0]]["position"])
                        if (graphe1.nodes[noeud1]["position"][0] < graphe1.nodes[noeud[0]]["position"][0] and graphe2.nodes[noeud2]["position"][0] > graphe2.nodes[noeud[1]]["position"][0]) or (graphe1.nodes[noeud1]["position"][0] > graphe1.nodes[noeud[0]]["position"][0] and graphe2.nodes[noeud2]["position"][0] < graphe2.nodes[noeud[1]]["position"][0]) :
        #                     print("test compa")
        #                     print(noeud)
        #                     print(chaine_noeud11)
        #                     print(chaine_noeud21)
        #                     print(distance_entre_noeuds_meme_chaine(graphe1, chaine_noeud11, noeud[0]))
        #                     print(distance_entre_noeuds_meme_chaine(graphe1, chaine_noeud21, noeud2[0]))
                            #if (distance_entre_noeuds_meme_chaine(graphe1, chaine_noeud11, noeud[0]) < distance_entre_noeuds_meme_chaine(graphe1, chaine_noeud21, noeud2[0]) and distance_entre_noeuds_meme_chaine(graphe2, chaine_noeud12, noeud[1]) > distance_entre_noeuds_meme_chaine(graphe2, chaine_noeud22, noeud2[1])) or  (distance_entre_noeuds_meme_chaine(graphe1, chaine_noeud11, noeud[0]) > distance_entre_noeuds_meme_chaine(graphe1, chaine_noeud21, noeud2[0]) and distance_entre_noeuds_meme_chaine(graphe2, chaine_noeud12, noeud[1]) < distance_entre_noeuds_meme_chaine(graphe2, chaine_noeud22, noeud2[1])) :
                            if (abs(noeud1-modele) < abs(noeud[0]-modele) and abs(noeud2 - modele) > abs(noeud[1]-modele)) or (abs(noeud1-modele) > abs(noeud[0]-modele) and abs(noeud2 - modele) < abs(noeud[1]-modele)) :
                                print(noeud)
                                return False
                            else : 
                                print("petit rat")
    
    return True

def dans_liste(liste_superposes, noeud1, noeud2):
    for elt in liste_superposes :
        if (noeud1 == elt[0] and noeud2 != elt[1]) or (noeud2 == elt[1] and noeud1 != elt[0]) :
            return True
    return False

def recup_sommets_possibles(noeud1, graphe1, dico_chaines_1, comp_2):
    sommets_possibles = []
    #print(comp_2)
    for elt in dico_chaines_1[graphe1.nodes[noeud1]["chaine"][0]-1][noeud1] :
        #print(elt)
        if elt in comp_2 :
            sommets_possibles.append(elt)
            
    return sommets_possibles
            

def new_heuristique(graphe1_deb, graphe2_deb):
    
    ### Initialisation ###
    
    ''' ajout de l'attribut nb_vnc aux 2 graphes '''
    ajout_nb_voisins_non_cov(graphe1_deb)
    ajout_nb_voisins_non_cov(graphe2_deb)
    
#     print(graphe1_deb.nodes.data())
#     print(graphe2_deb.nodes.data())
    
    ''' on recupere les noeuds de chaque chaine qu'on stocke dans des listes de listes '''
    chaines_1 = recup_chaines(graphe1_deb)
    chaines_2 = recup_chaines(graphe2_deb)
    
    ''' on enleve le premier element de chaque chaine car c'est un sommet du motif '''
    for chaine in chaines_1 :
        chaine.remove(chaine[0])
        
    for chaine in chaines_2 :
        chaine.remove(chaine[0])
    
    ''' on enleve les sommets du motif des graphes '''
    
    nx.set_node_attributes(graphe1_deb, 0, 'marque')
    nx.set_node_attributes(graphe2_deb, 0, 'marque')
    nx.set_edge_attributes(graphe1_deb, 0, 'marque')
    nx.set_edge_attributes(graphe2_deb, 0, 'marque')
    
    graphe1 = graphe1_deb.copy()
    graphe2 = graphe2_deb.copy()
        
    graphe1.remove_nodes_from([1,2,3,4,5])
    graphe2.remove_nodes_from([1,2,3,4,5])
    
    print(chaines_1)
    print(chaines_2)
    
    ''' on recupere la liste des noeuds du 2e graphe pouvant etre superposes avec chaque noeud du 1er graphe et inversement (stocke dans une liste de dico) '''
    dico_chaines_1 = liste_meme_type(chaines_1, graphe1, chaines_2, graphe2)
    dico_chaines_2 = liste_meme_type(chaines_2, graphe2, chaines_1, graphe1)
    
    print(dico_chaines_1)
    print(dico_chaines_2)
    
    ''' on recupere les liens entre chaines '''
    chaines_a_mettre_ensemble_1 = chaines_ensemble(graphe1)
    chaines_a_mettre_ensemble_2 = chaines_ensemble(graphe2)
    
    print(chaines_a_mettre_ensemble_1)
    print(chaines_a_mettre_ensemble_2)
    
    ''' on etiquette les graphes avec les numeros des composantes et on renvoie les noeuds contenus dans chaque composante (dico de listes) '''
    composantes_1 = reetiquetage(graphe1, chaines_1, chaines_a_mettre_ensemble_1)
    composantes_2 = reetiquetage(graphe2, chaines_2, chaines_a_mettre_ensemble_2)
    
#     for noeud, data in graphe1.nodes(data=True) : 
#         print(noeud)
#         print(data["num_composante"])
#     for noeud, data in graphe2.nodes(data=True) : 
#         print(noeud)
#         print(data["num_composante"])
    print(composantes_1)
    print(composantes_2)
    
    ### ###

    liste_superposes = set() ## pour garder en memoire les paires de noeuds qu'on veut mettre ensemble (pour faire les tests d'incompatibilite)
    indice_1 = 0
    indice_2 = 0
    
    ''' on traite composante par composante '''
    while indice_1 < len(composantes_1) :
        #print(composantes_1)
        comp1 = list(composantes_1.keys())[indice_1] ## nom de la composante courante du graphe1
        indice_2 = 0
        nb_aretes_max = 0 
        sommet_max_2 = -1     
        sommet_max_1 = -1
        comp2_max = -1
        liste_superposes_a_garder = [] 
        liste_superposes_a_garder_ar = []
        
        noeud_maxi = maximum_voisin_non_cov(composantes_1[comp1], graphe1_deb)  ## recherche du noeud qui a le plus grand nombre de voisins non covalents dans comp1
        if noeud_maxi == -1 and len(composantes_1[comp1]) > 0 : ## cas ou on ne trouve pas de noeud_maxi alors que la composante n'est pas vide
            print("pas de noeud max !!")
            exit()

        if noeud_maxi != -1 :
            print("raaaapaaaalaaaa")
            print(noeud_maxi)
            #sommets_possibles = dico_chaines_1[graphe1.nodes[noeud_maxi]["chaine"][0]-1][noeud_maxi] ## on recupere l'ensemble des noeuds du graphe2 compatibles avec noeud_maxi
            
            ''' on teste chaque composante du graphe 2 à la recherche d'une compatibilite avec comp1 du graphe1 '''
            while indice_2 < len(composantes_2) :
                
                comp2 = list(composantes_2.keys())[indice_2]
                print("composantes :")
                print(comp1)
                print(comp2)
                
                ''' comp2 est compatible avec comp1 s'il y a au moins un numero identique dans le nom des deux composantes '''
                ok = False
                for elt in comp1 :
                    if elt in comp2 :
                        ok = True
                for elt in comp2 :
                    if elt in comp1 :
                        ok = True
                        
                print(ok)
                if ok :
                    #print("gros tas")
                    sommets_possibles = recup_sommets_possibles(noeud_maxi, graphe1, dico_chaines_1, composantes_2[comp2]) ## on recupere l'ensemble des noeuds du graphe2 compatibles avec noeud_maxi et appartenant a comp2
            
                    #print(sommets_possibles)
                    ''' on teste chaque sommet du graphe2 de comp2 compatibles avec noeud_maxi pour trouver le sous-graphe commun le plus grand en étendant aux voisins '''
                    for sommet in sommets_possibles :
                        
                        ## si le sommet n'est pas deja dans le sous-graphe commun (marque) et si la paire est compatible avec le reste du sous-graphe commun, on peut essayer d'etendre aux voisins
                        if graphe2.nodes[sommet]["marque"] == 0 and test_compatibilite_new(noeud_maxi, sommet, graphe1_deb, graphe2_deb, liste_superposes) :
                            ## sommet compatible
                            
                            liste_superposes_temp = [(noeud_maxi, sommet)] ## les paires de sommets qu'on voudrait ajouter dans le sous-graphe commun en partant de cette premiere paire (sert a traiter les voisins les uns apres les autres)
                            liste_superposes_ar_temp = set() ## les paires d'aretes qu'on voudrait ajouter au sous-graphe commun en partant de cette paire de sommets
                            liste_deja_teste = []
                            nb_aretes = 0
                            compteur = 0
                            
                            while compteur < len(liste_superposes_temp) :
#                                 print(compteur)

                                ''' on cherche a etendre le sous-graphe commun temporaire aux voisins de la paire de noeuds courante '''
                                noeud1 = liste_superposes_temp[compteur][0]
                                noeud2 = liste_superposes_temp[compteur][1]
                                
                                ## pas besoin de chercher les voisins pour les noeuds de type artificiel, ils n'en ont qu'un normalement
                                if graphe1.nodes[noeud1]["type"] != -1 :
                                
                                    voisins_1 = set(graphe1[noeud1])
                                    voisins_1.update(graphe1.predecessors(noeud1))
                                    voisins_2 = set(graphe2[noeud2])
                                    voisins_2.update(graphe2.predecessors(noeud2))
                                    
    #                                 print(voisins_1)
    #                                 print(voisins_2)
                                    nb_voisins_trouves = 0
                                    for voisin_1 in voisins_1 :
    #                                     print("voisin 1")
    #                                     print(voisin_1)
                                        for voisin_2 in voisins_2 :
    #                                         print("voisin 2")
    #                                         print(voisin_2)
        #                                     print(test_compatibilite_new(voisin_1, voisin_2, graphe1_deb, graphe2_deb, liste_superposes))
        #                                     print(test_compatibilite_new(voisin_1, voisin_2, graphe1_deb, graphe2_deb, liste_superposes_temp))
        #                                     print(meme_type_meme_chaine(voisin_1, voisin_2, graphe1, graphe2))
        #                                     print(graphe1.nodes[voisin_1]["marque"])
        #                                     print(graphe2.nodes[voisin_2]["marque"])
        #                                     print(liste_superposes_temp)
        #                                     print(dans_liste(liste_superposes_temp, voisin_1,voisin_2))
        
                                            ''' les voisins de la paire courante doivent :
                                            - ne pas appartenir deja au sous-graphe commun (marque) en paire avec un autre noeud
                                            - ne pas avoir deja ete traites comme paire courante (s'il y a des boucles dans un des graphes)
                                            - etre de meme type et appartenir a la meme chaine 
                                            - ne pas appartenir deja au sous-graphe commun temporaire en paire avec un autre noeud
                                            - ne pas introduire d'incompatibilite de sequence ni dans le sous-graphe commun temporaire, ni dans le sous-graphe commun final
                                            et de plus, il faut que les aretes entre le noeud et son voisin soient de même type
                                            '''
                                            if graphe1.nodes[voisin_1]["marque"] in [0, voisin_2] and graphe2.nodes[voisin_2]["marque"] in [0, voisin_1] and (voisin_1, voisin_2) not in liste_superposes_temp[:compteur] and (voisin_1, voisin_2) not in liste_deja_teste and meme_type_meme_chaine(voisin_1, voisin_2, graphe1, graphe2) and not dans_liste(liste_superposes_temp, voisin_1, voisin_2) and test_compatibilite_new(voisin_1, voisin_2, graphe1_deb, graphe2_deb, liste_superposes) and test_compatibilite_new(voisin_1, voisin_2, graphe1_deb, graphe2_deb, liste_superposes_temp) :
    #                                             print("tadaa")
                                                ## cas des aretes ou des arcs orientes dans le sens noeud vers voisin
                                                if (noeud1, voisin_1) in graphe1.edges() and (noeud2, voisin_2) in graphe2.edges() :
                                                    for edge_1 in graphe1[noeud1][voisin_1] :
                                                        for edge_2 in graphe2[noeud2][voisin_2] :
                                                                label1 = graphe1[noeud1][voisin_1][edge_1]["label"]
                                                                label2 = graphe2[noeud2][voisin_2][edge_2]["label"]
                                                                
                                                                if label1 == label2 :
                                                                    
    #                                                                 print(noeud1, voisin_1)
    #                                                                 print(noeud2, voisin_2)
                                                                    ## la paire (voisin_1, voisin_2) peut deja exister dans le sous-graphe commun final ou temporaire, dans ce cas, pas besoin de le rajouter (mais on peut avoir des aretes a rajouter)
                                                                    if (voisin_1, voisin_2) not in liste_superposes_temp and graphe1.nodes[voisin_1]["marque"] != voisin_2 and graphe2.nodes[voisin_2]["marque"] != voisin_1 :
                                                                        liste_superposes_temp.append((voisin_1, voisin_2))
                                                                    
                                                                    if label1 == 'B53' : ## si c'est une liaison B53, on ne l'ajoute que dans un sens dans la liste des aretes
                                                                        liste_superposes_ar_temp.add(((noeud1, voisin_1, edge_1), (noeud2, voisin_2, edge_2)))
                                                                    else : ## si ce n'est pas une liaison B53, il faut ajouter l'arete dans les deux sens
    #                                                                     print("roupoulou")
                                                                        nb_voisins_trouves += 1
                                                                        nb_aretes += min(graphe1.nodes[noeud1]["poids"], graphe2.nodes[noeud2]["poids"])
                                                                        liste_superposes_ar_temp.add(((noeud1, voisin_1, edge_1), (noeud2, voisin_2, edge_2)))
                                                                        ## mais l'arete dans l'autre sens ne porte peut-etre pas le meme numero (s'il y a plusieurs aretes entre les deux memes sommets)
                                                                        for edge_1_test in graphe1[voisin_1][noeud1] :
                                                                            for edge_2_test in graphe2[voisin_2][noeud2] :
                                                                                if graphe1[voisin_1][noeud1][edge_1_test]["label"] == graphe2[voisin_2][noeud2][edge_2_test]["label"] :
                                                                                    if graphe1[voisin_1][noeud1][edge_1_test]["label"] == label1 or (len(label1) == 3 and graphe1[voisin_1][noeud1][edge_1_test]["label"] == label1[0]+label1[2]+label1[1]) :  
                                                                                            liste_superposes_ar_temp.add(((voisin_1, noeud1, edge_1_test), (voisin_2, noeud2, edge_2_test)))                              
                                                ## cas des arcs dans le sens voisin vers noeud                            
                                                elif (voisin_1, noeud1) in graphe1.edges() and (voisin_2, noeud2) in graphe2.edges():
                                                    for edge_1 in graphe1[voisin_1][noeud1] :
                                                        for edge_2 in graphe2[voisin_2][noeud2] :
                                                                label1 = graphe1[voisin_1][noeud1][edge_1]["label"]
                                                                label2 = graphe2[voisin_2][noeud2][edge_2]["label"]
                                                                
                                                                ## traitement idem qu'au-dessus
                                                                if label1 == label2 :
                                                                    if label1 != 'B53' :
                                                                        print("bizarre, pas b53...")
                                                                        exit()
                                                                        
                                                                    if (voisin_1, voisin_2) not in liste_superposes_temp and graphe1.nodes[voisin_1]["marque"] != voisin_2 and graphe2.nodes[voisin_2]["marque"] != voisin_1 :
                                                                        liste_superposes_temp.append((voisin_1, voisin_2))
            
                                                                    liste_superposes_ar_temp.add(((voisin_1, noeud1, edge_1), (voisin_2, noeud2, edge_2)))
                                    if nb_voisins_trouves == 0 and graphe1.nodes[noeud1]["type"] in [2,3] :
#                                         print("roupoulou")
#                                         print((noeud1, noeud2))
#                                         print(liste_superposes_temp)
#                                         print(liste_superposes_ar_temp)
                                        liste_deja_teste.append((noeud1, noeud2))
                                        liste_superposes_temp.remove((noeud1, noeud2))
                                        a_enlever = []
                                        for arete in liste_superposes_ar_temp :
                                            if (arete[0][0] == noeud1 or arete[0][1] == noeud1) and (arete[1][0] == noeud2 or arete[1][1] == noeud2) :
                                                a_enlever.append(arete)
                                        for elt in a_enlever :
                                            liste_superposes_ar_temp.remove(elt)
                                            
#                                         print(liste_superposes_temp)
#                                         print(liste_superposes_ar_temp)
                                        compteur -= 1  
                                    #exit()
                                compteur += 1    
#                             print("et la")
#                             print(nb_aretes) 

                            ''' on recherche le sous-graphe commun temporaire qui a le plus grand nombre d'aretes parmi les sous-graphes communs issus des differentes premieres paires de sommets (noeud_maxi, sommet) dans toutes les composantes du graphe2 compatibles avec comp1
                            , et c'est celui-la qu'on garde en memoire  '''    
                            if nb_aretes > nb_aretes_max :
#                                 print("tu rentres la")
                                sommet_max_2 = sommet
                                sommet_max_1 = noeud_maxi
                                comp2_max = comp2
#                                 print(sommet_max_1)
#                                 print(sommet_max_2)
#                                 print(liste_superposes_temp)
#                                 print(liste_superposes_ar_temp)
                                nb_aretes_max = nb_aretes
                                liste_superposes_a_garder = list(liste_superposes_temp)
                                liste_superposes_a_garder_ar = list(liste_superposes_ar_temp)  
                            
                indice_2 += 1
#             print(composantes_1)
#             print(composantes_2)        
            
            
            if sommet_max_2 != -1 and sommet_max_1 != -1 : 
                        ''' on a trouve un morceau de graphe compatible dans comp1 en partant de noeud_maxi'''
                        print("gros tas")
                        print(comp2_max)
                        
                        ## on marque pour de bon les sommets et les aretes du sous-graphe commun temporaire choisi, comme faisant partie du sous-graphe commun max
                        ## et on les exclut de la suite de l'analyse en modifiant leur attribut num_composante : ils n'appartiennent plus à aucune composante
                        for paire in liste_superposes_a_garder :
                            liste_superposes.add(paire)
                            graphe1.nodes[paire[0]]["marque"] = paire[1]
                            graphe2.nodes[paire[1]]["marque"] = paire[0]
                            
                            graphe1.nodes[paire[0]]["num_composante"] = -1
                            graphe2.nodes[paire[1]]["num_composante"] = -1
                            
                            
                        
                        for paire in liste_superposes_a_garder_ar :
#                             print(paire)
#                             print(graphe1.edges.data())
#                             print(graphe2.edges.data())
                            graphe1[paire[0][0]][paire[0][1]][paire[0][2]]["marque"] = paire[1]
                            graphe2[paire[1][0]][paire[1][1]][paire[1][2]]["marque"] = paire[0]
                            
                            
                        ## on veut determiner si comp1 et comp2_max sont composees d'une seule chaine ou non (car on ne les traitera pas de la meme facon)
                        
                        nb_nb_1 = 0
                        for elt in comp1 :
                            if isinstance(elt, int) :
                                nb_nb_1 += 1
                        
                        nb_nb_2 = 0
                        for elt in comp2_max :
                            if isinstance(elt, int) :
                                nb_nb_2 += 1
                               
                                
                        
                        if nb_nb_1 == 1 and nb_nb_2 == 1 :
                            ## cas ou comp1 et comp2_max sont constituees d'une seule chaine
                            
                            ## dans ce cas : on va pouvoir separer en 2 comp1 et comp2_max (en supprimant les sommets ajoutes au graphe commun)
                            ## et associer la premiere partie de comp1 avec la premiere partie de comp2_max et la deuxieme partie de comp1 avec la deuxieme partie de comp2_max 
                            ## pour cela, on ajoute au nom de chaque chaine de la composante a ou b 
                            #  print("eh oh")
                            new_groupes_1 = []
                            compteur_non_traite = -1
                            compteur_traite = -1
                            for i in range(len(composantes_1[comp1])) :
                                if graphe1.nodes[composantes_1[comp1][i]]["num_composante"] == -1 :
                                    new_groupes_1.append(composantes_1[comp1][compteur_traite+1:compteur_non_traite+1])
                                    compteur_traite = i
                                elif  graphe1.nodes[composantes_1[comp1][i]]["num_composante"] != -1 :
                                    compteur_non_traite = i
                            if compteur_traite < len(composantes_1[comp1]) - 1 :
                                new_groupes_1.append(composantes_1[comp1][compteur_traite+1:compteur_non_traite+1])
#                             print(new_groupes_1)
                            
                        
                            new_groupes_2 = []
                            compteur_non_traite = -1
                            compteur_traite = -1
                            for i in range(len(composantes_2[comp2_max])) :
                                if graphe2.nodes[composantes_2[comp2_max][i]]["num_composante"] == -1 :
                                    new_groupes_2.append(composantes_2[comp2_max][compteur_traite+1:compteur_non_traite+1])
                                    compteur_traite = i
                                elif  graphe2.nodes[composantes_2[comp2_max][i]]["num_composante"] != -1 :
                                    compteur_non_traite = i
                            if compteur_traite < len(composantes_2[comp2_max]) - 1 :
                                new_groupes_2.append(composantes_2[comp2_max][compteur_traite+1:compteur_non_traite+1])
    
#                             print(new_groupes_2)
                            
                            for i in range(min(len(new_groupes_1), len(new_groupes_2))) :
                            
                                if len(new_groupes_1[i]) != 0 and len(new_groupes_2[i]) != 0  :
                                    
                                    for elt in new_groupes_1[i] :
                                        num = chr(ord('a') + i)
                                        for j in range(len(graphe1.nodes[elt]["num_composante"])) :
                                            graphe1.nodes[elt]["num_composante"][j] = str(graphe1.nodes[elt]["num_composante"][j])+ num
                                    composantes_1.update({tuple(graphe1.nodes[elt]["num_composante"]) : list(new_groupes_1[i])})
     
                                    for elt in new_groupes_2[i] :
                                        num = chr(ord('a') + i)
                                        for j in range(len(graphe2.nodes[elt]["num_composante"])) :
                                            graphe2.nodes[elt]["num_composante"][j] = str(graphe2.nodes[elt]["num_composante"][j])+ num 
                                    composantes_2.update({tuple(graphe2.nodes[elt]["num_composante"]) : list(new_groupes_2[i])})
                                    
                                    #del(composantes_1[comp1])
                            ## on supprime la composante comp2_max mais pas comp1 car le traitement de comp1 est termine, on n'y reviendra plus contrairement a comp2_max
                            del(composantes_2[comp2_max])
                                    
    
                                    
    
                        else :
                            
                            if nb_nb_1 > 1 :
                                ## cas ou comp1 est constituee de plusieurs chaines
                                ## dans ce cas, on recherche les nouvelles composantes connexes dans comp1 en ayant supprime les sommets de comp1 qui ont ete ajoutes au sous-graphe commun
                                ## et on n'interdit pas de comparaison avec les composantes du graphe2 car les informations d'incompatibilite sont plus dures a determiner
                                ## pour cela, on ajoute a ou b ou c etc.. au nom de la nouvelle composante (mais pas à l'interieur du nom de chaque chaine dans la composante, comme ca, elles seront comparables)
                                new_groupes_1 = []
                                for i in range(len(composantes_1[comp1])) :
                                    if graphe1.nodes[composantes_1[comp1][i]]["num_composante"] != -1 : 
                                        new_groupes_1.append(composantes_1[comp1][i])
                                        
                                composantes = recherche_composante_connexe(graphe1.subgraph(new_groupes_1))
                                
                                compte = 0
                                for compo in composantes :
                                    for elt in compo :
                                        graphe1.nodes[elt]["num_composante"].append(chr(ord('a') + compte))
                                    compte += 1
                                    composantes_1.update({tuple(graphe1.nodes[elt]["num_composante"]) : list(compo)})
                            else :
                                ## cas ou comp1 est constituee d'une seule chaine
                                ## idem que le premier cas rencontre pour comp1
                                new_groupes_1 = []
                                compteur_non_traite = -1
                                compteur_traite = -1
                                for i in range(len(composantes_1[comp1])) :
                                    if graphe1.nodes[composantes_1[comp1][i]]["num_composante"] == -1 :
                                        new_groupes_1.append(composantes_1[comp1][compteur_traite+1:compteur_non_traite+1])
                                        compteur_traite = i
                                    elif  graphe1.nodes[composantes_1[comp1][i]]["num_composante"] != -1 :
                                        compteur_non_traite = i
                                if compteur_traite < len(composantes_1[comp1]) - 1 :
                                    new_groupes_1.append(composantes_1[comp1][compteur_traite+1:compteur_non_traite+1])
                                                                        
                                for i in range(len(new_groupes_1)) :
                                    if len(new_groupes_1[i]) >  0 :
                                        for elt in new_groupes_1[i] :
                                            graphe1.nodes[elt]["num_composante"].append(chr(ord('a') + i))
                                        composantes_1.update({tuple(graphe1.nodes[elt]["num_composante"]) : list(new_groupes_1[i])})
                                print("ramousnif")
                                print(composantes_1)
    
                            if nb_nb_2 > 1 :
                                ## cas ou comp2_max est constituee de plusieurs chaines
                                ## idem que pour comp1 dans ce cas
                                new_groupes_2 = []
                                for i in range(len(composantes_2[comp2_max])) :
                                    if graphe2.nodes[composantes_2[comp2_max][i]]["num_composante"] != -1 : 
                                        new_groupes_2.append(composantes_2[comp2_max][i])
                                        
                                composantes = recherche_composante_connexe(graphe2.subgraph(new_groupes_2))
                                
                                compte = 0
                                for compo in composantes :
                                    for elt in compo :
                                        graphe2.nodes[elt]["num_composante"].append(chr(ord('a') + compte))
                                    compte += 1
                                    composantes_2.update({tuple(graphe2.nodes[elt]["num_composante"]) : list(compo)})
                                print("roupoulou")
                                print(composantes_2)
                                del(composantes_2[comp2_max])
                            else :
                                ## cas ou comp2_max est constituee d'une seule chaine
                                ## idem que le premier cas rencontre pour comp2_max
                                new_groupes_2 = []
                                compteur_non_traite = -1
                                compteur_traite = -1
                                for i in range(len(composantes_2[comp2_max])) :
                                    if graphe2.nodes[composantes_2[comp2_max][i]]["num_composante"] == -1 :
                                        new_groupes_1.append(composantes_2[comp2][compteur_traite+1:compteur_non_traite+1])
                                        compteur_traite = i
                                    elif  graphe2.nodes[composantes_2[comp2_max][i]]["num_composante"] != -1 :
                                        
                                        compteur_non_traite = i
                                if compteur_traite < len(composantes_2[comp2_max]) - 1 :
                                    new_groupes_2.append(composantes_2[comp2_max][compteur_traite+1:compteur_non_traite+1])
                                                                        
                                for i in range(len(new_groupes_2)) :
                                    if len(new_groupes_2[i]) > 0 :
                                        for elt in new_groupes_2[i] :
                                            graphe2.nodes[elt]["num_composante"].append(chr(ord('a') + i))
                                        composantes_2.update({tuple(graphe2.nodes[elt]["num_composante"]) : list(new_groupes_2[i])})
                                del(composantes_2[comp2_max])
                            
    
                        
                    
                
            else :
                        ''' on n'a pas trouve de morceau de graphe compatible avec noeud_maxi dans comp1'''
                        graphe1.nodes[noeud_maxi]["marque"] = -1
                        
                        
                        ## idem que dans l'autre cas 
                        ## on veut determiner si comp1 est composee d'une seule chaine ou non (car on ne le traitera pas de la meme facon)

                        nb_nb_1 = 0
                        for elt in comp1 :
                            if isinstance(elt, int) :
                                nb_nb_1 += 1
                        
                        if nb_nb_1 == 1 :
                            ## cas ou comp1 est constitue d'une seule chaine
                            ## dans ce cas, on va creer deux nouvelles composantes : une avec les sommets se trouvant avant le noeud_maxi (sur la sequence) annote 'a'
                            ## et l'autre avec les sommets se trouvant apres le noeud_maxi (sur la sequence) annotee 'b'
                            new_groupes = []
                            for i in range(len(composantes_1[comp1])) :
                                if composantes_1[comp1][i] == noeud_maxi :
                                    new_groupes.append(list(composantes_1[comp1][:i]))
                                    new_groupes.append(list(composantes_1[comp1][i+1:]))
#                             print(new_groupes)        
                                    
                            if len(new_groupes) == 0 :
                                print("bizarre, pas trouve noeud_maxi ?")
                                exit()
                            else : ## len(new_groupes) == 2 en principe
                                if len(new_groupes) != 2 :
                                    print("plus de deux nouvelles composantes, bizarre")
                                    exit()
                                if len(new_groupes[0]) > 0 : 
                                    for elt in new_groupes[0] :
                                        graphe1.nodes[elt]["num_composante"].append('a')
                                    composantes_1.update({tuple(graphe1.nodes[elt]["num_composante"]) : list(new_groupes[0])})
                                if len(new_groupes[1]) > 0 : 
                                    for elt in new_groupes[1] :
                                        graphe1.nodes[elt]["num_composante"].append('b')
                                    composantes_1.update({tuple(graphe1.nodes[elt]["num_composante"]) : list(new_groupes[1])})
                        else :
                            ## cas ou comp1 est constitue de plusieurs chaines
                            ## dans ce cas, on recherche les composantes connexes de comp1 en supprimant le noeud_maxi
                            ## on les annote avec a ou b ou c etc.. en tant que nouvel element de la liste dans le nom de la composante
                            new_groupes_1 = []
                            for i in range(len(composantes_1[comp1])) :
                                    if graphe1.nodes[composantes_1[comp1][i]]["num_composante"] != -1 : 
                                        new_groupes_1.append(composantes_1[comp1][i])
                                        
                            composantes = recherche_composante_connexe(graphe1.subgraph(new_groupes_1))
                                
                            compte = 0
                            for compo in composantes :
                                for elt in compo :
                                    graphe1.nodes[elt]["num_composante"].append(chr(ord('a') + compte))
                                compte += 1
                                composantes_1.update({tuple(graphe1.nodes[elt]["num_composante"]) : list(compo)})
       
                                
                        graphe1.nodes[noeud_maxi]["num_composante"] = -1
#             print(graphe1.nodes.data())
            print(composantes_1)
            print(composantes_2)
        indice_1 += 1

    ''' on construit le sous-graphe commun a partir des attributs marque des noeuds et des aretes des graphes '''        
    sous_graphe_commun = nx.MultiDiGraph()

    for noeud, data in graphe1.nodes(data=True) :
#         print(data)
        if data["marque"] not in [-1,0] :
            sous_graphe_commun.add_node((noeud, data["marque"])) 
    
    for u,v,key,data in graphe1.edges(data=True, keys=True) :
#         print(data)
        if data["marque"] != 0 :
            sous_graphe_commun.add_edge((u, data["marque"][0]), (v, data["marque"][1]), label=data["label"]) 
    
#     for noeud, data in sous_graphe_commun.nodes(data=True) :
#                             print(noeud)
#                             print(data)
#     for u,v, data in sous_graphe_commun.edges(data=True) :
#                             print(u,v)
#                             print(data)
#     print("petit rat")

        
    ''' on ajoute les sommets et les aretes du motif au sous-graphe commun '''
    for i in range(1,6) :
        sous_graphe_commun.add_node((i,i))
        graphe1_deb.nodes[i]["marque"] = i
        graphe2_deb.nodes[i]["marque"] = i
        
    sous_graphe_commun.add_edge((1,1),(2,2), label="CSS")
    sous_graphe_commun.add_edge((2,2),(1,1), label="CSS")
    sous_graphe_commun.add_edge((1,1),(5,5), label="TSS")
    sous_graphe_commun.add_edge((5,5),(1,1), label="TSS")
    sous_graphe_commun.add_edge((2,2),(5,5), label="CWW")
    sous_graphe_commun.add_edge((5,5),(2,2), label="CWW")
    sous_graphe_commun.add_edge((3,3),(4,4), label="CSS")
    sous_graphe_commun.add_edge((4,4),(3,3), label="CSS")
    sous_graphe_commun.add_edge((3,3),(1,1), label="B53")
    sous_graphe_commun.add_edge((2,2),(4,4), label="B53")
    
    
    ''' il reste juste a traiter les sommets etant lie au motif directement par des liaisons non cov (qui sont alors isoles du graphe quand on supprime les sommets du motif, ou non s'ils sont lies aussi a d'autres sommets du graphe) '''
    ## on recherche dans les voisins directs des sommets du motif des correspondances possibles (on pourra tomber sur des liaisons B53 qu'on ajoutera au sous-graphe commun mais ça n'a pas d'importance pour la sim (en fait on ajoutera seulement ceux dans le sens de la sequence selon la facon dont c'est construit ici))
    for i in range(1,6) :
#         print(i)
        for voisin in graphe1_deb[i] :
            for edge in graphe1_deb[i][voisin] :
                #if voisin not in [1,2,3,4,5] or ((i,i),(voisin, voisin)) not in sous_graphe_commun.edges() or graphe1_deb[i][voisin][edge]["label"] != sous_graphe_commun[(i,i)][(voisin, voisin)][0]["label"] :    
                    for voisin_2 in graphe2_deb[i] :
                        for edge_2 in graphe2_deb[i][voisin_2] :
#                             print(liste_superposes)
#                             print(voisin)
#                             print(voisin_2)
#                             print(meme_type_meme_chaine(voisin, voisin_2, graphe1_deb, graphe2_deb))
#                             print(dans_liste(liste_superposes, voisin, voisin_2))
#                             print(test_compatibilite_new(voisin, voisin_2, graphe1_deb, graphe2_deb, liste_superposes))
                            ''' comme pour les autres sommets on recherche des voisins :
                            - de meme type et de meme chaine
                            - qui ne soit pas deja dans une paire avec un autre sommet
                            - qui n'introduisent pas d'incompatibilite de sequence dans le sous-graphe commun si on les met en paire
                            on va ici les ajouter directement au sous-graphe commun
                            '''
                            if meme_type_meme_chaine(voisin, voisin_2, graphe1_deb, graphe2_deb) and not dans_liste(liste_superposes, voisin, voisin_2) and test_compatibilite_new(voisin, voisin_2, graphe1_deb, graphe2_deb, liste_superposes) :
                                label1 = graphe1_deb[i][voisin][edge]["label"]
                                label2 = graphe2_deb[i][voisin_2][edge_2]["label"]
                                if label1 == label2 :
                                    ## la aussi on verifie que la paire (voisin_1, voisin_2) n'ait pas deja ete ajoutee au sous-graphe commun 
                                    ## si c'est le cas, on ne la rajoute pas mais on verifie quand meme les aretes
                                    if graphe1_deb.nodes[voisin]["marque"] in [0, voisin_2] and graphe2_deb.nodes[voisin_2]["marque"] in [0, voisin] :#and (voisin, voisin_2) not in sous_graphe_commun.nodes() :
                                        graphe1_deb.nodes[voisin]["marque"] = voisin_2
                                        graphe2_deb.nodes[voisin_2]["marque"] = voisin
                                        sous_graphe_commun.add_node((voisin, voisin_2))
                                        liste_superposes.add((voisin, voisin_2))
                                    
                                    ## ajout des aretes en commun si elle n'existe pas deja
                                    if ((i, i), (voisin, voisin_2)) not in sous_graphe_commun.edges() :
                                        sous_graphe_commun.add_edge((i, i), (voisin, voisin_2), label=label1)
                                        if label1 != 'B53' :
                                            if len(label1) == 3 :
                                                sous_graphe_commun.add_edge((voisin, voisin_2), (i, i), label=label1[0]+label1[2]+label1[1])
                                            else :
                                                sous_graphe_commun.add_edge((voisin, voisin_2), (i, i), label=label1)
                                    else :
                                        deja_vu = False
                                        for edge_c in sous_graphe_commun[(i, i)][(voisin, voisin_2)] :
                                            if sous_graphe_commun[(i, i)][(voisin, voisin_2)][edge_c]["label"] == label1 :
                                                deja_vu = True
                                        if not deja_vu :
                                            sous_graphe_commun.add_edge((i, i), (voisin, voisin_2), label=label1)
                                            if label1 != 'B53' :
                                                if len(label1) == 3 :
                                                    sous_graphe_commun.add_edge((voisin, voisin_2), (i, i), label=label1[0]+label1[2]+label1[1])
                                                else :
                                                    sous_graphe_commun.add_edge((voisin, voisin_2), (i, i), label=label1)
    
    return sous_graphe_commun
    
if __name__ == '__main__':
    with open("fichier_diff_de89e8b.pickle", 'rb') as fichier_diff :
            mon_depickler = pickle.Unpickler(fichier_diff)
            liste_pas_pareil_1 = mon_depickler.load()
    
    with open("fichier_diff_426d48b.pickle", 'rb') as fichier_diff :
            mon_depickler = pickle.Unpickler(fichier_diff)
            liste_pas_pareil_2 = mon_depickler.load()
            
    
    for elt in liste_pas_pareil_2 :
        if elt not in liste_pas_pareil_1 :
            print(elt)
            
    #exit()
    
    
#     with open("fichier_diff.pickle", 'rb') as fichier_diff :
#             mon_depickler = pickle.Unpickler(fichier_diff)
#             liste_pas_pareil = mon_depickler.load()
#             print(liste_pas_pareil)
#             print(len(liste_pas_pareil))
#             exit()
# #     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_5dm6_4_2.pickle", "rb") as fichier_graphe1 :
#         mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
#         graphe1 = mon_depickler_1.load()
# #                         
#     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_4ybb_18_2.pickle", "rb") as fichier_graphe2 :
#         mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
#         graphe2 = mon_depickler_2.load()  
#                
#         print(graphe2.nodes.data())
#         print(graphe2.edges.data())
#         sous_graphe_commun, sous_graphe_commun1, sous_graphe_commun2 = heuristique(graphe1, graphe2)
#         print(sous_graphe_commun.nodes.data())
#         print(sous_graphe_commun.edges.data())
#         sim = calcul_sim_aretes_avec_coeff(graphe1, graphe2, sous_graphe_commun, "petit rat", 1,1,1)
#                            
#         print("petit rat")
#         print(sim)
#                           
#         dico_graphe_sim = {(('5dm6', 4), ('4ybb', 18)) : {"graphe" : sous_graphe_commun, "sim" : sim}}
#                           
#         with open("dico_algo_heuristique.pickle", 'wb') as fichier_sortie :
#             mon_pickler = pickle.Pickler(fichier_sortie)
#             mon_pickler.dump(dico_graphe_sim)
        
#         with open("/media/coline/Maxtor/dico_new.pickle", 'rb') as fichier_sim :
#             mon_depickler = pickle.Unpickler(fichier_sim)
#             dico_sim = mon_depickler.load() 
#              
#         if (('4ybb', 35), ('4y4o', 13)) in dico_sim.keys() :
#             print(dico_sim[(('4ybb', 35), ('4y4o', 13))]) 
#         else :
#             print(dico_sim[(('4y4o', 13),('4ybb', 35))]) 
         
    with open("liste_sim_diff_algo_heuristique_avec_modif_encore_modif_distance_2_plus_rapide_v7.pickle", 'rb') as fichier_sortie :
        mon_depickler = pickle.Unpickler(fichier_sortie)
        liste = mon_depickler.load()     
             
        with open("liste_sim_diff_algo_heuristique_avec_modif_encore_modif_distance_1_plus_rapide.pickle", 'rb') as fichier_sortie2 :
            mon_depickler = pickle.Unpickler(fichier_sortie2)
            liste2 = mon_depickler.load()
             
        print(len(liste))
                      
        somme_diff = 0.0
        maxi_diff = 0.0
        elt_maxi_diff = -1
        compteur_en_dessous = 0
        compteur_eleve = 0
        compteur = 0
        for elt in liste :
            #print(elt)
#             for elt2 in liste2 :
#                 if elt[0] == elt2[0] and elt[1] == elt2[1] :
#                     print(elt)
#                     print(elt2)
#                     compteur += 1
            diff = elt[2]-elt[1]
            if diff < 0 :
                print(elt)
                print(diff)
                compteur_en_dessous += 1
            else :
                print(elt)
                print(diff)
                if elt[2] > 0.7 :
                    compteur_eleve += 1
            somme_diff += diff
            if diff > maxi_diff :
                print(elt[0])
                maxi_diff = diff
                elt_maxi_diff = elt[0]
            print(somme_diff)
        print(compteur)              
                       
#         print(maxi_diff)
#         print(elt_maxi_diff) 
#         print(somme_diff/len(liste))
#         print(compteur_en_dessous)
#         print(compteur_eleve)
#     liste = []
#     tps1 = time.time()
#     compteur_idem = 0
#     compteur = 0
#     
#     dico_graphes_heuristique = {}
#     
#     with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_612.pickle", "rb") as fichier_dico_graphe_algo_exact :
#         mon_depickler = pickle.Unpickler(fichier_dico_graphe_algo_exact)
#         dico_graphe_sim = mon_depickler.load()
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
#         compteur_p = 0
#                        
#         for i in range(len(liste_a_garder)) :
#             #if liste_a_garder[i] == ('3g78', 2) :
#                 for j in range(i+1, len(liste_a_garder)) :
#                     #if liste_a_garder[j] == ('4w2f', 40) :
#                         print(compteur)
#                         with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s_2.pickle"%(liste_a_garder[i][0], liste_a_garder[i][1]), "rb") as fichier_graphe1 :
#                             mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
#                             graphe1 = mon_depickler_1.load()
#                                          
#                         with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s_2.pickle"%(liste_a_garder[j][0], liste_a_garder[j][1]), "rb") as fichier_graphe2 :
#                             mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
#                             graphe2 = mon_depickler_2.load()  
#                                          
#                         print(liste_a_garder[i])
#                         print(liste_a_garder[j])
#                         sous_graphe_commun, sous_graphe_commun1, sous_graphe_commun2 = heuristique(graphe1, graphe2)
#                         print(sous_graphe_commun.nodes.data())
#                         print(sous_graphe_commun.edges.data())
#                         sim = calcul_sim_aretes_avec_coeff(graphe1, graphe2, sous_graphe_commun, "petit rat", 1,1,1)
#                         
#                         dico_graphes_heuristique.update({(liste_a_garder[i], liste_a_garder[j]) : {"graphe" : sous_graphe_commun, "sim" : sim}})
#                                  
#                         if (liste_a_garder[i], liste_a_garder[j]) in dico_graphe_sim.keys() :
#                             print(sim)
#                             print(dico_graphe_sim[(liste_a_garder[i], liste_a_garder[j])]["graphe"].nodes.data())
#                             print(dico_graphe_sim[(liste_a_garder[i], liste_a_garder[j])]["graphe"].edges.data())
#                             print(dico_graphe_sim[(liste_a_garder[i], liste_a_garder[j])]["sim"])
#                                              
#                             for u,v,data in dico_graphe_sim[(liste_a_garder[i], liste_a_garder[j])]["graphe"].edges(data=True) :
#                                 if (u,v) not in sous_graphe_commun.edges() :
#                                     print(u,v)
#                                                      
#                             if sim == dico_graphe_sim[(liste_a_garder[i], liste_a_garder[j])]["sim"] : 
#                                 compteur_idem += 1
#                             else :
#                                 liste.append(((liste_a_garder[i], liste_a_garder[j]), sim, dico_graphe_sim[(liste_a_garder[i], liste_a_garder[j])]["sim"]))
#                             compteur += 1
#                                              
#                         elif (liste_a_garder[j], liste_a_garder[i]) in dico_graphe_sim.keys() :
#                             print(sim)
#                             print(dico_graphe_sim[(liste_a_garder[j], liste_a_garder[i])]["graphe"].nodes.data())
#                             print(dico_graphe_sim[(liste_a_garder[j], liste_a_garder[i])]["graphe"].edges.data())
#                             print(dico_graphe_sim[(liste_a_garder[j], liste_a_garder[i])]["sim"])
#                                              
#                             for u,v,data in dico_graphe_sim[(liste_a_garder[j], liste_a_garder[i])]["graphe"].edges(data=True) :
#                                 if (u,v) not in sous_graphe_commun.edges() :
#                                     print(u,v)
#                                                      
#                             if sim == dico_graphe_sim[(liste_a_garder[j], liste_a_garder[i])]["sim"] :
#                                 compteur_idem += 1
#                             else :
#                                 liste.append(((liste_a_garder[j], liste_a_garder[i]), sim, dico_graphe_sim[(liste_a_garder[j], liste_a_garder[i])]["sim"]))
#                             compteur += 1
#                         else :
#                             print("dommage")
#                         compteur_p += 1            
#                         #break
#                 #break
#                    
#     with open("liste_sim_diff_algo_heuristique_avec_modif_encore_modif_distance_2_plus_rapide_v7.pickle", 'wb') as fichier_sortie :
#         mon_pickler = pickle.Pickler(fichier_sortie)
#         mon_pickler.dump(liste)        
#     
#     with open("dico_graphes_heuristique_v7.pickle", 'wb') as fichier_sortie_2 :
#         mon_pickler = pickle.Pickler(fichier_sortie_2)
#         mon_pickler.dump(dico_graphes_heuristique) 
#                    
#     print("meme sim : %d"%compteur_idem)
#     print("tot : %d"%compteur)    
#     tps2 = time.time()
#     print("temps d'execution : %f"%(tps2-tps1))
#     
    
    ''' new_heuristique '''
      
    with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_612.pickle", "rb") as fichier_dico_graphe_algo_exact :
        mon_depickler = pickle.Unpickler(fichier_dico_graphe_algo_exact)
        dico_graphe_sim = mon_depickler.load()
                      
    types_arn = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16S_arnm", "arnt_16S"]
                        
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
                        if elt == ('4w2f', 16) :
                            print(typ)
                        if elt not in liste_a_garder :
                            liste_a_garder.append(elt)
        
        tps1 = time.time()  
        
        dico_graphe_sim_heuri = {}
        liste_pas_pareil = []
        compter_diff = 0
        for i in range(len(liste_a_garder)) :
            #if liste_a_garder[i] == ('5dm6',3) :
                for j in range(i+1, len(liste_a_garder)) :
                    #if liste_a_garder[j] == ('4ybb', 36) :
                        print(compteur)
                        with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s_2.pickle"%(liste_a_garder[i][0], liste_a_garder[i][1]), "rb") as fichier_graphe1 :
                            mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
                            graphe1_vrai = mon_depickler_1.load()
                                          
                        with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s_2.pickle"%(liste_a_garder[j][0], liste_a_garder[j][1]), "rb") as fichier_graphe2 :
                            mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
                            graphe2_vrai = mon_depickler_2.load()  
                        
                        print(liste_a_garder[i])
                        print(liste_a_garder[j])
                        
                        sous_graphe_commun_1 = new_heuristique(graphe1_vrai, graphe2_vrai)
#                         for noeud, data in sous_graphe_commun.nodes(data=True) :
#                             print(noeud)
#                             print(data)
#                         for u,v, data in sous_graphe_commun.edges(data=True) :
#                             print(u,v)
#                             print(data)
#                         print(graphe1_vrai.nodes())
                        sim_1 = calcul_sim_aretes_avec_coeff(graphe1_vrai, graphe2_vrai, sous_graphe_commun_1, "rat", 1, 1, 1)
                        
                        sous_graphe_commun_2 = new_heuristique(graphe2_vrai, graphe1_vrai)
                        sim_2 = calcul_sim_aretes_avec_coeff(graphe2_vrai, graphe1_vrai, sous_graphe_commun_2, "rat", 1, 1, 1)
                        
                        print(sim_1)
                        print(sim_2)
                        
                        #dico_graphe_sim_heuri.update({(liste_a_garder[i], liste_a_garder[j]) : {"graphe" : sous_graphe_commun_1, "sim" : sim_1}})
                        #dico_graphe_sim_heuri.update({(liste_a_garder[j], liste_a_garder[i]) : {"graphe" : sous_graphe_commun_2, "sim" : sim_2}})

                        
                        if sim_1 >= sim_2 :
                           
                            dico_graphe_sim_heuri.update({(liste_a_garder[i], liste_a_garder[j]) : {"graphe" : sous_graphe_commun_1, "sim" : sim_1}})
                               
                            if (liste_a_garder[i], liste_a_garder[j]) in dico_graphe_sim.keys():
                                if sim_1 != dico_graphe_sim[(liste_a_garder[i], liste_a_garder[j])]["sim"] :
                                    liste_pas_pareil.append((liste_a_garder[i], liste_a_garder[j]))
                            else :
                                if sim_1 != dico_graphe_sim[(liste_a_garder[j], liste_a_garder[i])]["sim"] :
                                    liste_pas_pareil.append((liste_a_garder[i], liste_a_garder[j]))  
                        else :
                            dico_graphe_sim_heuri.update({(liste_a_garder[j], liste_a_garder[i]) : {"graphe" : sous_graphe_commun_2, "sim" : sim_2}})
                               
                            if (liste_a_garder[i], liste_a_garder[j]) in dico_graphe_sim.keys():
                                if sim_2 != dico_graphe_sim[(liste_a_garder[i], liste_a_garder[j])]["sim"] :
                                    liste_pas_pareil.append((liste_a_garder[j], liste_a_garder[i]))
                            else :
                                if sim_2 != dico_graphe_sim[(liste_a_garder[j], liste_a_garder[i])]["sim"] :
                                    liste_pas_pareil.append((liste_a_garder[j], liste_a_garder[i]))  
                          
                        if sim_1 != sim_2 :
                            compter_diff += 1
                                
                        compteur += 1
#                         
        print(len(liste_pas_pareil))
        print(compter_diff)
         
        with open("fichier_diff.pickle", 'wb') as fichier_diff :
            mon_pickler = pickle.Pickler(fichier_diff)
            mon_pickler.dump(liste_pas_pareil)
# #         
#         print(liste_pas_pareil)
        with open("dico_algo_heuristique_new_v.pickle", 'wb') as fichier_sortie :
            mon_pickler = pickle.Pickler(fichier_sortie)
            mon_pickler.dump(dico_graphe_sim_heuri)
            
        tps2 = time.time()
        print("temps d'execution : ")
        print(str(tps2-tps1))
    