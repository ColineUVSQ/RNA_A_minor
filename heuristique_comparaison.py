'''
Created on 4 nov. 2019

@author: coline
'''

import networkx as nx
import pickle
import time
import copy
from recup_data.new_algo_comparaison import calcul_sim_aretes_avec_coeff,\
    recup_chaines
from recup_data.constantes import NEW_EXTENSION_PATH_TAILLE

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
       

def calcul_sim_aretes_avec_coeff_heuristique(graphe1, graphe2, graphe_commun, cle, coeffc, coeffa, coeffn):
    aretes_1 = calcul_aretes_avec_coeff_heuristique(graphe1, cle, coeffc, coeffa, coeffn)
    aretes_2 = calcul_aretes_avec_coeff_heuristique(graphe2, cle, coeffc, coeffa, coeffn)
    aretes_commun = calcul_aretes_communes_avec_coeff_heuristique(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn)
    
#     print(aretes_1)
#     print(aretes_2)
#     print(aretes_commun)
    
    return aretes_commun/max(aretes_1, aretes_2)

def meme_type_meme_chaine(noeud1, noeud2, graphe1, graphe2):
    meme_chaine = False
    for elt in graphe1.nodes[noeud1]["chaine"] :
        if elt in graphe2.nodes[noeud2]["chaine"] :
            meme_chaine = True
    if (graphe1.nodes[noeud1]["type"] == graphe2.nodes[noeud2]["type"] and meme_chaine) :
        return True
    else :
        return False

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
    les sommets de la meme chaine qui sont de meme type (et aussi -1 pour pouvoir considerer le cas ou le sommet n'est pas superpose) '''
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


''' on veut essayer d'etendre la superposition de sommet1 et sommet2 aux sommets voisins '''
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


''' on veut essayer d'etendre la superposition de sommet1 et sommet2 aux sommets voisins '''
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
#             for voisin_2 in graphe2.predecessors(noeud[1]) :
#                 for voisin_voisin_2 in graphe2.predecessors(voisin_2) :
#                             if voisin_voisin_2 not in sous_graphe_commun_2.nodes() :
#                                 for voisin_1 in graphe1.predecessors(noeud[0]) :
#                                     #if voisin_1 not in [2,3,4] :
#                                         
#                                         for voisin_voisin_1 in graphe1[voisin_1] :
#                                             if voisin_voisin_1 not in sous_graphe_commun_1.nodes() :
#                                                 if (voisin_voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir  and (voisin_voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)] :
#                                                     if meme_type_meme_chaine(voisin_voisin_1, voisin_voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_voisin_1, voisin_voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_voisin_1, voisin_voisin_2)) :                                                 
#                                                         liste_noeuds_a_voir.append((voisin_voisin_1, voisin_voisin_2))
#                                         for voisin_voisin_1 in graphe1.predecessors(voisin_1) :
#                                             if voisin_voisin_1 not in sous_graphe_commun_1.nodes() :
#                                                 if (voisin_voisin_1, voisin_voisin_2) not in liste_noeuds_a_voir and (voisin_voisin_1, voisin_voisin_2) not in [(1,1), (2,2), (3,3), (4,4)]:
#                                                     if meme_type_meme_chaine(voisin_voisin_1, voisin_voisin_2, graphe1, graphe2) and test_compatibilite(sous_graphe_commun, (voisin_voisin_1, voisin_voisin_2), graphe1, graphe2) and not dans_graphe(sous_graphe_commun, (voisin_voisin_1, voisin_voisin_2)) :                                                 
#                                                         liste_noeuds_a_voir.append((voisin_voisin_1, voisin_voisin_2))
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
    
    
    
if __name__ == '__main__':
#     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_5dm6_4_2.pickle", "rb") as fichier_graphe1 :
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
         
#     with open("liste_sim_diff_algo_heuristique_avec_modif_encore_modif_distance_2_plus_rapide_v3.pickle", 'rb') as fichier_sortie :
#         mon_depickler = pickle.Unpickler(fichier_sortie)
#         liste = mon_depickler.load()     
#             
#         with open("liste_sim_diff_algo_heuristique_avec_modif_encore_modif_distance_1_plus_rapide.pickle", 'rb') as fichier_sortie2 :
#             mon_depickler = pickle.Unpickler(fichier_sortie2)
#             liste2 = mon_depickler.load()
#             
#         print(len(liste))
#                      
#         somme_diff = 0.0
#         maxi_diff = 0.0
#         elt_maxi_diff = -1
#         compteur_en_dessous = 0
#         compteur_eleve = 0
#         compteur = 0
#         for elt in liste :
#             #print(elt)
# #             for elt2 in liste2 :
# #                 if elt[0] == elt2[0] and elt[1] == elt2[1] :
# #                     print(elt)
# #                     print(elt2)
# #                     compteur += 1
#             diff = elt[2]-elt[1]
#             if diff < 0 :
#                 print(elt)
#                 print(diff)
#                 compteur_en_dessous += 1
#             else :
#                 print(elt)
#                 print(diff)
#                 if elt[2] > 0.7 :
#                     compteur_eleve += 1
#             somme_diff += diff
#             if diff > maxi_diff :
#                 print(elt[0])
#                 maxi_diff = diff
#                 elt_maxi_diff = elt[0]
#             print(somme_diff)
#         print(compteur)              
                       
#         print(maxi_diff)
#         print(elt_maxi_diff) 
#         print(somme_diff/len(liste))
#         print(compteur_en_dessous)
#         print(compteur_eleve)
    liste = []
    tps1 = time.time()
    compteur_idem = 0
    compteur = 0
    
    dico_graphes_heuristique = {}
    
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
                            
#     with open("occ_multi_chaine.pickle", 'rb') as fichier_multi_chaines :
#         mon_depickler = pickle.Unpickler(fichier_multi_chaines)
#         liste_plusieurs_chaines = mon_depickler.load() 
                       
        print(len(liste_a_garder))
        compteur_p = 0
                       
        for i in range(len(liste_a_garder)) :
            #if liste_a_garder[i] == ('3g78', 2) :
                for j in range(i+1, len(liste_a_garder)) :
                    #if liste_a_garder[j] == ('4w2f', 40) :
                        print(compteur)
                        with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s_2.pickle"%(liste_a_garder[i][0], liste_a_garder[i][1]), "rb") as fichier_graphe1 :
                            mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
                            graphe1 = mon_depickler_1.load()
                                         
                        with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s_2.pickle"%(liste_a_garder[j][0], liste_a_garder[j][1]), "rb") as fichier_graphe2 :
                            mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
                            graphe2 = mon_depickler_2.load()  
                                         
                        print(liste_a_garder[i])
                        print(liste_a_garder[j])
                        sous_graphe_commun, sous_graphe_commun1, sous_graphe_commun2 = heuristique(graphe1, graphe2)
                        print(sous_graphe_commun.nodes.data())
                        print(sous_graphe_commun.edges.data())
                        sim = calcul_sim_aretes_avec_coeff(graphe1, graphe2, sous_graphe_commun, "petit rat", 1,1,1)
                        
                        dico_graphes_heuristique.update({(liste_a_garder[i], liste_a_garder[j]) : {"graphe" : sous_graphe_commun, "sim" : sim}})
                                 
                        if (liste_a_garder[i], liste_a_garder[j]) in dico_graphe_sim.keys() :
                            print(sim)
                            print(dico_graphe_sim[(liste_a_garder[i], liste_a_garder[j])]["graphe"].nodes.data())
                            print(dico_graphe_sim[(liste_a_garder[i], liste_a_garder[j])]["graphe"].edges.data())
                            print(dico_graphe_sim[(liste_a_garder[i], liste_a_garder[j])]["sim"])
                                             
                            for u,v,data in dico_graphe_sim[(liste_a_garder[i], liste_a_garder[j])]["graphe"].edges(data=True) :
                                if (u,v) not in sous_graphe_commun.edges() :
                                    print(u,v)
                                                     
                            if sim == dico_graphe_sim[(liste_a_garder[i], liste_a_garder[j])]["sim"] : 
                                compteur_idem += 1
                            else :
                                liste.append(((liste_a_garder[i], liste_a_garder[j]), sim, dico_graphe_sim[(liste_a_garder[i], liste_a_garder[j])]["sim"]))
                            compteur += 1
                                             
                        elif (liste_a_garder[j], liste_a_garder[i]) in dico_graphe_sim.keys() :
                            print(sim)
                            print(dico_graphe_sim[(liste_a_garder[j], liste_a_garder[i])]["graphe"].nodes.data())
                            print(dico_graphe_sim[(liste_a_garder[j], liste_a_garder[i])]["graphe"].edges.data())
                            print(dico_graphe_sim[(liste_a_garder[j], liste_a_garder[i])]["sim"])
                                             
                            for u,v,data in dico_graphe_sim[(liste_a_garder[j], liste_a_garder[i])]["graphe"].edges(data=True) :
                                if (u,v) not in sous_graphe_commun.edges() :
                                    print(u,v)
                                                     
                            if sim == dico_graphe_sim[(liste_a_garder[j], liste_a_garder[i])]["sim"] :
                                compteur_idem += 1
                            else :
                                liste.append(((liste_a_garder[j], liste_a_garder[i]), sim, dico_graphe_sim[(liste_a_garder[j], liste_a_garder[i])]["sim"]))
                            compteur += 1
                        else :
                            print("dommage")
                        compteur_p += 1            
                        #break
                #break
                   
    with open("liste_sim_diff_algo_heuristique_avec_modif_encore_modif_distance_2_plus_rapide_v6.pickle", 'wb') as fichier_sortie :
        mon_pickler = pickle.Pickler(fichier_sortie)
        mon_pickler.dump(liste)        
    
    with open("dico_graphes_heuristique_v6.pickle", 'wb') as fichier_sortie_2 :
        mon_pickler = pickle.Pickler(fichier_sortie_2)
        mon_pickler.dump(dico_graphes_heuristique) 
                   
    print("meme sim : %d"%compteur_idem)
    print("tot : %d"%compteur)    
    tps2 = time.time()
    print("temps d'execution : %f"%(tps2-tps1))
#     
    
                
    
    