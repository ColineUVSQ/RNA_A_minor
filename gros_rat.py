'''
Created on 10 sept. 2019

@author: coline
'''

'''
Created on 29 août 2019

@author: coline
'''

import networkx as nx
import os
import pickle
import matplotlib.pyplot as plt
import multiprocessing
import tabnanny

def chaines_reliees(graphe):
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

def calcul_aretes_communes_avec_coeff(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
    for u,v,data in graphe_commun.edges(data=True) :
        if data["label"] != 'B53' :
            if data["label"] == '0' :
                if coeffa == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
            elif graphe1.nodes[u[0]]["type"] == 1 and graphe1.nodes[v[0]]["type"] == 1 :
                if coeffc == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
            else :
                if coeffn == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
                    
    somme_aretes = somme_aretes/2 - 4
    
    return somme_aretes
       

def calcul_sim_aretes_avec_coeff(graphe1, graphe2, graphe_commun, cle, coeffc, coeffa, coeffn):
    aretes_1 = calcul_aretes_avec_coeff(graphe1, cle, coeffc, coeffa, coeffn)
    aretes_2 = calcul_aretes_avec_coeff(graphe2, cle, coeffc, coeffa, coeffn)
    aretes_commun = calcul_aretes_communes_avec_coeff(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn)
    
#     print(aretes_1)
#     print(aretes_2)
#     print(aretes_commun)
    
    return aretes_commun/max(aretes_1, aretes_2)


def test_compatibilite(graphe_commun, noeud, graphe1, graphe2):
    
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

def recup_chaines(graphe):
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
    for noeud in graphe.nodes() :
        #print(noeud)
        if (noeud[0] == couple_a_chercher[0] and noeud[1] != couple_a_chercher[1]) or (noeud[1] == couple_a_chercher[1] and noeud[0] != couple_a_chercher[0])  :
            return True
    return False

def meme_type_meme_chaine(noeud1, noeud2, graphe1, graphe2):
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
        
        
        #print(liste_etapes)
#         print(etape)
        for i in range(etape, len(liste_etapes)) :
#             print(liste_etapes[i])
            if "noeuds" in liste_etapes[i].keys() :
                graphe_commun_temp.remove_nodes_from(liste_etapes[i]["noeuds"])
            if "aretes" in liste_etapes[i].keys() :
                graphe_commun_temp.remove_edges_from(liste_etapes[i]["aretes"])
        
        if y != -1 :
            
            ## ajout de la paire (x,y) au graphe commun
            taille_liste_etapes = len(liste_etapes)
            if not dans_graphe(graphe_commun_temp, (x,y)) and (x,y) not in graphe_commun_temp.nodes() :
                graphe_commun_temp.add_node((x,y))
                liste_etapes.append({"noeuds" :[(x,y)]})
            
            ## recherche de voisins non cov de x et y a ajouter aussi
            
            for voisin1 in graphe1[x] :
                
                for edge1 in graphe1[x][voisin1] :
                    if graphe1[x][voisin1][edge1]["label"] != 'B53' :
                        label1 = graphe1[x][voisin1][edge1]["label"]
                        #long_range1 = graphe1[x][voisin1][edge1]["long_range"]
                        for voisin2 in graphe2[y] :
                            for edge2 in graphe2[y][voisin2] :
                                if graphe2[y][voisin2][edge2]["label"] != 'B53' :
#                                     if x == 3 :
#                                         print('voisin')
#                                         print(voisin1)
                                    if label1 == graphe2[y][voisin2][edge2]["label"] : #and long_range1 == graphe2[y][voisin2][edge2]["long_range"] : ## meme label des aretes
                                        if not dans_graphe(graphe_commun_temp, (voisin1, voisin2))  : ## verif elements de la paire des voisins ne sont pas deja dans ue autre paire
                                            if meme_type_meme_chaine(voisin1, voisin2, graphe1, graphe2)  :
#                                                 print(liste_etapes)
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
                                                    graphe_commun_temp.add_edge((x,y), (voisin1, voisin2), label=label1)
                                                    if len(label1) == 3 :
                                                        graphe_commun_temp.add_edge((voisin1, voisin2), (x,y), label=label1[0]+label1[2]+label1[1])
                                                    elif len(label1) == 4 :
                                                        graphe_commun_temp.add_edge((voisin1, voisin2), (x,y), label=label1[0]+label1[2]+label1[1]+label1[3:])
                                                    else : ## aretes artificielles
                                                        graphe_commun_temp.add_edge((voisin1, voisin2), (x,y), label=label1)
                                                        
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
            sim = calcul_aretes_communes_avec_coeff(graphe_commun_temp, graphe1, graphe2, cle, 1, 1, 1)
            if sim > sim_max :
                sim_max = sim
                graphe_commun_max = graphe_commun_temp.copy() 
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
                suivant = chaines_1[num_chaine][num_x+1]
                for super_possible in dico_chaines_1[num_chaine][suivant] :
                    if super_possible == -1 or (test_compatibilite(graphe_commun_temp, (suivant, super_possible), graphe1, graphe2) and not dans_graphe(graphe_commun_temp, (suivant, super_possible)) and (suivant, super_possible) not in graphe_commun_temp.nodes() ):
                        pile.append((suivant, super_possible, num_chaine, num_x+1, len(liste_etapes)))  
                        
            else :
                ## on change de chaine
                suivant = chaines_1[num_chaine+1][0]
#                 print(suivant)
                for super_possible in dico_chaines_1[num_chaine+1][suivant] :
                    if super_possible == -1 or suivant in [1,2,3,4] or (test_compatibilite(graphe_commun_temp, (suivant, super_possible), graphe1, graphe2) and not dans_graphe(graphe_commun_temp, (suivant, super_possible)) and (suivant, super_possible) not in graphe_commun_temp.nodes() ):
                        pile.append((suivant, super_possible, num_chaine+1, 0, len(liste_etapes)))
    
    return graphe_commun_max
                        
#                 print("pile courante")
#                 print(pile)

def comparaison(graphe1, graphe2, cle):
    chaines_1 = recup_chaines(graphe1)
    chaines_2 = recup_chaines(graphe2)
    
    
    
    #print(chaines_1)
    #print(chaines_2)
    
    dico_chaines_1 = liste_meme_type(chaines_1, graphe1, chaines_2, graphe2)
    #dico_chaines_2 = liste_meme_type(chaines_2, graphe2, chaines_1, graphe1)
    #print(dico_chaines_1)
    
    chaines_reliees_1 = chaines_reliees(graphe1)
    chaines_reliees_2 = chaines_reliees(graphe2)
    #print(chaines_reliees_1)
    #print(chaines_reliees_2)
    
    chaines_reliees_tot = []
    chaines_reliees_tot.extend(list(chaines_reliees_1))
    #print(chaines_reliees_1)
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
                
    #print(chaines_a_mettre_ensemble)
        
                
    
    graphe_commun_temp = nx.MultiDiGraph()
    
    for i in range(1,6) :
        graphe_commun_temp.add_node((i,i))
    graphe_commun_temp.add_edge((1,1),(2,2), label="CSS")
    graphe_commun_temp.add_edge((2,2),(1,1), label="CSS")
    graphe_commun_temp.add_edge((1,1),(5,5), label="TSS")
    graphe_commun_temp.add_edge((5,5),(1,1), label="TSS")
    graphe_commun_temp.add_edge((2,2),(5,5), label="CWW")
    graphe_commun_temp.add_edge((5,5),(2,2), label="CWW")
    graphe_commun_temp.add_edge((3,3),(4,4), label="CSS")
    graphe_commun_temp.add_edge((4,4),(3,3), label="CSS")
           
    #print(chaines_reliees_tot)
    for elt in chaines_a_mettre_ensemble :
        chaines = []
        dico_chaines = []
        for chaine in elt :
            chaines.append(chaines_1[chaine-1])
            dico_chaines.append(dico_chaines_1[chaine-1])
        graphe_commun_temp = algo_principal(graphe1, graphe2, chaines, dico_chaines, graphe_commun_temp, cle)
        
    graphe_commun_max = graphe_commun_temp.copy()
    sim_max = calcul_sim_aretes_avec_coeff(graphe1, graphe2, graphe_commun_max, cle, 1,1,1)
    
        
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
            graphe_commun_max.add_edge(noeud, (num_voisin1_b53, num_voisin2_b53), label='B53')
            
            
    return graphe_commun_max, sim_max


def toutes_comparaison():
    
    liste_fichiers = []
    for fic in os.listdir("Nouvelles_donnees") :
        if "pickle" in fic : #and len(fic.split("_")) == 6 :
            liste_fichiers.append(fic)
    print(liste_fichiers)
    print(len(liste_fichiers))
    
    with open("Nouvelles_donnees/Resultats/dico_graphe_sim_new_algo.pickle", 'rb') as fichier_lecture :
        mon_depickler = pickle.Unpickler(fichier_lecture)
        dico_graphe = mon_depickler.load()
        
        compteur = 0
        for i in range(len(liste_fichiers)) :
            #print(liste_fichiers[i])
            
            for j in range(i+1, len(liste_fichiers)) :
                #print(liste_fichiers[j])
                
                #if liste_fichiers[i] == "fichier_4V9F_0_48_2_3.pickle" and liste_fichiers[j] == "fichier_3JCS_1_25_16_3.pickle" :
                if compteur < 10000 and (liste_fichiers[i][8:len(liste_fichiers[i])-7], liste_fichiers[j][8:len(liste_fichiers[j])-7]) not in dico_graphe.keys() :  
                
                    with open(NEW_EXTENSION_PATH_TAILLE+liste_fichiers[i], 'rb') as fichier1 :
                        mon_depickler = pickle.Unpickler(fichier1)
                        graphe1 = mon_depickler.load()
                        
                        
                        
                    with open(NEW_EXTENSION_PATH_TAILLE+liste_fichiers[j], 'rb') as fichier2 :
                        mon_depickler = pickle.Unpickler(fichier2)
                        graphe2 = mon_depickler.load()
        #                 tab = chaines_reliees(graphe2)
        #                 print(tab)
        #                 if len(tab) > 0 :
        #                     compteur +=1
                            
            
                        
                         
                        graphe_commun_max, sim_max = comparaison(graphe1, graphe2, "petit rat")   
                         
                        #dico_sim.update({(liste_fichiers[i][8:len(liste_fichiers[i])-7], liste_fichiers[j][8:len(liste_fichiers[j])-7]) : sim_max})
                        print("sim max")
                        print(sim_max)
                        print("graphe commun max")
                        print(graphe_commun_max.nodes.data())
                        print(graphe_commun_max.edges.data())
                        
                        dico_graphe.update({(liste_fichiers[i][8:len(liste_fichiers[i])-7], liste_fichiers[j][8:len(liste_fichiers[j])-7]) : {"sim" : sim_max, "graphe" : graphe_commun_max}})
                        
                        compteur += 1
            #break
    
    with open("Nouvelles_donnees/Resultats/dico_graphe_sim_new_algo.pickle", 'wb') as fichier_graphe :
        mon_pickler = pickle.Pickler(fichier_graphe)
        mon_pickler.dump(dico_graphe)
            
#     with open("Nouvelles_donnees/Resultats/dico_sim_new_algo.pickle", 'wb') as fichier_ecriture :
#         mon_pickler = pickle.Pickler(fichier_ecriture)
#         mon_pickler.dump(dico_sim)
 
'''09/09/2019
effectuer les comparaisons au sein des groupes d'homologues => fait jusquau groupe 3 des ARNr 23S'''  
def comparaison_homologues(num):     
    
    #compteur = 0
        
    with open("groupes_23S_homologues.pickle", 'rb') as fichier_homologues :
        mon_depickler_1 = pickle.Unpickler(fichier_homologues)
        groupes_homologues = mon_depickler_1.load()
        compteur = 0
        for groupe in groupes_homologues :
            if compteur == num :
                os.makedirs("Nouvelles_donnees/groupe_23S_homologues_%s"%compteur)
                #dico_graphe = {}
                for i in range(len(groupe)) :
                    print(groupe[i])
                    for j in range(i+1, len(groupe)) :
                        print(groupe[j])
                        with open("Nouvelles_donnees/fichier_"+str(groupe[i][0])+"_"+str(groupe[i][1])+".pickle", 'rb') as fichier1 :
                                mon_depickler = pickle.Unpickler(fichier1)
                                graphe1 = mon_depickler.load()
                                
                                
                                
                        with open("Nouvelles_donnees/fichier_"+str(groupe[j][0])+"_"+str(groupe[j][1])+".pickle", 'rb') as fichier2 :
                                mon_depickler = pickle.Unpickler(fichier2)
                                graphe2 = mon_depickler.load()
                #                 tab = chaines_reliees(graphe2)
                #                 print(tab)
                #                 if len(tab) > 0 :
                #                     compteur +=1
                                    
                    
                                
                                 
                                graphe_commun_max, sim_max = comparaison(graphe1, graphe2, "petit rat")   
                                 
                                #dico_sim.update({(liste_fichiers[i][8:len(liste_fichiers[i])-7], liste_fichiers[j][8:len(liste_fichiers[j])-7]) : sim_max})
                                #print("sim max")
                                #print(sim_max)
                                #print("graphe commun max")
                                #print(graphe_commun_max.nodes.data())
                                #print(graphe_commun_max.edges.data())
                                
                                with open(NEW_EXTENSION_PATH_TAILLE+"/Resultats/graphe_%s_%s.pickle"%(str(groupe[i][0])+"_"+str(groupe[i][1]), str(groupe[j][0])+"_"+str(groupe[j][1])), 'wb') as fichier_ecriture :
                                    mon_pickler = pickle.Pickler(fichier_ecriture)
                                    mon_pickler.dump(list(graphe_commun_max.edges()))
                                
                                #dico_graphe.update({(str(groupe[i][0])+"_"+str(groupe[i][1]), str(groupe[j][0])+"_"+str(groupe[j][1])) : {"sim" : sim_max, "graphe" : graphe_commun_max}})
                                
    #                             compteur += 1
                    #break
            
                #with open("Nouvelles_donnees/Resultats/dico_graphe_sim_new_algo_23s_homologues_groupe_%s.pickle"%compteur, 'wb') as fichier_graphe :
                #    mon_pickler = pickle.Pickler(fichier_graphe)
                #    mon_pickler.dump(dico_graphe)
            compteur += 1
            
#     with open("Nouvelles_donnees/Resultats/dico_sim_new_algo.pickle", 'wb') as fichier_ecriture :
#         mon_pickler = pickle.Pickler(fichier_ecriture)
#         mon_pickler.dump(dico_si

'''09/09/19
chercher les vrais identiques : d'un point de vue extension aussi '''
def obs_resultats():
#     with open("Nouvelles_donnees/Resultats/dico_graphe_sim_new_algo_23s_homologues_groupe_%s.pickle"%5, 'rb') as fichier_graphe :
#         mon_depickler = pickle.Unpickler(fichier_graphe)
#         dico_graphe_sim = mon_depickler.load()
#         
#         print(len(dico_graphe_sim))
#         #print(dico_graphe_sim.keys())
#         
#         compter_idem = 0
#         tab_sim = []
#         for cle in dico_graphe_sim.keys() :
#             tab_sim.append(dico_graphe_sim[cle]["sim"]) 
#             if dico_graphe_sim[cle]["sim"] ==  1.0 :
# #                 print(cle)
# #                 print(dico_graphe_sim[cle])
#                 #print(dico_graphe_sim[cle]["graphe"].edges.data())
#                 compter_idem += 1
#         plt.plot(tab_sim)
#         plt.show()
        #print(compter_idem)
        compteur_nb_groupes_vraiment_identiques = 0
        with open("groupes_23S_homologues.pickle", 'rb') as fichier_homologues :
            mon_depickler_1 = pickle.Unpickler(fichier_homologues)
            groupes_homologues = mon_depickler_1.load()
             
            with open("groupes_23S_identiques.pickle", 'rb') as fichier_identiques :
                mon_depickler_2 = pickle.Unpickler(fichier_identiques)
                groupes_identiques = mon_depickler_2.load()
                 
                 
                
                 
                compter_vraiment_idem = 0
                compteur = 0
                for groupes in groupes_identiques :
                    
                    
                    if compteur <= 5 :
                        with open("Nouvelles_donnees/Resultats/dico_graphe_sim_new_algo_23s_homologues_groupe_%s.pickle"%compteur, 'rb') as fichier_graphe :
                            mon_depickler = pickle.Unpickler(fichier_graphe)
                            dico_graphe_sim = mon_depickler.load()
                        
                            print(len(groupes))
                            for groupe in groupes :
                                groupes_vraiment_identiques = []
                                print(len(groupe))
                                print(groupe)
    #                             for i in range(len(groupe)) :
    #                                 for j in range(i+1, len(groupe)) :
    #                                     if (str(groupe[i][0])+"_"+str(groupe[i][1]), str(groupe[j][0])+"_"+str(groupe[j][1])) in dico_graphe_sim.keys() :
    #                                         print(dico_graphe_sim[(str(groupe[i][0])+"_"+str(groupe[i][1]), str(groupe[j][0])+"_"+str(groupe[j][1]))])
    #                                     else  :
    #                                         print(dico_graphe_sim[(str(groupe[j][0])+"_"+str(groupe[j][1]), str(groupe[i][0])+"_"+str(groupe[i][1]))])
                                for i in range(len(groupe)) :
                                    for j in range(i+1, len(groupe)) :
                                        #print((groupe[i], groupe[j]))
                                        if (str(groupe[i][0])+"_"+str(groupe[i][1]), str(groupe[j][0])+"_"+str(groupe[j][1])) in dico_graphe_sim.keys() :
                                            if dico_graphe_sim[(str(groupe[i][0])+"_"+str(groupe[i][1]), str(groupe[j][0])+"_"+str(groupe[j][1]))]["sim"] == 1.0 :
                                                #print((groupe[i], groupe[j]))
                                            #else :
                                                ok = False
                                                compter_vraiment_idem += 1
                                                for elt in groupes_vraiment_identiques :
                                                    if groupe[j] in elt and groupe[i] not in elt :
                                                        elt.append(groupe[i])
                                                        ok = True
                                                    elif groupe[i] in elt and groupe[j] not in elt :
                                                        elt.append(groupe[j])
                                                        ok = True
                                                    elif groupe[i] in elt and groupe[j] in elt :
                                                        ok = True
                                                        
                                                if not ok :
                                                    groupes_vraiment_identiques.append([groupe[i], groupe[j]])
                                        else  :
                                            if dico_graphe_sim[(str(groupe[j][0])+"_"+str(groupe[j][1]), str(groupe[i][0])+"_"+str(groupe[i][1]))]["sim"] == 1.0 :
    #                                             print((groupe[j], groupe[i]))
    #                                         else :
                                                ok = False
                                                compter_vraiment_idem += 1
                                                for elt in groupes_vraiment_identiques :
                                                    if groupe[j] in elt and groupe[i] not in elt :
                                                        elt.append(groupe[i])
                                                        ok = True
                                                    elif groupe[i] in elt and groupe[j] not in elt :
                                                        elt.append(groupe[j])
                                                        ok = True
                                                    elif groupe[i] in elt and groupe[j] in elt :
                                                        ok = True
                                                        
                                                if not ok :
                                                    groupes_vraiment_identiques.append([groupe[i], groupe[j]])
                                                    
                                print(groupes_vraiment_identiques)
                                print("nb de groupes vraiment identiques")
                                print(len(groupes_vraiment_identiques))
                                compteur_nb_groupes_vraiment_identiques += len(groupes_vraiment_identiques)
                                
                                compteur_seuls = 0
                                for elt in groupe :
                                    y_est = False
                                    for groupes in groupes_vraiment_identiques :
                                        for elt_id in groupes :
                                            if elt == elt_id :
                                                y_est = True
                                    if not y_est :
                                        compteur_seuls += 1
                                
                                print("tout seul")
                                print(compteur_seuls)
                                
                                compteur_nb_groupes_vraiment_identiques += compteur_seuls
                                for i in range(len(groupes_vraiment_identiques)) :
                                    for j in range(i+1, len(groupes_vraiment_identiques)) :
                                        for elt1 in groupes_vraiment_identiques[i] :
                                            for elt2 in groupes_vraiment_identiques[j] : 
                                                if elt1 == elt2 :
                                                    print("gros rat")
                                                    print(elt1)
                                
                             
                    compteur += 1
        print(compter_vraiment_idem)
        print(compteur_nb_groupes_vraiment_identiques)
#         plt.plot(tab_sim)
#         plt.show()
if __name__ == '__main__':
    
    with open("groupes_23S_homologues.pickle", "rb") as fichier_homologues :
        mon_depickler_1 = pickle.Unpickler(fichier_homologues)
        groupes_homologues = mon_depickler_1.load()
        print("petit rat")
        comparaison_homologues(6)
        #liste_a_faire = []
        
        #for i in range(6, len(groupes_homologues)) :
        #    if "dico_graphe_sim_new_algo_23s_homologues_groupe_%d.pickle"%i not in os.listdir("Nouvelles_donnees/Resultats") :
        #        liste_a_faire.append(i)
        #print(liste_a_faire)
        #p = multiprocessing.Pool(8)
        #result = p.map(comparaison_homologues, liste_a_faire)
        #p.close()
        #p.join()
        
   
    
    
