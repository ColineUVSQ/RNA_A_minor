'''
Created on 4 févr. 2019

@author: coline
'''
import os
import pickle
import multiprocessing
import itertools        
import networkx as nx
import math

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


def recherche_couples_grands_graphes(graphe1, graphe2, chaines_1, chaines_2):
   
    couples_possibles = []
    for num_chaine in range(0, 4) :
        print(chaines_1[num_chaine])
        print(chaines_2[num_chaine])
        
        with open("liste_combi_9.pickle", 'rb') as fichier_liste :
                mon_depickler_liste = pickle.Unpickler(fichier_liste)
                liste_combi_9 = mon_depickler_liste.load()
                #print(liste_combi_9)
                liste_couples_2 = []
                liste_couples = []
                for groupe in liste_combi_9 :
                    liste = []
                    for elt in groupe :
                        
                        if graphe1.nodes[chaines_1[num_chaine][elt[0]-1]]["type"] == graphe2.nodes[chaines_2[num_chaine][elt[1]-1]]["type"] : #and (abs(graphe1.nodes[elt[0]]["poids"] - graphe2.nodes[elt[1]]["poids"]) < 2) :  
                            if  (chaines_1[num_chaine][elt[0]-1], chaines_2[num_chaine][elt[1]-1]) not in liste :
                                liste.append((chaines_1[num_chaine][elt[0]-1], chaines_2[num_chaine][elt[1]-1]))
                                #print("ramou")
#                         print(elt)
#                         print(groupe_ok)
#                         time.sleep(2)
                      
                    if len(liste) != 0 :
                        liste_couples.append(liste)        
                 
                a_enlever = []
                for i in range(0,len(liste_couples)) :
                    for j in range(0, len(liste_couples)) :
                            if j != i :
                                nb_existe_deja = 0
                                for m in range(0, len(liste_couples[i])) :
                                    for l in range(0, len(liste_couples[j])) :
                                        if liste_couples[i][m][0] == liste_couples[j][l][0] and liste_couples[i][m][1] == liste_couples[j][l][1]  :
                                            nb_existe_deja += 1
                                if nb_existe_deja == len(liste_couples[i]) and i not in a_enlever :
                                    if len(liste_couples[i]) == len(liste_couples[j]) :
                                        a_enlever.append(j)
                                    else :
                                        a_enlever.append(i)

                for i in range(len(liste_couples)) :
                    if i not in a_enlever :
                        liste_couples_2.append(liste_couples[i])
                
                couples_possibles.append(liste_couples_2)
    return couples_possibles


def test_compatibilite_grands_graphes(graphe_commun, graphe1, graphe2):
    
    for i in range(len(list(graphe_commun.nodes()))) :
        for j in range(i+1, len(list(graphe_commun.nodes()))) :
            noeud1 = list(graphe_commun.nodes())[i]
            noeud2 = list(graphe_commun.nodes())[j]
            meme_chaine = True
            for elt in graphe1.nodes[noeud1[0]]["chaine"] :
                if elt not in graphe1.nodes[noeud2[0]]["chaine"] :
                    meme_chaine = False
            if noeud1 != noeud2 and meme_chaine and graphe1.nodes[noeud1[0]]["type"] != None and graphe1.nodes[noeud2[0]]["type"] != None and noeud1[0] not in [1,2,3,4,5] and noeud2[0] not in [1,2,3,4,5]:
                if (noeud1[0] < noeud2[0] and noeud1[1] > noeud2[1]) or (noeud1[0] > noeud2[0] and noeud1[1] < noeud2[1]) :
                    return False
    
    return True

def dans_graphe(graphe, couple_a_chercher):
    pas_bon = []
    for noeud in graphe.nodes() :
        #print(noeud)
        if (noeud[0] == couple_a_chercher[0] and noeud[1] != couple_a_chercher[1]) or (noeud[1] == couple_a_chercher[1] and noeud[0] != couple_a_chercher[0])  :
            pas_bon.append(noeud)
    if len(pas_bon) == 0 :
        return False
    else :
        return pas_bon

def comparaison(elt, graphe1, graphe2, elt_1, elt_2, elt_3, couples_possibles, num, graphe_commun, compt) :
    deja_vus_voisin_1 = []
    deja_vus_voisin_2 = []
    sommets_en_plus_1 = []
    sommets_en_plus_2 = []
    tab_temp = []
    
#     succ_1 = True
#     if len(elt) > 1 :
#         if elt[1][0] in graphe1.successors(elt[0][0]) :
#             succ_1 = True
#         else :
#             succ_1 = False
#             
#     succ_2 = True
#     if len(elt) > 1 :
#         if elt[1][1] in graphe1.successors(elt[0][1]) :
#             succ_2 = True
#         else :
#             succ_2 = False
    if num == 2 or num == 3 :
        succ_1 = True
        succ_2 = True
    else :
        succ_1 = False
        succ_2 = False

    #print()
    #print(succ_1)
    #print(succ_2)
        
    for i in range(0, len(elt)) :
        #print(elt_1[i])
    
        if dans_graphe(graphe_commun, elt[i]) == False :
            if elt[i] not in graphe_commun.nodes() :
                graphe_commun.add_node(elt[i])  
        
        if succ_1 == True :
            voisins_1 = graphe1.successors(elt[i][0])
        else :
            voisins_1 = graphe1.predecessors(elt[i][0])            
        for voisin_1 in voisins_1 :
                #print(voisin_1)
                if succ_2 == True :              
                    voisins_2 = graphe2.successors(elt[i][1])
                else :
                    voisins_2 = graphe2.predecessors(elt[i][1])
                
                
                for voisin_2 in voisins_2 :
                    #print("elt_1")
                    #print(elt_1[i])
                    #print(voisin_1)
                    #print(voisin_2)
                    
                    if succ_1 == True :
                        dict_voisin_1 = graphe1[elt[i][0]][voisin_1]
                    else :
                        dict_voisin_1 = graphe1[voisin_1][elt[i][0]]
                    
                    if succ_2 == True :
                        dict_voisin_2 = graphe2[elt[i][1]][voisin_2]
                    else :
                        dict_voisin_2 = graphe2[voisin_2][elt[i][1]]
                    #print(dict_voisin_1)
                    #print(dict_voisin_2)
                    
                    for edge_1 in dict_voisin_1 :
                        for edge_2 in dict_voisin_2 :
                            if dict_voisin_1[edge_1]["label"] == 'CWW' :
                                label_1 = 'CAN'
                            elif dict_voisin_1[edge_1]["label"] == 'B53' : 
                                label_1 = 'COV'
                            else :
                                label_1 = 'NON_CAN'
                            
                            if dict_voisin_2[edge_2]["label"] == 'CWW' :
                                label_2 = 'CAN'
                            elif dict_voisin_2[edge_2]["label"] == 'B53' : 
                                label_2 = 'COV'
                            else :
                                label_2 = 'NON_CAN'
                            
                            bonne_chaine = False
                            for ch in graphe1.nodes[voisin_1]["chaine"] :
                                if ch in graphe2.nodes[voisin_2]["chaine"] :
                                    bonne_chaine = True
                            
                            if bonne_chaine :
                                    
                                if label_1 == label_2 and dict_voisin_1[edge_1]["long_range"] == dict_voisin_2[edge_2]["long_range"] :
                                    if label_1 == 'COV' :
                                        #if (graphe1.nodes[elt[i][0]]["position"][0] - graphe1.nodes[voisin_1]["position"][0] < 0 and graphe1.nodes[1]["position"][0] - graphe1.nodes[3]["position"][0] < 0) or (graphe1.nodes[elt[i][0]]["position"][0] - graphe1.nodes[voisin_1]["position"][0] > 0 and graphe1.nodes[1]["position"][0] - graphe1.nodes[3]["position"][0] > 0) : 
                                        #    if (graphe2.nodes[elt[i][1]]["position"][0] - graphe2.nodes[voisin_2]["position"][0] < 0 and graphe2.nodes[1]["position"][0] - graphe2.nodes[3]["position"][0] < 0) or (graphe2.nodes[elt[i][1]]["position"][0] - graphe2.nodes[voisin_2]["position"][0] > 0 and graphe2.nodes[1]["position"][0] - graphe2.nodes[3]["position"][0] > 0) : 
                                        if (succ_1 == True and succ_2 == True) or (succ_1 == False and succ_2 == False) :# dans le meme sens dans les deux graphes
                                                if voisin_1 == num and voisin_2 == num : # elt courant est lie au motif
                                                    if dans_graphe(graphe_commun, elt[i]) == False :
                                                        if elt[i] not in graphe_commun.nodes() :
                                                            graphe_commun.add_node(elt[i])
                                                        deja_ajoute = False
                                                        for voisin in graphe_commun[elt[i]] :
                                                            if voisin == (voisin_1, voisin_2) :
                                                                for edge in graphe_commun[elt[i]][voisin] :
                                                                    if graphe_commun[elt[i]][voisin][edge]["type"] == label_1 :
                                                                        deja_ajoute = True
                                                        if deja_ajoute == False :
                                                            graphe_commun.add_edge((num,num),elt[i], type='COV', long_range=False)
                                                    #nb_aretes += 1
                                                    #print("ramou")
                                                else : #elt courant non lie au motif
                                                    for couple in elt : # on cherche si un couple deja dans la liste ne correspond pas aux voisins
                                                        if dans_graphe(graphe_commun, elt[i]) == False and voisin_1 == couple[0] and voisin_2 == couple[1]  and elt[i][0] not in deja_vus_voisin_1 and elt[i][1] not in deja_vus_voisin_2 :
                                                            deja_vus_voisin_1.append(voisin_1)
                                                            deja_vus_voisin_2.append(voisin_2)
    #                                                         nb_aretes += 1
                                                            if elt[i] not in graphe_commun.nodes() :
                                                                graphe_commun.add_node(elt[i])
                                                            deja_ajoute = False
                                                            for voisin in graphe_commun[elt[i]] :
                                                                if voisin == (voisin_1, voisin_2) :
                                                                    for edge in graphe_commun[elt[i]][voisin] :
                                                                        if graphe_commun[elt[i]][voisin][edge]["type"] == label_1 :
                                                                            deja_ajoute = True
                                                            if deja_ajoute == False :
                                                                graphe_commun.add_edge(elt[i], (voisin_1, voisin_2), type='COV', long_range=False)
                                                            #print("petit rat1")
                                                            #print("ramou2")   
                                    else : ## lies par une liaison non covalente
                                        ## a ce moment-la on regarde si les sommets voisins pourraient se superposer
                                        if graphe1.nodes[voisin_1]["type"] == graphe2.nodes[voisin_2]["type"] :
                                            #print("ramousnif")
                                            pas_vu_1_motif = True
                                            pas_vu_2_motif = True
                                            for j in range(1,6) : ## on regarde si l'un des voisins n'est pas dans le motif
                                                if voisin_1 == j and voisin_2 == j : #si les deux sont dans le motif c'est bon
                                                    #nb_aretes += 1
                                                    if dans_graphe(graphe_commun, elt[i]) == False :
                                                        if elt[i] not in graphe_commun.nodes() :
                                                            graphe_commun.add_node(elt[i])
                                                        deja_ajoute = False
                                                        for voisin in graphe_commun[elt[i]] :
                                                            if voisin == (voisin_1, voisin_2) :
                                                                for edge in graphe_commun[elt[i]][voisin] :
                                                                    if graphe_commun[elt[i]][voisin][edge]["type"] == label_1 :
                                                                        deja_ajoute = True
                                                        if deja_ajoute == False :
                                                            graphe_commun.add_edge(elt[i], (voisin_1, voisin_2), type=label_1, long_range = dict_voisin_1[edge_1]["long_range"])
                                                elif voisin_1 == j :
                                                    pas_vu_1_motif = False
                                                elif voisin_2 == j :
                                                    pas_vu_2_motif = False
                #                             if pas_vu_1 and pas_vu_2 and elt[i][0] not in deja_vus_voisin_1 and elt[i][1] not in deja_vus_voisin_2 and voisin_1 not in sommets_en_plus_1 and voisin_2 not in sommets_en_plus_2 :
                #                                 print("ajout arete")
                #                                 print(voisin_1)
                #                                 print(voisin_2)
                #                                 deja_vus_voisin_1.append(voisin_1)
                #                                 deja_vus_voisin_2.append(voisin_2)
                #                                 deja_vus_voisin_1.append(elt[i][0])
                #                                 deja_vus_voisin_2.append(elt[i][1])
                #                                 nb_aretes += 1
                                            pas_vu_1 = True
                                            pas_vu_2 = True
                                            for couple in elt :
                                                if voisin_1 == couple[0] :
                                                    pas_vu_1 = False
                                                if voisin_2 == couple[1] : 
                                                    pas_vu_2 = False 
                                            for couple in elt_1 :
                                                if voisin_1 == couple[0] :
                                                    pas_vu_1 = False
                                                if voisin_2 == couple[1] : 
                                                    pas_vu_2 = False 
                                            for couple in elt_2 :
                                                if voisin_1 == couple[0] :
                                                    pas_vu_1 = False
                                                if voisin_2 == couple[1] : 
                                                    pas_vu_2 = False 
                                            for couple in elt_3 :
                                                if voisin_1 == couple[0] :
                                                    pas_vu_1 = False
                                                if voisin_2 == couple[1] : 
                                                    pas_vu_2 = False 
    #                                         if voisin_1 == 9 and voisin_2 == 9 :
    #                                             print(graphe_commun.nodes.data())
    #                                             print(dans_graphe(graphe_commun, elt[i]))
                                            if dans_graphe(graphe_commun, elt[i]) == False and dans_graphe(graphe_commun, (voisin_1, voisin_2)) == False :
                                                
                                                if pas_vu_1_motif and pas_vu_2_motif and elt[i][0] not in deja_vus_voisin_1 and elt[i][1] not in deja_vus_voisin_2 and voisin_1 not in sommets_en_plus_1 and voisin_2 not in sommets_en_plus_2 : ## on ajoute un sommet si le sommet n'existe nulle part
                                                    if pas_vu_1 and pas_vu_2 and test_compatibilite_grands_graphes(graphe_commun, graphe1, graphe2) : #les sommets n'existent pas, on les rajoute
                                                        #nb_sommets += 1
                        #                                                     print("petit ramousnif")
                        #                                                     print(voisin_1)
                        #                                                     print(voisin_2)
                                                        #print("petit rat")
                                                        #print(voisin_1)
                                                        #print(voisin_2)
                                                        
                                                        sommets_en_plus_1.append(voisin_1)
                                                        sommets_en_plus_2.append(voisin_2)
                                                        
                                                        graphe_commun.add_node((voisin_1, voisin_2))
                                                        if (voisin_1, voisin_2) != (5,5) :
                                                            tab_temp.append((voisin_1, voisin_2))
                                                    if elt[i] not in graphe_commun.nodes() :
                                                        graphe_commun.add_node(elt[i])
                                                    deja_ajoute = False
                                                    for voisin in graphe_commun[elt[i]] :
                                                        if voisin == (voisin_1, voisin_2) :
                                                            for edge in graphe_commun[elt[i]][voisin] :
                                                                if graphe_commun[elt[i]][voisin][edge]["type"] == label_1 :
                                                                    deja_ajoute = True
                                                    if deja_ajoute == False :
                                                        graphe_commun.add_edge(elt[i], (voisin_1, voisin_2), type=label_1, long_range = dict_voisin_1[edge_1]["long_range"])   
                                        #couples_possibles[0].append(tab_temp) 
                #print(nb_aretes)
                #print(nb_sommets)
    couples_en_plus = []
    couples_en_plus.extend(elt)
    couples_en_plus.extend(tab_temp)
    #print(couples_en_plus)
    #print(couples_possibles[0])
    max_deja = 0
    for groupe in couples_possibles[num-1] :
        deja = 0
        for couple in groupe :
            for couple_2 in couples_en_plus :
                if couple[0] == couple_2[0] and couple[1] == couple_2[1] :
                    deja += 1
        if max_deja < deja :
            max_deja = deja
        #print(deja)
        
    if max_deja < len(couples_en_plus) :
        couples_possibles[num-1].append(couples_en_plus)
        del(couples_possibles[num-1][compt])
    
#     deja_vus_voisin_1 = list(set(deja_vus_voisin_1))
#     tab_sommets_aretes =[]
#     tab_sommets_aretes.append(nb_sommets + len(deja_vus_voisin_1))
#     tab_sommets_aretes.append(nb_aretes)
    #print("ramousnif")
    #print(graphe_commun.nodes.data())
    #print(graphe_commun.edges.data())
    #print(graphe_commun.nodes.data())
    return graphe_commun  


def comparaison_chaine_par_chaine(elt, graphe1, graphe2, num, graphe_commun, chaines_1, chaines_2, occ_a_minor_1, occ_a_minor_2, graphe_commun_max) :
    deja_vus_voisin_1 = []
    deja_vus_voisin_2 = []
    sommets_en_plus_1 = []
    sommets_en_plus_2 = []
    
#     succ_1 = True
#     if len(elt) > 1 :
#         if elt[1][0] in graphe1.successors(elt[0][0]) :
#             succ_1 = True
#         else :
#             succ_1 = False
#             
#     succ_2 = True
#     if len(elt) > 1 :
#         if elt[1][1] in graphe1.successors(elt[0][1]) :
#             succ_2 = True
#         else :
#             succ_2 = False
    if num == 2 or num == 3 :
        succ_1 = True
        succ_2 = True
    else :
        succ_1 = False
        succ_2 = False

    #print()
    #print(succ_1)
    #print(succ_2)
        
    for i in range(0, len(elt)) :
        #print(elt_1[i])
          
        if succ_1 == True :
            voisins_1 = graphe1.successors(elt[i][0])
        else :
            voisins_1 = graphe1.predecessors(elt[i][0])
             
        
#                    
        #print("elt : " + str(elt[i]))
        for voisin_1 in voisins_1 :
            #print(voisin_1)
            if succ_2 == True :              
                voisins_2 = graphe2.successors(elt[i][1])
            else :
                voisins_2 = graphe2.predecessors(elt[i][1]) 
            for voisin_2 in voisins_2 :
                #print("elt_1")
                #print(elt_1[i])
                #print(voisin_1)
                #print(voisin_2)
                
                if succ_1 == True :
                    dict_voisin_1 = graphe1[elt[i][0]][voisin_1]
                else :
                    dict_voisin_1 = graphe1[voisin_1][elt[i][0]]
                
                if succ_2 == True :
                    dict_voisin_2 = graphe2[elt[i][1]][voisin_2]
                else :
                    dict_voisin_2 = graphe2[voisin_2][elt[i][1]]
                #print(dict_voisin_1)
                #print(dict_voisin_2)
                
                for edge_1 in dict_voisin_1 :
                    for edge_2 in dict_voisin_2 :
                        if dict_voisin_1[edge_1]["label"] == 'CWW' :
                            label_1 = 'CAN'
                        elif dict_voisin_1[edge_1]["label"] == 'B53' : 
                            label_1 = 'COV'
                        else :
                            label_1 = 'NON_CAN'
                        
                        if dict_voisin_2[edge_2]["label"] == 'CWW' :
                            label_2 = 'CAN'
                        elif dict_voisin_2[edge_2]["label"] == 'B53' : 
                            label_2 = 'COV'
                        else :
                            label_2 = 'NON_CAN'
                        
                        bonne_chaine = False
                        for ch in graphe1.nodes[voisin_1]["chaine"] :
                            if ch in graphe2.nodes[voisin_2]["chaine"] :
                                bonne_chaine = True
                        
                        if bonne_chaine :
                                
                            if label_1 == label_2 :
                                if label_1 == 'COV' :
                                    
                                    #if (graphe1.nodes[elt[i][0]]["position"][0] - graphe1.nodes[voisin_1]["position"][0] < 0 and graphe1.nodes[1]["position"][0] - graphe1.nodes[3]["position"][0] < 0) or (graphe1.nodes[elt[i][0]]["position"][0] - graphe1.nodes[voisin_1]["position"][0] > 0 and graphe1.nodes[1]["position"][0] - graphe1.nodes[3]["position"][0] > 0) : 
                                    #    if (graphe2.nodes[elt[i][1]]["position"][0] - graphe2.nodes[voisin_2]["position"][0] < 0 and graphe2.nodes[1]["position"][0] - graphe2.nodes[3]["position"][0] < 0) or (graphe2.nodes[elt[i][1]]["position"][0] - graphe2.nodes[voisin_2]["position"][0] > 0 and graphe2.nodes[1]["position"][0] - graphe2.nodes[3]["position"][0] > 0) : 
                                    if (succ_1 == True and succ_2 == True) or (succ_1 == False and succ_2 == False) :# dans le meme sens dans les deux graphes
                                        if dans_graphe(graphe_commun, elt[i]) == False and dans_graphe(graphe_commun, (voisin_1, voisin_2)) == False and dans_graphe(graphe_commun_max, elt[i]) == False and dans_graphe(graphe_commun_max, (voisin_1, voisin_2)) == False:
                                            
                                            if voisin_1 == occ_a_minor_1[num-1] and voisin_2 == occ_a_minor_2[num-1] : # elt courant est lie au motif
                                                    
                                                    if elt[i] not in graphe_commun.nodes() :
                                                        graphe_commun.add_node(elt[i])
                                                    deja_ajoute = False
                                                    for voisin in graphe_commun[elt[i]] :
                                                        if voisin == (voisin_1, voisin_2) :
                                                            for edge in graphe_commun[elt[i]][voisin] :
                                                                if graphe_commun[elt[i]][voisin][edge]["type"] == label_1 :
                                                                    deja_ajoute = True
                                                    if deja_ajoute == False :
                                                        graphe_commun.add_edge((occ_a_minor_1[num-1],occ_a_minor_2[num-1]),elt[i], type='COV')
                                                #nb_aretes += 1
                                                #print("ramou")
                                            else : #elt courant non lie au motif
                                                for couple in elt : # on cherche si un couple deja dans la liste ne correspond pas aux voisins
                                                    if voisin_1 == couple[0] and voisin_2 == couple[1]  and elt[i][0] not in deja_vus_voisin_1 and elt[i][1] not in deja_vus_voisin_2 :
                                                        deja_vus_voisin_1.append(voisin_1)
                                                        deja_vus_voisin_2.append(voisin_2)
#                                                         nb_aretes += 1
                                                        if elt[i] not in graphe_commun.nodes() :
                                                            graphe_commun.add_node(elt[i])
                                                        deja_ajoute = False
                                                        for voisin in graphe_commun[elt[i]] :
                                                            if voisin == (voisin_1, voisin_2) :
                                                                for edge in graphe_commun[elt[i]][voisin] :
                                                                    if graphe_commun[elt[i]][voisin][edge]["type"] == label_1 :
                                                                        deja_ajoute = True
                                                        if deja_ajoute == False :
                                                            graphe_commun.add_edge(elt[i], (voisin_1, voisin_2), type='COV')
                                                        #print("petit rat1")
                                                        #print("ramou2")   
                                else : ## lies par une liaison non covalente
                                    if graphe1.nodes[voisin_1]["type"] == graphe2.nodes[voisin_2]["type"] :
                                        pas_vu_1_motif = True
                                        pas_vu_2_motif = True

                                        for j in range(5) : ## on regarde si l'un des voisins n'est pas dans le motif
                                            if voisin_1 == occ_a_minor_1[j] and voisin_2 == occ_a_minor_2[j] : #si les deux sont dans le motif c'est bon
                                                #nb_aretes += 1
                                                if dans_graphe(graphe_commun, elt[i]) == False and dans_graphe(graphe_commun, (voisin_1, voisin_2)) and dans_graphe(graphe_commun_max, elt[i]) == False and dans_graphe(graphe_commun_max, (voisin_1, voisin_2)) == False :
                                                    if elt[i] not in graphe_commun.nodes() :
                                                        graphe_commun.add_node(elt[i])
                                                    deja_ajoute = False
                                                    for voisin in graphe_commun[elt[i]] :
                                                        if voisin == (voisin_1, voisin_2) :
                                                            for edge in graphe_commun[elt[i]][voisin] :
                                                                if graphe_commun[elt[i]][voisin][edge]["type"] == label_1 :
                                                                    deja_ajoute = True
                                                    if deja_ajoute == False :
                                                        graphe_commun.add_edge(elt[i], (voisin_1, voisin_2), type=label_1)
                                            elif voisin_1 == occ_a_minor_1[j] :
                                                pas_vu_1_motif = False
                                            elif voisin_2 == occ_a_minor_2[j] :
                                                pas_vu_2_motif = False
            #                             if pas_vu_1 and pas_vu_2 and elt[i][0] not in deja_vus_voisin_1 and elt[i][1] not in deja_vus_voisin_2 and voisin_1 not in sommets_en_plus_1 and voisin_2 not in sommets_en_plus_2 :
            #                                 print("ajout arete")
            #                                 print(voisin_1)
            #                                 print(voisin_2)
            #                                 deja_vus_voisin_1.append(voisin_1)
            #                                 deja_vus_voisin_2.append(voisin_2)
            #                                 deja_vus_voisin_1.append(elt[i][0])
            #                                 deja_vus_voisin_2.append(elt[i][1])
            #                                 nb_aretes += 1
#                                         pas_vu_1 = True
#                                         pas_vu_2 = True
#                                         for k in range(4) :
#                                             if k != num :
#                                                 for e in chaines_1[k] :
#                                                     if voisin_1 == e :
#                                                         pas_vu_1 = False
#                                                 for e in chaines_2[k] :
#                                                     if voisin_2 == e :
#                                                         pas_vu_2 = False
#                                         if voisin_1 == 9 and voisin_2 == 9 :
#                                             print(graphe_commun.nodes.data())
#                                             print(dans_graphe(graphe_commun, elt[i]))
                                        if dans_graphe(graphe_commun, elt[i]) == False and dans_graphe(graphe_commun, (voisin_1, voisin_2)) == False and dans_graphe(graphe_commun_max, elt[i]) == False and dans_graphe(graphe_commun_max, (voisin_1, voisin_2)) == False :
                                            
                                            
                                            if pas_vu_1_motif and pas_vu_2_motif and elt[i][0] not in deja_vus_voisin_1 and elt[i][1] not in deja_vus_voisin_2 and voisin_1 not in sommets_en_plus_1 and voisin_2 not in sommets_en_plus_2 and test_compatibilite_grands_graphes(graphe_commun, graphe1, graphe2): ## on ajoute un sommet si le sommet n'existe nulle part
                                                if (voisin_1 not in chaines_1[0] and voisin_1 not in chaines_1[1] and voisin_1 not in chaines_1[2] and voisin_1 not in chaines_1[3] and voisin_2 not in chaines_2[0] and voisin_2 not in chaines_2[1] and voisin_2 not in chaines_2[2] and voisin_2 not in chaines_2[3]) or (voisin_1 in chaines_1[0] and voisin_2 in chaines_2[0] or voisin_1 in chaines_1[1] and voisin_2 in chaines_2[1] or voisin_1 in chaines_1[2] and voisin_2 in chaines_2[2] or voisin_1 in chaines_1[3] and voisin_2 in chaines_2[3]) :
#                                                 if pas_vu_1 and pas_vu_2 : #les sommets n'existent pas, on les rajoute
                                                    #nb_sommets += 1
                    #                                                     print("petit ramousnif")
                    #                                                     print(voisin_1)
                    #                                                     print(voisin_2)
                                                    #print("petit rat")
                                                    #print(voisin_1)
                                                    #print(voisin_2)
                                                    
                                                    sommets_en_plus_1.append(voisin_1)
                                                    sommets_en_plus_2.append(voisin_2)
    #                                                 if (voisin_1, voisin_2) != (chaines_1[num][4],chaines_2[num][4]) :
    #                                                     tab_temp.append((voisin_1, voisin_2))
                                                    graphe_commun.add_node((voisin_1, voisin_2))
                                                    if elt[i] not in graphe_commun.nodes() :
                                                        graphe_commun.add_node(elt[i])
                                                    deja_ajoute = False
                                                    for voisin in graphe_commun[elt[i]] :
                                                        if voisin == (voisin_1, voisin_2) :
                                                            for edge in graphe_commun[elt[i]][voisin] :
                                                                if graphe_commun[elt[i]][voisin][edge]["type"] == label_1 :
                                                                    deja_ajoute = True
                                                    if deja_ajoute == False :
                                                        graphe_commun.add_edge(elt[i], (voisin_1, voisin_2), type=label_1)   
                                    #couples_possibles[0].append(tab_temp) 
                #print(nb_aretes)
                #print(nb_sommets)
#     couples_en_plus = []
#     couples_en_plus.extend(elt)
#     couples_en_plus.extend(tab_temp)
    #print(couples_en_plus)
    #print(couples_possibles[0])
#     max_deja = 0
#     for groupe in couples_possibles[num-1] :
#         deja = 0
#         for couple in groupe :
#             for couple_2 in couples_en_plus :
#                 if couple[0] == couple_2[0] and couple[1] == couple_2[1] :
#                     deja += 1
#         if max_deja < deja :
#             max_deja = deja
#         #print(deja)
#         
#     if max_deja < len(couples_en_plus) :
#         couples_possibles[num-1].append(couples_en_plus) 
    
#     deja_vus_voisin_1 = list(set(deja_vus_voisin_1))
#     tab_sommets_aretes =[]
#     tab_sommets_aretes.append(nb_sommets + len(deja_vus_voisin_1))
#     tab_sommets_aretes.append(nb_aretes)
    #print("ramousnif")
    #print(graphe_commun.nodes.data())
    #print(graphe_commun.edges.data())
    #print(graphe_commun.nodes.data())
    return graphe_commun  

def nombre_aretes_arcs(graphe1, graphe2):
    compteur_arc = 0
    for (u, v, keys, t) in graphe1.edges(data="label", keys = True) :
        if t == "B53" :
            compteur_arc += 1
    compteur_arc_arete_1 = compteur_arc + (graphe1.number_of_edges() - compteur_arc)/2
    
    compteur_arc = 0
    for (u, v, keys, t) in graphe2.edges(data="label", keys = True) :
        if t == "B53" :
            compteur_arc += 1
    compteur_arc_arete_2 = compteur_arc + (graphe2.number_of_edges() - compteur_arc)/2
    return compteur_arc_arete_1, compteur_arc_arete_2
    
def recup_chaines(graphe, occ_a_minor):               
    chaines = [[occ_a_minor[0]]]
    i = 1
    for elt in occ_a_minor :
        if i <= 4 :
            compteur = elt
            if i != 1 : chaines.append([elt])
            liaison_B53 = True
            while liaison_B53 :
                liaison_B53 = False
                temp = compteur
                for voisin in graphe.successors(compteur) :
                    for arc in graphe[compteur][voisin] :
                        if voisin not in occ_a_minor[:4] and voisin not in chaines[len(chaines)-1] and graphe[compteur][voisin][arc]["label"] == 'B53' :
                            liaison_B53 = True
                            temp = voisin
                            chaines[len(chaines)-1].append(voisin)
                            
                for voisin in graphe.predecessors(compteur) :
                    for arc in graphe[voisin][compteur] :
                        if voisin not in occ_a_minor[:4] and voisin not in chaines[len(chaines)-1] and graphe[voisin][compteur][arc]["label"] == 'B53' :
                            liaison_B53 = True
                            temp = voisin
                            chaines[len(chaines)-1].append(voisin)
                compteur = temp
        i = i +1
    return chaines

def recup_couples_possibles(liste_epuree, chaines_1, chaines_2):
    couples_possibles = [[], [], [], []]
    for i in range(4) :
        for groupe in liste_epuree :
            couples_temp = [(chaines_1[i][0], chaines_2[i][0])]
            for elt in groupe :
                if len(chaines_1[i]) > elt[0]-1 and len(chaines_2[i]) > elt[1]-1 :
                    couples_temp.append((chaines_1[i][elt[0]-1], chaines_2[i][elt[1]-1]))
            couples_possibles[i].append(list(couples_temp))
            
            couples_temp = [(chaines_1[i][0], chaines_2[i][0])]
            for elt in groupe :
                if len(chaines_1[i]) > elt[1]-1 and len(chaines_2[i]) > elt[0]-1 :
                    couples_temp.append((chaines_1[i][elt[1]-1], chaines_2[i][elt[0]-1]))
            couples_possibles[i].append(list(couples_temp))
    return couples_possibles

def donne_graphe_motif(occ_a_minor_1, occ_a_minor_2):   
    graphe_motif = nx.MultiGraph()
    for i in range(5) :
        graphe_motif.add_node((occ_a_minor_1[i],occ_a_minor_2[i]))
        
    graphe_motif.add_edge((occ_a_minor_1[0],occ_a_minor_2[0]),(occ_a_minor_1[1],occ_a_minor_2[1]), type="NON_CAN")
    graphe_motif.add_edge((occ_a_minor_1[0],occ_a_minor_2[0]),(occ_a_minor_1[2],occ_a_minor_2[2]), type="COV")
    graphe_motif.add_edge((occ_a_minor_1[0],occ_a_minor_2[0]),(occ_a_minor_1[4],occ_a_minor_2[4]), type="NON_CAN")
    graphe_motif.add_edge((occ_a_minor_1[1],occ_a_minor_2[1]),(occ_a_minor_1[3],occ_a_minor_2[3]), type="COV")
    graphe_motif.add_edge((occ_a_minor_1[1],occ_a_minor_2[1]),(occ_a_minor_1[4],occ_a_minor_2[4]), type="CAN")
    graphe_motif.add_edge((occ_a_minor_1[2],occ_a_minor_2[2]),(occ_a_minor_1[3],occ_a_minor_2[3]), type="NON_CAN")
    
    return graphe_motif

def ajout_attribut_chaine(graphe, chaines, graphe_motif, num):
#     print(graphe.nodes.data())
    nx.set_node_attributes(graphe, -1, "chaine")
#     print(graphe.nodes.data())
#     print(chaines)
#     print(str([chaines[0][0], chaines[1][0], chaines[2][0], chaines[3][0], occ_a_minor[4]]))
    edges = []
    for edge in graphe_motif.edges() :
        edges.append((edge[0][num], edge[1][num])) 
      
    for i in range(4) :
        for e in chaines[i] :
#             print(e)
            if graphe.nodes[e]["chaine"] == -1 :
                graphe.nodes[e]["chaine"] = []
            if i+1 not in graphe.nodes[e]["chaine"] :

                graphe.nodes[e]["chaine"].append(i+1)

#             print("voisins")
            for voisin in graphe[e] :
#                 print(voisin)
                if graphe.nodes[voisin]["chaine"] == -1 :
                    graphe.nodes[voisin]["chaine"] = []
                if i+1 not in graphe.nodes[voisin]["chaine"] and (e,voisin) not in edges and (voisin,e) not in edges :
                    graphe.nodes[voisin]["chaine"].append(i+1)
    
    for noeud in graphe.nodes() :
        if graphe.nodes[noeud]["chaine"] == -1 :
            print(noeud)
#         print(graphe.nodes.data())
                
    return graphe

def donne_arcs_aretes_par_chaine(graphe, i):
    compteur_arc = 0
    compteur_arete = 0
    for (u, v, keys, t) in graphe.edges(data="label", keys = True) :
        if i in graphe.nodes[u]["chaine"] and i in graphe.nodes[v]["chaine"] :
            if t == "B53" :
                compteur_arc += 1
            else :
                compteur_arete += 1
    compteur_arc_arete = compteur_arc + compteur_arete/2 - 6

    return compteur_arc_arete

def donne_aretes_par_chaine(graphe, i):
    compteur_arete = 0
    for (u, v, keys, t) in graphe.edges(data="label", keys = True) :
        if i in graphe.nodes[u]["chaine"] and i in graphe.nodes[v]["chaine"] :
            if t != 'B53' :
                compteur_arete += 1
    compteur_arete = compteur_arete/2
    
    return compteur_arete

def donne_sommets_par_chaine(graphe, occ_a_minor, i):
    compteur_sommets = 0
#     print(graphe1.nodes.data())
    for noeud, data in graphe.nodes(data=True) :
        if i in data["chaine"] and noeud not in occ_a_minor :
            compteur_sommets += 1
    
    return compteur_sommets


def donne_sommets_par_chaine_commun(graphe1, graphe2, graphe_commun, occ_a_minor_1, occ_a_minor_2,  i):
    compteur_sommets = 0
#     print(graphe1.nodes.data())
    for noeud in graphe_commun.nodes() :
        if i in graphe1.nodes[noeud[0]]["chaine"] and i in graphe2.nodes[noeud[1]]["chaine"] and noeud[0] not in occ_a_minor_1 and noeud[1] not in occ_a_minor_2  :
            compteur_sommets += 1
    
    return compteur_sommets

def donne_aretes_par_chaine_commun(graphe1, graphe2, graphe_commun, i):
    compteur_arete = 0
    for (u, v, t) in graphe_commun.edges(data="type") :
        if i in graphe1.nodes[u[0]]["chaine"] and i in graphe2.nodes[u[1]]["chaine"] and i in graphe1.nodes[v[0]]["chaine"] and i in graphe2.nodes[v[1]]["chaine"] :
            if t != 'COV' :
                compteur_arete += 1
    
    return compteur_arete

def sous_graphe_commun_max_grands_graphes(couple):
    dico_comp =  {}
#     dico_sim = {}
    with open("liste_combi_9.pickle", 'rb') as fichier_liste :
        mon_depickler_liste = pickle.Unpickler(fichier_liste)
        liste_epuree = mon_depickler_liste.load()
          
        with open("grands_graphes.pickle", 'rb') as fichier :
            mon_depickler = pickle.Unpickler(fichier)
            dico_graphes = mon_depickler.load()
              
            with open("a-minor_test2.pickle", 'rb') as fichier_pickle :
                mon_depickler_aminor = pickle.Unpickler(fichier_pickle)
                tab_aminor = mon_depickler_aminor.load()
#                         dico_comp = {}
#                         dico_sim = {}
#                 tab_deja_fait = []
                
                if couple in dico_graphes.keys() :
                    element_1 = couple[0]
                    element_2 = couple[1]
                else :
                    element_2 = couple[0]
                    element_1 = couple[1]   
                print(element_1)
                print(element_2)

#                         for cle_1 in dico_graphes.keys() :
#                             element_1 = cle_1
                    #liste_a_faire = [('4V9F', '0',48,21), ('4V9F', '0',30,23), ('1FJG', 'A', 109, 6), ('5J7L', 'DA', 197, 4), ('5FDU', '1A', 197, 3), ('5DM6', 'X', 48, 9)]
#                     element_1 = ('1FJG', 'A', 48,8)
                if "graphe_comp_"+str(element_1) + "_" + str(element_2)+".pickle" not in os.listdir("result_graphes_comp_grands_graphes") : 
                    #print(element_1)
                    print(element_1)
                    print(element_2)
                    for occ in tab_aminor :
                        if occ["num_PDB"] == element_1[0] and occ["num_ch"] == element_1[1] and occ["num_motif"] == element_1[2] and occ["num_occ"] == element_1[3] :
                            graphe1 = dico_graphes[element_1]
                            chaines_1 = recup_chaines(graphe1, occ["a_minor"])

                                
#                                     
#                                 print(graphe1.nodes.data())
#                                     for cle_2 in dico_graphes.keys() :
#                                         element_2 = cle_2
                            #element_2  = ('5J5B', 'BA', 48,23) 
                            #print(element_2)
#                                 if element_1 != element_2 and (element_1, element_2) not in tab_deja_fait and (element_2,element_1) not in tab_deja_fait :
#                                 tab_deja_fait.append((element_1, element_2))
                            for occ_2 in tab_aminor : 
                                if occ_2["num_PDB"] == element_2[0] and occ_2["num_ch"] == element_2[1] and occ_2["num_motif"] == element_2[2] and occ_2["num_occ"] == element_2[3] :
                                    graphe2 = dico_graphes[element_2]
                                    chaines_2 = recup_chaines(graphe2, occ_2["a_minor"])
                                    print(chaines_1)
                                    print(chaines_2)
                                    
                                    graphe_motif = donne_graphe_motif(occ["a_minor"], occ_2["a_minor"])
                                    
                                    if len(nx.get_node_attributes(graphe2, "chaine")) == 0 :
                                        graphe2 = ajout_attribut_chaine(graphe2, chaines_2, graphe_motif, 1)
                                    
                                    if len(nx.get_node_attributes(graphe1, "chaine")) == 0 :
                                        graphe1 = ajout_attribut_chaine(graphe1, chaines_1, graphe_motif, 0)
                                     
                                    couples_possibles = recherche_couples_grands_graphes(graphe1, graphe2, chaines_1, chaines_2)
#                                         for k in range(len(couples_possibles[3])) :
#                                             print(couples_possibles[3][k])

                                    if couples_possibles != None :
                                        nombre = 0
                                        for i in range(4) : 
                                            for elt in couples_possibles[i] :
                                                nombre += 1
                                                
                                        if nombre < 365 :
                                        #print(couples_possibles)
                                            memory_error = False
                                            for i in range(4) :
                                                if 'memory error' in couples_possibles[i] :
                                                    memory_error = True
                                                    break
           
                                                   
                                            if memory_error == False :
                                                    #print("ramousnif")
                                                   
                                                graphe_commun_max = nx.MultiGraph()

                                                maxi = 0
                                                couples_max = [[]]
                                                
      
#                                                         
                                                for compt_1 in range(max(len(couples_possibles[0]), 1)) :
                                                    if len(couples_possibles[0]) > 0 :
                                                        elt_1 = couples_possibles[0][compt_1]
                                                    else :
                                                        elt_1 = []
                                                    for compt_2 in range(max(len(couples_possibles[1]), 1)) :
                                                        if len(couples_possibles[1]) > 0 :
                                                            elt_2 = couples_possibles[1][compt_2]
                                                        else :
                                                            elt_2 = []
                                                        for compt_3 in range(max(len(couples_possibles[2]), 1)) :
                                                            if len(couples_possibles[2]) > 0 :
                                                                elt_3 = couples_possibles[2][compt_3]
                                                            else :
                                                                elt_3 = [] 
                                                            for compt_4 in range(max(len(couples_possibles[3]), 1)) :
                                                                if len(couples_possibles[3]) > 0 :
                                                                    elt_4 = couples_possibles[3][compt_4]
                                                                else :
                                                                    elt_4 = []

                                                                graphe_commun = graphe_motif.copy()
                                                                 
                                                                graphe_commun = comparaison(elt_1, graphe1, graphe2, elt_2, elt_3, elt_4, couples_possibles, 1, graphe_commun, compt_1)
                                                                graphe_commun = comparaison(elt_2, graphe1, graphe2, elt_1, elt_3, elt_4, couples_possibles, 2, graphe_commun, compt_2)
                                                                
                                                                graphe_commun = comparaison(elt_3, graphe1, graphe2, elt_1, elt_2, elt_4, couples_possibles, 3, graphe_commun, compt_3)                        
                                                                graphe_commun = comparaison(elt_4, graphe1, graphe2, elt_1, elt_2, elt_3, couples_possibles, 4, graphe_commun, compt_4)
                                                                     
                                                                if test_compatibilite_grands_graphes(graphe_commun, graphe1, graphe2):
                                                                    sim = calcul_sommets_aretes_grands_graphes_commun(graphe_commun_max)/max(calcul_sommets_aretes_grands_graphes(graphe1), calcul_sommets_aretes_grands_graphes(graphe2))
                                                                    if maxi <= sim and sim <= 1.0 :
                                                                        maxi = sim
                                                                        del(couples_max[:])
                                                                        couples_max.append(elt_1)
                                                                        couples_max.append(elt_2)
                                                                        couples_max.append(elt_3)
                                                                        couples_max.append(elt_4)
                                                                        graphe_commun_max = graphe_commun.copy()
                                        else :
                                            memory_error = False
                                            for i in range(4) :
                                                if 'memory error' in couples_possibles[i] :
                                                    memory_error = True
                                                    break
           
                                                   
                                            if memory_error == False :
                                                    #print("ramousnif")
                                                   
                                                graphe_commun_max = nx.MultiGraph()
                                                #result = exectimeout(50, toutes_comparaisons, args=(couples_possibles, graphe1, graphe2, graphe_motif, graphe_commun_max))
                                                maxi = 0
                                                couples_max = [[]]
                                                
    #                                                         compteur_arc_arete_1, compteur_arc_arete_2 = nombre_aretes_arcs(graphe1, graphe2)        
    # #                                                         
    #                                                         for compt_1 in range(max(len(couples_possibles[0]), 1)) :
    #                                                             if len(couples_possibles[0]) > 0 :
    #                                                                 elt_1 = couples_possibles[0][compt_1]
    #                                                             else :
    #                                                                 elt_1 = []
    #                                                             for compt_2 in range(max(len(couples_possibles[1]), 1)) :
    #                                                                 if len(couples_possibles[1]) > 0 :
    #                                                                     elt_2 = couples_possibles[1][compt_2]
    #                                                                 else :
    #                                                                     elt_2 = []
    #                                                                 for compt_3 in range(max(len(couples_possibles[2]), 1)) :
    #                                                                     if len(couples_possibles[2]) > 0 :
    #                                                                         elt_3 = couples_possibles[2][compt_3]
    #                                                                     else :
    #                                                                         elt_3 = [] 
    #                                                                     for compt_4 in range(max(len(couples_possibles[3]), 1)) :
    #                                                                         if len(couples_possibles[3]) > 0 :
    #                                                                             elt_4 = couples_possibles[3][compt_4]
    #                                                                         else :
    #                                                                             elt_4 = []
    #                                                                         print(compt_1)
    #                                                                         print(compt_2)
    #                                                                         print(compt_3)
    #                                                                         print(compt_4)
    #                                                                         graphe_commun = graphe_motif.copy()
    #                                                                         
    #                                                                         graphe_commun = comparaison(elt_1, graphe1, graphe2, 1, graphe_commun)
    #                                                                         graphe_commun = comparaison(elt_2, graphe1, graphe2, 2, graphe_commun)
    #                                                                         #print("petit rat")
    #                                                                         #print(graphe_commun)
    #                                                                         graphe_commun = comparaison(elt_3, graphe1, graphe2, 3, graphe_commun)                            
    #                                                                         graphe_commun = comparaison(elt_4, graphe1, graphe2, 4, graphe_commun)
    #                                                                             
    #                                                                         
    #                                                                         sim = ((graphe_commun.number_of_nodes() + graphe_commun.number_of_edges())*(graphe_commun.number_of_nodes() + graphe_commun.number_of_edges()))/((graphe1.number_of_nodes()+compteur_arc_arete_1)*(graphe2.number_of_nodes()+compteur_arc_arete_2))
    #                                                                         if maxi <= sim and sim <= 1.0 :
    #                                                                             maxi = sim
    #                                                                             del(couples_max[:])
    #                                                                             couples_max.append(elt_1)
    #                                                                             couples_max.append(elt_2)
    #                                                                             couples_max.append(elt_3)
    #                                                                             couples_max.append(elt_4)
    #                                                                             graphe_commun_max = graphe_commun.copy()    

                                                for i in range(4) :

                                                    if i != 0 : couples_max.append([])
                                                    maxi = 0.0
                                                    graphe_commun_max_temp = nx.MultiGraph()
                                                      
                                                    compteur = 0
                                                    graphe_commun = nx.MultiGraph()
                                                    graphe_commun = graphe_motif.copy()
    
                                                    for elt in couples_possibles[i] :
    #                                                                 print(i)
    #                                                                 print(compteur)
                                                        compteur += 1
                                                        graphe_commun = comparaison_chaine_par_chaine(elt, graphe1, graphe2, i+1, graphe_commun, chaines_1, chaines_2, occ["a_minor"], occ_2["a_minor"], graphe_commun_max)
            #                                                                    graphe_commun.remove_edges_from([((1,1), (2,2)),((1,1), (3,3)),((1,1), (5,5)),((2,2), (4,4)),((2,2),(5,5)),((3,3),(4,4))])
            #                                                                     graphe_commun.remove_nodes_from([(1,1),(2,2),(3,3),(4,4),(5,5)])
            #                                                                     if i == 3 :
            #                                                                         print(compteur)
            #                                                                         print(graphe_commun.nodes.data())
                                                                
                                                            
                                                        if test_compatibilite_grands_graphes(graphe_commun, graphe1, graphe2) :
                                                            sim = (donne_sommets_par_chaine_commun(graphe1, graphe2, graphe_commun, occ["a_minor"], occ_2["a_minor"], i+1) + donne_aretes_par_chaine_commun(graphe1, graphe2, graphe_commun, i+1))/max(donne_sommets_par_chaine(graphe1, occ["a_minor"], i+1) + donne_aretes_par_chaine(graphe1, i+1), donne_sommets_par_chaine(graphe2, occ_2["a_minor"], i+1) + donne_aretes_par_chaine(graphe2, i+1))
                                                            if maxi <= sim and sim <= 1.0 :
                                                                maxi = sim
                                                                del(couples_max[len(couples_max)-1][:])
                                                                couples_max[len(couples_max)-1].append(elt)
                                                                graphe_commun_max_temp = graphe_commun.copy()
                                                             
    #                                                             print(graphe_commun_max_temp.number_of_nodes())
    #                                                             print(graphe_commun_max_temp.nodes())
    #                                                             print(graphe_commun_max_temp.number_of_edges())
    #                                                             print(compteur_sommets_1)
    #                                                             print(compteur_arc_arete_1_par_chaine)
    #                                                             print(compteur_sommets_2)
    #                                                             print(compteur_arc_arete_2_par_chaine)
    #                                                             print(i)     
    #                                                             print(maxi)    
                                                    print(graphe_commun_max_temp.nodes.data())   
    #                                                             a_enlever = []    
    #                                                             for noeud, attr in graphe_commun_max_temp.nodes(data=True) :
    #                                                                 pas_bon = dans_graphe(graphe_commun_max_temp, noeud)
    #                                                                 if pas_bon != False :
    #                                                                     print("probleme")
    #                                                                     pas_bon.append(noeud)
    #                                                                     max_pas_bon = -1
    #                                                                     pos_max  = -1
    #                                                                     compte = 0
    #                                                                     for e in pas_bon :
    #                                                                         if len(graphe_commun_max_temp[e]) > max_pas_bon :
    #                                                                             max_pas_bon = len(graphe_commun_max_temp[e])
    #                                                                             pos_mas = compte
    #                                                                         
    #                                                                         if e[0] in occ["a_minor"] or e[1] in occ_2["a_minor"] :
    #                                                                             pos_max = compte
    #                                                                             break
    #                                                                         
    #                                                                         compte += 1
    #                                                                         
    #                                                                     compte = 0
    #                                                                     
    #                                                                     for e in pas_bon :
    #                                                                         if compte != pos_max :
    # #                                                                             for edge in graphe_commun_max_temp.edges() :
    # #                                                                                 if edge[0] == e or edge[1] == e :
    # #                                                                                     if edge not in a_enlever :
    # #                                                                                         a_enlever.append(edge)
    #                                                                             if e not in a_enlever :
    #                                                                                 a_enlever.append(e)
    #                                                                             
    #                                                                         compte += 1
    # #                                                                     for edge in a_enlever :
    # #                                                                         graphe_commun_max_temp.remove_edge(*edge)
    #                                                                 
    #                                                             print(a_enlever)
    #                                                             for e in a_enlever :
    #                                                                 graphe_commun_max_temp.remove_node(e)
                                                    
                                                    for noeud, attr in graphe_commun_max_temp.nodes(data=True) :
                                                        if noeud not in graphe_commun_max.nodes() :
                                                            graphe_commun_max.add_node(noeud, **attr)
                                                                
                                                    for u, v, attr_arete in graphe_commun_max_temp.edges(data=True) :
                                                        if (u,v) not in graphe_commun_max.edges() :
                                                            graphe_commun_max.add_edge(u,v, **attr_arete)
                                                               
                                                maxi = calcul_sommets_aretes_grands_graphes_commun(graphe_commun_max)/max(calcul_sommets_aretes_grands_graphes(graphe1), calcul_sommets_aretes_grands_graphes(graphe2))
                                        
                                        dico_comp.update({(element_1, element_2) : graphe_commun_max.copy()})
#                                                                 dico_sim.update({(element_1, element_2) : maxi})  
                                            
                                        with open("result_graphes_comp_grands_graphes/graphe_comp_"+str(element_1) + "_" + str(element_2)+".pickle", 'wb') as fichier_graphes :
                                            mon_pickler_3 = pickle.Pickler(fichier_graphes)
                                            mon_pickler_3.dump(dico_comp) 
    
    
                                        #with open("fichier_similarite_grands_graphes_nouvelle_metrique.txt", 'a') as fichier_ecriture :    
                                            #fichier_ecriture.write(str(element_1) + '\n' + str(element_2) + '\n')
                                            #fichier_ecriture.write("Graphe : ")
                                            #fichier_ecriture.write(str(graphe_commun_max.nodes.data())+ '\n')
                                            #fichier_ecriture.write(str(graphe_commun_max.edges.data())+ '\n')
                                            #fichier_ecriture.write("Similarite :" + str(maxi))
                                            #fichier_ecriture.write('\n\n')                                                 


    
if __name__ == '__main__':
    
    liste_a_faire = []
    
    with open("grands_graphes.pickle", 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        dico_graphes = mon_depickler.load()
    
        for cle_1 in dico_graphes.keys() :
            for cle_2 in dico_graphes.keys() :
                if cle_1 != cle_2 and (cle_2, cle_1) not in liste_a_faire and (cle_1, cle_2) not in liste_a_faire :
                    liste_a_faire.append((cle_1, cle_2))

    
    print(len(liste_a_faire))
    
    p = multiprocessing.Pool(8)
    result = p.map(sous_graphe_commun_max_grands_graphes, liste_a_faire)
    p.close()
    p.join()
    
    
            
                                         
                                        
                                
                                    
                                    
                            
                                
                                
