'''
Created on 19 déc. 2018

@author: coline
'''

import pickle
import networkx as nx
import os
import time
import multiprocessing

def exectimeout(timeout, fonc, *args, **kwargs):
    """exécute fonc(args, kwargs), avec le temps maxi timeout"""
    pool = multiprocessing.Pool(processes=1) # crée 1 processus
    poolexec = pool.apply_async(func=fonc, *args, **kwargs) # eval. asynchrone
    try:
        poolexec.get(timeout=timeout)
        result = poolexec.get()
        pool.terminate()
        return result
    except multiprocessing.TimeoutError:
        pool.terminate()
        raise # renvoi de l'exception à l'appelant

def dans_graphe(graphe, couple_a_chercher):
    for noeud in graphe.nodes() :
        #print(noeud)
        if (noeud[0] == couple_a_chercher[0] and noeud[1] != couple_a_chercher[1]) or (noeud[1] == couple_a_chercher[1] and noeud[0] != couple_a_chercher[0])  :
            return True
    return False

def comparaison(elt, graphe1, graphe2, elt_1, elt_2, elt_3, couples_possibles, num, graphe_commun) :
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
                                                        graphe_commun.add_edge((num,num),elt[i], type='COV')
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
                                                            graphe_commun.add_edge(elt[i], (voisin_1, voisin_2), type='COV')
                                                        #print("petit rat1")
                                                        #print("ramou2")   
                                else : ## lies par une liaison non covalente
                                    ## a ce moment-la on regarde si les sommets voisins pourraient se superposer
                                    if graphe1.nodes[voisin_1]["type"] == graphe2.nodes[voisin_2]["type"] and abs(graphe1.nodes[voisin_1]["poids"] - graphe2.nodes[voisin_2]["poids"]) < 2 :
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
                                                        graphe_commun.add_edge(elt[i], (voisin_1, voisin_2), type=label_1)
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
                                                if pas_vu_1 and pas_vu_2 : #les sommets n'existent pas, on les rajoute
                                                    #nb_sommets += 1
                    #                                                     print("petit ramousnif")
                    #                                                     print(voisin_1)
                    #                                                     print(voisin_2)
                                                    #print("petit rat")
                                                    #print(voisin_1)
                                                    #print(voisin_2)
                                                    
                                                    sommets_en_plus_1.append(voisin_1)
                                                    sommets_en_plus_2.append(voisin_2)
                                                    if (voisin_1, voisin_2) != (5,5) :
                                                        tab_temp.append((voisin_1, voisin_2))
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
    
#     deja_vus_voisin_1 = list(set(deja_vus_voisin_1))
#     tab_sommets_aretes =[]
#     tab_sommets_aretes.append(nb_sommets + len(deja_vus_voisin_1))
#     tab_sommets_aretes.append(nb_aretes)
    #print("ramousnif")
    #print(graphe_commun.nodes.data())
    #print(graphe_commun.edges.data())
    #print(graphe_commun.nodes.data())
    return graphe_commun  


                                                                        
# for i in range(len(os.listdir("graphes_extension/fichiers_pickle/"))) :
#     for j in range(i+1, len(os.listdir("graphes_extension/fichiers_pickle/"))) :
#         element1 = os.listdir("graphes_extension/fichiers_pickle/")[i]
#         element2 = os.listdir("graphes_extension/fichiers_pickle/")[j]

def sous_graphe_commun_max(nom_extension):
    c = 0
    
    #if 'tab_calcules_sim_'+nom_extension+'.pickle' in os.listdir() :
    if 'tab_calcules_sim_qui_manquent.pickle' in os.listdir() :
        #with open("tab_calcules_sim_"+nom_extension+".pickle", 'rb') as fichier_deja_fait :
        with open("tab_calcules_sim_qui_manquent.pickle", 'rb') as fichier_deja_fait :
            mon_depickler_tab = pickle.Unpickler(fichier_deja_fait)
            tab_deja_fait = mon_depickler_tab.load()
    else :
        tab_deja_fait = []
    
    #if 'fichier_max_'+nom_extension+'.pickle' in os.listdir() : 
    if 'fichier_max_qui_manquent.pickle' in os.listdir() :
        #with open("fichier_max_"+nom_extension+".pickle", 'rb') as fichier_max_pickle :
        with open("fichier_max_qui_manquent.pickle", 'rb') as fichier_max_pickle :
            mon_depickler_sim = pickle.Unpickler(fichier_max_pickle)
            sim_max = mon_depickler_sim.load()
    else :
        sim_max = []
    
    #tab_deja_fait.remove(("fichier_5J7L_DA_30_36","fichier_4V88_A5_137_2"))
    #with open("fichier_similarite_"+nom_extension+".txt", 'w') as fichier_ecriture :
    with open("fichier_similarite_soft_manquent.txt", 'a') as fichier_ecriture :
    #     with open("fichier_similarite_test_2XD0_V_36_21_2.txt", 'a') as fichier_ecriture_2 :
    #         with open("fichier_similarite_4V9F_0_62_12_2.txt", 'a') as fichier_ecriture_3 :
                for fic in os.listdir("graphes_extension/fichiers_couples_qui_manquent/") :
                        tps1 = time.time()
                    #print(fic)
                    #if nom_extension in fic :#or "fichier_2XD0_V_36_21" in fic or "fichier_4V9F_0_62_12" in fic:
                        #c = c+1
                        #print(fic)
                        element1 = fic.split('_')[2] + '_' + fic.split('_')[3] + '_' + fic.split('_')[4] + '_' + fic.split('_')[5] + '_' + fic.split('_')[6]
                        element2 = fic.split('_')[7] + '_' + fic.split('_')[8] + '_' + fic.split('_')[9] + '_' + fic.split('_')[10] + '_' + fic.split('_')[11][:len(fic.split('_')[11])-7]
                        #element1 = "fichier_5J7L_DA_48_30"
                        #element2 = "fichier_5DM6_X_25_15"
                        #fic = "couples_possibles_fichier_1FJG_A_48_11_fichier_1FJG_A_48_8.pickle"
                        if (element1, element2) not in tab_deja_fait :
                            tab_deja_fait.append((element1,element2))
                            c = c+1
                            print(c)
                            with open("graphes_extension/"+element1+".pickle", 'rb') as fichier1 :
                                mon_depickler1 = pickle.Unpickler(fichier1)
                                graphe1 = mon_depickler1.load()     
                                with open("graphes_extension/"+element2+".pickle", 'rb') as fichier2 :
                                    mon_depickler2 = pickle.Unpickler(fichier2)
                                    graphe2 = mon_depickler2.load()
                     
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
                                                                                
                                    
                                    with open("graphes_extension/fichiers_couples_qui_manquent/"+fic, 'rb') as fichier_pickle :
                                                    memory_error = False
                                                    graphe_motif = nx.MultiGraph()
                                                    for i in range(1,6) :
                                                        graphe_motif.add_node((i,i))
                                                    graphe_motif.add_edge((1,1),(2,2), type="NON_CAN")
                                                    graphe_motif.add_edge((1,1),(3,3), type="COV")
                                                    graphe_motif.add_edge((1,1),(5,5), type="NON_CAN")
                                                    graphe_motif.add_edge((2,2),(4,4), type="COV")
                                                    graphe_motif.add_edge((2,2),(5,5), type="CAN")
                                                    graphe_motif.add_edge((3,3),(4,4), type="NON_CAN")
                                                    
                                                                            
                                                    mon_depickler = pickle.Unpickler(fichier_pickle)
                                                    couples_possibles = mon_depickler.load()
                                                    #print(fic)
                                                    print(couples_possibles)
                                                    if couples_possibles != None :
                                                    #print(couples_possibles)
                                                        for i in range(4) :
                                                            if 'memory error' in couples_possibles[i] :
                                                                memory_error = True
                                                                break
                                                            for chaine in couples_possibles[i] :
                                                                chaine.insert(0, (i+1, i+1))
                                                        #print(couples_possibles)
                                                        
                                                        
                                                        
                                                        if memory_error == False :
                                                            #print("ramousnif")
                                                            t = time.time()
                                                            
                                                            graphe_commun_max = nx.MultiGraph()
                                                            #result = exectimeout(50, toutes_comparaisons, args=(couples_possibles, graphe1, graphe2, graphe_motif, graphe_commun_max))
                                                            maxi = 0
                                                            couples_max = []
                                                            
                                                            
#                                                             for compt_1 in range(max(len(couples_possibles[0]), 1)) :
#                                                                 if len(couples_possibles[0]) > 0 :
#                                                                     elt_1 = couples_possibles[0][compt_1]
#                                                                 else :
#                                                                     elt_1 = []
#                                                                 for compt_2 in range(max(len(couples_possibles[1]), 1)) :
#                                                                     if len(couples_possibles[1]) > 0 :
#                                                                         elt_2 = couples_possibles[1][compt_2]
#                                                                     else :
#                                                                         elt_2 = []
#                                                                     for compt_3 in range(max(len(couples_possibles[2]), 1)) :
#                                                                         if len(couples_possibles[2]) > 0 :
#                                                                             elt_3 = couples_possibles[2][compt_3]
#                                                                         else :
#                                                                             elt_3 = [] 
#                                                                         for compt_4 in range(max(len(couples_possibles[3]), 1)) :
#                                                                             if len(couples_possibles[3]) > 0 :
#                                                                                 elt_4 = couples_possibles[3][compt_4]
#                                                                             else :
#                                                                                 elt_4 = []
                                                            elt = []
                                                            for i in range(4) : 
                                                                max_taille_couples = 0
                                                                pos_max_taille_couples = -1
                                                                compteur = 0
                                                                for couple in couples_possibles[i] :
                                                                    if len(couple) > max_taille_couples :
                                                                        max_taille_couples = len(couple)
                                                                        pos_max_taille_couples = compteur
                                                                    compteur += 1 
                                                                print("pos_max_taille_couples")
                                                                print(pos_max_taille_couples)
                                                                if pos_max_taille_couples != -1 :
                                                                    elt.append(couples_possibles[i][pos_max_taille_couples])
                                                                else :
                                                                    elt.append([])  
 
                                                            graphe_commun = graphe_motif.copy()
                                                            
                                                            graphe_commun = comparaison(elt[0], graphe1, graphe2, elt[1], elt[2], elt[3], couples_possibles, 1, graphe_commun)
                                                            graphe_commun = comparaison(elt[1], graphe1, graphe2, elt[0], elt[2], elt[3], couples_possibles, 2, graphe_commun)
                                                            #print("petit rat")
                                                            #print(graphe_commun)
                                                            graphe_commun = comparaison(elt[2], graphe1, graphe2, elt[0], elt[1], elt[3], couples_possibles, 3, graphe_commun)                            
                                                            graphe_commun = comparaison(elt[3], graphe1, graphe2, elt[0], elt[1], elt[2], couples_possibles, 4, graphe_commun)
                                                                
                                                            
                                                            sim = ((graphe_commun.number_of_nodes() + graphe_commun.number_of_edges())*(graphe_commun.number_of_nodes() + graphe_commun.number_of_edges()))/((graphe1.number_of_nodes()+compteur_arc_arete_1)*(graphe2.number_of_nodes()+compteur_arc_arete_2))
                                                            if maxi <= sim and sim <= 1.0 :
                                                                maxi = sim
                                                                del(couples_max[:])
                                                                couples_max.append(elt[0])
                                                                couples_max.append(elt[1])
                                                                couples_max.append(elt[2])
                                                                couples_max.append(elt[3])
                                                                graphe_commun_max = graphe_commun.copy()   
                                                           
                                                            print("maxi")
                                                            print(maxi)
                                                            print("couples max")
                                                            print(couples_max)
                                                            print("graphe")
                                                            print(graphe_commun_max.nodes.data())
                                                            print(graphe_commun_max.edges.data())
                                                            print("couples possibles")
                                                            print(couples_possibles)
                                                            
                                                             
                                                            sim_max.append(maxi)
                                                            
                                                            
                                                            print(graphe_commun_max)
                                                            with open("fichier_max_petit.pickle", 'wb') as fichier_max_pickle :
                                                                mon_pickler = pickle.Pickler(fichier_max_pickle)
                                                                mon_pickler.dump(sim_max)
                                                                
                                                            #if element1 == "fichier_4V9F_0_62_12" or element2 == "fichier_4V9F_0_62_12" :
                                                                fichier_ecriture.write(element1 + '\n' + element2 + '\n')
                                                                fichier_ecriture.write("Graphe : ")
                                                                fichier_ecriture.write(str(graphe_commun_max.nodes.data())+ '\n')
                                                                fichier_ecriture.write(str(graphe_commun_max.edges.data())+ '\n')
                                                                fichier_ecriture.write("Similarite :" + str(maxi))
                                                                fichier_ecriture.write('\n\n')  
                                                            
                                                            
                                
                                
                        print("nombre")
                        print(c)
                        print("temps")
                        print(time.time() - tps1)
                        #with open("tab_calcules_sim_"+nom_extension+".pickle", 'wb') as fichier_deja_fait :
                        with open("tab_calcules_sim_qui_manquent.pickle", 'wb') as fichier_deja_fait :
                            mon_pickler_2 = pickle.Pickler(fichier_deja_fait)
                            mon_pickler_2.dump(tab_deja_fait)

if __name__ == '__main__':
    sous_graphe_commun_max("fichier_5DM6_X_127_7")             
