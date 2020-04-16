'''
Created on 28 janv. 2019

@author: coline
'''
import pickle
import networkx as nx

def calcul_sim_non_cov_sans_motif(graphe1, graphe2, graphe_commun, chaines_1, chaines_2, chaines_commun, i, le_bon):
    
    compteur_arete_1 = 0
    for (u, v, t) in graphe1.edges(data="label") :
        if (u in chaines_1 or v in chaines_1) :
            if t != 'B53' :
                compteur_arete_1 += 1
                
    if i == 2 or i == 3 :
        compteur_arete_1 = compteur_arete_1 - 2
    else :
        compteur_arete_1 = compteur_arete_1 - 4
    compteur_arete_1 = compteur_arete_1/2
#     if element1 == "fichier_5DM6_X_197_1" and element2 == "fichier_5DM6_X_48_9" :

    compteur_arete_2 = 0
    for (u, v, t) in graphe2.edges(data="label") :
        if (u in chaines_2 or v in chaines_2) :
            if t != 'B53' :
                compteur_arete_2 += 1
    if i == 2 or i == 3 :
        compteur_arete_2 = compteur_arete_2 - 2
    else :
        compteur_arete_2 = compteur_arete_2 - 4
    compteur_arete_2 = compteur_arete_2/2
#     if element1 == "fichier_5DM6_X_197_1" and element2 == "fichier_5DM6_X_48_9" :
    
    compteur_arete_commun = 0
    for (u, v, t) in graphe_commun.edges(data="type") :
        if (u in chaines_commun or v in chaines_commun) :
            if t != 'COV' :
                compteur_arete_commun += 1
                
    if i == 2 or i == 3 :
        compteur_arete_commun = compteur_arete_commun - 1
    else :
        compteur_arete_commun = compteur_arete_commun - 2
        
    #compteur_arete_commun = compteur_arete_commun/2
    
#     if element1 == "fichier_5DM6_X_197_1" and element2 == "fichier_5DM6_X_48_9" :
#         print(compteur_arete_commun)
#         print(compteur_arete_1)
#         print(compteur_arete_2)
    
#     print(compteur_arete_1)
#     print(compteur_arete_2)
#     print(compteur_arete_commun)

    if le_bon :
        print(compteur_arete_1)
        print(compteur_arete_2)
        print(compteur_arete_commun)
    
    if max(compteur_arete_1, compteur_arete_2) > 0 :
        sim = compteur_arete_commun/max(compteur_arete_1, compteur_arete_2)
    else :
        sim = 0.0
    return sim

def sim_par_chaine() :
    with open("graphe_complet_pondere_sim.pickle", 'rb') as fichier_graphe_complet :
        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
        graphe_complet = mon_depickler_complet.load()
        
        #nx.set_edge_attributes(graphe_complet, [], "sim_par_chaine")
        
        for u,v,data in graphe_complet.edges(data=True) :
            graphe_complet.edges[u,v]["sim_par_chaine"] = []
            for i in range(4) :
                chaine_1 = graphe_complet.nodes[u]["chaines"][i].copy()
                chaine_2 = graphe_complet.nodes[v]["chaines"][i].copy()
                
                chaine_1.append(i+1)
                chaine_2.append(i+1)
                
                with open("graphes_extension/fichier_"+graphe_complet.nodes[u]["nom"]+".pickle", 'rb') as fichier_graphe1 :
                    mon_depickler_graphe1 = pickle.Unpickler(fichier_graphe1)
                    graphe1 = mon_depickler_graphe1.load()
                    with open("graphes_extension/fichier_"+graphe_complet.nodes[v]["nom"]+".pickle", 'rb') as fichier_graphe2:
                        mon_depickler_graphe2 = pickle.Unpickler(fichier_graphe2)
                        graphe2 = mon_depickler_graphe2.load()
                        
                        with open("dico_graphe_epure_en_tout.pickle", 'rb') as fichier_commun :
                            mon_depickler_commun = pickle.Unpickler(fichier_commun)
                            dico_graphe = mon_depickler_commun.load()
                        
                        if ("fichier_"+graphe_complet.nodes[u]["nom"],"fichier_"+graphe_complet.nodes[v]["nom"]) in dico_graphe.keys() :
                            graphe_commun = dico_graphe[("fichier_"+graphe_complet.nodes[u]["nom"],"fichier_"+graphe_complet.nodes[v]["nom"])]
                            
                            chaine_commun = []
                            for noeud in graphe_commun.nodes() :
                                if noeud[0] in chaine_1 and noeud[1] in chaine_2 :
                                        chaine_commun.append(noeud)
                            chaine_commun = list(set(chaine_commun + [(i+1)]))
                        else :
                            graphe_commun = dico_graphe[("fichier_"+graphe_complet.nodes[v]["nom"], "fichier_"+graphe_complet.nodes[u]["nom"])]
                            
                            chaine_commun = []
                            for noeud in graphe_commun.nodes() :
                                if noeud[1] in chaine_1 and noeud[0] in chaine_2 :
                                        chaine_commun.append(noeud)
                            chaine_commun = list(set(chaine_commun + [(i+1)]))
                        
                        
                        sim = calcul_sim_non_cov_sans_motif(graphe1, graphe2, graphe_commun, chaine_1, chaine_2, chaine_commun, i, 0)
                        
                        graphe_complet.edges[u,v]["sim_par_chaine"].append(sim) 
                        
    with open("graphe_complet_pondere_sim.pickle", 'wb') as fichier_graphe_complet_2 :
        mon_pickler_complet = pickle.Pickler(fichier_graphe_complet_2)
        mon_pickler_complet.dump(graphe_complet)
        
if __name__ == '__main__':
    sim_par_chaine()
        

