'''
Created on 24 janv. 2019

@author: coline
'''

import pickle
import networkx as nx
import numpy as np
import csv

from recup_data.calcul_sim import calcul_sim_non_cov_sans_motif


def recherche_chaines(graphe_complet, sommet) :
    with open("graphes_extension/fichier_"+graphe_complet.nodes[sommet]["nom"]+".pickle", 'rb') as fichier_graphe1 :
        mon_depickler_graphe1 = pickle.Unpickler(fichier_graphe1)
        graphe = mon_depickler_graphe1.load()
        chaines = [[],[],[],[]]
            
        for noeud,ch in graphe.nodes(data="chaine") :
            if 1 in ch :
                    chaines[0].append(noeud)
            if 2 in ch :
                    chaines[1].append(noeud)
            if 3 in ch :
                    chaines[2].append(noeud)
            if 4 in ch :
                    chaines[3].append(noeud)
        
        chaines_temp = [[],[],[],[]]
        for u,v in graphe.edges() :
            for i in range(4) :
                if u in chaines[i] and v not in chaines[i] and v not in [1,2,3,4,5]:
                    chaines_temp[i].append(v)
                elif v in chaines[i] and u not in chaines[i] and u not in [1,2,3,4,5] :
                    chaines_temp[i].append(u)
        
        for i in range(4) :
            chaines[i] = chaines[i] + chaines_temp[i]         
            chaines[i] = list(set(chaines[i]))
            
        return chaines

def calcul_sim_chaines_2_4():
    with open("graphe_complet_pondere_sim.pickle", 'rb') as fichier_graphe_complet :
        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
        graphe_complet = mon_depickler_complet.load()
        print(graphe_complet.edges.data())
        
        
        nx.set_node_attributes(graphe_complet, [], "chaines")
        for u in graphe_complet.nodes() :
            graphe_complet.nodes[u]["chaines"] = recherche_chaines(graphe_complet, u)
        print(graphe_complet.nodes.data())
        
        nx.set_edge_attributes(graphe_complet, -1.0, "sim_chaines_2_4")
        for u,v,data in graphe_complet.edges(data=True) :
            chaine_1 = list(set(graphe_complet.nodes[u]["chaines"][1] + graphe_complet.nodes[u]["chaines"][3] + [1,2,3,4,5]))
            chaine_2 = list(set(graphe_complet.nodes[v]["chaines"][1] + graphe_complet.nodes[v]["chaines"][3] + [1,2,3,4,5]))
            
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
                        chaine_commun = list(set(chaine_commun + [(1,1),(2,2),(3,3),(4,4),(5,5)]))
                    else :
                        graphe_commun = dico_graphe[("fichier_"+graphe_complet.nodes[v]["nom"], "fichier_"+graphe_complet.nodes[u]["nom"])]
                    
                        chaine_commun = []
                        for noeud in graphe_commun.nodes() :
                            if noeud[1] in chaine_1 and noeud[0] in chaine_2 :
                                    chaine_commun.append(noeud)
                        chaine_commun = list(set(chaine_commun + [(1,1),(2,2),(3,3),(4,4),(5,5)]))
                    
                    sim_2_4 = calcul_sim_non_cov_sans_motif(graphe1.subgraph(chaine_1), graphe2.subgraph(chaine_2), graphe_commun.subgraph(chaine_commun))
                    print(sim_2_4)
                    
                    graphe_complet.edges[u,v]["sim_chaines_2_4"] = sim_2_4
                    
    
    with open("graphe_complet_pondere_sim.pickle", 'wb') as fichier_graphe_complet_2 :
        mon_pickler_complet = pickle.Pickler(fichier_graphe_complet_2)
        mon_pickler_complet.dump(graphe_complet)
            
                        
            
        
        #for a, b, data in sorted(graphe_complet.edges(data=True), key=lambda x: x[2]['poids']):
            
            
if __name__ == '__main__':
    
    #calcul_sim_chaines_2_4()
    
    with open("graphe_complet_pondere_sim.pickle", 'rb') as fichier_graphe_complet :
        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
        graphe_complet = mon_depickler_complet.load()
          
        print(graphe_complet.edges[0,2])
        print(graphe_complet.nodes[0])
        print(graphe_complet.nodes[2])
        
        with open("stats_chaines_2_4.csv", 'w', newline='') as fichier :
            csvwriter = csv.writer(fichier)
            i = 0.1
            while i < 1.0 :
                tab = []
                for a, b, data in graphe_complet.edges(data=True) :
                    if data['sim_chaines_2_4'][0] < i and data['sim_chaines_2_4'][0] >= i-0.1:
                        tab.append(data['sim_chaines_2_4'][0]-data['poids'])
                csvwriter.writerow([i, np.mean(tab), np.std(tab)])
                print(tab)
                i = i+0.1
        
        
        #print(np.mean(np.array(tab)))
        
        
#         print(graphe_complet.edges[28,64])
#         
#         chaine_1 = list(set(graphe_complet.nodes[28]["chaines"][1] + graphe_complet.nodes[28]["chaines"][3] + [1,2,3,4,5]))
#         chaine_2 = list(set(graphe_complet.nodes[64]["chaines"][1] + graphe_complet.nodes[64]["chaines"][3] + [1,2,3,4,5]))
#         
#         with open("graphes_extension/fichier_"+graphe_complet.nodes[28]["nom"]+".pickle", 'rb') as fichier_graphe1 :
#             mon_depickler_graphe1 = pickle.Unpickler(fichier_graphe1)
#             graphe1 = mon_depickler_graphe1.load()
#             with open("graphes_extension/fichier_"+graphe_complet.nodes[64]["nom"]+".pickle", 'rb') as fichier_graphe2:
#                 mon_depickler_graphe2 = pickle.Unpickler(fichier_graphe2)
#                 graphe2 = mon_depickler_graphe2.load()
#                 
#                 with open("dico_graphe_epure_en_tout.pickle", 'rb') as fichier_commun :
#                     mon_depickler_commun = pickle.Unpickler(fichier_commun)
#                     dico_graphe = mon_depickler_commun.load()
#                 
#                 if ("fichier_"+graphe_complet.nodes[28]["nom"],"fichier_"+graphe_complet.nodes[64]["nom"]) in dico_graphe.keys() :
#                     graphe_commun = dico_graphe[("fichier_"+graphe_complet.nodes[28]["nom"],"fichier_"+graphe_complet.nodes[64]["nom"])]
#                     
#                     chaine_commun = []
#                     for noeud in graphe_commun.nodes() :
#                         if noeud[0] in chaine_1 and noeud[1] in chaine_2 :
#                                 chaine_commun.append(noeud)
#                     chaine_commun = list(set(chaine_commun + [(1,1),(2,2),(3,3),(4,4),(5,5)]))
#                 else :
#                     graphe_commun = dico_graphe[("fichier_"+graphe_complet.nodes[64]["nom"], "fichier_"+graphe_complet.nodes[28]["nom"])]
#                     
#                     chaine_commun = []
#                     for noeud in graphe_commun.nodes() :
#                         if noeud[1] in chaine_1 and noeud[0] in chaine_2 :
#                                 chaine_commun.append(noeud)
#                     chaine_commun = list(set(chaine_commun + [(1,1),(2,2),(3,3),(4,4),(5,5)]))
#                 print(graphe_commun.nodes.data())
#                 print(chaine_1)
#                 print(chaine_2)
#                 
#                 
#                 sim_2_4 = calcul_sim_non_cov_sans_motif(graphe1.subgraph(chaine_1), graphe2.subgraph(chaine_2), graphe_commun.subgraph(chaine_commun))
#                 print(sim_2_4)
#         with open("graphes_extension/fichier_"+graphe_complet.nodes[64]["nom"]+".pickle", 'rb') as fichier_graphe1 :
#             mon_depickler_graphe1 = pickle.Unpickler(fichier_graphe1)
#             graphe1 = mon_depickler_graphe1.load()
#             
#             print(graphe1.nodes.data())
#             print(graphe1[10])
            
#             with open("graphs_2.92.pickle", 'rb') as fichier_graphes :
#                 mon_depickler_graphes = pickle.Unpickler(fichier_graphes)
#                 graphes = mon_depickler_graphes.load()
#                 
#                 print(graphes[('5J5B','BA')][1529])
        #print(graphe_complet.edges.data())
            