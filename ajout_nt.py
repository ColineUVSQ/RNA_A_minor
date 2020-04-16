'''
Created on 12 oct. 2018

@author: Coline Gi
'''
import pickle
import os
import networkx as nx

with open("graphs_2.92.pickle", 'rb') as fichier_tout :
        mon_depickler_graphes = pickle.Unpickler(fichier_tout)
        graphes = mon_depickler_graphes.load()
        
        with open("fichier_graphes_avec_nt.txt", 'w') as fichier_sortie :
            for element in os.listdir('graphes_extension/'):
                if "pickle" in element :
                    with open("graphes_extension/"+element, 'rb') as fichier_entree :
                        mon_depickler = pickle.Unpickler(fichier_entree)
                        G = mon_depickler.load()
                        
                        nom_cle = (element.split("_")[1], element.split("_")[2]) 
                        
                        nx.set_node_attributes(G, None, "nt")
                        print(nom_cle)
                        for noeud in G.nodes() :
                            print(G.nodes[noeud]["position"])
                            G.nodes[noeud]["nt"] = graphes[nom_cle].nodes[G.nodes[noeud]["position"]]["nt"]
                        
                        fichier_sortie.write(str(nom_cle) + " " + element.split("_")[3] + " " + element.split("_")[4] + "\n")
                        fichier_sortie.write(str(G.nodes.data())+"\n")
                        fichier_sortie.write(str(G.edges.data())+"\n")
                    with open("graphes_extension/"+element, 'wb') as fichier_pickle_sortie :
                        mon_pickler = pickle.Pickler(fichier_pickle_sortie)
                        mon_pickler.dump(G)
            
                        
                        