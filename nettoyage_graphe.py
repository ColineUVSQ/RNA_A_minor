'''
Created on 11 déc. 2018

@author: coline
'''
import os
import pickle

for fichier in os.listdir("graphes_extension/") :
    if "pickle" in fichier :
        with open("graphes_extension/"+fichier, 'rb') as fichier_graphe :
            mon_depickler_1 = pickle.Unpickler(fichier_graphe)
            graphe = mon_depickler_1.load()
            #print(graphe.edges.data())
        
            liste_edges = []
            a_enlever = []
            for (u, v, keys, t) in graphe.edges(data="label", keys = True) :
                if (u,v,t) not in liste_edges :
                    liste_edges.append((u,v,t))
                else : 
                    a_enlever.append((u,v,t, keys))
            
            avant = graphe.number_of_edges()
            for elt in a_enlever :
                graphe.remove_edge(elt[0], elt[1], key = elt[3])
            
            apres = graphe.number_of_edges()
            
        if avant != apres :
            print(fichier)
            print(avant)
            print(apres)
            
            with open("graphes_extension/"+fichier, 'wb') as fichier_graphe_new :
                mon_pickler = pickle.Pickler(fichier_graphe_new)
                mon_pickler.dump(graphe)
            
                
                
               
                    