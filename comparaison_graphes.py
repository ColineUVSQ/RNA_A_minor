'''
Created on 12 déc. 2018

@author: coline

comparaison de graphes : noeuds et arêtes
'''
import pickle

        
with open("dico_graphe_4V9F_0_25_4.pickle", 'rb') as fichier_graphe :
    mon_depickler = pickle.Unpickler(fichier_graphe)
    dico_graphe_4v9f = mon_depickler.load()
    
    with open("dico_graphe.pickle", 'rb') as fichier_graphe_tout :
        mon_depickler_2 = pickle.Unpickler(fichier_graphe_tout)
        dico_graphe = mon_depickler_2.load()
        
        for cle in dico_graphe_4v9f.keys() :
            
            diff = 0
            for noeud in dico_graphe_4v9f[cle].nodes() :
                if noeud not in dico_graphe[cle].nodes() :
                    diff += 1
            
            for edge in dico_graphe_4v9f[cle].edges() :
                if edge not in dico_graphe[cle].edges() :
                    diff += 1 
            
            print(cle)
            print(diff)
        
        
        