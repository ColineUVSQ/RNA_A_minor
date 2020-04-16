'''
Created on 24 juil. 2019

@author: coline

Choix d'un des homologues

'''
from recup_data.constantes import EXTENSION_PATH, HOMOLOGUES
import pickle
import numpy as np
from random import randint

''' On choisit l'homologue qui perd le moins d'aretes incidentes en augmentant les seuils '''
def calcul_nombre_aretes_chaque_seuil():
    with open(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes"+"graphe_complet_pondere_sim_toutes_aretes_coeff_all1_max_taille_taille_max.pickle", 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        graphe_complet = mon_depickler.load()
        
        print(graphe_complet.number_of_nodes())
        print(graphe_complet.number_of_edges())
       
        maxi_somme = 0
        maxi_groupe = []
        for _ in range(2000) :
            graphe_test = graphe_complet.copy()
            a_enlever = []
            for groupe in HOMOLOGUES :
                choix = randint(0, len(groupe)-1) 
                compteur = 0
                for elt in groupe :
                    if compteur != choix :
                        a_enlever.append(elt)
                    compteur += 1
            
            a_enlever_noeuds = []
            for noeud, data in graphe_test.nodes(data=True) :
                if data["nom"] in a_enlever :
                    a_enlever_noeuds.append(noeud)
                    
            for elt in a_enlever_noeuds :
                graphe_test.remove_node(elt)
             
            print(graphe_test.number_of_nodes())
            print(graphe_test.number_of_edges())  
             
            print(graphe_test.nodes.data()) 
            somme_aretes = 0
            for i in np.arange(0, 1.1, 0.1) :
                a_enlever_aretes = []
                #print(a_enlever_aretes)
                #print(graphe_test.number_of_edges())
                for u,v,data in graphe_test.edges(data=True) :
                    #print(u)
                    #print(v)
                    if data["poids"] < i :
                        a_enlever_aretes.append((u,v))
                
                #print(a_enlever_aretes)         
                for elt in a_enlever_aretes :
                    graphe_test.remove_edge(elt[0], elt[1])
                
                   
                print("%1.1f : %d"%(i,graphe_test.number_of_edges()))
                somme_aretes += graphe_test.number_of_edges()
            print(somme_aretes)
            if somme_aretes > maxi_somme :
                maxi_somme = somme_aretes
                maxi_groupe = list(a_enlever)
                
        print(maxi_somme)
        print(maxi_groupe)
            
            
if __name__ == '__main__':
    calcul_nombre_aretes_chaque_seuil()