'''
Created on 9 oct. 2018

@author: Coline Gi
'''
import json
import pickle

fichier_json = open('dataset.json', 'r', encoding="utf-8") 
fichier_interactions_longue_distance = open('dico_liaisons_entier.txt', 'r', encoding="utf8")

with fichier_json as json_data:
    data_dict = json.load(json_data)
    with open('fichiers_pickle/dico_liaisons_entier_plus_types.pickle', 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        tab_temp = mon_depickler.load()
        tab_new = []
        dic_new = {}
        for occ in tab_temp :
            noeuds = data_dict[occ["num_motif"]-1]["l_graphs"][occ["num_occ"]-1]["nodes"]
            liens = data_dict[occ["num_motif"]-1]["l_graphs"][occ["num_occ"]-1]["links"]
            dico_adj = {}
            for noeud in noeuds :
                dico_adj.update({noeud["id"] : []})## on initialise chaque tableau d'adjacence
            #print(dico_adj)
            for lien in liens :
                #print(dico_adj[noeuds[lien["source"]]["id"]])
                source = noeuds[lien["source"]]["id"]
                target = noeuds[lien["target"]]["id"]
                
                deja_vu = False
                for elt in dico_adj[source] :
                    if elt["voisin"] == target :
                        deja_vu = True
                if deja_vu == False :
                    dico_adj[source].append({"voisin" : target, "type" : lien["label"]})
                
                deja_vu = False
                for elt in dico_adj[target] :
                    if elt["voisin"] == source :
                        deja_vu = True
                if deja_vu == False :
                    dico_adj[target].append({"voisin" : source, "type" : lien["label"]})
            dict_temp = {"num_PDB" : occ["num_PDB"], "num_ch" : occ["num_ch"], "num_motif" : occ["num_motif"], "num_occ" : occ["num_occ"], "liaisons" : dico_adj}
            tab_new.append(dict_temp)
            
            dic_new.update({(occ["num_PDB"], occ["num_ch"], occ["num_motif"], occ["num_occ"]) : dico_adj})
        
        with open("fichiers_pickle/tab_liaisons_entier_liste_adj.pickle", 'wb') as fichier_sortie :
            mon_pickler = pickle.Pickler(fichier_sortie)
            mon_pickler.dump(tab_new)
        
        with open("fichiers_pickle/dico_liaisons_entier_liste_adj.pickle", 'wb') as fichier2 :
            mon_pickler2 = pickle.Pickler(fichier2)
            mon_pickler2.dump(dic_new)
    
        print(tab_new)    