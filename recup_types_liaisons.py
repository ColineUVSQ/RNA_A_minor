'''
Created on 8 oct. 2018

@author: Coline Gi
'''
import json
import pickle

fichier_json = open('dataset.json', 'r', encoding="utf-8") 
fichier_interactions_longue_distance = open('dico_liaisons_entier.txt', 'r', encoding="utf8")

with fichier_json as json_data:
    data_dict = json.load(json_data)
    with fichier_interactions_longue_distance as data_tab :
        tab_temp = json.load(data_tab)
        tab_temp_plus_types_liaisons = []
        with open('fichiers_pickle/dico_liaisons_entier_plus_types.pickle', 'wb') as fichier_sortie :
            mon_pickler = pickle.Pickler(fichier_sortie)
            for occ in tab_temp :
                liens = data_dict[occ["num_motif"]-1]["l_graphs"][occ["num_occ"]-1]["links"]
                noeuds = data_dict[occ["num_motif"]-1]["l_graphs"][occ["num_occ"]-1]["nodes"]
                tab_liens = []
                for lien in liens :
                    noeud_source = lien["source"]
                    noeud_target = lien["target"]
                    
                    tab_liens.append({"noeud_source" : noeuds[noeud_source]["id"], "noeud_target" : noeuds[noeud_target]["id"], "type" : lien["label"]})
                
                dict_temp = {"num_PDB" : occ["num_PDB"], "num_ch" : occ["num_ch"], "num_motif" : occ["num_motif"], "num_occ" : occ["num_occ"], "liaisons" : tab_liens}
                tab_temp_plus_types_liaisons.append(dict_temp)
            #print(tab_tot)
            mon_pickler.dump(tab_temp_plus_types_liaisons)
            print(tab_temp_plus_types_liaisons)
            
                
            
            