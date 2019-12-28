'''
Created on 5 oct. 2018

@author: Coline Gi

Recuperation des motifs A-minor de CaRNAval
'''
import json
import pickle

fichier_json = open('dataset.json', 'r', encoding="utf-8") 

with fichier_json as json_data:
    data_dict = json.load(json_data)
    with open('fichiers_pickle/dico_liaisons_entier_liste_adj.pickle', 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        tab_temp = mon_depickler.load()
        #print(tab_temp)
        tab_new = []
        for occ in tab_temp :
            noeuds = data_dict[occ["num_motif"]-1]["l_graphs"][occ["num_occ"]-1]["nodes"]
            liens = data_dict[occ["num_motif"]-1]["l_graphs"][occ["num_occ"]-1]["links"]
            #print(occ)
            stem = False
            for noeud in noeuds :
                if noeud["part"] == "Stem" :
                    stem = True
            if stem == True :
                print(occ)
                sommets = [-1]*5
                for noeud in occ["liaisons"].keys() :
                    if -1 in sommets :
                        sommets = [-1]*5
                        liaison_TSS = False
                        liaison_CSS = False
                        liaison_B53 = False
                        for liaison in occ["liaisons"][noeud] :
                            if liaison["type"] == 'CSS' :
                                liaison_CSS = True
                            if liaison["type"] == 'TSS' :
                                liaison_TSS = True
                            if liaison["type"] == 'B53' :
                                liaison_B53 = True
                        if liaison_CSS and liaison_TSS and liaison_B53 :
                            sommets[0] = noeud
                            for liaison in occ["liaisons"][noeud] :
                                if liaison["type"] == 'CSS' :
                                    for liaison_voisin in occ["liaisons"][liaison["voisin"]] :
                                        if liaison_voisin["type"] == 'CWW' :
                                            if sommets[4] != -1 :
                                                if sommets[4] == liaison_voisin["voisin"] :
                                                        sommets[1] = liaison["voisin"]
                                            else :
                                                sommets[1] = liaison["voisin"]
                                if liaison["type"] == 'TSS' :
                                    for liaison_voisin in occ["liaisons"][liaison["voisin"]] :
                                        if liaison_voisin["type"] == 'CWW' :
                                            if sommets[1] != -1 :
                                                if sommets[1] == liaison_voisin["voisin"] :
                                                        sommets[4] = liaison["voisin"]
                                            else :
                                                sommets[4] = liaison["voisin"] 
                                if liaison["type"] == 'B53' :
                                    for liaison_voisin in occ["liaisons"][liaison["voisin"]] :
                                        if liaison_voisin["type"] == 'CSS' :
                                            if sommets[3] != -1 :
                                                if sommets[3] == liaison_voisin["voisin"] :
                                                        sommets[2] = liaison["voisin"]
                                            else :
                                                sommets[2] = liaison["voisin"]
                                                sommets[3] = liaison_voisin["voisin"]
                print(sommets)
                if -1 not in sommets :
                    tab_new.append({"num_PDB": occ["num_PDB"], "num_ch": occ["num_ch"], "num_motif": occ["num_motif"], "num_occ": occ["num_occ"], "liaisons": occ["liaisons"],"a_minor": sommets })
        print(len(tab_new))
        with open("fichiers_pickle/a-minor_test2.pickle", 'wb') as fichier_pickle :
            mon_pickler = pickle.Pickler(fichier_pickle)
            mon_pickler.dump(tab_new)          
        with open("a-minor_test2.txt", 'w') as fichier_a_minor :
            fichier_a_minor.write(str(tab_new))

        #        if len(noeud_cible) != 0 and nb_liaisons_CSS >= 4 and nb_liaisons_TSS >= 2 :
        #            tab_new.append(occ)
        #with open("a-minor.txt", 'w') as fichier_a_minor :
        #    fichier_a_minor.write(str(tab_new))
        #    print(tab_new)  
        #    print(len(tab_new))
        
        #with open("fichiers_pickle/a-minor.pickle", 'wb') as fichier_pickle :
        #    mon_pickler = pickle.Pickler(fichier_pickle)
        #    mon_pickler.dump(tab_new)                
                        
                    #for liaison in occ["liaisons"] :
                    #    if 
               # garde_liaison = []
#                 for liaison_1 in occ["liaisons"] :
#                     for liaison_2 in occ["liaisons"] :    
#                         if liaison_2[0] != liaison_2[1] - 1 :
#                             if len(garde_liaison) == 0 :
#                                 garde_liaison.append(liaison_2)
#                             elif len(garde_liaison) == 1:
#                                 if (liaison_2[0] == garde_liaison[0]+1 and liaison_2[1] == garde_liaison[1]-1)  or (liaison_2[0] == garde_liaison[0]-1 and liaison_2[1] == garde_liaison[1]+1) :
#                                     garde_liaison.append(liaison_2)
                        
                #tab_new.append(occ)
            
#print(len(tab_new))
#print(tab_new)
