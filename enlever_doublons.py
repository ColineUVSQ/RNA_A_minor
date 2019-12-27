'''
Created on 5 oct. 2018

@author: Coline Gi

Enlever les doublons de CaRNAval pour avoir chaque occurrence de pseudonoeuds une seule fois
'''
import json
import time

fichier_json = open('dataset.json', 'r', encoding="utf-8") 

tab_temp = []
tps1 = time.clock()
##Construction du tableau sans les doublons (version plus rapide, mais moins facile � comprendre)
with fichier_json as json_data:
    data_dict = json.load(json_data)
    fichier_doublons = open("data_elimines.txt", "w")
    compteur_motif = 0
    for motif in data_dict :
        compteur_motif = compteur_motif+1
        compteur_occ = 0
        for occurrence in motif["names"] :
            compteur_occ = compteur_occ+1;
            #print(occurrence[1])
            #print(occurrence[0][0])
            chaine_num_PDB = occurrence[0][0][0]
            chaine_num_ch = occurrence[0][0][1]
            #print(chaine)
            present1 = False
            present2 = False
            for dico in tab_temp :
                if dico["num_PDB"] == chaine_num_PDB and dico["num_ch"] == chaine_num_ch:
                    present1 = True
                    for liaison in dico["liaisons"] :
                        if liaison not in occurrence[1] :
                            present1 = False
                    if present1 == True :
                        compte = 0
                        for dic in tab_temp :  
                            if dic == dico :
                                fichier_doublons.write(str(tab_temp[compte]))
                                del(tab_temp[compte])
                            compte = compte+1
                        dict_temp = {"num_PDB" : chaine_num_PDB, "num_ch" : chaine_num_ch, "num_motif" : compteur_motif, "num_occ" : compteur_occ, "liaisons" : occurrence[1]}
                        
                    present2 = True
                    for liaison in occurrence[1] :
                        if liaison not in dico["liaisons"] :
                            present2 = False
            if present1 == True and present2 == False:
                tab_temp.append(dict_temp)
            if present1 == False and present2 == False :
                dict_temp = {"num_PDB" : chaine_num_PDB, "num_ch" : chaine_num_ch, "num_motif" : compteur_motif, "num_occ" : compteur_occ, "liaisons" : occurrence[1]}
                tab_temp.append(dict_temp)
            
            #print(dico["liaisons"][0])
    compteur1 = 0
    tab_elimine = []
    for dico1 in tab_temp :
        for dico2 in tab_temp :
            if dico2 != dico1 and dico1["num_PDB"] == dico2["num_PDB"] and dico1["num_ch"] == dico2["num_ch"] :
                present = True
                for liaison in dico1["liaisons"] :
                    if liaison not in dico2["liaisons"] :
                        present = False
                if present == True :
                    tab_elimine.append(compteur1)
                    #print("dico1")
                    #print(dico1)
                    #print("dico2")
                    #print(dico2)
        compteur1 = compteur1+1
    print(tab_elimine)
    print(len(tab_elimine))
    for i in tab_elimine :
        del(tab_temp[i])
   
    
    with open('dico_liaisons_entier.txt', 'w') as file:
        json.dump(tab_temp, file)
    #print(tab_temp)
    print(len(tab_temp))
    tps2 = time.clock()
    print(tps2-tps1)
    fichier_doublons.close()