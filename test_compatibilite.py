'''
Created on 4 déc. 2018

@author: coline
'''
import pickle
import os
# def incomp(group1, group2, incompatibles):
#     if len(group1) == 1 and len(group2) == 1 :
#         if group1 in incompatibles.keys() and group2 in incompatibles[group1] or group2 in incompatibles.keys() and group1 in incompatibles[group2] :
#             return True
#         else :
#             return False
#     elif len(group2) == 1:
#         return incomp(group1[1:len(group1)], group2, incompatibles)
#     elif len(group1) == 1:
#         return incomp(group1, group2[1:len(group2)], incompatibles)
#     else :
        
#max = 0
for fic in os.listdir("graphes_extension/fichiers_couples_qui_manquent/") :
    with open("graphes_extension/fichiers_couples_qui_manquent/"+fic, 'rb') as fichier_pickle :
            if fic not in os.listdir("graphes_extension/fichiers_couples_qui_manquent_incomp/") :
                with open("graphes_extension/fichiers_couples_qui_manquent_incomp/"+fic, 'wb') as fichier_sortie :
                    mon_depickler = pickle.Unpickler(fichier_pickle)
                    couples_possibles = mon_depickler.load()
                    print(couples_possibles)
                    element1 = fic.split('_')[2] + '_' + fic.split('_')[3] + '_' + fic.split('_')[4] + '_' + fic.split('_')[5] + '_' + fic.split('_')[6]
                    element2 = fic.split('_')[7] + '_' + fic.split('_')[8] + '_' + fic.split('_')[9] + '_' + fic.split('_')[10] + '_' + fic.split('_')[11][:len(fic.split('_')[11])-7]
                      
                    with open("graphes_extension/"+element1+".pickle", 'rb') as fichier1 :
                            mon_depickler1 = pickle.Unpickler(fichier1)
                            graphe1 = mon_depickler1.load()     
                            with open("graphes_extension/"+element2+".pickle", 'rb') as fichier2 :
                                mon_depickler2 = pickle.Unpickler(fichier2)
                                graphe2 = mon_depickler2.load()
                                
                                couples_possibles_new = []
                                for i in range(4) :
                                    new_couples = []
                                    if 'memory error' in couples_possibles[i] :
                                        memory_error = True
                                        break
                                    for chaine in couples_possibles[i] :
                                        new_chaine = []
                                        for couple in chaine :
                                            new_chaine.append(couple)
                                        new_couples.append(new_chaine)
                                    couples_possibles_new.append(new_couples)
                                
                                print(couples_possibles_new)
                                new_couples_possibles = []
                                 
                                for chaine in couples_possibles_new :
                                    new_chaine = []
                                    for possib in chaine :
                                        incompatibles = {}
                                        for i in range(len(possib)) :
                                            for j in range(i+1, len(possib)) :
                                                if (graphe1.nodes[possib[i][0]]["position"] < graphe1.nodes[possib[j][0]]["position"] and graphe2.nodes[possib[i][1]]["position"] > graphe2.nodes[possib[j][1]]["position"]) or (graphe1.nodes[possib[i][0]]["position"] > graphe1.nodes[possib[j][0]]["position"] and graphe2.nodes[possib[i][1]]["position"] < graphe2.nodes[possib[j][1]]["position"]) :
                                                    if possib[i] not in incompatibles.keys() :
                                                        incompatibles.update({possib[i]  : [possib[j]]})
                                                    else :
                                                        incompatibles[possib[i]].append(possib[j])
            #                             print(incompatibles)
            #                             for elt in incompatibles.keys() : 
            #                                 print(len(incompatibles[elt]))
            #                                 if max < len(incompatibles[elt]) :
            #                                     max = len(incompatibles[elt])
            #print(max)
                                        num_chaine = len(new_chaine)
                                        new_chaine.append([])
                                        print("incomp")
                                        print(incompatibles)
                                        print("debut")
                                        print(new_chaine)
                                        traites = []
                                        for i in range(len(possib)) :
                                            print("i")
                                            print(possib[i])
                                            if possib[i] not in traites :
                                                if possib[i] in incompatibles.keys():
                                                    if len(incompatibles[possib[i]]) == 1 : 
                                                        print("ramousnif")
                                                        new_chaine[len(new_chaine)-1].append(possib[i])
                                                        new_chaine.append([])
                                                        new_chaine[len(new_chaine)-1].append(incompatibles[possib[i]][0])
                                                        ### une chaine avec le incomp, et une chaine avec son incompatible
                                                    else :
                                                        new_chaine.append([])
                                                        new_chaine[len(new_chaine)-1].append(possib[i])
                                                        for k in range(len(incompatibles[possib[i]])) :
                                                            print("k")
                                                            
                                                            incomp1 = incompatibles[possib[i]][k]
                                                            print(incomp1)
                                                            peut_ajouter = False
                                                            new_possib = num_chaine
                                                            while new_possib < len(new_chaine) and peut_ajouter == False :
                                                                peut_ajouter = True
                                                                for elt in new_chaine[new_possib] :
                                                                    if (incomp1 in incompatibles.keys() and elt in incompatibles[incomp1]) or (elt in incompatibles.keys() and incomp1 in incompatibles[elt]) :
                                                                        peut_ajouter = False
                                                                new_possib += 1
                                                            if new_possib == len(new_chaine) :
                                                                new_chaine.append([])
                                                                new_chaine[len(new_chaine)-1].append(incomp1)
                                                            else :
                                                                new_chaine[new_possib-1].append(incomp1)
                                                            print(new_chaine)
                                                            traites.append(incomp1)
                                        for i in range(len(possib)) :
                                            incomp = False
                                            for cle in incompatibles.keys() :
                                                if possib[i] in incompatibles[cle] :
                                                    incomp = True
                                            if possib[i] not in incompatibles.keys() and incomp == False :
                                                for new_possib in range(num_chaine, len(new_chaine)) :
                                                    new_chaine[new_possib].append(possib[i])
                                    new_couples_possibles.append(new_chaine)
                                 
                                print(new_couples_possibles)  
                                
                                mon_pickler = pickle.Pickler(fichier_sortie)
                                mon_pickler.dump(new_couples_possibles)
                                         
                        
                                        
                    
                                
                    