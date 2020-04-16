'''
Created on 5 déc. 2018

@author: coline
'''

import pickle
import time

# with open("fichiers_tries.pickle", "rb") as fichier_tri :
#         mon_depickler_tri = pickle.Unpickler(fichier_tri)
#         tri = mon_depickler_tri.load()
    #for compte_tri_1 in range(124) :
        #element1 = tri[compte_tri_1]
element1 = "fichier_5DM6_X_328_2.pickle"
print(element1)
with open("graphes_extension/"+element1, 'rb') as fichier1 :
                    mon_depickler1 = pickle.Unpickler(fichier1)
                    graphe1 = mon_depickler1.load()
                    chaines_1 = [[1]]
                    for i in range(1,5) :
                        compteur = i
                        if i != 1 : chaines_1.append([i])
                        liaison_B53 = True
                        while liaison_B53 :
                            liaison_B53 = False
                            temp = compteur
                            for voisin in graphe1.successors(compteur) :
                                for arc in graphe1[compteur][voisin] :
                                    if voisin not in [1,2,3,4] and voisin not in chaines_1[len(chaines_1)-1] and graphe1[compteur][voisin][arc]["label"] == 'B53' :
                                        liaison_B53 = True
                                        temp = voisin
                                        chaines_1[len(chaines_1)-1].append(voisin)
                                        
                            for voisin in graphe1.predecessors(compteur) :
                                for arc in graphe1[voisin][compteur] :
                                    if voisin not in [1,2,3,4] and voisin not in chaines_1[len(chaines_1)-1] and graphe1[voisin][compteur][arc]["label"] == 'B53' :
                                        liaison_B53 = True
                                        temp = voisin
                                        chaines_1[len(chaines_1)-1].append(voisin)
                            compteur = temp
                    somme_temps = 0.0
                #compte_tri = 0 
                #for compte_tri in range(124) :
                    #element2 = tri[compte_tri]
                    tps1 = time.time()
                    element2 = "fichier_5FDU_1A_62_14.pickle"
                    #element2 = os.listdir("Coline/graphes_extension/")[fic]
                    #element2 = input("Entrer element 2 : ")
                    if element2 != element1 and "pickle" in element2 :# and "couples_possibles_"+ element1[:len(element1)-7] + "_" + element2 not in os.listdir("fichiers_couples") and "couples_possibles_"+ element2[:len(element2)-7] + "_" + element1 not in os.listdir("fichiers_couples_2") :     
                        print(element2)
                        #with open("fichiers_tot_couples_possibles.txt", 'a') as fichier_tot :
                            #fichier_tot.write(element1 + "\n")
                            #fichier_tot.write("Chaines : " + str(chaines_1) + "\n")
                        with open("graphes_extension/"+element2, 'rb') as fichier2 :
                            mon_depickler2 = pickle.Unpickler(fichier2)
                            graphe2 = mon_depickler2.load()   
                            chaines_2 = [[1]]
                            for i in range(1,5) :
                                compteur = i
                                if i != 1 : chaines_2.append([i])
                                liaison_B53 = True
                                while liaison_B53 :
                                    liaison_B53 = False
                                    temp = compteur
                                    for voisin in graphe2.successors(compteur) :
                                        for arc in graphe2[compteur][voisin] :
                                            if voisin not in [1,2,3,4] and voisin not in chaines_2[len(chaines_2)-1] and graphe2[compteur][voisin][arc]["label"] == 'B53' :
                                                liaison_B53 = True
                                                temp = voisin
                                                chaines_2[len(chaines_2)-1].append(voisin)
                                                
                                    for voisin in graphe2.predecessors(compteur) :
                                        for arc in graphe2[voisin][compteur] :
                                            if voisin not in [1,2,3,4] and voisin not in chaines_2[len(chaines_2)-1] and graphe2[voisin][compteur][arc]["label"] == 'B53' :
                                                liaison_B53 = True
                                                temp = voisin
                                                chaines_2[len(chaines_2)-1].append(voisin)
                                    compteur = temp
                                    
                            print(chaines_1)
                            print(chaines_2)
                            
                            couples_possibles = []
                            i = 0
                            j = 0
                            elt_1 = chaines_1[1][i]
                            elt_2 = chaines_2[1][j]
                            while i  < len(chaines_1[1]) :
                                elt_1 = chaines_1[1][i] 
                                while j < len(chaines_2[1]) :
                                    elt_1 = chaines_1[1][i]
                                    elt_2 = chaines_2[1][j]
                                    if graphe1.nodes[elt_1]["type"] == graphe2.nodes[elt_2]["type"] and (abs(graphe1.nodes[elt_1]["poids"] - graphe2.nodes[elt_2]["poids"]) < 2) :
                                        couples_possibles.append((elt_1, elt_2))
                                        i = i+1
                                        j = j+1
                                        
                                    else :
                                        j = j+1
                                i = i+1
                                                           
                            print(couples_possibles)           
                                    
                                    
                                    
                            