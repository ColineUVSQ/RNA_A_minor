'''
Created on 13 nov. 2018

@author: coline
'''
import os
import pickle
import time
import itertools
from itertools import combinations

import multiprocessing
 
 
#############################################################################
def exectimeout(timeout, fonc, *args, **kwargs):
    """exécute fonc(args, kwargs), avec le temps maxi timeout"""
    pool = multiprocessing.Pool(processes=1) # crée 1 processus
    poolexec = pool.apply_async(func=fonc, *args, **kwargs) # eval. asynchrone
    try:
        poolexec.get(timeout=timeout)
        result = poolexec.get()
        pool.terminate()
        return result
    except multiprocessing.TimeoutError:
        pool.terminate()
        raise # renvoi de l'exception à l'appelant
 

def recherche_couples(couples_possibles, chaines_1, chaines_2):
   
    for num_chaine in range(0, 4) :
        print(chaines_1[num_chaine])
        print(chaines_2[num_chaine])
        paires = []
        try :
            for element in itertools.product(chaines_1[num_chaine][1:],chaines_2[num_chaine][1:]):
                paires.append(element)
        except MemoryError as error :
            couples_possibles.append("memory error 1")
            break

        paires_epurees = []
        ## voir si les deux sommets qu'on superpose ont les bonnes proprietes
        for elt in paires :
            if graphe1.nodes[elt[0]]["type"] == graphe2.nodes[elt[1]]["type"] and (abs(graphe1.nodes[elt[0]]["poids"] - graphe2.nodes[elt[1]]["poids"]) < 2) :  
                paires_epurees.append(elt)
        print(len(paires_epurees))
        print(paires_epurees)

        liste = []
        for i in range(1, min(len(chaines_1[num_chaine])-1, len(chaines_2[num_chaine])-1 ) +1) :
            try :
                combinaisons = list(combinations(paires_epurees,i))
            except MemoryError as error  :
                couples_possibles.append("memory error 2")
                break
            print(len(combinaisons))
            print(combinaisons)
              
            for elt in combinaisons : ##voir si un sommet n est pas superpose avec deux autres sommets distincts
                pas_bon = False
                i = 0
                while i < len(elt) :
                    j = i + 1
                    while j < len(elt) :
                        if elt[i][0] == elt[j][0] or elt[i][1] == elt[j][1] :
                            pas_bon = True
                        j = j+1
                    i = i+1
                if pas_bon == False :
                    liste.append(elt)
              
        #print(len(liste))
        #print(liste)
         
        liste_epuree = []
        a_enlever = []
        for i in range(0,len(liste)) :
            for j in range(0, len(liste)) :
                    if j != i :
                        nb_existe_deja = 0
                        for k in range(0, len(liste[i])) :
                            for l in range(0, len(liste[j])) :
                                if liste[i][k][0] == liste[j][l][0] and liste[i][k][1] == liste[j][l][1]  :
                                    nb_existe_deja += 1
                        if nb_existe_deja == len(liste[i]) and i not in a_enlever :
                            if len(liste[i]) == len(liste[j]) :
                                a_enlever.append(j)
                            else :
                                a_enlever.append(i)
        for i in range(0, len(liste)) :
            if i not in a_enlever :
                liste_epuree.append(liste[i])
              
          
        couples_possibles.append(liste_epuree)
        #print(couples_possibles)
    return couples_possibles
        

with open("fichiers_tries.pickle", "rb") as fichier_tri :
        mon_depickler_tri = pickle.Unpickler(fichier_tri)
        tri = mon_depickler_tri.load()
    #for compte_tri_1 in range(98) :
        #element1 = tri[compte_tri_1]
        element1 = "fichier_1FJG_A_48_8.pickle"
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
                    compte_tri = 0 
                #for compte_tri in range(98) :
                    #element2 = tri[compte_tri]
                    tps1 = time.time()
                    element2 = "fichier_1FJG_A_48_11.pickle"
                    #element2 = os.listdir("Coline/graphes_extension/")[fic]
                    #element2 = input("Entrer element 2 : ")
                    if element2 != element1 and "pickle" in element2 and "couples_possibles_"+ element1[:len(element1)-7] + "_" + element2 not in os.listdir("graphes_extension/fichiers_couples_qui_manquent_epures") and "couples_possibles_"+ element2[:len(element2)-7] + "_" + element1 not in os.listdir("graphes_extension/fichiers_couples_qui_manquent_epures") :     
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
                            #fichier_tot.write(element2 + "\n")
                            #fichier_tot.write("Chaines : " + str(chaines_2) + "\n")
                        
        #                    faire = True
        #                    for elt_1 in chaines_1 :
        #                        if len(elt_1) >= 10 :
        #                            faire = False
        #                    for elt_2 in chaines_2 :
        #                        if len(elt_2) >= 10 :
        #                            faire = False

        #                    if faire :
                            couples_possibles = [] 
                            
                            t = time.time()
                            try:
                                result = exectimeout(180, recherche_couples, args=(couples_possibles, chaines_1, chaines_2))
                                
                                #print("Résultat:" + str(result))
                                
                                #fichier_tot.write("Couples possibles : " + str(result) + "\n") 
                                with open("graphes_extension/fichiers_couples_qui_manquent/couples_possibles_"+element1[:len(element1)-7]+"_"+element2[:len(element2)-7]+".pickle", "wb") as fichier_sauvegarde :
                                    mon_pickler = pickle.Pickler(fichier_sauvegarde)
                                    mon_pickler.dump(result)
                                tps2 = time.time()
                            
                                #fichier_tot.write("Temps : " + str(tps2 - tps1) + "\n")
                            except multiprocessing.TimeoutError:
                                result = None
                                #fichier_tot.write("Timeout expire : Arret au bout de %.3f sec.\n".format(time.time()-tps1))
                                print("Arrêt au bout de %.3f sec." % (time.time()-t,))
                                with open("graphes_extension/fichiers_couples_qui_manquent/couples_possibles_"+element1[:len(element1)-7]+"_"+element2[:len(element2)-7]+".pickle", "wb") as fichier_sauvegarde :
                                    mon_pickler = pickle.Pickler(fichier_sauvegarde)
                                    mon_pickler.dump(result)
                                

        
        
        
        
