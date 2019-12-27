'''
Created on 26 mars 2019

@author: coline

Processus pour lancer les superpositions (version avril 2019)

cad sans prise en compte des types de liaisons non canoniques
'''
import pickle
import os
import networkx as nx
from recup_data.extension import obtenir_extension
import multiprocessing
from itertools import product
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict

''' Renvoi de tous les ensembles de couples possibles de sommets entre les deux graphes 
avec comme critere sommets de meme type et appartenant a la meme chaine et tous les couples compatibles deux à deux par la sequence '''
def recherche_couples(couples_possibles, graphe1, graphe2, chaines_1, chaines_2):
   
    for num_chaine in range(0, 4) :
        print(chaines_1[num_chaine])
        print(chaines_2[num_chaine])
        
        if len(chaines_1[num_chaine]) == 10 and len(chaines_2[num_chaine]) == 10 :
            with open("liste_combi_9.pickle", 'rb') as fichier_liste :
                mon_depickler_liste = pickle.Unpickler(fichier_liste)
                liste_combi_9 = mon_depickler_liste.load()
                
                liste_couples_2 = []
                liste_couples = []
                for groupe in liste_combi_9 :
                    liste = []
                    for elt in groupe :
                        
                        if graphe1.nodes[chaines_1[num_chaine][elt[0]-1]]["type"] == graphe2.nodes[chaines_2[num_chaine][elt[1]-1]]["type"] : #and (abs(graphe1.nodes[elt[0]]["poids"] - graphe2.nodes[elt[1]]["poids"]) < 2) :  
                            if  (chaines_1[num_chaine][elt[0]-1], chaines_2[num_chaine][elt[1]-1]) not in liste :
                                liste.append((chaines_1[num_chaine][elt[0]-1], chaines_2[num_chaine][elt[1]-1]))
#                         print(elt)
#                         print(groupe_ok)
#                         time.sleep(2) 
                        
                    if len(liste) != 0 :
                        liste_couples.append(liste)
                        
                a_enlever = []
                for i in range(0,len(liste_couples)) :
                    for j in range(0, len(liste_couples)) :
                            if j != i :
                                nb_existe_deja = 0
                                for m in range(0, len(liste_couples[i])) :
                                    for l in range(0, len(liste_couples[j])) :
                                        if liste_couples[i][m][0] == liste_couples[j][l][0] and liste_couples[i][m][1] == liste_couples[j][l][1]  :
                                            nb_existe_deja += 1
                                if nb_existe_deja == len(liste_couples[i]) and i not in a_enlever :
                                    if len(liste_couples[i]) == len(liste_couples[j]) :
                                        a_enlever.append(j)
                                    else :
                                        a_enlever.append(i)
                for i in range(len(liste_couples)) :
                    if i not in a_enlever :
                        liste_couples_2.append(liste_couples[i])
                
                couples_possibles.append(liste_couples_2)
                    
                
        else :
            paires = []
            try :
                for element in product(chaines_1[num_chaine][1:],chaines_2[num_chaine][1:]):
                    paires.append(element)
            except MemoryError as error :
                couples_possibles.append("memory error 1")
                break
    
            paires_epurees = []
            ## voir si les deux sommets qu'on superpose ont les bonnes proprietes
            for elt in paires :
                if graphe1.nodes[elt[0]]["type"] == graphe2.nodes[elt[1]]["type"] : #and (abs(graphe1.nodes[elt[0]]["poids"] - graphe2.nodes[elt[1]]["poids"]) < 2) :  
                    paires_epurees.append(elt)
    #         print(len(paires_epurees))
    #         print(paires_epurees)
    
            liste = []
            new_couples = []
            new_chaine = []
            liste_epuree = []
            k = 0
            #for i in range(1, min(len(chaines_1[num_chaine])-1, len(chaines_2[num_chaine])-1 ) +1) :
            while min(len(chaines_1[num_chaine])-1, len(chaines_2[num_chaine])-1 )-k >= 1 :
    #             print(min(len(chaines_1[num_chaine])-1, len(chaines_2[num_chaine])-1 )-k)
    #             print(num_chaine)
                num_l = len(liste)
                try :
                    combinaisons = list(combinations(paires_epurees,min(len(chaines_1[num_chaine])-1, len(chaines_2[num_chaine])-1 )-k))
                except MemoryError as error  :
                    couples_possibles.append("memory error 2")
                    break
    #             print(len(combinaisons))
    #             print(combinaisons)
                  
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
    
    #             print(len(liste))
    #             print(liste)
                num_co = len(new_couples)
                
                for chaine in liste[num_l:] :
                    #print(chaine)
                    new_ch = []
                    for couple in chaine :
                        new_ch.append(couple)
                    new_couples.append(new_ch)
    #             print("new couples")
    #             print(len(new_couples))
    #             print(new_couples)
    
                num_ca = len(new_chaine) 
                for possib in new_couples[num_co:] :
                        #print(possib)
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
                        if len(incompatibles) == 0 :
                            new_chaine.append(possib)
        #print(max)
    #                     num_ch = len(new_chaine)
    #                     new_chaine.append([])
    #     #                 print("incomp")
    #     #                 print(incompatibles)
    #     #                 print("debut")
    #     #                 print(new_chaine)
    #                     traites = []
    #                     for i in range(len(possib)) :
    #                         #print("i")
    #                         #print(possib[i])
    #                         if possib[i] not in traites :
    #                             if possib[i] in incompatibles.keys():
    #                                 if len(incompatibles[possib[i]]) == 1 : 
    #                                     #print("ramousnif")
    #                                     new_chaine[len(new_chaine)-1].append(possib[i])
    #                                     new_chaine.append([])
    #                                     new_chaine[len(new_chaine)-1].append(incompatibles[possib[i]][0])
    #                                     ### une chaine avec le incomp, et une chaine avec son incompatible
    #                                 else :
    #                                     new_chaine.append([])
    #                                     new_chaine[len(new_chaine)-1].append(possib[i])
    #                                     for k in range(len(incompatibles[possib[i]])) :
    #                                         print("k")
    #                                         
    #                                         incomp1 = incompatibles[possib[i]][k]
    #                                         #print(incomp1)
    #                                         peut_ajouter = False
    #                                         new_possib = num_ch
    #                                         while new_possib < len(new_chaine) and peut_ajouter == False :
    #                                             peut_ajouter = True
    #                                             for elt in new_chaine[new_possib] :
    #                                                 if (incomp1 in incompatibles.keys() and elt in incompatibles[incomp1]) or (elt in incompatibles.keys() and incomp1 in incompatibles[elt]) :
    #                                                     peut_ajouter = False
    #                                             new_possib += 1
    #                                         if new_possib == len(new_chaine) :
    #                                             new_chaine.append([])
    #                                             new_chaine[len(new_chaine)-1].append(incomp1)
    #                                         else :
    #                                             new_chaine[new_possib-1].append(incomp1)
    #                                         #print(new_chaine)
    #                                         traites.append(incomp1)
    #                     for i in range(len(possib)) :
    #                         incomp = False
    #                         for cle in incompatibles.keys() :
    #                             if possib[i] in incompatibles[cle] :
    #                                 incomp = True
    #                         if possib[i] not in incompatibles.keys() and incomp == False :
    #                             for new_possib in range(num_ch, len(new_chaine)) :
    #                                 new_chaine[new_possib].append(possib[i])
    #             print("new_chaine")
    #             print(new_chaine) 
                
                
                a_enlever = []
                for i in range(0,len(new_chaine)) :
                    for j in range(0, len(new_chaine)) :
                            if j != i :
                                nb_existe_deja = 0
                                for m in range(0, len(new_chaine[i])) :
                                    for l in range(0, len(new_chaine[j])) :
                                        if new_chaine[i][m][0] == new_chaine[j][l][0] and new_chaine[i][m][1] == new_chaine[j][l][1]  :
                                            nb_existe_deja += 1
                                if nb_existe_deja == len(new_chaine[i]) and i not in a_enlever :
                                    if len(new_chaine[i]) == len(new_chaine[j]) :
                                        a_enlever.append(j)
                                    else :
                                        a_enlever.append(i)
                for i in range(num_ca, len(new_chaine)) :
                    if i not in a_enlever :
                        liste_epuree.append(new_chaine[i])
                
                k = k+1  
                
            couples_possibles.append(liste_epuree)
    for i in range(4) :
        if len(couples_possibles[i]) == 0 :
            couples_possibles[i].append([])
    return couples_possibles
        
''' Pour une paire d'extensions (notee element), recherche des numeros des sommets appartenant a chaque chaine 1 a 4  
et appel de la recherche de l'ensemble des couples pour les stocker dans un fichier pickle'''        
def recup_couples(element, taille_ext):
    print(element)
    element1 = element[0]
    element2 = element[1]
    print(element1)
    print(element2)
    # element1 = "fichier_4V9F_0_30_23.pickle"
    # element2 = "fichier_5DM6_X_48_9.pickle"
    # 
    # for fic1 in range(len(os.listdir("graphes_extension"))) :
    #         element1 = os.listdir("graphes_extension")[fic1]
    #         
    #         for fic2 in range(fic1+1, len(os.listdir("graphes_extension"))) :
    #             element2 = os.listdir("graphes_extension")[fic2]      
    if "pickle" in element1 and "pickle" in element2 and element2 != element1 and  "couples_possibles_test_"+ element1[:len(element1)-7] + "_" + element2 not in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s"%(taille_ext)) and "couples_possibles_test_"+ element2[:len(element2)-7] + "_" + element1 not in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s"%(taille_ext)) : 
                   # if "couples_possibles_"+ element1[:len(element1)-7] + "_" + element2 not in os.listdir("nouvelle_metrique") and "couples_possibles_"+ element2[:len(element2)-7] + "_" + element1 not in os.listdir("nouvelle_metrique") and "couples_possibles_"+ element1[:len(element1)-7] + "_" + element2 not in os.listdir("nouvelle_metrique") and "couples_possibles_"+ element2[:len(element2)-7] + "_" + element1 not in os.listdir("nouvelle_metrique")  and "couples_possibles_"+ element1[:len(element1)-7] + "_" + element2 not in os.listdir("nouvelle_metrique") and "couples_possibles_"+ element2[:len(element2)-7] + "_" + element1 not in os.listdir("nouvelle_metrique") :
        
        
        #element1 = "fichier_5J7L_DA_272_2.pickle"
        # with open("fichiers_tries.pickle", "rb") as fichier_tri :
        #     mon_depickler_tri = pickle.Unpickler(fichier_tri)
        #     tri = mon_depickler_tri.load()
        #     for compte_tri_1 in range(98) :
        #         element1 = tri[compte_tri_1]
        #liste_a_faire = ["fichier_1FJG_A_58_23.pickle", "fichier_1FJG_A_138_3.pickle", "fichier_4PRF_B_25_69.pickle", "fichier_4V9F_0_44_3.pickle", "fichier_4V88_A6_17_55.pickle", "fichier_4V88_A6_48_12.pickle", "fichier_5J5B_BA_48_7.pickle", "fichier_5J5B_BA_138_2.pickle", "fichier_5J7L_DA_50_21.pickle"]
        #     for k in range(len(liste_a_faire)) :
        #
                        print(element1)

                        
                        print(element2)
                        

                        with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/%s"%(taille_ext,element1), 'rb') as fichier1 :
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
                                        
            
            #                             compte_tri = 0 
            #                             for compte_tri in range(98) :
            #                                 element2 = tri[compte_tri]
            
            #                                 for j in range(k+1, len(liste_a_faire)) :
            #                                     element2 = liste_a_faire[j]
                                        #element2 = os.listdir("Coline/graphes_extension/")[fic]
                                        #element2 = input("Entrer element 2 : ")
                                        
                                        
                                                        
            #                                     if pas_bon == False and element2 != element1 and "pickle" in element2 and "couples_possibles_"+ element1[:len(element1)-7] + "_" + element2 not in os.listdir("nouvelle_metrique") and "couples_possibles_"+ element2[:len(element2)-7] + "_" + element1 not in os.listdir("nouvelle_metrique") and "couples_possibles_"+ element1[:len(element1)-7] + "_" + element2 not in os.listdir("nouvelle_metrique") and "couples_possibles_"+ element2[:len(element2)-7] + "_" + element1 not in os.listdir("nouvelle_metrique")  and "couples_possibles_"+ element1[:len(element1)-7] + "_" + element2 not in os.listdir("nouvelle_metrique") and "couples_possibles_"+ element2[:len(element2)-7] + "_" + element1 not in os.listdir("nouvelle_metrique") :     
                                                    #with open("fichiers_tot_couples_possibles.txt", 'a') as fichier_tot :
                                                        #fichier_tot.write(element1 + "\n")
                                                        #fichier_tot.write("Chaines : " + str(chaines_1) + "\n")
                                        with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/%s"%(taille_ext, element2), 'rb') as fichier2 :
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
                                        
                                            
                                            print(len(chaines_1))
                                            print(len(chaines_2))
                                            
                                            couples_possibles = recherche_couples(couples_possibles, graphe1, graphe2, chaines_1, chaines_2)
                                            print(couples_possibles)   
                                                
                                            with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/couples_possibles_test_%s_%s.pickle"%(taille_ext,element1[:len(element1)-7],element2[:len(element2)-7]), "wb") as fichier_sauvegarde :
                                                    mon_pickler = pickle.Pickler(fichier_sauvegarde)
                                                    mon_pickler.dump(couples_possibles)
                                            
#                                             try:
#                                                 if len(chaines_1) == 10 and len(chaines_2) == 10 :
#                                                     result = exectimeout(1800, recherche_couples_plus_rapide, args=(couples_possibles, chaines_1, chaines_2, graphe1, graphe2))
#                                                 else :
#                                                     result = exectimeout(1800, recherche_couples, args=(couples_possibles, chaines_1, chaines_2, graphe1, graphe2))
#                                                     
#                                                 #print("Résultat:" + str(result))
#                                                 
#                                                 #fichier_tot.write("Couples possibles : " + str(result) + "\n") 
#                                                 with open("nouvelle_metrique_reste_a_faire/couples_possibles_"+element1[:len(element1)-7]+"_"+element2[:len(element2)-7]+".pickle", "wb") as fichier_sauvegarde :
#                                                     mon_pickler = pickle.Pickler(fichier_sauvegarde)
#                                                     mon_pickler.dump(result)
#                                                 tps2 = time.time()
#                                             
#                                                 #fichier_tot.write("Temps : " + str(tps2 - tps1) + "\n")
#                                             except multiprocessing.TimeoutError:
#                                                 result = None
#                                                 #fichier_tot.write("Timeout expire : Arret au bout de %.3f sec.\n".format(time.time()-tps1))
#                                                 print("Arrêt au bout de %.3f sec." % (time.time()-t,))
#                                                 with open("nouvelle_metrique_reste_a_faire/couples_possibles_"+element1[:len(element1)-7]+"_"+element2[:len(element2)-7]+".pickle", "wb") as fichier_sauvegarde :
#                                                     mon_pickler = pickle.Pickler(fichier_sauvegarde)
#                                                     mon_pickler.dump(result)


''' idem que dans calcul_sim'''
def calcul_aretes_avec_coeff_avec_motif(graphe, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
    
    for u,v,data in graphe.edges(data=True) :
        if data["label"] != 'B53' :
            if data["label"] == '0' :
                if coeffa == 1 :
                    somme_aretes += graphe.nodes[u]["poids"]
            elif graphe.nodes[u]["type"] == 1 and graphe.nodes[v]["type"] == 1 :
                if coeffc == 1 :
                    somme_aretes += graphe.nodes[u]["poids"]
            else :
                if coeffn == 1 :
                    somme_aretes += graphe.nodes[u]["poids"]
    somme_aretes = somme_aretes/2
    
    return somme_aretes

def calcul_aretes_communes_avec_coeff_avec_motif(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
    for u,v,data in graphe_commun.edges(data=True) :
        if data["type"] != 'COV' :
            if data["type"] == '0' :
                if coeffa == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
            elif graphe1.nodes[u[0]]["type"] == 1 and graphe1.nodes[v[0]]["type"] == 1 :
                if coeffc == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
            else :
                if coeffn == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
                    
    somme_aretes = somme_aretes
    
    return somme_aretes
       

def calcul_sim_aretes_avec_coeff_avec_motif(graphe1, graphe2, graphe_commun, cle, coeffc, coeffa, coeffn):
    aretes_1 = calcul_aretes_avec_coeff_avec_motif(graphe1, cle, coeffc, coeffa, coeffn)
    aretes_2 = calcul_aretes_avec_coeff_avec_motif(graphe2, cle, coeffc, coeffa, coeffn)
    aretes_commun = calcul_aretes_communes_avec_coeff_avec_motif(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn)
    
    return aretes_commun/max(aretes_1, aretes_2)

def calcul_aretes_avec_coeff(graphe, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
    
    for u,v,data in graphe.edges(data=True) :
        if data["label"] != 'B53' :
            if data["label"] == '0' :
                if coeffa == 1 :
                    somme_aretes += graphe.nodes[u]["poids"]
            elif graphe.nodes[u]["type"] == 1 and graphe.nodes[v]["type"] == 1 :
                if coeffc == 1 :
                    somme_aretes += graphe.nodes[u]["poids"]
            else :
                if coeffn == 1 :
                    somme_aretes += graphe.nodes[u]["poids"]
    somme_aretes = somme_aretes/2 - 4
    
    return somme_aretes

def calcul_aretes_communes_avec_coeff(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
    for u,v,data in graphe_commun.edges(data=True) :
        if data["type"] != 'COV' :
            if data["type"] == '0' :
                if coeffa == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
            elif graphe1.nodes[u[0]]["type"] == 1 and graphe1.nodes[v[0]]["type"] == 1 :
                if coeffc == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
            else :
                if coeffn == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
                    
    somme_aretes = somme_aretes - 4
    
    return somme_aretes
       

def calcul_sim_aretes_avec_coeff(graphe1, graphe2, graphe_commun, cle, coeffc, coeffa, coeffn):
    aretes_1 = calcul_aretes_avec_coeff(graphe1, cle, coeffc, coeffa, coeffn)
    aretes_2 = calcul_aretes_avec_coeff(graphe2, cle, coeffc, coeffa, coeffn)
    aretes_commun = calcul_aretes_communes_avec_coeff(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn)
    
    return aretes_commun/max(aretes_1, aretes_2)


''' Recherche si le graphe_commun obtenu ne contient pas d'incompatibilite sur la sequence
cad deux paires de sommets dans un ordre dans l'un des graphes et dans l'autre dans l'autre '''
def test_compatibilite(graphe_commun, graphe1, graphe2):
    
    for i in range(len(list(graphe_commun.nodes()))) :
        for j in range(i+1, len(list(graphe_commun.nodes()))) :
            noeud1 = list(graphe_commun.nodes())[i]
            noeud2 = list(graphe_commun.nodes())[j]
            meme_chaine = True
            for elt in graphe1.nodes[noeud1[0]]["chaine"] :
                if elt not in graphe1.nodes[noeud2[0]]["chaine"] :
                    meme_chaine = False
            if noeud1 != noeud2 and meme_chaine and graphe1.nodes[noeud1[0]]["type"] != None and graphe1.nodes[noeud2[0]]["type"] != None and noeud1[0] not in [1,2,3,4,5] and noeud2[0] not in [1,2,3,4,5]:
                if (graphe1.nodes[noeud1[0]]["position"] < graphe1.nodes[noeud2[0]]["position"] and graphe2.nodes[noeud1[1]]["position"] > graphe2.nodes[noeud2[1]]["position"]) or (graphe1.nodes[noeud1[0]]["position"] > graphe1.nodes[noeud2[0]]["position"] and graphe2.nodes[noeud1[1]]["position"] < graphe2.nodes[noeud2[1]]["position"]) :
                    return False
    
    return True

        
''' Recherche si l'un des sommets du couple a chercher ne se trouve pas deja dans un autre couple du graphe (commun) '''
def dans_graphe(graphe, couple_a_chercher):
    for noeud in graphe.nodes() :
        #print(noeud)
        if (noeud[0] == couple_a_chercher[0] and noeud[1] != couple_a_chercher[1]) or (noeud[1] == couple_a_chercher[1] and noeud[0] != couple_a_chercher[0])  :
            return True
    return False

''' Ajout ou non d'un nouvel element au graphe commun et de ses potentiels voisins non covalents
criteres : aucun des elements de la paire de sommets n'a encore ete ajoute, 
si des aretes a ajouter, doivent etre de meme type '''
def comparaison(elt, graphe1, graphe2, elt_1, elt_2, elt_3, couples_possibles, num, graphe_commun, compt) :

    deja_vus_voisin_1 = []
    deja_vus_voisin_2 = []
    sommets_en_plus_1 = []
    sommets_en_plus_2 = []
    tab_temp = []
    
#     succ_1 = True
#     if len(elt) > 1 :
#         if elt[1][0] in graphe1.successors(elt[0][0]) :
#             succ_1 = True
#         else :
#             succ_1 = False
#             
#     succ_2 = True
#     if len(elt) > 1 :
#         if elt[1][1] in graphe1.successors(elt[0][1]) :
#             succ_2 = True
#         else :
#             succ_2 = False
    if num == 2 or num == 3 :
        succ_1 = True
        succ_2 = True
    else :
        succ_1 = False
        succ_2 = False

    #print()
    #print(succ_1)
    #print(succ_2)
        
    for i in range(0, len(elt)) :
        #print(elt_1[i])
        if elt[i][0] in graphe1.nodes() and elt[i][1] in graphe2.nodes() :
        
            if dans_graphe(graphe_commun, elt[i]) == False :
                if elt[i] not in graphe_commun.nodes() :
                    graphe_commun.add_node(elt[i])  
            
            if succ_1 == True :
                voisins_1 = graphe1.successors(elt[i][0])
            else :
                voisins_1 = graphe1.predecessors(elt[i][0])            
            for voisin_1 in voisins_1 :
                    #print(voisin_1)
                    if succ_2 == True :              
                        voisins_2 = graphe2.successors(elt[i][1])
                    else :
                        voisins_2 = graphe2.predecessors(elt[i][1])
                    
                    
                    for voisin_2 in voisins_2 :
                        #print("elt_1")
                        #print(elt_1[i])
                        #print(voisin_1)
                        #print(voisin_2)
                        
                        if succ_1 == True :
                            dict_voisin_1 = graphe1[elt[i][0]][voisin_1]
                        else :
                            dict_voisin_1 = graphe1[voisin_1][elt[i][0]]
                        
                        if succ_2 == True :
                            dict_voisin_2 = graphe2[elt[i][1]][voisin_2]
                        else :
                            dict_voisin_2 = graphe2[voisin_2][elt[i][1]]
                        #print(dict_voisin_1)
                        #print(dict_voisin_2)
                        
                        for edge_1 in dict_voisin_1 :
                            for edge_2 in dict_voisin_2 :
                                if dict_voisin_1[edge_1]["label"] == 'CWW' :
                                    label_1 = 'CAN'
                                elif dict_voisin_1[edge_1]["label"] == 'B53' : 
                                    label_1 = 'COV'
                                elif dict_voisin_1[edge_1]["label"] == '0' :
                                    label_1 = 'ART'
                                else :
                                    label_1 = 'NON_CAN'  
                
                                
                                if dict_voisin_2[edge_2]["label"] == 'CWW' :
                                    label_2 = 'CAN'
                                elif dict_voisin_2[edge_2]["label"] == 'B53' : 
                                    label_2 = 'COV'
                                elif dict_voisin_2[edge_2]["label"] == '0' :
                                    label_2 = 'ART'
                                else :
                                    label_2 = 'NON_CAN'
                                
                                bonne_chaine = False
                                for ch in graphe1.nodes[voisin_1]["chaine"] :
                                    if ch in graphe2.nodes[voisin_2]["chaine"] :
                                        bonne_chaine = True
                                
                                if bonne_chaine :
                                    if label_1 == label_2 and dict_voisin_1[edge_1]["long_range"] == dict_voisin_2[edge_2]["long_range"] :
                                        if label_1 == 'COV' :
                                            #if (graphe1.nodes[elt[i][0]]["position"][0] - graphe1.nodes[voisin_1]["position"][0] < 0 and graphe1.nodes[1]["position"][0] - graphe1.nodes[3]["position"][0] < 0) or (graphe1.nodes[elt[i][0]]["position"][0] - graphe1.nodes[voisin_1]["position"][0] > 0 and graphe1.nodes[1]["position"][0] - graphe1.nodes[3]["position"][0] > 0) : 
                                            #    if (graphe2.nodes[elt[i][1]]["position"][0] - graphe2.nodes[voisin_2]["position"][0] < 0 and graphe2.nodes[1]["position"][0] - graphe2.nodes[3]["position"][0] < 0) or (graphe2.nodes[elt[i][1]]["position"][0] - graphe2.nodes[voisin_2]["position"][0] > 0 and graphe2.nodes[1]["position"][0] - graphe2.nodes[3]["position"][0] > 0) : 
                                            if (succ_1 == True and succ_2 == True) or (succ_1 == False and succ_2 == False) :# dans le meme sens dans les deux graphes
                                                    if voisin_1 == num and voisin_2 == num : # elt courant est lie au motif
                                                        if dans_graphe(graphe_commun, elt[i]) == False :
                                                            if elt[i] not in graphe_commun.nodes() :
                                                                graphe_commun.add_node(elt[i])
                                                            deja_ajoute = False
                                                            for voisin in graphe_commun[elt[i]] :
                                                                if voisin == (voisin_1, voisin_2) :
                                                                    for edge in graphe_commun[elt[i]][voisin] :
                                                                        if graphe_commun[elt[i]][voisin][edge]["type"] == label_1 :
                                                                            deja_ajoute = True
                                                            if deja_ajoute == False :
                                                                graphe_commun.add_edge((num,num),elt[i], type='COV', long_range=False)
                                                        #nb_aretes += 1
                                                        #print("ramou")
                                                    else : #elt courant non lie au motif
                                                        for couple in elt : # on cherche si un couple deja dans la liste ne correspond pas aux voisins
                                                            if dans_graphe(graphe_commun, elt[i]) == False and voisin_1 == couple[0] and voisin_2 == couple[1]  and elt[i][0] not in deja_vus_voisin_1 and elt[i][1] not in deja_vus_voisin_2 :
                                                                deja_vus_voisin_1.append(voisin_1)
                                                                deja_vus_voisin_2.append(voisin_2)
        #                                                         nb_aretes += 1
                                                                if elt[i] not in graphe_commun.nodes() :
                                                                    graphe_commun.add_node(elt[i])
                                                                deja_ajoute = False
                                                                for voisin in graphe_commun[elt[i]] :
                                                                    if voisin == (voisin_1, voisin_2) :
                                                                        for edge in graphe_commun[elt[i]][voisin] :
                                                                            if graphe_commun[elt[i]][voisin][edge]["type"] == label_1 :
                                                                                deja_ajoute = True
                                                                if deja_ajoute == False :
                                                                    graphe_commun.add_edge(elt[i], (voisin_1, voisin_2), type='COV', long_range=False)
                                                                #print("petit rat1")
                                                                #print("ramou2")   
                                        else : ## lies par une liaison non covalente
                                            ## a ce moment-la on regarde si les sommets voisins pourraient se superposer
                                            if graphe1.nodes[voisin_1]["type"] == graphe2.nodes[voisin_2]["type"] :
                                                #print("ramousnif")
                                                pas_vu_1_motif = True
                                                pas_vu_2_motif = True
                                                for j in range(1,6) : ## on regarde si l'un des voisins n'est pas dans le motif
                                                    if voisin_1 == j and voisin_2 == j : #si les deux sont dans le motif c'est bon
                                                        #nb_aretes += 1
                                                        if dans_graphe(graphe_commun, elt[i]) == False :
                                                            if elt[i] not in graphe_commun.nodes() :
                                                                graphe_commun.add_node(elt[i])
                                                            deja_ajoute = False
                                                            for voisin in graphe_commun[elt[i]] :
                                                                if voisin == (voisin_1, voisin_2) :
                                                                    for edge in graphe_commun[elt[i]][voisin] :
                                                                        if graphe_commun[elt[i]][voisin][edge]["type"] == label_1 :
                                                                            deja_ajoute = True
                                                            if deja_ajoute == False :
                                                                graphe_commun.add_edge(elt[i], (voisin_1, voisin_2), type=label_1, long_range = dict_voisin_1[edge_1]["long_range"])
                                                    elif voisin_1 == j :
                                                        pas_vu_1_motif = False
                                                    elif voisin_2 == j :
                                                        pas_vu_2_motif = False
                    #                             if pas_vu_1 and pas_vu_2 and elt[i][0] not in deja_vus_voisin_1 and elt[i][1] not in deja_vus_voisin_2 and voisin_1 not in sommets_en_plus_1 and voisin_2 not in sommets_en_plus_2 :
                    #                                 print("ajout arete")
                    #                                 print(voisin_1)
                    #                                 print(voisin_2)
                    #                                 deja_vus_voisin_1.append(voisin_1)
                    #                                 deja_vus_voisin_2.append(voisin_2)
                    #                                 deja_vus_voisin_1.append(elt[i][0])
                    #                                 deja_vus_voisin_2.append(elt[i][1])
                    #                                 nb_aretes += 1
                                                pas_vu_1 = True
                                                pas_vu_2 = True
                                                for couple in elt :
                                                    if voisin_1 == couple[0] :
                                                        pas_vu_1 = False
                                                    if voisin_2 == couple[1] : 
                                                        pas_vu_2 = False 
                                                for couple in elt_1 :
                                                    if voisin_1 == couple[0] :
                                                        pas_vu_1 = False
                                                    if voisin_2 == couple[1] : 
                                                        pas_vu_2 = False 
                                                for couple in elt_2 :
                                                    if voisin_1 == couple[0] :
                                                        pas_vu_1 = False
                                                    if voisin_2 == couple[1] : 
                                                        pas_vu_2 = False 
                                                for couple in elt_3 :
                                                    if voisin_1 == couple[0] :
                                                        pas_vu_1 = False
                                                    if voisin_2 == couple[1] : 
                                                        pas_vu_2 = False 
        #                                         if voisin_1 == 9 and voisin_2 == 9 :
        #                                             print(graphe_commun.nodes.data())
        #                                             print(dans_graphe(graphe_commun, elt[i]))
                                                if dans_graphe(graphe_commun, elt[i]) == False and dans_graphe(graphe_commun, (voisin_1, voisin_2)) == False :
                                                    
                                                    if pas_vu_1_motif and pas_vu_2_motif and elt[i][0] not in deja_vus_voisin_1 and elt[i][1] not in deja_vus_voisin_2 and voisin_1 not in sommets_en_plus_1 and voisin_2 not in sommets_en_plus_2 : ## on ajoute un sommet si le sommet n'existe nulle part
                                                        if pas_vu_1 and pas_vu_2 and test_compatibilite(graphe_commun, graphe1, graphe2) : #les sommets n'existent pas, on les rajoute
                                                            #nb_sommets += 1
                            #                                                     print("petit ramousnif")
                            #                                                     print(voisin_1)
                            #                                                     print(voisin_2)
                                                            #print("petit rat")
                                                            #print(voisin_1)
                                                            #print(voisin_2)
                                                            
                                                            sommets_en_plus_1.append(voisin_1)
                                                            sommets_en_plus_2.append(voisin_2)
                                                            
                                                            graphe_commun.add_node((voisin_1, voisin_2))
                                                            if (voisin_1, voisin_2) != (5,5) :
                                                                tab_temp.append((voisin_1, voisin_2))
                                                        if elt[i] not in graphe_commun.nodes() :
                                                            graphe_commun.add_node(elt[i])
                                                        deja_ajoute = False
                                                        for voisin in graphe_commun[elt[i]] :
                                                            if voisin == (voisin_1, voisin_2) :
                                                                for edge in graphe_commun[elt[i]][voisin] :
                                                                    if graphe_commun[elt[i]][voisin][edge]["type"] == label_1 :
                                                                        deja_ajoute = True
                                                        if deja_ajoute == False :
                                                            graphe_commun.add_edge(elt[i], (voisin_1, voisin_2), type=label_1, long_range = dict_voisin_1[edge_1]["long_range"])   
                                        #couples_possibles[0].append(tab_temp) 
                #print(nb_aretes)
                #print(nb_sommets)
    couples_en_plus = []
    couples_en_plus.extend(elt)
    couples_en_plus.extend(tab_temp)
    #print(couples_en_plus)
    #print(couples_possibles[0])
    max_deja = 0
    for groupe in couples_possibles[num-1] :
        deja = 0
        for couple in groupe :
            for couple_2 in couples_en_plus :
                if couple[0] == couple_2[0] and couple[1] == couple_2[1] :
                    deja += 1
        if max_deja < deja :
            max_deja = deja
        #print(deja)
        
    if max_deja < len(couples_en_plus) :
        couples_possibles[num-1][compt] = list(couples_en_plus)
        #del(couples_possibles[num-1][compt])
        #print(compt)
    
#     deja_vus_voisin_1 = list(set(deja_vus_voisin_1))
#     tab_sommets_aretes =[]
#     tab_sommets_aretes.append(nb_sommets + len(deja_vus_voisin_1))
#     tab_sommets_aretes.append(nb_aretes)
    #print("ramousnif")
    #print(graphe_commun.nodes.data())
    #print(graphe_commun.edges.data())
    #print(graphe_commun.nodes.data())
    return graphe_commun 

''' Construction du graphe commun a partir de deux graphes d'extensions et 
des ensembles de couples de sommets possibles calcules par ailleurs '''
def sous_graphe_commun_max(fichier, graphe1, graphe2, taille_ext):
    c = 0
    

    dico_graphes = {}
#                         if "pickle" in fic :
    fic = fichier
    if "pickle" in fic :
        element1 = fic.split('_')[3] + '_' + fic.split('_')[4] + '_' + fic.split('_')[5] + '_' + fic.split('_')[6] + '_' + fic.split('_')[7]
        element2 = fic.split('_')[8] + '_' + fic.split('_')[9] + '_' + fic.split('_')[10] + '_' + fic.split('_')[11] + '_' + fic.split('_')[12][:len(fic.split('_')[12])-7]
       #element1 = "fichier_4V9F_0_134_5"
       #element2 = "fichier_5FDU_1A_272_1"
        print(element1)
        print(element2)
        if "graphe_comp_test_"+fic not in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s"%(str(taille_ext))) : 
        
                    c = c+1
                    print(c)
#             with open("graphes_extension/"+element1+".pickle", 'rb') as fichier1 :
#                 mon_depickler1 = pickle.Unpickler(fichier1)
#                 graphe1 = mon_depickler1.load()     
#                 with open("graphes_extension/"+element2+".pickle", 'rb') as fichier2 :
#                     mon_depickler2 = pickle.Unpickler(fichier2)
#                     graphe2 = mon_depickler2.load()
    
    #                                     compteur_arc = 0
    #                                     for (u, v, keys, t) in graphe1.edges(data="label", keys = True) :
    #                                         if t == "B53" :
    #                                             compteur_arc += 1
    #                                     compteur_arc_arete_1 = compteur_arc + (graphe1.number_of_edges() - compteur_arc)/2
    #                                     
    #                                     compteur_arc = 0
    #                                     for (u, v, keys, t) in graphe2.edges(data="label", keys = True) :
    #                                         if t == "B53" :
    #                                             compteur_arc += 1
    #                                     compteur_arc_arete_2 = compteur_arc + (graphe2.number_of_edges() - compteur_arc)/2
                                                               
                   
                   #with open("graphes_extension/fichiers_couples_a_faire/"+fic, 'rb') as fichier_pickle :
                   
                    memory_error = False
                    
                    graphe_motif = nx.MultiGraph()
                    for i in range(1,6) :
                        graphe_motif.add_node((i,i))
                    graphe_motif.add_edge((1,1),(2,2), type="NON_CAN", long_range=True)
                    graphe_motif.add_edge((1,1),(3,3), type="COV", long_range=False)
                    graphe_motif.add_edge((1,1),(5,5), type="NON_CAN", long_range=True)
                    graphe_motif.add_edge((2,2),(4,4), type="COV", long_range=False)
                    graphe_motif.add_edge((2,2),(5,5), type="CAN", long_range=False)
                    graphe_motif.add_edge((3,3),(4,4), type="NON_CAN", long_range=True)
                    
                 
                    with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/%s"%(taille_ext, fic), 'rb') as fichier_pickle :                        
                        mon_depickler = pickle.Unpickler(fichier_pickle)
                        couples_possibles = mon_depickler.load()   
                   #couples_possibles = [[[(31,35),(32,36),(40,44),(42,46),(43,40)]],[[(9,8),(12,12),(18,16),(20,18)]],[[(44,48)]],[[(27,20),(28,21),(29,22)]]]
                        print(fic)
                        print(couples_possibles)
                        
                        #couples_possibles[0] = [[(25,23),(29,26),(31,29)]]
                        print(couples_possibles)
                        if couples_possibles != None :
                            print(couples_possibles)
                            
                            couples_possibles_new = []
                            for i in range(4) :
                                if 'memory error' in couples_possibles[i] :
                                    memory_error = True
                                    break
                                try :
                                    for chaine in couples_possibles[i] :
                                        chaine.insert(0, (i+1, i+1))
                                    couples_possibles_new.append(couples_possibles[i])
                                except AttributeError :
                                    new_couples = []
                                    for chaine in couples_possibles[i] :
                                        new_chaine = []
                                        for couple in chaine :
                                            new_chaine.append(couple)
                                        new_couples.append(new_chaine)
                                    couples_possibles_new.append(new_couples)
                            print(couples_possibles_new)
                            
                            #couples_possibles_new = [[[(1, 1), (25, 23), (29, 26), (31, 29)]], [[(2, 2), (6, 6), (8, 7), (11, 11)]], [[(3, 3), (33, 30), (34, 31), (28, 25)]], [[(4, 4), (15, 13), (16, 14), (18, 18), (20, 20), (23, 22)]]]
                            
#                                
                            if memory_error == False :
                                   #print("ramousnif")
                                   
                                    graphe_commun_max = nx.MultiGraph()
                                    graphe_commun_max = graphe_motif.copy()
                                   #result = exectimeout(50, toutes_comparaisons, args=(couples_possibles, graphe1, graphe2, graphe_motif, graphe_commun_max))
                                    maxi = 0
                                    couples_max = []
                                   
                                    print(couples_possibles_new)
                                    for compt_1 in range(max(len(couples_possibles_new[0]), 1)) :
                                        if len(couples_possibles_new[0]) > 0 :
                                            elt_1 = couples_possibles_new[0][compt_1]
                                        else :
                                            elt_1 = []
                                        for compt_2 in range(max(len(couples_possibles_new[1]), 1)) :
                                            if len(couples_possibles_new[1]) > 0 :
                                                elt_2 = couples_possibles_new[1][compt_2]
                                            else :
                                                elt_2 = []
                                            for compt_3 in range(max(len(couples_possibles_new[2]), 1)) :
                                                if len(couples_possibles_new[2]) > 0 :
                                                    elt_3 = couples_possibles_new[2][compt_3]
                                                else :
                                                    elt_3 = [] 
                                                for compt_4 in range(max(len(couples_possibles_new[3]), 1)) :
                                                    if len(couples_possibles_new[3]) > 0 :
                                                        elt_4 = couples_possibles_new[3][compt_4]
                                                    else :
                                                        elt_4 = []
                                                    
                                                    graphe_commun = graphe_motif.copy()
                                                   
                                                    graphe_commun = comparaison(elt_1, graphe1, graphe2, elt_2, elt_3, elt_4, couples_possibles_new, 1, graphe_commun, compt_1)
                                                    graphe_commun = comparaison(elt_2, graphe1, graphe2, elt_1, elt_3, elt_4, couples_possibles_new, 2, graphe_commun, compt_2)
                                                   #print("petit rat")
                                                   #print(graphe_commun)
                                                    graphe_commun = comparaison(elt_3, graphe1, graphe2, elt_1, elt_2, elt_4, couples_possibles_new, 3, graphe_commun, compt_3)                            
                                                    graphe_commun = comparaison(elt_4, graphe1, graphe2, elt_1, elt_2, elt_3, couples_possibles_new, 4, graphe_commun, compt_4)
                                                    print(graphe_commun[(4,4)])
                                                    if test_compatibilite(graphe_commun, graphe1, graphe2) :
                                                        sim = calcul_sim_aretes_avec_coeff(graphe1, graphe2, graphe_commun, fic, 1, 1, 1)
                                                        #sim = calcul_sim_avec_poids(graphe1, graphe2, graphe_commun, fic)
                                                    
                                                    
                                                    
        #                                                                 if round(sim,2) == 0.22 :
        #                                                                     print(elt_1)
        #                                                                     print(elt_2)
        #                                                                     print(elt_3)
        #                                                                     print(elt_4)
                                                   #time.sleep(2)
                                                   
                                                   #sim = ((graphe_commun.number_of_nodes() + graphe_commun.number_of_edges())*(graphe_commun.number_of_nodes() + graphe_commun.number_of_edges()))/((graphe1.number_of_nodes()+compteur_arc_arete_1)*(graphe2.number_of_nodes()+compteur_arc_arete_2))
                                                        if maxi < sim and sim <= 1.0 :
                                                            print(sim)
                                                            print(elt_1)
                                                            maxi = sim
                                                            del(couples_max[:])
                                                            couples_max.append(elt_1)
                                                            couples_max.append(elt_2)
                                                            couples_max.append(elt_3)
                                                            couples_max.append(elt_4)
                                                            graphe_commun_max = graphe_commun.copy()   
                                  
                                    print("maxi")
                                    print(maxi)
                                    print("couples max")
                                    print(couples_max)
                                    print("graphe")
                                    print(graphe_commun_max.nodes.data())
                                    print(graphe_commun_max.edges.data())
                                    print("couples possibles")
                                    print(couples_possibles)
                                   
                                    
#                                     sim_max.append(maxi)
                                    dico_graphes.update({(element1, element2) : graphe_commun_max})
                                    
                                    print(couples_possibles)
                                    
                                    with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/graphe_comp_test_%s"%(str(taille_ext),fic), 'wb') as fichier_graphes :
                                        mon_pickler_3 = pickle.Pickler(fichier_graphes)
                                        mon_pickler_3.dump(dico_graphes) 
                                    
                                    print("nombre")
                                    print(c)
                                    print("temps")
                                    
                                    print(graphe_commun_max)
                                    
                                    
                                    
#                                     with open("fichier_max_nouvelle_metrique.pickle", 'wb') as fichier_max_pickle :
#                                         mon_pickler = pickle.Pickler(fichier_max_pickle)
#                                         mon_pickler.dump(sim_max)
                                       
                                   #if element1 == "fichier_4V9F_0_62_12" or element2 == "fichier_4V9F_0_62_12" :
#                                     with open("fichier_similarite_nouvelle_metrique_%s.txt"%(str(taille_ext)), 'a') as fichier_ecriture :    
#                                         fichier_ecriture.write(element1 + '\n' + element2 + '\n')
#                                         fichier_ecriture.write("Graphe : ")
#                                         fichier_ecriture.write(str(graphe_commun_max.nodes.data())+ '\n')
#                                         fichier_ecriture.write(str(graphe_commun_max.edges.data())+ '\n')
#                                         fichier_ecriture.write("Similarite :" + str(maxi))
#                                         fichier_ecriture.write('\n\n')  
                                    
                                    return graphe_commun_max, maxi
                                   #with open("tab_calcules_sim_"+nom_extension+".pickle", 'wb') as fichier_deja_fait :
#                                     with open("tab_calcules_sim_nouvelle_metrique.pickle", 'wb') as fichier_deja_fait :
#                                         mon_pickler_2 = pickle.Pickler(fichier_deja_fait)
#                                         mon_pickler_2.dump(tab_deja_fait)
#                                         with open("dico_graphes_communs_max_nouvelle_metrique.pickle", 'wb') as fichier_graphes :
#                                             mon_pickler_3 = pickle.Pickler(fichier_graphes)
#                                             mon_pickler_3.dump(dico_graphes)

''' Appel de la recherche des ensembles de couples puis de la recherche du graphe commun max 
pour une paire d'extensions, sur toutes les tailles d'extension de 3 a 10 '''
def calcul_par_taille(couple):
    
    liste_a_faire = [[], [], [],[], [], [], ['5FJC_A_138_1', '5DM6_X_48_9', '1FJG_A_271_1', '5J5B_BA_48_14', '4V9F_0_62_12', '3JCS_1_25_16', '5FDU_1A_137_6', '1GID_B_25_20', '5FDU_1A_272_1', '5J7L_DA_30_15', '1FJG_A_48_11', '5DM6_X_48_10', '5J7L_DA_25_10', '5FDU_1A_134_3', '5J7L_DA_272_2', '4V9F_0_30_4', '4V88_A6_48_12', '5J7L_DA_197_4', '4V88_A6_17_55', '4V9F_0_48_13', '3UCZ_R_62_15', '4V9F_0_207_3', '3JCS_1_48_18', '4PRF_B_25_69', '5FDU_1A_62_14', '4YAZ_R_36_25', '5DM6_X_134_2', '5DM6_X_227_2', '5FDU_1A_25_78', '4V9F_0_48_26', '1FJG_A_48_17', '5FDU_1A_25_68', '5FDU_1A_48_25', '1FJG_A_48_8', '5J5B_BA_138_2', '5J7L_DA_30_36', '4V88_A5_48_3', '4V9F_0_25_4', '5FDU_1A_197_3', '4V88_A5_287_1', '3JCS_1_282_1', '3JCS_2_30_56', '5J7L_DA_170_5', '4V9F_0_44_3', '4V88_A5_25_30', '5J7L_DA_134_1', '1U9S_A_25_74', '5J7L_DA_48_15', '5DM6_X_127_7', '5J5B_BA_48_23', '1FJG_A_294_1', '5J7L_DA_62_5', '5J5B_BA_294_2', '4V9F_0_48_2', '4V9F_0_48_16', '5J7L_DA_50_21', '5FDU_1A_48_19', '5J7L_DA_48_20', '4V9F_0_25_56', '5J5B_BA_58_3', '5DM6_X_25_15', '5J7L_DA_36_3', '1FJG_A_58_23', '5J5B_BA_48_5', '4V9F_0_137_5', '4V9F_0_134_5', '1FJG_A_138_3', '4V88_A5_25_47', '2XD0_V_36_21', '5J5B_BA_48_7', '4V88_A5_48_6', '5DM6_X_25_34', '5DM6_X_328_2', '5FDU_1A_290_2', '5J7L_DA_48_1', '4V9F_0_48_21', '5J7L_DA_25_12', '4FAU_A_207_4', '1FJG_A_109_6', '5FDU_1A_74_7', '4V88_A5_30_5', '3JCS_1_25_46', '5DM6_X_48_28', '1U9S_A_58_11', '5FDU_1A_30_17', '4V9F_0_127_6', '4L81_A_25_77', '4V9F_0_287_2'], [], []]
    u = couple[0]
    v = couple[1]
#     for u,v, data in graphe_comp.edges(data=True) :
    with open("Extensions/Metrique_toutes_aretes/graphe_complet_pondere_sim_coeff_n1_a1_c1.pickle", 'rb') as fichier_graphe_complet :
        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
        graphe_comp = mon_depickler_complet.load()
        
        dico_elts = {}
        #if graphe_comp.edges[u,v]["poids"] > 0.5 :
        for i in range(3,10) :
                for u,v in graphe_comp.edges() :
                    if graphe_comp.nodes[u]["nom"] in liste_a_faire[i-1] or graphe_comp.nodes[v]["nom"] in liste_a_faire[i-1] :

                        print("i : %d"%(i))
                        with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/fichier_%s.pickle"%(str(i), graphe_comp.nodes[u]["nom"]), 'rb') as fichier1 :
                            mon_depickler_graphe1 = pickle.Unpickler(fichier1)
                            graphe1 = mon_depickler_graphe1.load()
                             
                            with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/fichier_%s.pickle"%(str(i), graphe_comp.nodes[v]["nom"]), 'rb') as fichier2 :
                                mon_depickler_graphe2 = pickle.Unpickler(fichier2)
                                graphe2 = mon_depickler_graphe2.load()
                                element = ("fichier_"+graphe_comp.nodes[u]["nom"]+".pickle", "fichier_"+graphe_comp.nodes[v]["nom"]+".pickle")
                                recup_couples(element, i)
                                
                                fichier = "couples_possibles_test_%s_%s.pickle"%(element[0][:len(element[0])-7], element[1][:len(element[1])-7])
        #                         for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/couples_possibles_%s_%s"%(i, element[0], element[1])) :
        #                             if graphe_comp.nodes[v]["nom"] in fic and graphe_comp.nodes[u]["nom"] in fic :
        #                                 fichier = fic
        #                         if fichier == None :
        #                             print("probleme")  
        #                             break
                                print("ramousnif")
                                if "graphe_comp_test_"+fichier not in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s"%(str(i))) : 
                                    graphe_commun, sim = sous_graphe_commun_max(fichier, graphe1, graphe2, i)
                                else :
                                    with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/graphe_comp_test_%s"%(i, fichier), 'rb') as fichier_commun:
                                        mon_depickler_graphe_commun = pickle.Unpickler(fichier_commun)
                                        graphe_commun = mon_depickler_graphe_commun.load()
                                        for cle in graphe_commun.keys() :
                                            sim = calcul_sim_aretes_avec_coeff(graphe1, graphe2, graphe_commun[cle], cle, 1,1,1) 
                                #if sim > data["poids"] - data["poids"]/10 and sim < data["poids"] + data["poids"]/10 :
                                if (graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]) not in dico_elts.keys() :
                                    dico_elts.update({(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]) : [sim]})
                                else :
                                    dico_elts[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])].append(sim)     
                                print(sim)
     
#         with open("crible_taille_extensions_%s_%s.pickle"%(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]), 'wb') as fichier_pickle :
#             mon_pickler = pickle.Pickler(fichier_pickle)
#             mon_pickler.dump(dico_elts)

''' Etude des resultats par cluster '''         
def resultats():    
    petit_cluster = ['4V9F_0_30_23', '4V9F_0_48_16', '1FJG_A_271_1', '1FJG_A_109_6']     
    
    with open("Extensions/Metrique_toutes_aretes/couples_possibles_bonne_version/couples_possibles_fichier_1FJG_A_109_6_fichier_4V9F_0_30_23.pickle", 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        couples_possibles = mon_depickler.load()
        
        print(couples_possibles)
    
    for i in range(1,10) :
        graphe_cluster = nx.Graph()   
        for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s"%(i)) :
            if "comp" in fic :
                elt1 = fic.split("_")[5] + "_" + fic.split("_")[6] + "_" + fic.split("_")[7] + "_" + fic.split("_")[8]
                elt2 = fic.split("_")[10] + "_" + fic.split("_")[11] + "_" + fic.split("_")[12] + "_" + fic.split("_")[13][:len(fic.split("_")[13])-7]
                compteur1 = -1
                compteur2 = -1
                
                if elt1 in petit_cluster and elt2 in petit_cluster :
                    deja_1 = False
                    deja_2 = False
                    for noeud, nom in graphe_cluster.nodes(data="nom") :
                        if nom == elt1 :
                            deja_1 = noeud
                        if nom == elt2 :
                            deja_2 = noeud 
                    if deja_1 == False : 
                        compteur1 = graphe_cluster.number_of_nodes()
                        graphe_cluster.add_node(compteur1, nom=elt1)
                    else :
                        compteur1 = deja_1
                    if deja_2 == False : 
                        compteur2 = graphe_cluster.number_of_nodes()
                        graphe_cluster.add_node(compteur2, nom=elt2)
                    else :
                        compteur2 = deja_2
                    
                    deja_edge = False
                    if deja_1 and deja_2 :
                        for u,v in graphe_cluster.edges() :
                            if (graphe_cluster.nodes[u]["nom"] == elt1 and graphe_cluster.nodes[v]["nom"] == elt2) or (graphe_cluster.nodes[v]["nom"] == elt1 and graphe_cluster.nodes[u]["nom"] == elt2) :
                                deja_edge = True
     
                    if deja_edge == False :
                        with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/%s"%(i, fic), 'rb') as fichier_commun :
                            mon_depickler = pickle.Unpickler(fichier_commun)
                            graphe_commun = mon_depickler.load()
                            
                            
                            with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/fichier_%s.pickle"%(i, elt1), 'rb') as fichier_graphe1 :
                                mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
                                graphe1 = mon_depickler_1.load()
                                
                                with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/fichier_%s.pickle"%(i, elt2), 'rb') as fichier_graphe2 :
                                    mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
                                    graphe2 = mon_depickler_2.load()
                                    
                                    for cle in graphe_commun.keys() :
                                        print(cle)
                                        print(elt1)
                                        print(elt2)
                                        print(graphe_commun[cle].nodes.data())
                                        print(graphe1.nodes.data())
                                        print(graphe2.nodes.data())
                                        print(i)
                                        sim = calcul_sim_aretes_avec_coeff(graphe1, graphe2, graphe_commun[cle], fic, 1,1,1)
    
                                    graphe_cluster.add_edge(compteur1, compteur2, poids = sim)
                                
                    
        with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/graphe_cluster.pickle"%(i), 'wb') as fichier_graphe:
            mon_pickler = pickle.Pickler(fichier_graphe) 
            mon_pickler.dump(graphe_cluster)                      

''' Etude des resultats des tailles d'extension inferieures a 10 par rapport a ceux de k=10 '''
def resultats_entier():
    liste_10pourcent = []
    for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles") :
        if "crible" in fic :
            elt1 = fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6]
            elt2 = fic.split("_")[7] + "_" + fic.split("_")[8] + "_" + fic.split("_")[9] + "_" + fic.split("_")[10][:len(fic.split("_")[10])-7]
            
            with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/%s"%(fic), 'rb') as fichier :
                mon_depickler = pickle.Unpickler(fichier)
                crible = mon_depickler.load()
                
                with open("Extensions/Metrique_toutes_aretes/graphe_complet_pondere_sim_coeff_n1_a1_c1.pickle", 'rb') as fichier_graphe_complet :
                        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
                        graphe_comp = mon_depickler_complet.load()
                         
                        for u,v, sim in graphe_comp.edges(data="poids") :
                            if (graphe_comp.nodes[u]["nom"] == elt1 and graphe_comp.nodes[v]["nom"] == elt2) or (graphe_comp.nodes[u]["nom"] == elt2 and graphe_comp.nodes[v]["nom"] == elt1) :
                 
                                for cle in crible :
                                    crible[cle].append(sim)
                
                for cle in crible :
                    print(fic)
                    print(len(crible[cle]))
                    print(crible)
                    compteur = 8
                    while compteur >= 0 and (crible[cle][compteur] > crible[cle][9] - crible[cle][9]/10 and crible[cle][compteur] < crible[cle][9] + crible[cle][9]/10) :
                        compteur -= 1
                    liste_10pourcent.append(compteur)
    print(liste_10pourcent)
    
    nombre_par_valeur = [0]*10
    for i in range (1,9) :
        for val in liste_10pourcent :
            if val == i :
                nombre_par_valeur[i] += 1
    
    print(nombre_par_valeur)
    
    
    plt.hist(liste_10pourcent, range = (1, 9), bins = 9, color = 'yellow', edgecolor = 'red')
    plt.xlabel("Taille d'extension")
    plt.ylabel("Nombre de paires d'extension")
    plt.title("Nombre de paires de k-extension pour lesquelles la similarité se trouve \n à plus ou moins de 10% de la similarité à k = 10")
    plt.show()
    plt.savefig("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/histogramme_10pourcent.png")
                
                    
''' Etude des differences de sim en fonction de la taille d'extension au sein d'un cluster '''                  
def resultats_difference_taille(petit_cluster, nom):
    figure = plt.figure(figsize=(10,10))

    plt.gcf().subplots_adjust(left = 0.1, bottom = 0.1,
                       right = 0.7, top = 0.9, wspace = 0, hspace = 0.5)

    for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles") :
        if "crible" in fic :
            elt1 = fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6]
            elt2 = fic.split("_")[7] + "_" + fic.split("_")[8] + "_" + fic.split("_")[9] + "_" + fic.split("_")[10][:len(fic.split("_")[10])-7]
            
            if elt1 in petit_cluster and elt2 in petit_cluster :
            
                with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/%s"%(fic), 'rb') as fichier :
                    mon_depickler = pickle.Unpickler(fichier)
                    crible = mon_depickler.load()
                    
                    crible_avec_motif = []
                    
                                    
                    with open("Extensions/Metrique_toutes_aretes/graphe_complet_pondere_sim_coeff_n1_a1_c1.pickle", 'rb') as fichier_graphe_complet :
                        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
                        graphe_comp = mon_depickler_complet.load()
                         
                        for u,v, sim in graphe_comp.edges(data="poids") :
                            if (graphe_comp.nodes[u]["nom"] == elt1 and graphe_comp.nodes[v]["nom"] == elt2) or (graphe_comp.nodes[u]["nom"] == elt2 and graphe_comp.nodes[v]["nom"] == elt1) :
                                for compteur in range(1,10) :
                                    with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/graphe_comp_couples_possibles_fichier_%s_fichier_%s.pickle"%(compteur, elt1, elt2), 'rb') as fichier_commun :
                                        mon_depickler = pickle.Unpickler(fichier_commun)
                                        graphe_commun = mon_depickler.load()
                                    
                                    
                                        with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/fichier_%s.pickle"%(compteur, elt1), 'rb') as fichier_graphe1 :
                                            mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
                                            graphe1 = mon_depickler_1.load()
                                            
                                            with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/fichier_%s.pickle"%(compteur, elt2), 'rb') as fichier_graphe2 :
                                                mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
                                                graphe2 = mon_depickler_2.load()
                                                
                                                for cle in graphe_commun.keys() :
                                                    crible_avec_motif.append(calcul_sim_aretes_avec_coeff_avec_motif(graphe1, graphe2, graphe_commun[cle], fic, 1, 1, 1))
                               
                                if "graphe_comp_test_couples_possibles_fichier_%s_fichier_%s.pickle"%(elt1, elt2) in os.listdir("Extensions/Metrique_toutes_aretes/graphes_comparaison_pickle") :
                                    with open("Extensions/Metrique_toutes_aretes/graphes_comparaison_pickle/graphe_comp_test_couples_possibles_fichier_%s_fichier_%s.pickle"%(elt1, elt2), 'rb') as fichier_commun :
                                        mon_depickler = pickle.Unpickler(fichier_commun)
                                        graphe_commun_10 = mon_depickler.load()
                                        
                                        with open("Extensions/Metrique_toutes_aretes/graphes_extension_avec_liaisons_b53_qui_manquent/fichier_%s.pickle"%(elt1), 'rb') as fichier_graphe1 :
                                            mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
                                            graphe1 = mon_depickler_1.load()
                                            
                                            with open("Extensions/Metrique_toutes_aretes/graphes_extension_avec_liaisons_b53_qui_manquent/fichier_%s.pickle"%(elt2), 'rb') as fichier_graphe2 :
                                                mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
                                                graphe2 = mon_depickler_2.load()
                                                
                                                for cle in graphe_commun_10.keys() :
                                                    crible_avec_motif.append(calcul_sim_aretes_avec_coeff_avec_motif(graphe1, graphe2, graphe_commun_10[cle], fic, 1, 1, 1))
    
                                else :
                                    with open("Extensions/Metrique_toutes_aretes/graphes_comparaison_pickle/graphe_comp_test_couples_possibles_fichier_%s_fichier_%s.pickle"%(elt2, elt1), 'rb') as fichier_commun :
                                        mon_depickler = pickle.Unpickler(fichier_commun)
                                        graphe_commun_10 = mon_depickler.load()
                                
                                
                                        with open("Extensions/Metrique_toutes_aretes/graphes_extension_avec_liaisons_b53_qui_manquent/fichier_%s.pickle"%(elt2), 'rb') as fichier_graphe1 :
                                            mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
                                            graphe1 = mon_depickler_1.load()
                                            
                                            with open("Extensions/Metrique_toutes_aretes/graphes_extension_avec_liaisons_b53_qui_manquent/fichier_%s.pickle"%(elt1), 'rb') as fichier_graphe2 :
                                                mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
                                                graphe2 = mon_depickler_2.load()
                                                
                                                for cle in graphe_commun_10.keys() :
                                                    crible_avec_motif.append(calcul_sim_aretes_avec_coeff_avec_motif(graphe1, graphe2, graphe_commun_10[cle], fic, 1, 1, 1))
    
                                
                                for cle in crible :
                                    crible[cle].append(sim)
                                    
                        for cle in crible :
                            y1 = crible[cle]
                            compteur = 9
                            while compteur >= 0 and  crible[cle][compteur] > crible[cle][9] - crible[cle][9]/10 and crible[cle][compteur] < crible[cle][9] + crible[cle][9]/10 :
                                compteur -= 1
                            print(compteur)
                        
                        compteur_avec_motif = 9
                        while compteur >= 0 and crible_avec_motif[compteur_avec_motif] > crible_avec_motif[9] - crible_avec_motif[9]/10 and crible_avec_motif[compteur_avec_motif] < crible_avec_motif[9] + crible_avec_motif[9]/10 :
                            compteur_avec_motif -= 1
                        print(compteur_avec_motif)       
    
                        print(crible)
                        print(len(crible_avec_motif))
                        print(crible_avec_motif)
                        
                        x1 = range(0,10)
                        #ax1 = figure.add_subplot(211) 
                        #ax2 = figure.add_subplot(212) 
                        
                        
                        ax1 = plt.subplot(211)
                        plt.ylim(0.0, 1.1)
                        plt.xlim(0, 10)
                        ax1.xaxis.set_ticks(np.arange(0,10,1))
                        ax1.yaxis.set_ticks(np.arange(0.0,1.1,0.1))
                        ax1.set_title("Valeur de similarité en fonction de la taille de l'extension (sans compter le motif)")
                        ax1.set_xlabel("Taille de l'extension")
                        ax1.set_ylabel("Valeur de similarité")
                        plt.plot(x1,y1, label=cle)
                        plt.legend(loc='lower right', bbox_to_anchor= (1.5, -0.25), ncol=1, 
                borderaxespad=0.1, frameon=False)
                        
                        x2 = range(0,10)
                        y2 = crible_avec_motif
                        ax2 = plt.subplot(212)
                        plt.ylim(0.0, 1.1)
                        plt.xlim(0, 10)
                        #ax2.set_xlabel(np.arange(0,10,1))
                        ax2.xaxis.set_ticks(np.arange(0,10,1))
                        #ax2.set_ylabel(np.arange(0.0,1.0,0.1))
                        ax2.yaxis.set_ticks(np.arange(0.0,1.1,0.1))
                        
                        ax2.set_title("Valeur de similarité en fonction de la taille de l'extension (en comptant le motif)")
                        ax2.set_xlabel("Taille de l'extension")
                        ax2.set_ylabel("Valeur de similarité")
                        plt.plot(x2,y2, label=cle)
                        plt.legend(loc='lower right', bbox_to_anchor= (1.5, -0.25), ncol=1, 
                borderaxespad=0.1, frameon=False)
                        
    plt.savefig("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/%s_similarites_par_taille_extension.png"%(nom))
    plt.show()    
    

''' Comparaison des rangs des paires d'extension entre k=5 et k=10 '''
def comparaisons_rang_taille_5_10(petit_cluster):
    liste_rangs = []
    for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles") :
        if "crible" in fic :
            elt1 = fic.split("_")[3] + "_" + fic.split("_")[4] + "_" + fic.split("_")[5] + "_" + fic.split("_")[6]
            elt2 = fic.split("_")[7] + "_" + fic.split("_")[8] + "_" + fic.split("_")[9] + "_" + fic.split("_")[10][:len(fic.split("_")[10])-7]
            
            if elt1 in petit_cluster and elt2 in petit_cluster :
            
                with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/sim_extensions_toutes_aretestaille_5.pickle", 'rb') as fichier_sim : 
                        mon_depickler_sim = pickle.Unpickler(fichier_sim)
                        dico_sim_5 = mon_depickler_sim.load()  
                        
                        with open("Extensions/Metrique_toutes_aretes/sim_extensions_toutes_aretes_coeffn1_a1_c1.pickle", 'rb') as fichier_sim_10 : 
                            mon_depickler_sim_10 = pickle.Unpickler(fichier_sim_10)
                            dico_sim_10 = mon_depickler_sim_10.load()  
                            
                            dico_sim_sorted_5 = OrderedDict(sorted(dico_sim_5.items(), key= lambda t: t[1], reverse=True))
                            dico_sim_sorted_10 = OrderedDict(sorted(dico_sim_10.items(), key= lambda t: t[1], reverse=True))
                            
                            etape = 1
                            rang_5 = 0
                            for elt_sim in dico_sim_sorted_5 :
                                if elt_sim == (elt1, elt2) : 
                                    rang_5 = etape
                                etape += 1   
                            
                            etape = 1
                            rang_10 = 0
                            for elt_sim in dico_sim_sorted_10 :
                                if elt_sim == (elt1, elt2) or elt_sim == (elt2, elt1) : 
                                    rang_10 = etape
                                etape += 1      
                            liste_rangs.append((rang_5, rang_10))               
                            
    print(liste_rangs)                    
                    
if __name__ == '__main__':
    #comparaisons_rang_taille_5_10(['1U9S_A_58_11', '3UCZ_R_62_15', '4V9F_0_44_3', '1FJG_A_48_11', '1FJG_A_58_23', '4V9F_0_48_26', '5J7L_DA_50_21'])
    #resultats_entier()
    #resultats_difference_taille(['1U9S_A_58_11', '3UCZ_R_62_15', '4V9F_0_44_3', '1FJG_A_48_11', '1FJG_A_58_23', '4V9F_0_48_26', '5J7L_DA_50_21'], "motif_gnra")
#     for i in range(1,10) :
#         obtenir_extension(i)


    with open("Extensions/Metrique_toutes_aretes/graphe_complet_pondere_sim_coeff_n1_a1_c1.pickle", 'rb') as fichier_graphe_complet :
            mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
            graphe_comp = mon_depickler_complet.load()
            print(graphe_comp.number_of_nodes())
            p = multiprocessing.Pool(1)
              
            #petit_cluster = ['4V9F_0_30_23', '4V9F_0_48_16', '1FJG_A_271_1', '1FJG_A_109_6']
#             liste_a_faire = []
#             for u,v in graphe_comp.edges() :
#                 if graphe_comp.nodes[u]["nom"] in ['1FJG_A_48_8', '5J7L_DA_48_15'] and graphe_comp.nodes[v]["nom"] in ['5J7L_DA_48_15','1FJG_A_48_8'] :
#                     liste_a_faire.append((u,v))
#             print(list(graphe_comp.edges()) )
#             print(liste_a_faire)
            result = p.map(calcul_par_taille, list(graphe_comp.edges()) )
            p.close()
            p.join()
    #resultats()
    
    ##valeurs de sim avec differentes tailles d'extension (resultats pour 0.6 2)
      
                         
    
                                    
                                    
                                    