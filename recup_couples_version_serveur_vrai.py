'''
Created on 4 juin 2019

@author: coline

Implementation de la methode 1 de comparaison de graphes (recherche de tous les couples possibles de noeuds + recherche du sous-graphe commun max a partir de ces ensembles de couples possibles)
version CaRNAval 2.0 (version lancee sur le serveur)
'''
from recup_data.constantes import EXTENSION_PATH_TAILLE, EXTENSION_PATH,\
    CLUSTERING_PEREZ_VERSION_NON_CAN_2
import pickle
import os
import networkx as nx
import multiprocessing
from itertools import product
from itertools import combinations

''' renvoie l'ensemble des couples possibles de superposition de noeuds entre graphe1 et graphe2 (idem que dans recup_couples_V2)'''
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
        
''' fait la recherche de couples possibles pour une paire de graphes 
pour cela, recherche les noeuds appartenant a chacune des chaines
appelle la fonction recherche_couples et stocke le resultat dans un fichier pickle (idem que dans recup_couples_V2) '''
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
    if "pickle" in element1 and "pickle" in element2 and element2 != element1 and  "couples_possibles_"+ element1[:len(element1)-7] + "_" + element2 not in os.listdir(EXTENSION_PATH_TAILLE%taille_ext) and "couples_possibles_"+ element2[:len(element2)-7] + "_" + element1 not in os.listdir(EXTENSION_PATH_TAILLE%taille_ext) : 
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
                        

                        with open(EXTENSION_PATH_TAILLE%taille_ext+element1, 'rb') as fichier1 :
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
                                        with open(EXTENSION_PATH_TAILLE%taille_ext+element2, 'rb') as fichier2 :
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
                                                
                                            with open(EXTENSION_PATH_TAILLE%taille_ext+"couples_possibles_%s_%s.pickle"%(element1[:len(element1)-7],element2[:len(element2)-7]), "wb") as fichier_sauvegarde :
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

''' meme fonction que dans calcul_sim '''
def calcul_aretes_avec_coeff(graphe, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
    
    for u,v,data in graphe.edges(data=True) :
        if data["label"] != 'B53' :
            if data["label"] == '0' :
                if coeffa == 1 :
                    somme_aretes += graphe.nodes[u]["poids"]
            elif graphe.nodes[u]["type"] == 1 and graphe.nodes[v]["type"] == 1 or (data["type"] == 'CWW' and data["long_range"] == False) :
                if coeffc == 1 :
                    somme_aretes += graphe.nodes[u]["poids"]
            else :
                if coeffn == 1 :
                    somme_aretes += graphe.nodes[u]["poids"]
    somme_aretes = somme_aretes/2 - 4
    
    return somme_aretes

''' meme fonction que dans calcul_sim '''
def calcul_aretes_communes_avec_coeff(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn):
    somme_aretes = 0
    for u,v,data in graphe_commun.edges(data=True) :
        if data["type"] != 'B53' and (u,v) not in liste_aretes :
            if data["type"] == '0' :
                if coeffa == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
            elif graphe1.nodes[u[0]]["type"] == 1 and graphe1.nodes[v[0]]["type"] == 1 or (data["type"] == 'CWW' and data["long_range"] == False) :
                if coeffc == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
            else :
                if coeffn == 1 :
                    somme_aretes += min(graphe1.nodes[u[0]]["poids"], graphe2.nodes[u[1]]["poids"]) 
            liste_aretes.append((u,v))       
    somme_aretes = somme_aretes/2 - 4
    
    return somme_aretes
       
''' meme fonction que dans calcul_sim '''
def calcul_sim_aretes_avec_coeff(graphe1, graphe2, graphe_commun, cle, coeffc, coeffa, coeffn):
    aretes_1 = calcul_aretes_avec_coeff(graphe1, cle, coeffc, coeffa, coeffn)
    aretes_2 = calcul_aretes_avec_coeff(graphe2, cle, coeffc, coeffa, coeffn)
    aretes_commun = calcul_aretes_communes_avec_coeff(graphe_commun, graphe1, graphe2, cle, coeffc, coeffa, coeffn)
    
#     print(aretes_1)
#     print(aretes_2)
#     print(aretes_commun)
    
    return aretes_commun/max(aretes_1, aretes_2)

''' renvoie vrai s'il n'y a pas d'incompatibilite de sequence dans graphe_commun par rapport a graphe1 et graphe2 (idem que dans sous_graphe_commun_max)'''
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

        
''' renvoie vrai si l'un des noeuds du couple_a_chercher est deja present dans une paire dans le graphe et pas l'autre
    faux sinon (idem que dans sous_graphe_commun_max)'''
def dans_graphe(graphe, couple_a_chercher):
    for noeud in graphe.nodes() :
        #print(noeud)
        if (noeud[0] == couple_a_chercher[0] and noeud[1] != couple_a_chercher[1]) or (noeud[1] == couple_a_chercher[1] and noeud[0] != couple_a_chercher[0])  :
            return True
    return False

''' ajout des noeuds stockes dans elt au graphe commun si possible (et leurs voisins non covalents s'il en a et que c'est possible aussi) (tres proche de la version de sous_graphe_commun_max)
elt : ensemble des couples de noeuds de graphe1 et graphe2 a ajouter potentiellement au graphe_commun (appartiennent tous a la meme chaine)
elt_1, elt_2, elt_3 :  ensemble des couples de noeuds des autres chaines pour verifier la compatibilite
couples_possibles : ensemble de tous les couples de noeuds possibles pour chaque chaine
num : numero de la chaine
compt : indice de elt dans couples_possibles '''
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
                                #if dict_voisin_1[edge_1]["label"] == 'CWW' :
                                #    label_1 = 'CAN'
                                #elif dict_voisin_1[edge_1]["label"] == 'B53' : 
                                   #    label_1 = 'COV'
                                #elif dict_voisin_1[edge_1]["label"] == '0' :
                                #    label_1 = 'ART'
                                #else :
                                #    label_1 = 'NON_CAN'  
                
                                
                                #if dict_voisin_2[edge_2]["label"] == 'CWW' :
                                #    label_2 = 'CAN'
                                #elif dict_voisin_2[edge_2]["label"] == 'B53' : 
                                #    label_2 = 'COV'
                                #elif dict_voisin_2[edge_2]["label"] == '0' :
                                #    label_2 = 'ART'
                                #else :
                                #    label_2 = 'NON_CAN'
                                
                                label_1 = dict_voisin_1[edge_1]["label"]
                                label_2 = dict_voisin_2[edge_2]["label"]


                                bonne_chaine = False
                                for ch in graphe1.nodes[voisin_1]["chaine"] :
                                    if ch in graphe2.nodes[voisin_2]["chaine"] :
                                        bonne_chaine = True
                                
                                
      
                                
                                if bonne_chaine :
                                        
                                    if label_1 == label_2 and dict_voisin_1[edge_1]["long_range"] == dict_voisin_2[edge_2]["long_range"] :
                                        if label_1 == 'B53' :
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
                                                                if succ_1 and succ_2 :
                                                                    graphe_commun.add_edge(elt[i], (num,num), type='B53', long_range=False)
                                                                else :
                                                                    graphe_commun.add_edge((num, num), elt[i], type='B53', long_range=False)
                                                        #nb_aretes += 1
                                                        #print("ramou")
                                                    else : #elt courant non lie au motif
                                                        for couple in elt : # on cherche si un couple deja dans la liste ne correspond pas aux voisins
                                                            if dans_graphe(graphe_commun, elt[i]) == False and voisin_1 == couple[0] and voisin_2 == couple[1]  : #and elt[i][0] not in deja_vus_voisin_1 and elt[i][1] not in deja_vus_voisin_2 :
                                                                #deja_vus_voisin_1.append(voisin_1)
                                                                #deja_vus_voisin_2.append(voisin_2)
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
                                                                    if succ_1 and succ_2 :
                                                                        graphe_commun.add_edge(elt[i], (voisin_1, voisin_2), type='B53', long_range=False)
                                                                    else :
                                                                        graphe_commun.add_edge((voisin_1, voisin_2), elt[i], type='B53', long_range=False)
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
                                                                ## on va ajouter une liaison de chaque cote si elles existent bien dans les graphes de depart  
                                                                if succ_1 == True and succ_2 == True : 
                                                                    graphe_commun.add_edge(elt[i], (voisin_1, voisin_2), type=label_1, long_range = dict_voisin_1[edge_1]["long_range"])
                                                                    ordre = ((voisin_1, voisin_2), elt[i])
                                                                elif succ_1 == False and succ_2 == False :
                                                                    graphe_commun.add_edge((voisin_1, voisin_2), elt[i], type=label_1, long_range = dict_voisin_1[edge_1]["long_range"])
                                                                    ordre = (elt[i], (voisin_1, voisin_2))
                                                                else : 
                                                                    print("bizarre 2")
                                                                if (ordre[0][0], ordre[1][0]) in graphe1.edges() and (ordre[0][1], ordre[1][1]) in graphe2.edges() :
                                                                    
                                                                    label = "manque"
                                                                    for edge in graphe1[ordre[0][0]][ordre[1][0]] :
                                                                        if (label_1 == '0' and graphe1[ordre[0][0]][ordre[1][0]][edge]["label"] == '0') or graphe1[ordre[0][0]][ordre[1][0]][edge]["label"][:3] == label_1[0] + label_1[2] + label_1[1] :
                                                                            label = graphe1[ordre[0][0]][ordre[1][0]][edge]["label"]
#                                                                         if elt[i] == (20,20) :
#                                                                             print("tout petit rat")
#                                                                             print(graphe1[ordre[0][0]][ordre[1][0]][edge]["label"][3:])
#                                                                             print()   
#                                                                     if elt[i] == (20,20) :
#                                                                         print("gros rat")
#                                                                         print(voisin_1)
#                                                                         print(voisin_2)
#                                                                         print(label)
#                                                                         print(label_1)
#                                                                         print(ordre)
                                                                    graphe_commun.add_edge(ordre[0], ordre[1], type=label, long_range = dict_voisin_1[edge_1]["long_range"])
                                                                else : 
                                                                    print("bizarre")
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
                                                    if (3, 3) in elt and (37, 27) in elt and (39, 29) in elt and (40, 30) in elt and (41, 31) in elt and \
                                                            (1, 1) in elt_1 and (2, 2) in elt_2 and (6, 6) in elt_2 and (9, 7) in elt_2 and (12, 8) in elt_2 and (13, 9) in elt_2 and (15, 11) in elt_2 and\
                                                            (4, 4) in elt_3 and (17, 13) in elt_3 and (18, 14) in elt_3 and (20, 16) in elt_3 and (23, 20) in elt_3 and (25, 22) in elt_3 and (27, 24) in elt_3 and \
                                                            elt[i] == (40, 30) and (voisin_1, voisin_2) == (38, 28) :
                                                                print("petit rat")
                                                                print(voisin_1)
                                                                print(voisin_2)
                                                                print(label_1)
                                                                print(label_2)
                                                                print(graphe_commun[elt[i]])
                                                                print(succ_1)
                                                                print(succ_2)
                                                    if pas_vu_1_motif and pas_vu_2_motif  : #and elt[i][0] not in deja_vus_voisin_1 and elt[i][1] not in deja_vus_voisin_2 and voisin_1 not in sommets_en_plus_1 and voisin_2 not in sommets_en_plus_2 : ## on ajoute un sommet si le sommet n'existe nulle part
                                                        
                                                        graphe_commun_temp = graphe_commun.copy()
                                                        if elt[i] not in graphe_commun.nodes() :
                                                            graphe_commun_temp.add_node(elt[i])
                                                            
                                                        if (voisin_1, voisin_2) not in graphe_commun.nodes() :
                                                            graphe_commun_temp.add_node((voisin_1, voisin_2))
                                                        
                                                        if test_compatibilite(graphe_commun_temp, graphe1, graphe2) :
                                                        
                                                            if pas_vu_1 and pas_vu_2 : #les sommets n'existent pas, on les rajoute
                                                                    graphe_commun.add_node((voisin_1, voisin_2))
                                                                    if (voisin_1, voisin_2) != (5,5) :
                                                                        tab_temp.append((voisin_1, voisin_2))
                                                                        
                                                            if elt[i] not in graphe_commun.nodes() :
                                                                graphe_commun.add_node(elt[i])
                                                            
                                                            
                                                            deja_ajoute = False
                                                            
                                                            
                                                            
                                                            if succ_1 and succ_2 :
                                                                
                                                                for voisin in graphe_commun[elt[i]] :
                                                                    if voisin == (voisin_1, voisin_2) :
                                                                        for edge in graphe_commun[elt[i]][voisin] :
                                                                            if graphe_commun[elt[i]][voisin][edge]["type"] == label_1 :
                                                                                deja_ajoute = True
                                                            elif not succ_1 and not succ_2 :
                                                                for voisin in graphe_commun[elt[i]] :
                                                                    if voisin == (voisin_1, voisin_2) :
                                                                        for edge in graphe_commun[voisin][elt[i]] :
                                                                            if graphe_commun[voisin][elt[i]][edge]["type"] == label_1 :
                                                                                deja_ajoute = True
                                                                
                                                            if deja_ajoute == False and (voisin_1, voisin_2) in graphe_commun.nodes() :
                                                                
#                                                                 if elt[i] == (7,7) :
#                                                                     print(graphe_commun.edges.data())
#                                                                     print("ratmou1")
                                                                if succ_1 == True and succ_2 == True : 
                                                                        #print("petit petit rat")
                                                                        graphe_commun.add_edge(elt[i], (voisin_1, voisin_2), type=label_1, long_range = dict_voisin_1[edge_1]["long_range"])
                                                                        ordre = ((voisin_1, voisin_2), elt[i])
                                                                elif succ_1 == False and succ_2 == False :
                                                                        graphe_commun.add_edge((voisin_1, voisin_2), elt[i], type=label_1, long_range = dict_voisin_1[edge_1]["long_range"])
                                                                        ordre = (elt[i], (voisin_1, voisin_2))
                                                                else : 
                                                                        print("bizarre 2")
    #                                                             print(elt[i])
    #                                                             print(ordre)
    #                                                             print(ordre[0][0], ordre[1][0])
                                                                if (ordre[0][0], ordre[1][0]) in graphe1.edges() and (ordre[0][1], ordre[1][1]) in graphe2.edges() :
                                                                        label = "manque"
                                                                        
                                                                        for edge in graphe1[ordre[0][0]][ordre[1][0]] :
                                                                            if (label_1 == '0' and graphe1[ordre[0][0]][ordre[1][0]][edge]["label"] == '0') or graphe1[ordre[0][0]][ordre[1][0]][edge]["label"][:3] == label_1[0] + label_1[2] + label_1[1] :
                                                                                label = graphe1[ordre[0][0]][ordre[1][0]][edge]["label"]
#                                                                         if elt[i] == (4,4) :
#                                                                             print(voisin_1)
#                                                                             print(voisin_2)
#                                                                             print(label)
#                                                                             print("petit petit rat")
                                                                        graphe_commun.add_edge(ordre[0], ordre[1], type=label, long_range = dict_voisin_1[edge_1]["long_range"])
#                                                                         print(label_1)
#                                                                         print(label)
#                                                                         print("tamousnif")
                                                                else : 
                                                                        print("bizarre") 
                                                                        
#                                                                 if elt[i] == (7,7) :
#                                                                     print(graphe_commun.edges.data())
#                                                                     print("ratmou")
#                                                                     print(voisin_1)
#                                                                     print(voisin_2)  
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

''' effectue la recherche de sous-graphe commun maximum pour une paire d'occurrences
en testant toutes les possibilites de sous_graphes qu'on peut obtenir a partir des ensembles de couples de noeuds stockes dans couples_possibles (obtenus dans recup_couples_V2)
(appel a comparaison pour chaque couple de noeuds) (proche de la version dans sous_graphe_commun_max.py)
stockage du graphe commun resultat dans un fichier pickle (1 noeud dans le graphe commun = (1 noeud du premier graphe, 1 noeud du deuxieme graphe))
fichier : noms des deux occurrences 
graphe1 : graphe d'extension de l'occurrence 1
graphe2 : graphe d'extension de l'occurrence 2
taille_ext : taille de l'extension'''
def sous_graphe_commun_max(fichier, graphe1, graphe2, taille_ext):
    c = 0
    

    dico_graphes = {}
#                         if "pickle" in fic :
    fic = fichier
    if "pickle" in fic :
        element1 = fic.split('_')[2] + '_' + fic.split('_')[3] + '_' + fic.split('_')[4] + '_' + fic.split('_')[5] + '_' + fic.split('_')[6]
        element2 = fic.split('_')[7] + '_' + fic.split('_')[8] + '_' + fic.split('_')[9] + '_' + fic.split('_')[10] + '_' + fic.split('_')[11][:len(fic.split('_')[11])-7]
       #element1 = "fichier_4V9F_0_134_5"
       #element2 = "fichier_5FDU_1A_272_1"
        print(element1)
        print(element2)
        print(fic)
        if "graphe_comp_test_"+fic not in os.listdir(EXTENSION_PATH_TAILLE%taille_ext) : 
        #if True : 
            c = c+1
            print(c)
            with open("graphes_extension/"+element1+".pickle", 'rb') as fichier1 :
                mon_depickler1 = pickle.Unpickler(fichier1)
                graphe1 = mon_depickler1.load()     
                with open("graphes_extension/"+element2+".pickle", 'rb') as fichier2 :
                    mon_depickler2 = pickle.Unpickler(fichier2)
                    graphe2 = mon_depickler2.load()
    
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
                    
                    graphe_motif = nx.MultiDiGraph()
                    for i in range(1,6) :
                        graphe_motif.add_node((i,i))
                    graphe_motif.add_edge((1,1),(2,2), type="CSS", long_range=True)
                    graphe_motif.add_edge((2,2),(1,1), type="CSS", long_range=True)
                    graphe_motif.add_edge((3,3),(1,1), type="B53", long_range=False)
                    graphe_motif.add_edge((1,1),(5,5), type="TSS", long_range=True)
                    graphe_motif.add_edge((5,5),(1,1), type="TSS", long_range=True)
                    graphe_motif.add_edge((2,2),(4,4), type="B53", long_range=False)
                    graphe_motif.add_edge((2,2),(5,5), type="CWW", long_range=False)
                    graphe_motif.add_edge((5,5),(2,2), type="CWW", long_range=False)
                    graphe_motif.add_edge((3,3),(4,4), type="CSS", long_range=True)
                    graphe_motif.add_edge((4,4),(3,3), type="CSS", long_range=True)
                    
                    print("ramousnif")
                    with open(EXTENSION_PATH_TAILLE%taille_ext+fic, 'rb') as fichier_pickle :        
                        print("ramou")                
                        mon_depickler = pickle.Unpickler(fichier_pickle)
                        couples_possibles = mon_depickler.load()   
                   #couples_possibles = [[[(31,35),(32,36),(40,44),(42,46),(43,40)]],[[(9,8),(12,12),(18,16),(20,18)]],[[(44,48)]],[[(27,20),(28,21),(29,22)]]]
                        print(fic)
                        print(couples_possibles)
                        
                        #couples_possibles[0] = [[(25,23),(29,26),(31,29)]]
                        if couples_possibles != None :
                            
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
                                   
                                    graphe_commun_max = nx.MultiDiGraph()
                                    graphe_commun_max = graphe_motif.copy()
                                   #result = exectimeout(50, toutes_comparaisons, args=(couples_possibles, graphe1, graphe2, graphe_motif, graphe_commun_max))
                                    maxi = 0
                                    couples_max = []
                                   
                                   
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
                                                    
                                                    
                                                    if (3, 3) in elt_3 and (37, 27) in elt_3 and (39, 29) in elt_3 and (40, 30) in elt_3 and (41, 31) in elt_3 and \
                                                        (1, 1) in elt_1 and (2, 2) in elt_2 and (6, 6) in elt_2 and (9, 7) in elt_2 and (12, 8) in elt_2 and (13, 9) in elt_2 and (15, 11) in elt_2 and\
                                                        (4, 4) in elt_4 and (17, 13) in elt_4 and (18, 14) in elt_4 and (20, 16) in elt_4 and (23, 20) in elt_4 and (25, 22) in elt_4 and (27, 24) in elt_4 :
                                                        print("ramousn")
                                                    
                                                    
                                                    if test_compatibilite(graphe_commun, graphe1, graphe2) :
                                                        sim = calcul_sim_aretes_avec_coeff(graphe1, graphe2, graphe_commun, fic, 1, 1, 1)
                                                        
                                                        if (3, 3) in elt_3 and (37, 27) in elt_3 and (39, 29) in elt_3 and (40, 30) in elt_3 and (41, 31) in elt_3 and \
                                                        (1, 1) in elt_1 and (2, 2) in elt_2 and (6, 6) in elt_2 and (9, 7) in elt_2 and (12, 8) in elt_2 and (13, 9) in elt_2 and (15, 11) in elt_2 and\
                                                        (4, 4) in elt_4 and (17, 13) in elt_4 and (18, 14) in elt_4 and (20, 16) in elt_4 and (23, 20) in elt_4 and (25, 22) in elt_4 and (27, 24) in elt_4 :
                                                            print(sim)
                                                            print(graphe_commun.nodes.data())
                                                            print(graphe_commun.edges.data())
 
                                                        
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
                                    
                                    with open(EXTENSION_PATH_TAILLE%taille_ext+"graphe_comp_test_%s"%(fic), 'wb') as fichier_graphes :
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
                                    with open("fichier_similarite_nouvelle_metrique_3_%s.txt"%(str(taille_ext)), 'a') as fichier_ecriture :    
                                        fichier_ecriture.write(element1 + '\n' + element2 + '\n')
                                        fichier_ecriture.write("Graphe : ")
                                        fichier_ecriture.write(str(graphe_commun_max.nodes.data())+ '\n')
                                        fichier_ecriture.write(str(graphe_commun_max.edges.data())+ '\n')
                                        fichier_ecriture.write("Similarite :" + str(maxi))
                                        fichier_ecriture.write('\n\n')  
                                    
                                    print(couples_possibles[2])
                                    
                                    return graphe_commun_max, maxi
                                   #with open("tab_calcules_sim_"+nom_extension+".pickle", 'wb') as fichier_deja_fait :
#                                     with open("tab_calcules_sim_nouvelle_metrique.pickle", 'wb') as fichier_deja_fait :
#                                         mon_pickler_2 = pickle.Pickler(fichier_deja_fait)
#                                         mon_pickler_2.dump(tab_deja_fait)
#                                         with open("dico_graphes_communs_max_nouvelle_metrique.pickle", 'wb') as fichier_graphes :
#                                             mon_pickler_3 = pickle.Pickler(fichier_graphes)
#                                             mon_pickler_3.dump(dico_graphes)
        return False

''' effectue la recherche de couples possibles de noeuds et de sous-graphes communs max pour une paire d'occurrences
pour plusieurs tailles d'extension 
peut potentiellement stocker tous les resultats dans un fichier crible'''
def calcul_par_taille(couple):

    u = couple[0]
    v = couple[1]
#     for u,v, data in graphe_comp.edges(data=True) :
    with open(EXTENSION_PATH%1 + "graphe_complet_pondere_sim_toutes_aretes_coeff_all1_taille_%d.pickle"%1, 'rb') as fichier_graphe_complet :
        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
        graphe_comp = mon_depickler_complet.load()
        
        dico_elts = {}
        #if graphe_comp.edges[u,v]["poids"] > 0.5 :
        for i in range(4,5) :
                #for u,v in graphe_comp.edges() :
                    #if graphe_comp.nodes[u]["nom"] in liste_a_faire[i-1] or graphe_comp.nodes[v]["nom"] in liste_a_faire[i-1] :
                        print("i : %d"%(i))
                        print("u : %d v : %d"%(u,v))
                        with open(EXTENSION_PATH_TAILLE%i+"fichier_%s.pickle"%(graphe_comp.nodes[u]["nom"]), 'rb') as fichier1 :
                            mon_depickler_graphe1 = pickle.Unpickler(fichier1)
                            graphe1 = mon_depickler_graphe1.load()
                            print(graphe1)
                            with open(EXTENSION_PATH_TAILLE%i+"fichier_%s.pickle"%(graphe_comp.nodes[v]["nom"]), 'rb') as fichier2 :
                                mon_depickler_graphe2 = pickle.Unpickler(fichier2)
                                graphe2 = mon_depickler_graphe2.load()
                                print(graphe2)
                                element = ("fichier_"+graphe_comp.nodes[u]["nom"]+".pickle", "fichier_"+graphe_comp.nodes[v]["nom"]+".pickle")
                                #recup_couples(element, i)
                                
                                fichier = "couples_possibles_%s_%s.pickle"%(element[0][:len(element[0])-7], element[1][:len(element[1])-7])
            #                         for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s/couples_possibles_%s_%s"%(i, element[0], element[1])) :
            #                             if graphe_comp.nodes[v]["nom"] in fic and graphe_comp.nodes[u]["nom"] in fic :
            #                                 fichier = fic
            #                         if fichier == None :
            #                             print("probleme")  
            #                             break
                                print(type(fichier))
                                print(type(graphe1))
                                print(type(graphe2))
                                print(type(i))
                                res = sous_graphe_commun_max(fichier, graphe1, graphe2, i)
                                if res != False :
                                    graphe_commun = res[0]
                                    sim = res[1]
                                 
                                #if sim > data["poids"] - data["poids"]/10 and sim < data["poids"] + data["poids"]/10 :
                                    if (graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]) not in dico_elts.keys() :
                                        dico_elts.update({(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]) : [sim]})
                                    else :
                                        dico_elts[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])].append(sim)     
            
                                   
     
#     with open("graphes_extension_autres_tailles_new/crible_taille_extensions_%s_%s_3.pickle"%(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]), 'wb') as fichier_pickle :
#             mon_pickler = pickle.Pickler(fichier_pickle)
#             mon_pickler.dump(dico_elts)
if __name__ == '__main__':
    
#     for i in range(1,10) :
#         obtenir_extension(i)

        liste_aretes = []
        with open(EXTENSION_PATH%1+"graphe_complet_pondere_sim_toutes_aretes_coeff_all1_taille_%d.pickle"%1, 'rb') as fichier_graphe_complet :
            mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
            graphe_comp = mon_depickler_complet.load()
            for u,v in graphe_comp.edges() :
                if graphe_comp.nodes[u]["nom"] == "4PRF_B_25_69" or graphe_comp.nodes[v]["nom"] == "4PRF_B_25_69" :
                    liste_aretes.append((u,v))
        
        for elt in liste_aretes :
            calcul_par_taille(elt)
            
        print(len(liste_aretes))

#             liste_aretes = []
#             for u,v in graphe_comp.edges() :
#                 deja_fait = False
#                 for fic in os.listdir(EXTENSION_PATH_TAILLE%1) :
#                     if "pickle" in fic and "graphe_comp" in fic :
#                         cle1 = fic.split("_")[5] + "_" + fic.split("_")[6] + "_" + fic.split("_")[7] + "_" + fic.split("_")[8]
#                         cle2 = fic.split("_")[11] + "_" + fic.split("_")[12] + "_" + fic.split("_")[13] + "_" + fic.split("_")[14]
#                         if graphe_comp.nodes[u]["nom"] in [cle1, cle2] and graphe_comp.nodes[v]["nom"] in [cle1, cle2] :
#                             deja_fait = True
#                 if not deja_fait :
#                     liste_aretes.append((u,v))
#             print(len(liste_aretes))
 
            #for i in range(3,10) :
            #p = multiprocessing.Pool(1)
            #liste_a_faire_edges = []
#             petit_cluster = ['4V9F_0_30_23', '4V9F_0_48_16', '1FJG_A_271_1', '1FJG_A_109_6']
#             liste_a_faire = []
#            for u,v in graphe_comp.edges() :
#                if graphe_comp.nodes[u]["nom"] in liste_a_faire[i-1] or graphe_comp.nodes[v]["nom"] in liste_a_faire[i-1] :
#                    liste_a_faire_edges.append((u,v))
            #result = p.map(calcul_par_taille, liste_aretes)
            #print(liste_a_faire_edges)
                #result = p.apply_async(calcul_par_taille(liste_a_faire_edges, i), liste_a_faire_edges, i)
            #p.close()
            #p.join()
            
            
#             print((u,v))
#             calcul_par_taille((u,v))
