'''
Created on 30 janv. 2019

@author: coline
'''

import pickle
import itertools        

''' calcule l'ensemble des combinaisons de couples possibles entre deux listes d'entiers de 2 a 10
a partir de la liste des combinaisons de couples possibles entre deux listes d'entiers de 2 a 9 '''
def recherche_paires_9():
    with open("liste_combi_8.pickle", 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        liste_epuree = mon_depickler.load()
        
        print(liste_epuree)
        
        paires = [[]]
        try :
            for element in itertools.product([2,3,4,5,6,7,8,9,10], [10]):
                paires[0].append(element)
            del(paires[0][len(paires[0])-1])
            paires.append([])
            for element in itertools.product([10], [2,3,4,5,6,7,8,9,10]):
                if element not in paires :
                    paires[1].append(element)
            del(paires[1][len(paires[1])-1])
                
        except MemoryError as error :
            print("memory error 1")
        
        print(paires)
        
        liste_epuree_2 = []
        liste_ajoute = []
        
        for groupe in liste_epuree :
            pas_bon = False
            for elt in groupe :
                if 10 == elt[0] or 10 == elt[1] :
                    pas_bon = True 
            if pas_bon == False :
                liste_temp = []
                liste_temp = list(groupe)
                liste_temp.append((10,10))
                liste_epuree_2.append(list(liste_temp))
                if (10,10) not in liste_ajoute :
                    liste_ajoute.append((10,10))
            
            ok_1 = []
            for elt_1 in paires[0] :
                pas_bon = False
                for elt in groupe :
                    if elt_1[0] == elt[0] or elt_1[1] == elt[1] or (elt_1[0] < elt[0] and elt_1[1] > elt[1]) or (elt_1[0] > elt[0] and elt_1[1] < elt[1]):
                        pas_bon = True
                if pas_bon == False :
                    ok_1.append(elt_1)
            
            ok_2 = []     
            for elt_2 in paires[1] :
                pas_bon = False
                for elt in groupe :
                    if elt_2[0] == elt[0] or elt_2[1] == elt[1] or (elt_2[0] < elt[0] and elt_2[1] > elt[1]) or (elt_2[0] > elt[0] and elt_2[1] < elt[1]):
                        pas_bon = True
                if pas_bon == False :
                    ok_2.append(elt_2)
            
            print(ok_1)        
            for e in range(max(1, len(ok_1))) :
                for f in range(max(1,len(ok_2))) :
                    liste_temp = []
                    
                    liste_temp = list(groupe)
                    taille_liste = len(liste_temp)
                    
                    if len(ok_1) > 0 :
                        liste_temp.append(ok_1[e])
                        if ok_1[e] not in liste_ajoute :
                            liste_ajoute.append(ok_1[e])
                    if len(ok_2) > 0 :
                        liste_temp.append(ok_2[f])
                        if ok_2[f] not in liste_ajoute :
                            liste_ajoute.append(ok_2[f])
                    if taille_liste < len(liste_temp) :
                        liste_epuree_2.append(list(liste_temp))
        
        for elt in paires[0] :
            if elt not in liste_ajoute :
                liste_epuree_2.append([elt])
                
        for elt in paires[1] :
            if elt not in liste_ajoute :
                liste_epuree_2.append([elt])            

        
        print(len(liste_epuree_2))
#         print(liste_epuree_2)    
                
        with open("liste_combi_9.pickle", 'wb') as fichier :
            mon_pickler = pickle.Pickler(fichier)
            mon_pickler.dump(liste_epuree_2)
            
        with open("liste_combi_9.txt", 'w') as fichier_2 :
            for elt in liste_epuree_2 :
                fichier_2.write(str(elt) + '\n')

''' recherche l'ensemble des combinaisons de couples possibles entre deux listes d'entiers de 2 a 9  '''
def recherche_toutes_paires():

    chaines_1 = [2,3,4,5,6,7,8,9]
    chaines_2 = [2,3,4,5,6,7,8,9]
    paires = []
    try :
        for element in itertools.product(chaines_1, chaines_2):
            paires.append(element)
    except MemoryError as error :
        print("memory error 1")
        
    print(paires)
    
    liste = []
    new_couples = []
    new_chaine = []
    liste_epuree = [[(2,2), (3,3), (4,4), (5,5), (6,6), (7,7), (8,8), (9,9)]]
    k = 0
    #for i in range(1, min(len(chaines_1[num_chaine])-1, len(chaines_2[num_chaine])-1 ) +1) :
    while 7-k >= 1 :
#             print(min(len(chaines_1[num_chaine])-1, len(chaines_2[num_chaine])-1 )-k)
#             print(num_chaine)
        num_l = len(liste)
        try :
            combinaisons = list(itertools.combinations(paires,7-k))
        except MemoryError as error  :
            print("memory_error_2")
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
                        if (possib[i][0] < possib[j][0] and possib[i][1] > possib[j][1]) or (possib[i][0] > possib[j][0] and possib[i][1] < possib[j][1]) :
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

        
        ## ne garder que les combinaisons uniques
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
        
        print(liste_epuree)
        print(len(liste_epuree))
        
    with open("liste_combi_10.pickle", 'wb') as fichier :
        mon_pickler = pickle.Pickler(fichier)
        mon_pickler.dump(liste_epuree)
        
    #print(combinaisons)
#     paires = [[(1,1), (2,2), (3,3), (4,4), (5,5), (6,6), (7,7), (8,8), (9,9), (10, 10)]]
#     i = 9
#     compteur = 0
#     anciennes_paires = []
#     taille_anciennes_paires = 0
#     while i > 6 :
#         print(i)  
# 
#         l = 0
#         m = 0
#         ancienne_taille_ancienne_paire = taille_anciennes_paires
#         taille_anciennes_paires = len(anciennes_paires)
#             #print(taille_anciennes_paires)
#         while l < 10 :    
#  
#             for n in range(ancienne_taille_ancienne_paire, max(1,taille_anciennes_paires)) :
#                 paire_temp = []
#                 for k in range(1,i) :
#                     paire_temp.append((k,k))
#                     
#                 if i+m+1 < 11 :
#                     l = i+m+1
#                     paire_temp.append((i,l))
#                 
#                 if len(anciennes_paires) > 0 :
#                     prob = False
#                     for elt in anciennes_paires[n] :
#                         for elt_1 in paire_temp :
#                             if elt[0] == elt_1[0] or elt[1] == elt_1[1] :
#                                 prob = True
#                     if prob == False :
#                         for elt in anciennes_paires[n] :
#                             paire_temp.append(elt)
#                     #del(anciennes_paires[n])
#                 print(paire_temp)    
#                 existe_deja = False
#                 for j in range(len(anciennes_paires)) :
#                     groupe = anciennes_paires[j]
#                     
#                     nb_pareil = 0
#                     for elt in groupe :
#                         for elt_1 in paire_temp[i-1:] :
#                             if elt[0] == elt_1[0] and elt[1] == elt_1[1] :
#                                 nb_pareil += 1
#                     
#                     if nb_pareil == len(groupe) and nb_pareil == len(paire_temp[i-1:]) :
#                         existe_deja = True
#                
#                 
#                 if existe_deja == False :
#                     anciennes_paires.append(paire_temp[i-1:])
#                 #print(paire_temp)
#                 paires.append(paire_temp)    
#             m = m+1
#    
#             if taille_anciennes_paires == 0 :
#                 l = 11
#         print(anciennes_paires)
# #         for n in range(taille_anciennes_paires) :
# #             del(anciennes_paires[n])
# #         print(anciennes_paires) 
#                 
#         
#         i = i-1
#     print(anciennes_paires)
#     print(paires)        

''' pas tres utile je pense '''
def recherche_couples(couples_possibles, chaines_1, chaines_2):
   
    for num_chaine in range(0, 4) :
        print(chaines_1[num_chaine])
        print(chaines_2[num_chaine])
        
        
            
        
        
if __name__ == '__main__':
#     with open("grands_graphes.pickle", 'rb') as fichier :
#         mon_depickler = pickle.Unpickler(fichier)
#         dico_graphes = mon_depickler.load()
#         
#         with open("fichiers_pickle/a-minor_test2.pickle", 'rb') as fichier_pickle :
#             mon_depickler = pickle.Unpickler(fichier_pickle)
#             tab_aminor = mon_depickler.load()
#             
#             for occ in tab_aminor :
#                 if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 48 and occ["num_occ"] == 9 :
#             
#                     graphe1 = dico_graphes[('5DM6', 'X', 48, 9)]
#                     
#                     chaines_1 = [[occ["a_minor"][0]]]
#                     i = 1
#                     for elt in occ["a_minor"] :
#                         compteur = elt
#                         if i != 1 : chaines_1.append([elt])
#                         liaison_B53 = True
#                         while liaison_B53 :
#                             liaison_B53 = False
#                             temp = compteur
#                             for voisin in graphe1.successors(compteur) :
#                                 for arc in graphe1[compteur][voisin] :
#                                     if voisin not in occ["a_minor"] and voisin not in chaines_1[len(chaines_1)-1] and graphe1[compteur][voisin][arc]["label"] == 'B53' :
#                                         liaison_B53 = True
#                                         temp = voisin
#                                         chaines_1[len(chaines_1)-1].append(voisin)
#                                         
#                             for voisin in graphe1.predecessors(compteur) :
#                                 for arc in graphe1[voisin][compteur] :
#                                     if voisin not in occ["a_minor"] and voisin not in chaines_1[len(chaines_1)-1] and graphe1[voisin][compteur][arc]["label"] == 'B53' :
#                                         liaison_B53 = True
#                                         temp = voisin
#                                         chaines_1[len(chaines_1)-1].append(voisin)
#                             compteur = temp
#                         i = i +1
#             
#                             #with open("fichiers_tot_couples_possibles.txt", 'a') as fichier_tot :
#                                 #fichier_tot.write(element1 + "\n")
#                                 #fichier_tot.write("Chaines : " + str(chaines_1) + "\n")
#                     graphe2 = dico_graphes[('5J7L', 'DA', 197, 4)]
#                     
#                     for occ_2 in tab_aminor : 
#                         if occ_2["num_PDB"] == '5J7L' and occ_2["num_ch"] == 'DA' and occ_2["num_motif"] == 197 and occ_2["num_occ"] == 4 :
#                             chaines_2 = [[occ_2["a_minor"][0]]]
#                             
#                             i = 1
#                             for elt in occ_2["a_minor"] :
#                                     compteur = elt
#                                     if i != 1 : chaines_2.append([elt])
#                                     liaison_B53 = True
#                                     while liaison_B53 :
#                                         liaison_B53 = False
#                                         temp = compteur
#                                         for voisin in graphe2.successors(compteur) :
#                                             for arc in graphe2[compteur][voisin] :
#                                                 if voisin not in occ_2["a_minor"] and voisin not in chaines_2[len(chaines_2)-1] and graphe2[compteur][voisin][arc]["label"] == 'B53' :
#                                                     liaison_B53 = True
#                                                     temp = voisin
#                                                     chaines_2[len(chaines_2)-1].append(voisin)
#                                                     
#                                         for voisin in graphe2.predecessors(compteur) :
#                                             for arc in graphe2[voisin][compteur] :
#                                                 if voisin not in occ_2["a_minor"] and voisin not in chaines_2[len(chaines_2)-1] and graphe2[voisin][compteur][arc]["label"] == 'B53' :
#                                                     liaison_B53 = True
#                                                     temp = voisin
#                                                     chaines_2[len(chaines_2)-1].append(voisin)
#                                         compteur = temp
#                                     i = i+1
#                     print(chaines_1)
#                     print(chaines_2)
                    
                    
    #recherche_toutes_paires()
    recherche_paires_9()
                    
                    
                    