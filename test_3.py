'''
Created on 4 déc. 2018

@author: coline
'''

import pickle
import os
import time

# for fic in os.listdir("graphes_extension/fichiers_couples_qui_manquent_incomp/") :
#         tps1 = time.time()
#         if fic not in os.listdir("graphes_extension/fichiers_couples_qui_manquent_epures/")  :
#             with open("graphes_extension/fichiers_couples_qui_manquent_incomp/"+fic, 'rb') as fichier_sortie :
#                 mon_depickler = pickle.Unpickler(fichier_sortie)
#                 new_couples_possibles = mon_depickler.load()
#             
#                 new_couples_possibles_epuree = [[],[],[],[]]
#                 
#                 #print(len(new_couples_possibles))
#                 if len(new_couples_possibles) == 4 :
#                     for m in range(4) :
#                         #print(m)
#                         print(new_couples_possibles)
#                         if 'memory_error' not in new_couples_possibles[m] :
#                             a_enlever = []
#                             for i in range(0,len(new_couples_possibles[m])) :
#                                 for j in range(0, len(new_couples_possibles[m])) :
#                                         if j != i :
#                                             nb_existe_deja = 0
#                                             for k in range(0, len(new_couples_possibles[m][i])) :
#                                                 for l in range(0, len(new_couples_possibles[m][j])) :
#                                                     if new_couples_possibles[m][i][k][0] == new_couples_possibles[m][j][l][0] and new_couples_possibles[m][i][k][1] == new_couples_possibles[m][j][l][1]  :
#                                                         nb_existe_deja += 1
#                                             #print(nb_existe_deja)
#                                             if nb_existe_deja == len(new_couples_possibles[m][i]) and i not in a_enlever :
#                                                 if len(new_couples_possibles[m][i]) == len(new_couples_possibles[m][j]) :
#                                                     a_enlever.append(j)
#                                                 else :
#                                                     a_enlever.append(i)
#                             #print(a_enlever)
#                             for i in range(0, len(new_couples_possibles[m])) :
#                                 if i not in a_enlever :
#                                     new_couples_possibles_epuree[m].append(new_couples_possibles[m][i])
#                                 else :
#                                     a_enlever.remove(i)
#                         else :
#                             break
#                 print(fic)
#                 
#                 print(new_couples_possibles)
#                 print(new_couples_possibles_epuree)
#                 print("\n")
#                 
#                 with open("graphes_extension/fichiers_couples_qui_manquent_epures/"+fic, 'wb') as fichier_epure :
#                     mon_pickler = pickle.Pickler(fichier_epure)
#                     mon_pickler.dump(new_couples_possibles_epuree)


liste = ['5J7L_DA_191_3', '5J7L_DA_191_4', '5FDU_1A_301_1', '5J7L_DA_301_2', '5DM6_X_334_1', '5FDU_1A_334_2', '4V9F_0_335_1', '5J7L_DA_335_2', '3JCS_1_137_4', '4V88_A5_290_1', '4V88_A6_314_2', '5J7L_DA_218_3', '4V9F_0_251_2', '1FJG_A_62_8', '5J7L_DA_137_1', '4V9F_0_118_1', '4V9F_0_62_2', '5J7L_DA_271_2', '4V9F_0_224_1', '5DM6_X_197_1', '3GX5_A_138_6', '1FJG_A_317_2', '5J5B_BA_317_1', '1FJG_A_326_1', '5DM6_X_137_3', '5J5B_BA_314_1', '4V9F_0_134_6', '4V9F_0_328_1', '4V9F_0_197_2', '4V9F_0_62_16', '5J7L_DA_282_2', '4V88_A5_137_2', '5FDU_1A_224_3', '5J7L_DA_326_2']

for i in range(len(os.listdir("Extensions/graphes_extension/fichiers_couples_epures"))) :
    for j in range(i+1, len(os.listdir("Extensions/graphes_extension/fichiers_couples_epures"))) :
        fic_1 = os.listdir("Extensions/graphes_extension/fichiers_couples_epures")[i]
        fic_2 = os.listdir("Extensions/graphes_extension/fichiers_couples_epures")[j]
        
        if "pickle" in fic_1 and "pickle" in fic_2 :
        
            elt11 = fic_1.split("_")[3] + "_" + fic_1.split("_")[4] + "_" + fic_1.split("_")[5] + "_" + fic_1.split("_")[6]
            elt12 = fic_1.split("_")[8] + "_" + fic_1.split("_")[9] + "_" + fic_1.split("_")[10] + "_" + fic_1.split("_")[11][:len(fic_1.split("_")[11])-7]
            
            elt21 = fic_2.split("_")[3] + "_" + fic_2.split("_")[4] + "_" + fic_2.split("_")[5] + "_" + fic_2.split("_")[6]
            elt22 = fic_2.split("_")[8] + "_" + fic_2.split("_")[9] + "_" + fic_2.split("_")[10] + "_" + fic_2.split("_")[11][:len(fic_2.split("_")[11])-7]
            
            if (elt11 == elt21 and elt12 == elt22) or (elt12 == elt21 and elt11 == elt22) :
                print("pb")
                print(fic_1)
                print(fic_2)
                
            if elt11 in liste :
                print(elt11)
                print("pb")
            
            if elt12 in liste :
                print(elt12)
                print("pb")
