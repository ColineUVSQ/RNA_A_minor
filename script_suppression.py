'''
Created on 10 déc. 2018

@author: coline
'''

import os

liste_doublons = ['5J7L_DA_191_3', '5J7L_DA_191_4', '5FDU_1A_301_1', '5J7L_DA_301_2', '5DM6_X_334_1', '5FDU_1A_334_2', '4V9F_0_335_1', '5J7L_DA_335_2', '3JCS_1_137_4', '4V88_A5_290_1', '4V88_A6_314_2', '5J7L_DA_218_3', '4V9F_0_251_2', '1FJG_A_62_8', '5J7L_DA_137_1', '4V9F_0_118_1', '4V9F_0_62_2', '5J7L_DA_271_2', '4V9F_0_224_1', '5DM6_X_197_1', '3GX5_A_138_6', '1FJG_A_317_2', '5J5B_BA_317_1', '1FJG_A_326_1', '5DM6_X_137_3', '5J5B_BA_314_1', '4V9F_0_134_6', '4V9F_0_328_1', '4V9F_0_197_2', '4V9F_0_62_16', '5J7L_DA_282_2', '4V88_A5_137_2', '5FDU_1A_224_3', '5J7L_DA_326_2']
liste_graphes_a_refaire = []

for fic in os.listdir("graphes_extension/fichiers_couples_qui_manquent") :
    if "pickle" in fic :
        print(fic[17:len(fic)-7])
        liste_graphes_a_refaire.append(fic[17:len(fic)-7])

for fic in os.listdir("graphes_extension/fichiers_couples_qui_manquent/sans_boucle") :
    if "pickle" in fic :
        print(fic[17:len(fic)-7])
        liste_graphes_a_refaire.append(fic[17:len(fic)-7])

with open("script_suppression_graphes_a_refaire.sh", "w") as fichier :
    for elt in liste_graphes_a_refaire :
        fichier.write("find *" +elt+"* | xargs rm\n")