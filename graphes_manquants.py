'''
Created on 14 déc. 2018

@author: coline
'''
import pickle
import os

liste = ['5J7L_DA_191_3', '5J7L_DA_191_4', '5FDU_1A_301_1', '5J7L_DA_301_2', '5DM6_X_334_1', '5FDU_1A_334_2', '4V9F_0_335_1', '5J7L_DA_335_2', '3JCS_1_137_4', '4V88_A5_290_1', '4V88_A6_314_2', '5J7L_DA_218_3', '4V9F_0_251_2', '1FJG_A_62_8', '5J7L_DA_137_1', '4V9F_0_118_1', '4V9F_0_62_2', '5J7L_DA_271_2', '4V9F_0_224_1', '5DM6_X_197_1', '3GX5_A_138_6', '1FJG_A_317_2', '5J5B_BA_317_1', '1FJG_A_326_1', '5DM6_X_137_3', '5J5B_BA_314_1', '4V9F_0_134_6', '4V9F_0_328_1', '4V9F_0_197_2', '4V9F_0_62_16', '5J7L_DA_282_2', '4V88_A5_137_2', '5FDU_1A_224_3']

# tri_new = []
# with open("fichiers_tries.pickle", "rb") as fichier_tri :
#     mon_depickler_tri = pickle.Unpickler(fichier_tri)
#     tri = mon_depickler_tri.load()
#     
#     for elt in tri :
#         print(elt[8:len(elt)-7])
#         if elt[8:len(elt)-7] not in liste :
#             tri_new.append(elt)
#     print(len(tri_new))
# print(tri_new)
#     
# with open("fichiers_tries.pickle", "wb") as fichier_tri_new :
#     mon_pickler = pickle.Pickler(fichier_tri_new)
#     mon_pickler.dump(tri_new)

# with open("dico_graphe.pickle", 'rb') as fichier_graphe :
#     mon_depickler = pickle.Unpickler(fichier_graphe)
#     dico_graphe = mon_depickler.load()
#      
#     print(dico_graphe.keys())
#      
#     dico = {}
#  
#     for cle in dico_graphe.keys() :
#         if cle[0][8:] not in liste and cle[1][8:] not in liste :
#             if cle[0] not in dico.keys() :
#                 dico.update({cle[0] : [cle[1]] })
#             else :
#                 dico[cle[0]].append(cle[1])
#             if cle[1] not in dico.keys() :
#                 dico.update({cle[1] : [cle[0]] })
#             else :
#                 dico[cle[1]].append(cle[0])
#              
#     print(dico)
#      
#     for elt in dico.keys() :
#         if len(dico[elt]) < 97 :
#             print(elt)
#             print(len(dico[elt]))

# with open("enlever_doubles.sh", 'w') as script :
#     script.write("rm ")
with open("dico_graphe_epure.pickle", 'rb') as fichier_graphe :
        mon_depickler = pickle.Unpickler(fichier_graphe)
        dico_graphe = mon_depickler.load()
    
        with open("tab_nexiste_pas.pickle", 'wb') as fichier :
            mon_pickler = pickle.Pickler(fichier)
            
            nexiste_pas = []
        
            print(len(dico_graphe))
            deja_vu = []
            compteur = 0
            for fic_1 in os.listdir("graphes_extension") :
                for fic_2 in os.listdir("graphes_extension") :
                    if fic_1 != fic_2 and "pickle" in fic_1 and "pickle" in fic_2 and (fic_2,fic_1) not in deja_vu and (fic_1,fic_2) not in deja_vu:
                        print(fic_1[:len(fic_1)-7])
                        print(fic_2[:len(fic_2)-7])
                        if (fic_1[:len(fic_1)-7], fic_2[:len(fic_2)-7]) not in dico_graphe.keys() and (fic_2[:len(fic_2)-7], fic_1[:len(fic_1)-7]) not in dico_graphe.keys() :
                            deja_vu.append((fic_1,fic_2))
                            if "couples_possibles_"+ fic_1[:len(fic_1)-7] + "_" + fic_2 not in os.listdir("graphes_extension/fichiers_couples_qui_manquent") and "couples_possibles_"+ fic_2[:len(fic_2)-7] + "_" + fic_1 not in os.listdir("graphes_extension/fichiers_couples_qui_manquent") :
                                print(fic_1 + " " + fic_2)
                                compteur += 1
                                nexiste_pas.append((fic_1, fic_2))
            
            mon_pickler.dump(nexiste_pas)                    
            print(len(deja_vu))
            print(compteur)
         
            compte = 0
            deja_vu = []
            for fic in os.listdir("graphes_extension/fichiers_couples_qui_manquent") :
                element1 = fic.split('_')[2] + '_' + fic.split('_')[3] + '_' + fic.split('_')[4] + '_' + fic.split('_')[5] + '_' + fic.split('_')[6]
                element2 = fic.split('_')[7] + '_' + fic.split('_')[8] + '_' + fic.split('_')[9] + '_' + fic.split('_')[10] + '_' + fic.split('_')[11][:len(fic.split('_')[11])-7]
                 
                pas_bon = False
                for elt in liste :
                    if elt in element1 or elt in element2 :
                        pas_bon = True
                 
             
    #             if (element1, element2) in dico_graphe.keys() or (element2, element1) in dico_graphe.keys() :
    #                 print("deja_fait_groupe1")
    #             elif pas_bon == False :
    #                 print("deja_fait_groupe2")
    #             else :
    #                 print("doublon") 
                  
                
                if (element1, element2) in deja_vu or (element2, element1) in deja_vu :
                    print(fic)
                    print("en double")
                    #script.write(fic + " ")
                    compte += 1
                 
                     
                deja_vu.append((element1, element2))
            print(compte)


# with open("graphes_extension/fichiers_couples_qui_manquent/couples_possibles_fichier_4V9F_0_48_26_fichier_3JCS_1_25_46.pickle", 'rb') as fichier_1 :
#     mon_depickler = pickle.Unpickler(fichier_1)
#     fichier1 = mon_depickler.load()
#     with open("graphes_extension/fichiers_couples_qui_manquent_epures/couples_possibles_fichier_4V9F_0_48_26_fichier_3JCS_1_25_46.pickle", 'rb') as fichier_2 :
#         mon_depickler_2 = pickle.Unpickler(fichier_2)
#         fichier2 = mon_depickler_2.load()
#         
#         print(fichier1)
#         print(fichier2)
                        