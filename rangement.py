'''
Created on 28 nov. 2018

@author: coline
'''

import os
import pickle

## Regarder les couples deja formes pour voir ceux qui manquent encore

liste = ['5J7L_DA_191_3', '5J7L_DA_191_4', '5FDU_1A_301_1', '5J7L_DA_301_2', '5DM6_X_334_1', '5FDU_1A_334_2', '4V9F_0_335_1', '5J7L_DA_335_2', '3JCS_1_137_4', '4V88_A5_290_1', '4V88_A6_314_2', '5J7L_DA_218_3', '4V9F_0_251_2', '1FJG_A_62_8', '5J7L_DA_137_1', '4V9F_0_118_1', '4V9F_0_62_2', '5J7L_DA_271_2', '4V9F_0_224_1', '5DM6_X_197_1', '3GX5_A_138_6', '1FJG_A_317_2', '5J5B_BA_317_1', '1FJG_A_326_1', '5DM6_X_137_3', '5J5B_BA_314_1', '4V9F_0_134_6', '4V9F_0_328_1', '4V9F_0_197_2', '4V9F_0_62_16', '5J7L_DA_282_2', '4V88_A5_137_2', '5FDU_1A_224_3']


def ajout_dico(dico, fic, fichier):
    element1 = fic.split('_')[3] + '_' + fic.split('_')[4] + '_' + fic.split('_')[5] + '_' + fic.split('_')[6]
    element2 = fic.split('_')[8] + '_' + fic.split('_')[9] + '_' + fic.split('_')[10] + '_' + fic.split('_')[11][:len(fic.split('_')[11])-7]
    
    if element1 not in liste and element2 not in liste :
    
        if element1 not in dico.keys() :
            dico.update({element1 : [element2]})
        else :
            if element2 not in dico[element1] :
                dico[element1].append(element2)
            else :
                print(fic)
                print("probleme")
            
        if element2 not in dico.keys() :
            dico.update({element2 : [element1]})
        else :
            if element1 not in dico[element2] :
                dico[element2].append(element1)
            else :
                print(fic)
                print("probleme")
    else :
        fichier.write(str(fic) + " ")

    return dico

def ajout_dico_graphe(dico_graphe, dico) :
    
    for cle in dico_graphe.keys() :
        if cle[0][8:] not in liste and cle[1][8:] not in liste :
    
            if cle[0] not in dico.keys() :
                dico.update({cle[0] : [cle[1]]})
            else :
                if cle[1] not in dico[cle[0]] :
                    dico[cle[0]].append(cle[1])
                else :
                    print(cle)
                    print("probleme")
                
            if cle[1] not in dico.keys() :
                dico.update({cle[1]: [cle[0]]})
            else :
                if cle[0] not in dico[cle[1]] :
                    dico[cle[1]].append(cle[0])
                else :
                    print(cle)
                    print("probleme")
    return dico
    
#dico = {}
#  
# with open("fichier_suppression.sh", 'w') as fichier :
#     fichier.write("rm ")
#      
#     print("epures")
#     for fic in os.listdir("graphes_extension/fichiers_couples_epures/") :
#         #if "pickle" in fic :
#         dico = ajout_dico(dico, fic, fichier)
#              
#     print("manquent")
#     for fic in os.listdir("graphes_extension/fichiers_couples_qui_manquent") :
#          
#         #if "pickle" in fic :
#         dico = ajout_dico(dico, fic, fichier)
#      
#     print("petit")
#     for fic in os.listdir("graphes_extension/fichiers_couples_qui_manquent_petit") :
#          
#         #if "pickle" in fic :
#         dico = ajout_dico(dico, fic, fichier)

# print(len(dico))
# print(dico)
  
# for elt in dico.keys() :
#     print(elt)
#     print(len(dico[elt]))
#     
#     for graphe in os.listdir("graphes_extension") :
#         if "pickle" in graphe :
#             element = graphe.split('_')[1] + "_" + graphe.split('_')[2] + "_" + graphe.split('_')[3] + "_" + graphe.split('_')[4][:len(graphe.split('_')[4])-7]
#             
#             if elt != element and element not in dico[elt] :
#                 print(str(elt) + " " + str(element))

# for elt in dict.keys() :
#     if '44' in elt or '58' in elt : print(elt)
#     print(len(dict[elt]))


with open("dico_graphe_epure_en_tout.pickle", 'rb') as fichier_graphe :
        mon_depickler = pickle.Unpickler(fichier_graphe)
        dico_graphe = mon_depickler.load()
         
        print(len(dico_graphe))
#         print(dico_graphe.keys())
        
        dico = {}
        
        ajout_dico_graphe(dico_graphe, dico)
        #print(dico)
        
        for cle in dico.keys() :
            print(cle)
            if len(dico[cle]) < 89 :
                print(len(dico[cle]))
        
        
#         compteur = 0
#         for fic in os.listdir("graphes_extension/fichiers_couples_epures") :
#             if 'pickle' in fic :
#                 element1 = fic.split('_')[2]+ "_" + fic.split('_')[3] + '_' + fic.split('_')[4] + '_' + fic.split('_')[5] + '_' + fic.split('_')[6]
#                 element2 = fic.split('_')[7] + "_" + fic.split('_')[8] + '_' + fic.split('_')[9] + '_' + fic.split('_')[10] + '_' + fic.split('_')[11][:len(fic.split('_')[11])-7]
#                 compteur = compteur + 1
#                 #print(compteur)
#                  
#                 if (element1, element2) not in dico_graphe.keys() :
#                     print(str(element1) + " " + str(element2))
#                     #print(dico_graphe[(element1, element2)].nodes.data())
#                     print("la")
            
                


