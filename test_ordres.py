'''
Created on 22 oct. 2018

@author: Coline Gi
'''
import itertools
from itertools import combinations
import pickle
import os

# paires = []
# for element in itertools.product([1,2,3,4],[1,2,3,4]):
#    paires.append(element)
#    
# print(paires)
# 
# 
# print(list(combinations(paires, 4)))
# 
# liste = []
# 
# for elt in list(combinations(paires,4)) :
#     #print(elt)
#     pas_bon = False
#     i = 0
#     while i < len(elt) :
#         j = i + 1
#         while j < len(elt) :
#             if elt[i][0] == elt[j][0] or elt[i][1] == elt[j][1] :
#                 pas_bon = True
#             j = j+1
#         i = i+1
#     if pas_bon == False :
#        liste.append(elt) 
# 
# print(len(liste))        
# print(liste)


with open("script_renommage.py", 'w') as fichier_renommage :
    with open("Extensions/Metrique_toutes_aretes/noms_dans_lordre.pickle", 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        liste_liens = mon_depickler.load()
        
        for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_comparaison_figures") :
            if "png" in fic :
                elt1 =  fic.split("_")[2] + "_" + fic.split("_")[3] + "_" +fic.split("_")[4] + "_" +fic.split("_")[5]
                elt2 =  fic.split("_")[7] + "_" + fic.split("_")[8] + "_" +fic.split("_")[9] + "_" +fic.split("_")[10][:len(fic.split("_")[10])-4]
                
                for elt in liste_liens :
                    if elt[0] == elt2 and elt[1] == elt1 :
                        fic_remplace = "_fichier_%s_fichier_%s.png"%(elt2, elt1)
                        fichier_renommage.write("mv %s %s\n"%(fic, fic_remplace))
                
                
        