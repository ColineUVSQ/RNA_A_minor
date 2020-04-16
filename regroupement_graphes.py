'''
Created on 12 mars 2019

@author: coline
'''

import os
import pickle
import recup_data.constantes
from recup_data.constantes import EXTENSION_PATH_TAILLE
from recup_data.constantes import EXTENSION_PATH

## Creation d'un dictionnaire avec tous les graphes de comparaison à partir des fichiers isolés
def creation_dico_comp(taille_ext):
    compteur_pas_dico = 0
    dico_graphes = {}

    for fic in os.listdir(EXTENSION_PATH_TAILLE%taille_ext) :
        if "graphe_comp" in fic and "test" not in fic :
            with open(EXTENSION_PATH_TAILLE%taille_ext+fic, 'rb') as fichier :
                mon_depickler = pickle.Unpickler(fichier)
                dico_graphe = mon_depickler.load()
#                 print(fic)
#                 print(dico_graphe)
                if type(dico_graphe) is dict :
                    for cle in dico_graphe.keys() :
#                         print(dico_graphe[cle].edges.data())
#                         print(cle)
                        if len(cle[0].split("_")) != 5 or len(cle[1].split("_")) != 5 :
                            if len(cle[0].split("_")) != 5 :
                                cle_new_0 = cle[0].split("_")[0] + "_" + cle[0].split("_")[1] + "_" + cle[0].split("_")[2] + "_" + cle[0].split("_")[3] + "_" + cle[0].split("_")[4]
                            if len(cle[1].split("_")) != 5 :
                                cle_new_1 = cle[1].split("_")[0] + "_" + cle[1].split("_")[1] + "_" + cle[1].split("_")[2] + "_" + cle[1].split("_")[3] + "_" + cle[1].split("_")[4]
                            cle_new = (cle_new_0, cle_new_1)
                        else :
                            print("petit rat")
                            if "pickle" in cle[1] :
                                cle_1 = cle[1][:len(cle[1])-7] 
                            else :
                                cle_1 = cle[1]
                            cle_0 = cle[0]
                            print(cle)
                            cle_new = (cle_0, cle_1)
                        if cle_new in dico_graphes.keys() :
                            print(cle)
                            print(fichier)
                            print("bizarre")
                        dico_graphes.update({cle_new : dico_graphe[cle]})
                else :
                    print("bizarre")
                    print(fic)
                    compteur_pas_dico += 1
                    
    print(compteur_pas_dico)
    print(len(dico_graphes))        
    
#     if len(dico_graphes) < 4005 :
#         with open("/home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/result_graphes_extension_%s/dico_comp_complet_metrique_toutes_aretes_coeff_all1_taille_%s.pickle"%(taille_ext, taille_ext), 'rb') as fichier_ajout :
#             mon_depickler = pickle.Unpickler(fichier_ajout)
#             ancien_dico_graphes = mon_depickler.load()
#             
#             for cle in ancien_dico_graphes.keys() :
#                 if (cle[0], cle[1]) not in dico_graphes.keys() and (cle[1], cle[0]) not in dico_graphes.keys() :
#                     dico_graphes.update({cle : ancien_dico_graphes[cle]})
    
    print(len(dico_graphes))
    with open(EXTENSION_PATH%taille_ext+"dico_comp_complet_metrique_toutes_aretes_coeff_all1_taille_%s.pickle"%taille_ext, 'wb') as fichier_ecriture :
        mon_pickler = pickle.Pickler(fichier_ecriture)
        mon_pickler.dump(dico_graphes)
    
if __name__ == '__main__':
    creation_dico_comp(6)    

