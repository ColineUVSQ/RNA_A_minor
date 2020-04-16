'''
Created on 8 févr. 2019

@author: coline
'''
import pickle
from recup_data.sous_graphe_commun_max import dans_graphe
from recup_data.calcul_sim import calcul_sim_non_cov_sans_motif

liste = ['5J7L_DA_191_3', '5J7L_DA_191_4', '5FDU_1A_301_1', '5J7L_DA_301_2', '5DM6_X_334_1', '5FDU_1A_334_2', '4V9F_0_335_1', '5J7L_DA_335_2', '3JCS_1_137_4', '4V88_A5_290_1', '4V88_A6_314_2', '5J7L_DA_218_3', '4V9F_0_251_2', '1FJG_A_62_8', '5J7L_DA_137_1', '4V9F_0_118_1', '4V9F_0_62_2', '5J7L_DA_271_2', '4V9F_0_224_1', '5DM6_X_197_1', '3GX5_A_138_6', '1FJG_A_317_2', '5J5B_BA_317_1', '1FJG_A_326_1', '5DM6_X_137_3', '5J5B_BA_314_1', '4V9F_0_134_6', '4V9F_0_328_1', '4V9F_0_197_2', '4V9F_0_62_16', '5J7L_DA_282_2', '4V88_A5_137_2', '5FDU_1A_224_3', '5J7L_DA_326_2']

## pour ajouter aux comparaisons les nucléotides liés au motif => on en a plus besoin
with open("dico_graphe_epure_en_tout.pickle", 'rb') as fichier_graphe :
    mon_depickler = pickle.Unpickler(fichier_graphe)
    dico_graphe = mon_depickler.load()
    
    dico_sim = {}
    for cle in dico_graphe.keys() :
#         print(cle)
        
        pas_bon = False
        for elt in liste :
            if elt in cle[0] or elt in cle[1] :
                pas_bon = True
        
        if pas_bon == False :
                with open("graphes_extension/"+cle[0]+".pickle", 'rb') as fichier_graphe1 :
                    mon_depickler_graphe1 = pickle.Unpickler(fichier_graphe1)
                    graphe1 = mon_depickler_graphe1.load()
                    
                    with open("graphes_extension/"+cle[1]+".pickle", 'rb') as fichier_graphe2 :
                        mon_depickler_graphe2 = pickle.Unpickler(fichier_graphe2)
                        graphe2 = mon_depickler_graphe2.load()
                        
                        for i in range(1,5) :
                            for u1,v1,data1 in graphe1.edges(data=True) :
#                                 if ('fichier_1GID_B_25_20', 'fichier_4V9F_0_287_2') == cle :
#                                         print("ramousnif")
#                                         print(i)
#                                         print(u1)
#                                         print(v1)
#                                         print(type(v1))
#                                         print(data1["label"])
                                if (u1 == i and data1["label"] != 'B53') or (v1 == i and data1["label"] != 'B53') : 
#                                     if ('fichier_1GID_B_25_20', 'fichier_4V9F_0_287_2') == cle :
#                                         print("ramousnif")
#                                         print(u1)
#                                         print(v1)
                                    for u2,v2,data2 in graphe2.edges(data=True) :
                                                             
                                        if (u2 == i and data2["label"] != 'B53') or (v2 == i and data2["label"] != 'B53'): 

                                            if (data1["label"] == data2["label"]) or (data1["label"] != 'CWW' and data2["label"] != 'CWW'):
                                                
                                                if u1 == i :
                                                    voisin_1 = v1
                                                else :
                                                    voisin_1 = u1
                                                    
                                                if u2 == i :
                                                    voisin_2 = v2
                                                else :
                                                    voisin_2 = u2
                                                
                                                if voisin_1 not in [1,2,3,4,5] and voisin_2 not in [1,2,3,4,5] and graphe1.nodes[voisin_1]["type"] == graphe2.nodes[voisin_2]["type"] :
                                                    if (voisin_1,voisin_2) not in dico_graphe[cle].nodes() and dans_graphe(dico_graphe[cle],(voisin_1, voisin_2)) == False :     
                                                        dico_graphe[cle].add_node((voisin_1,voisin_2))
                                                        
                                                        if (voisin_1, voisin_2) not in dico_graphe[cle][(i,i)] :
                                                            if data1["label"] == 'CWW' :
                                                                typ = "CAN"
                                                            else :
                                                                typ = "NON_CAN"
                                                            #dico_graphe[cle].add_edge((voisin_1,voisin_2), (i,i), type=typ)
                                                            dico_graphe[cle].add_edge((i,i), (voisin_1,voisin_2), type=typ)
                                                    
                        sim =  calcul_sim_non_cov_sans_motif(graphe1, graphe2, dico_graphe[cle]) 
                        if ('fichier_5FDU_1A_197_3', 'fichier_5J7L_DA_197_4') == cle :  
                            print(sim)
                        dico_sim.update({cle : sim[0]})                                      
    with open("dico_graphe_epure_en_tout_test.pickle", 'wb') as fichier_ecriture :
        mon_pickler = pickle.Pickler(fichier_ecriture)
        mon_pickler.dump(dico_graphe) 
        print(len(dico_graphe))      
    with open("dico_graphe_epure_en_tout_test_sim.pickle", 'wb') as fichier_ecriture_sim :
        mon_pickler_sim= pickle.Pickler(fichier_ecriture_sim)
        mon_pickler_sim.dump(dico_sim)