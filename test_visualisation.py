'''
Created on 31 janv. 2019

@author: coline
'''
import pickle
import numpy as np

liste = ['5J7L_DA_191_3', '5J7L_DA_191_4', '5FDU_1A_301_1', '5J7L_DA_301_2', '5DM6_X_334_1', '5FDU_1A_334_2', '4V9F_0_335_1', '5J7L_DA_335_2', '3JCS_1_137_4', '4V88_A5_290_1', '4V88_A6_314_2', '5J7L_DA_218_3', '4V9F_0_251_2', '1FJG_A_62_8', '5J7L_DA_137_1', '4V9F_0_118_1', '4V9F_0_62_2', '5J7L_DA_271_2', '4V9F_0_224_1', '5DM6_X_197_1', '3GX5_A_138_6', '1FJG_A_317_2', '5J5B_BA_317_1', '1FJG_A_326_1', '5DM6_X_137_3', '5J5B_BA_314_1', '4V9F_0_134_6', '4V9F_0_328_1', '4V9F_0_197_2', '4V9F_0_62_16', '5J7L_DA_282_2', '4V88_A5_137_2', '5FDU_1A_224_3', '5J7L_DA_326_2']

with open("grands_graphes.pickle", 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        dico_graphes = mon_depickler.load()
        
        with open("fichiers_pickle/a-minor_test2.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            tab_aminor = mon_depickler.load()
            
            for occ in tab_aminor :
                
                est_la = False
                for elt in liste :
                    if occ["num_PDB"] in elt and occ["num_ch"] in elt and str(occ["num_motif"]) in elt and str(occ["num_occ"]) in elt  :
                        est_la = True
                        
                if est_la == False :
                    num = (occ["num_PDB"], occ["num_ch"], occ["num_motif"],occ["num_occ"])
                    
                    mini = np.min(dico_graphes[num].nodes())
                    maxi = np.max(dico_graphes[num].nodes())
                    
                    
                