'''
Created on 26 oct. 2018

@author: Coline Gi
'''
import pickle
import numpy as np

liste = ['5J7L_DA_191_3', '5J7L_DA_191_4', '5FDU_1A_301_1', '5J7L_DA_301_2', '5DM6_X_334_1', '5FDU_1A_334_2', '4V9F_0_335_1', '5J7L_DA_335_2', '3JCS_1_137_4', '4V88_A5_290_1', '4V88_A6_314_2', '5J7L_DA_218_3', '4V9F_0_251_2', '1FJG_A_62_8', '5J7L_DA_137_1', '4V9F_0_118_1', '4V9F_0_62_2', '5J7L_DA_271_2', '4V9F_0_224_1', '5DM6_X_197_1', '3GX5_A_138_6', '1FJG_A_317_2', '5J5B_BA_317_1', '1FJG_A_326_1', '5DM6_X_137_3', '5J5B_BA_314_1', '4V9F_0_134_6', '4V9F_0_328_1', '4V9F_0_197_2', '4V9F_0_62_16', '5J7L_DA_282_2', '4V88_A5_137_2', '5FDU_1A_224_3', '5J7L_DA_326_2']

with open("fichiers_pickle/a-minor_test2.pickle", 'rb') as fichier_pickle :
    mon_depickler = pickle.Unpickler(fichier_pickle)
    tab_aminor = mon_depickler.load()
    print(tab_aminor)
    
    tab_aminor_2 = []
    for elt in tab_aminor :
        est_dans_liste = False
        for elt_liste in liste :
            noms = elt_liste.split("_")
            print(noms)
            if elt["num_PDB"] == noms[0] and elt["num_ch"] == noms[1] and elt["num_motif"] == int(noms[2]) and elt["num_occ"] == int(noms[3]) :
                est_dans_liste = True
        if est_dans_liste == False and elt["num_motif"] != 301 and elt["num_motif"] != 191 and elt["num_motif"] != 334 and elt["num_motif"] != 335 :
            tab_aminor_2.append(elt)
            
    print(len(tab_aminor_2))
    
    with open("fichier_stats.txt", "w") as fichier_stats :
        ##Nombre d'occ du motif dans la meme molecule et calcul moyenne ##
        fichier_stats.write("Nombre d'occ du motif dans la meme molecule et calcul moyenne\n")
        dico_mol = {}
        for elt in tab_aminor_2 :
            if (elt["num_PDB"],elt["num_ch"]) not in dico_mol :
                dico_mol.update({(elt["num_PDB"],elt["num_ch"]) : 1})
            else :
                dico_mol[(elt["num_PDB"],elt["num_ch"])] += 1
        
        print(dico_mol)
        somme = 0
        tab = []
        for elt in dico_mol.keys() :
            fichier_stats.write(str(elt)+ ": " + str(dico_mol[elt]) + "\n")
            somme += dico_mol[elt]
            tab.append(dico_mol[elt])
        fichier_stats.write("Moyenne : " + str((somme*1.0)/len(dico_mol)) + "\n")
        fichier_stats.write("Ecart-type : " +str(np.std(np.array(tab))) + "\n \n")
        
        ##Proportion des AA dans les occurrences##
        with open("graphs_2.92.pickle", 'rb') as fichier_tout :
            mon_depickler_graphes = pickle.Unpickler(fichier_tout)
            graphes = mon_depickler_graphes.load()
        
            fichier_stats.write("Proportion des AA dans les occurrences\n")
            occ_AA = []
            for elt in tab_aminor_2 :
                num = (elt["num_PDB"], elt["num_ch"])
                if graphes[num].nodes[elt["a_minor"][0]]["nt"] == 'A' and graphes[num].nodes[elt["a_minor"][2]]["nt"] == 'A' :
                    occ_AA.append(elt)
            fichier_stats.write(str(len(occ_AA)) + "\n")
            fichier_stats.write(str((len(occ_AA)*1.0)/len(tab_aminor_2)*100) + "%\n")
            fichier_stats.write(str(occ_AA)+"\n")

            
            dico_types_sse = {}
            for elt in tab_aminor_2 :
                num = (elt["num_PDB"], elt["num_ch"])
                types_sse = ()
                for i in range(5) :
                    if graphes[num].nodes[elt["a_minor"][i]]["part"] not in types_sse :
                        types_sse = types_sse + (graphes[num].nodes[elt["a_minor"][i]]["part"],)
                print(types_sse)
                
                deja_trouvee = False
                for cle in dico_types_sse.keys() :
                    pareil = 0
                    for typ in types_sse :
                        if typ in cle :
                            pareil += 1
                    if pareil == len(cle) and pareil == len(types_sse) :
                        dico_types_sse[cle] += 1
                        deja_trouvee = True
                        print("ramousnif")
                        print(cle)
                
                if deja_trouvee == False :
                    dico_types_sse.update({types_sse : 1})
                
            
            print(dico_types_sse)
            somme = 0
            for cle in dico_types_sse :
                fichier_stats.write(str(cle) + " : " +  str(dico_types_sse[cle]) +  "\n")
                somme += dico_types_sse[cle]
                
            print(somme)
                
                    
                
                
                    
                    
                    