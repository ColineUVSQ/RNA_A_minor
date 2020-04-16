'''
Created on 6 déc. 2018

@author: coline
'''
import networkx as nx
import pickle

liste = ['5J7L_DA_191_3', '5J7L_DA_191_4', '5FDU_1A_301_1', '5J7L_DA_301_2', '5DM6_X_334_1', '5FDU_1A_334_2', '4V9F_0_335_1', '5J7L_DA_335_2', '3JCS_1_137_4', '4V88_A5_290_1', '4V88_A6_314_2', '5J7L_DA_218_3', '4V9F_0_251_2', '1FJG_A_62_8', '5J7L_DA_137_1', '4V9F_0_118_1', '4V9F_0_62_2', '5J7L_DA_271_2', '4V9F_0_224_1', '5DM6_X_197_1', '3GX5_A_138_6', '1FJG_A_317_2', '5J5B_BA_317_1', '1FJG_A_326_1', '5DM6_X_137_3', '5J5B_BA_314_1', '4V9F_0_134_6', '4V9F_0_328_1', '4V9F_0_197_2', '4V9F_0_62_16', '5J7L_DA_282_2', '4V88_A5_137_2', '5FDU_1A_224_3', '5J7L_DA_326_2']

with open("fichier_similarite_qui_manquent.txt", 'r') as fichier :

    ligne = " "
    dico_sim = {}
    dico_graphe = {}
    while len(ligne) != 0 :
        ligne1 = fichier.readline()
        ligne2 = fichier.readline()
        element1 = ligne1[:len(ligne1)-1]
        element2 = ligne2[:len(ligne2)-1]
        print(element1)
        print(element2)
        
        graphe_commun = nx.MultiGraph()
        #for i in range(2) :
        sommets = fichier.readline()
        aretes = fichier.readline()
        if len(sommets.split(":")) > 1:
            sommets = sommets.split(":")[1][2:len(sommets.split(":")[1])-2]
            #print(sommets)
            sommets = sommets.split("((")
            #print(sommets)
            for sommet in sommets :
                if len(sommet) > 0 :
                    #print(sommet)
                    sommet_1 = int(sommet.split(",")[0])
                    sommet_2 = int(sommet.split(",")[1][:len(sommet.split(",")[1])-1])
                    graphe_commun.add_node((sommet_1, sommet_2))
            print(graphe_commun.nodes.data())
            

            aretes = aretes[1:len(aretes)-2]
            aretes = aretes.split("((")
            #print(aretes)
            
            for arete in aretes :
                if len(arete) > 0 :
                    arete_morceau = arete.split(",")
                    sommet_1 = int(arete_morceau[0])
                    sommet_2 = int(arete_morceau[1][:len(arete_morceau[1])-1])
                    sommet_3 = int(arete_morceau[2][2:])
                    sommet_4 = int(arete_morceau[3][1:len(arete_morceau[3])-1])
                    
                    typ = arete_morceau[4].split(":")[1][2:len(arete_morceau[4].split(":")[1])-3]
                    
                    graphe_commun.add_edge((sommet_1, sommet_2), (sommet_3, sommet_4), type=typ)
            print(graphe_commun.edges.data())       
                    
        ligne = fichier.readline()
        sim = ligne[:len(ligne)-1].split(":")[1]
        dico_sim.update({(element1, element2) : sim})
        dico_graphe.update({(element1, element2) : graphe_commun})
        ligne = fichier.readline()
        
#     with open("dico_sim_petit.pickle", 'wb') as fichier_sim :
#         mon_pickler = pickle.Pickler(fichier_sim)
#         mon_pickler.dump(dico_sim)
        
    with open("dico_graphe_epure.pickle", 'rb') as fichier_graphe :
        mon_depickler = pickle.Unpickler(fichier_graphe)
        dico_graphe_grand = mon_depickler.load()
        
        print(len(dico_graphe_grand))
                
        dico_graphe_grand_epure = {}
        for cle in dico_graphe_grand.keys() :
            doublon = False
            print(cle)
            for elt in liste :
                if elt in cle[0] or elt in cle[1] :
                    doublon = True
            if doublon == False :
                dico_graphe_grand_epure.update({(cle[0], cle[1]) : dico_graphe_grand[cle] })
        
        print(len(dico_graphe_grand_epure))
        
        dico_tot = {}
        dico_tot.update(dico_graphe_grand_epure)
        
        dico_graphe_epure = {}
        for cle in dico_graphe.keys() :
            doublon = False
            for elt in liste :
                if elt in cle[0] or elt in cle[1] :
                    doublon = True
            if doublon == False :
                dico_graphe_epure.update({(cle[0], cle[1]) : dico_graphe[cle] })
        
        dico_tot.update(dico_graphe_epure)
        
        print(len(dico_tot))
        
        print("ramou")
        print(dico_tot[('fichier_5FDU_1A_74_7','fichier_4V9F_0_207_3')].nodes.data())
        
    with open("dico_graphe_epure_en_tout.pickle", 'wb') as fichier_graphe_epure :
        mon_pickler = pickle.Pickler(fichier_graphe_epure)
        mon_pickler.dump(dico_tot)
    