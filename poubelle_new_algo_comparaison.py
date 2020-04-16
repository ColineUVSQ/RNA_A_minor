'''
Created on 15 nov. 2019

@author: coline
'''
import pickle
import networkx as nx
import csv
import os
import matplotlib.pyplot as plt
from recup_data.constantes import NEW_EXTENSION_PATH_TAILLE, PATH_MMCIF
from recup_data.new_algo_comparaison import comparaison

def graphe_des_sim_inter_groupes(liste_num_ARN):
    liste_tout = []
    for elt in liste_num_ARN :
        if elt not in ["arnt_16s", "arnt_16s_arnm"] :
            with open("/media/coline/Maxtor/noeuds_a_garder_selon_centrality_%s.pickle"%elt, 'rb') as fichier_sortie :
                mon_depickler = pickle.Unpickler(fichier_sortie)
                liste_a_garder_noms = mon_depickler.load()    
                  
                for element in liste_a_garder_noms :
                    liste_tout.append((elt, element))
        else :
            with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
                mon_depickler = pickle.Unpickler(fichier_sortie)
                liste_a_garder_noms = mon_depickler.load()    
                  
                for element in liste_a_garder_noms :
                    liste_tout.append((elt, element))
  
    print(liste_tout)  
    graphe_complet = nx.Graph()
                 
    compteur = 1
    #for groupe in groupes_homologues :
    for elt in liste_tout :
        graphe_complet.add_node(compteur, type=elt[0], nom=elt[1])
        compteur += 1
    print(graphe_complet.nodes.data())
    print(graphe_complet.number_of_nodes())
# # 
    with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_avec_manque_['23S', '18S', '16S', 'Ribozyme', 'Riboswitch', 'SRP', '28S', '25S', 'Intron'].pickle", 'rb') as fichier_sim_sans_doublons_avec_manque :
            mon_depickler_2 = pickle.Unpickler(fichier_sim_sans_doublons_avec_manque)
            graphe_manque = mon_depickler_2.load()
            dico_manque_new = {}
            
            for u,v,data in graphe_manque.edges(data=True) :
                graphe_complet.add_edge(u,v,**data)
            
            #graphe_complet.add_edges_from(graphe_manque.edges(), **graphe_manque.edges.data())
            
            compter = 0
            for i in range(1, graphe_complet.number_of_nodes()+1) :
                noeud1 = i
                data1 = graphe_complet.nodes[i]
                for j in range(i+1, graphe_complet.number_of_nodes()+1) : 
                    noeud2 = j
                    data2 = graphe_complet.nodes[j]
                    print(compter)
                    
                    if (noeud1, noeud2) not in graphe_complet.edges() :
#                     if (data1["nom"][0]+"_"+str(data1["nom"][1]), data2["nom"][0]+"_"+str(data2["nom"][1])) in dico_manque.keys() :
#                         graphe_complet.add_edge(noeud1, noeud2, sim=dico_manque[(data1["nom"][0]+"_"+str(data1["nom"][1]), data2["nom"][0]+"_"+str(data2["nom"][1]))]["sim"])
#                     elif (data2["nom"][0]+"_"+str(data2["nom"][1]), data1["nom"][0]+"_"+str(data1["nom"][1])) in dico_manque.keys() :
#                         graphe_complet.add_edge(noeud1, noeud2, sim=dico_manque[(data2["nom"][0]+"_"+str(data2["nom"][1]), data1["nom"][0]+"_"+str(data1["nom"][1]))]["sim"])
#                     else :
                        with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+data1["nom"][0]+"_"+str(data1["nom"][1])+".pickle", 'rb') as fichier1 :
                                mon_depickler_1 = pickle.Unpickler(fichier1)
                                graphe1 = mon_depickler_1.load()
                                           
                                with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+data2["nom"][0]+"_"+str(data2["nom"][1])+".pickle", 'rb') as fichier2 :
                                        mon_depickler = pickle.Unpickler(fichier2)
                                        graphe2 = mon_depickler.load()
                                        graphe_commun_max, sim_max = comparaison(graphe1, graphe2, "petit rat") 
                                        graphe_complet.add_edge(noeud1, noeud2, sim=sim_max)
                                        dico_manque_new.update({(data1["nom"][0]+"_"+str(data1["nom"][1]), data2["nom"][0]+"_"+str(data2["nom"][1])) : {"graphe" : graphe_commun_max, "sim" : sim_max} })
                        compter += 1
        #                     print(len(graphe_complet[noeud]))
        #                     print(graphe_complet[noeud])
#                                   
#                     print(compter)
#                     print(graphe_complet.number_of_nodes())
#                     print(graphe_complet.number_of_edges())
#             
# #             compteur =0
# #             for groupe in groupes_homologues :
# #                 if ('4v7s',20) in groupe and ('4tud', 16) in groupe :
# #                     print("meme_groupe")
# #                 if (('4v7s',20) in groupe and not ('4tud', 16) in groupe) or (not ('4v7s',20) in groupe and ('4tud', 16) in groupe) :
# #                     print(compteur)
# #                     print("pas meme groupe")
# #                 compteur += 1
# #                 
# #             if ('4v7s_20', '4tud_16') in big_dico.keys() or  ('4tud_16', '4v7s_20') in big_dico.keys() :
# #                 print("trouve")
    #print(graphe_complet.number_of_nodes())
    #print(graphe_complet.number_of_edges())         
    with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_avec_manque_%s.pickle"%liste_num_ARN, 'wb') as fichier_graphe_sim_sans_doublons_avec_manque :
            mon_pickler = pickle.Pickler(fichier_graphe_sim_sans_doublons_avec_manque)
            mon_pickler.dump(graphe_complet) 
#           

            
            #print(graphe_complet.edges.data())      
            a_enlever = []
            for u,v,data in graphe_complet.edges(data=True) :
                if not isinstance(data["sim"], dict) :
                    if data["sim"] < 0.6 :
                        a_enlever.append((u,v))
                else :
                    if data["sim"]["sim"] < 0.6 :
                        a_enlever.append((u,v))
    #     print(a_enlever)
            for elt in a_enlever :
                graphe_complet.remove_edge(elt[0], elt[1])
            print(graphe_complet.number_of_nodes())
            print(graphe_complet.number_of_edges())
#                 
#                 
#                 
            with open("/media/coline/Maxtor/fichier_csv_new_data_seuil_0.6_inter_groupe_tout_avec_groupe_special.csv",'w') as fichier_csv:
                    csvwriter = csv.writer(fichier_csv)
                    csvwriter.writerow(["source", "target", "label"])
                     
                    seuil = 0.6
                    for u,v,data in graphe_complet.edges(data=True) :
                         #if data["sim"] > seuil :
                        if not isinstance(data["sim"], dict) :
                            csvwriter.writerow([u,v,round(data["sim"],2)])
                        
                        else :
                            csvwriter.writerow([u,v,round(data["sim"]["sim"],2)])
                        with open("/media/coline/Maxtor/fichier_csv_new_data_seuil_0.6_inter_groupe_noeud_tout_avec_groupe_special.csv",'w') as fichier_csv:
                            csvwriter2 = csv.writer(fichier_csv)
                            csvwriter2.writerow(["id", "label", "type"])
                                 
                            for noeud,data in graphe_complet.nodes(data=True) :
                                
                                csvwriter2.writerow([noeud, data["nom"][0] + "_"+str(data["nom"][1]), graphe_complet.nodes[noeud]["type"]])
#                   #         with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_25s.pickle", 'rb') as fichier_graphe_sim_sans_doublons :
#             mon_depickler = pickle.Unpickler(fichier_graphe_sim_sans_doublons)
#             graphe_complet = mon_depickler.load() 
    with open("/media/coline/Maxtor/dico_sim_manque_%s.pickle"%liste_num_ARN, 'wb') as fichier_sim_sans_doublons_avec_manque :
                mon_pickler_2 = pickle.Pickler(fichier_sim_sans_doublons_avec_manque)
                mon_pickler_2.dump(dico_manque_new)                                     

def ecriture_csv_gephi_graphe(num_ARN):
    with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_avec_manque_%s_new_noms_res_3a.pickle"%num_ARN, 'rb') as fichier_graphe_sim_sans_doublons_avec_manque :
        mon_depickler = pickle.Unpickler(fichier_graphe_sim_sans_doublons_avec_manque)
        graphe_complet = mon_depickler.load() 
#           
        with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_1_new_data.pickle"%4, 'rb') as fichier_rmsd :
            mon_depickler = pickle.Unpickler(fichier_rmsd)
            rmsd = mon_depickler.load()
            
            #print(graphe_complet.edges.data())      
            a_enlever = []
            for u,v,data in graphe_complet.edges(data=True) :
                if not isinstance(data["sim"], dict) :
                    if data["sim"] < 0.6 :
                        a_enlever.append((u,v))
                else :
                    if data["sim"]["sim"] < 0.6 :
                        a_enlever.append((u,v))
    #     print(a_enlever)
            for elt in a_enlever :
                graphe_complet.remove_edge(elt[0], elt[1])
            print(graphe_complet.number_of_nodes())
            print(graphe_complet.number_of_edges())
#                 
#           
            if isinstance(num_ARN, list) :    
#                 
                with open("/media/coline/Maxtor/fichier_csv_new_data_seuil_0.6_inter_groupe_complet.csv",'w') as fichier_csv:
                        csvwriter = csv.writer(fichier_csv)
                        csvwriter.writerow(["source", "target", "label"])
                         
                        seuil = 0.6
                        for u,v,data in graphe_complet.edges(data=True) :
                             #if data["sim"] > seuil :
                            if not isinstance(data["sim"], dict) :
                                u_nom = "fichier_%s_%d_taille_4.pdb"%(u[0], u[1])
                                v_nom = "fichier_%s_%d_taille_4.pdb"%(v[0], v[1])
                                if (u_nom, v_nom) in rmsd.keys() :
                                    val_rmsd = rmsd[(u_nom, v_nom)] 
                                else :
                                    val_rmsd = rmsd[(v_nom, u_nom)] 
                                csvwriter.writerow([u,v,(round(data["sim"],2), round(val_rmsd,2))])
                            
                            else :
                                csvwriter.writerow([u,v,round(data["sim"]["sim"],2)])
                            with open("/media/coline/Maxtor/fichier_csv_new_data_seuil_0.6_inter_groupe_noeud_tout_complet.csv",'w') as fichier_csv:
                                csvwriter2 = csv.writer(fichier_csv)
                                csvwriter2.writerow(["id", "label", "type"])
                                     
                                for noeud,data in graphe_complet.nodes(data=True) :
                                    
                                    csvwriter2.writerow([noeud, data["nom"][0] + "_"+str(data["nom"][1]), graphe_complet.nodes[noeud]["type"]])
            else :
                if(num_ARN != "arnt_16s_arnm" and num_ARN != "arnt_16s") :
                    with open("/home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/groupes_%s_homologues.pickle"%num_ARN, 'rb') as fichier_homologues :
                        mon_depickler = pickle.Unpickler(fichier_homologues)
                        groupes_homologues = mon_depickler.load()  
#                    

                     
                    with open("/media/coline/Maxtor/fichier_csv_new_data_seuil_0.6_%s_res_3a.csv"%num_ARN,'w') as fichier_csv:
                        csvwriter = csv.writer(fichier_csv)
                        csvwriter.writerow(["source", "target", "label"])
    #                     
                        seuil = 0.6
                        for u,v,data in graphe_complet.edges(data=True) :
                            if data["sim"] > seuil :
                                u_nom = "fichier_%s_%d_taille_4.pdb"%(u[0], u[1])
                                v_nom = "fichier_%s_%d_taille_4.pdb"%(v[0], v[1])
                                if (u_nom, v_nom) in rmsd.keys() :
                                    if rmsd[(u_nom, v_nom)] != None :
                                        val_rmsd = round(rmsd[(u_nom, v_nom)],2)
                                    else :
                                        val_rmsd = None
                                else :
                                    if rmsd[(v_nom, u_nom)] != None :
                                        val_rmsd = round(rmsd[(v_nom, u_nom)],2)
                                    else :
                                        val_rmsd = None
                                csvwriter.writerow([u[0]+ "_" + str(u[1]),v[0]+"_"+str(v[1]),(round(data["sim"],2), val_rmsd)])
    #                         
                        with open("/media/coline/Maxtor/fichier_csv_new_data_noeuds_seuil_0.6_%s_res_3a.csv"%num_ARN,'w') as fichier_csv:
                            csvwriter2 = csv.writer(fichier_csv)
                            csvwriter2.writerow(["id", "label", "homologues"])
                            print(groupes_homologues)
                            compter_noeuds = 1
                            for noeud,data in graphe_complet.nodes(data=True) :
                                num_homologues = -1
                                compteur = 1
                                for groupe in groupes_homologues :
                                    #print(data["nom"])
                                    if noeud in groupe :
                                        num_homologues = compteur
                                    compteur += 1
                                if num_homologues == -1 :
                                    print("gros rat")
                                csvwriter2.writerow([noeud[0] + "_"+str(noeud[1]), noeud[0] + "_"+str(noeud[1]), num_homologues])
                                compter_noeuds += 1
                else :
                    with open("/media/coline/Maxtor/fichier_csv_new_data_seuil_0.6_%s_res_3a.csv"%num_ARN,'w') as fichier_csv:
                        csvwriter = csv.writer(fichier_csv)
                        csvwriter.writerow(["source", "target", "label"])
    #                     
                        seuil = 0.6
                        for u,v,data in graphe_complet.edges(data=True) :
                            if data["sim"] > seuil :
                                u_nom = "fichier_%s_%d_taille_4.pdb"%(u[0], u[1])
                                v_nom = "fichier_%s_%d_taille_4.pdb"%(v[0], v[1])
                                if (u_nom, v_nom) in rmsd.keys() :
                                    if rmsd[(u_nom, v_nom)] != None :
                                        val_rmsd = round(rmsd[(u_nom, v_nom)],2)
                                    else :
                                        val_rmsd = None
                                else :
                                    if rmsd[(v_nom, u_nom)] != None :
                                        val_rmsd = round(rmsd[(v_nom, u_nom)],2)
                                    else :
                                        val_rmsd = None
                                csvwriter.writerow([u[0]+ "_" + str(u[1]),v[0]+"_"+str(v[1]),(round(data["sim"],2), val_rmsd)])
    #                         
                        with open("/media/coline/Maxtor/fichier_csv_new_data_noeuds_seuil_0.6_%s_res_3a.csv"%num_ARN,'w') as fichier_csv:
                            csvwriter2 = csv.writer(fichier_csv)
                            csvwriter2.writerow(["id", "label"])
                            compter_noeuds = 1
                            for noeud,data in graphe_complet.nodes(data=True) :
                                csvwriter2.writerow([noeud[0] + "_"+str(noeud[1]), noeud[0] + "_"+str(noeud[1])])
                                compter_noeuds += 1
                                
                           
#                                  
#                                num_clusters = -1   
#                                 for groupe in tab_clustering :
#                                     if data["nom"] in groupe :
#                                         num_clusters = compteur
#                                     compteur += 1
#                                 if num_clusters == -1 :
#                                     print("gros rat")
#                                 if data["nom"] in groupe :
#                                     numero_homologues = 0
#                                     for j in range(len(HOMOLOGUES)) :
#                                         
#                                         if data["nom"] in HOMOLOGUES[j] :
#                                             numero_homologues = j+1
#                                             
#                                             
#                                     if data["nom"] in groupe_base :
#                                         dans_groupe = 1
#                                     else :
#                                         dans_groupe = 0
                            
#                   #         with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_25s.pickle", 'rb') as fichier_graphe_sim_sans_doublons :
#             mon_depickler = pickle.Unpickler(fichier_graphe_sim_sans_doublons)
#             graphe_complet = mon_depickler.load() 

''' ancien '''
def toutes_comparaison():
    
    liste_fichiers = []
    for fic in os.listdir(NEW_EXTENSION_PATH_TAILLE) :
        if  "couples_possibles" not in fic and "pickle" in fic and "avec_coord" not in fic : #and len(fic.split("_")) == 6 :
            liste_fichiers.append(fic)
    print(liste_fichiers)
    print(len(liste_fichiers))
    
    with open("Nouvelles_donnees/Resultats/dico_graphe_sim_new_algo.pickle", 'rb') as fichier_lecture :
        mon_depickler = pickle.Unpickler(fichier_lecture)
        dico_graphe = mon_depickler.load()
        
        compteur = 0
        for i in range(len(liste_fichiers)) :
            #print(liste_fichiers[i])
            
            for j in range(i+1, len(liste_fichiers)) :
                #print(liste_fichiers[j])
                
                #if liste_fichiers[i] == "fichier_4V9F_0_48_2_3.pickle" and liste_fichiers[j] == "fichier_3JCS_1_25_16_3.pickle" :
                if compteur < 10000 and (liste_fichiers[i][8:len(liste_fichiers[i])-7], liste_fichiers[j][8:len(liste_fichiers[j])-7]) not in dico_graphe.keys() :  
                
                    with open(NEW_EXTENSION_PATH_TAILLE+liste_fichiers[i], 'rb') as fichier1 :
                        mon_depickler = pickle.Unpickler(fichier1)
                        graphe1 = mon_depickler.load()
                        
                        
                        
                    with open(NEW_EXTENSION_PATH_TAILLE+liste_fichiers[j], 'rb') as fichier2 :
                        mon_depickler = pickle.Unpickler(fichier2)
                        graphe2 = mon_depickler.load()
        #                 tab = chaines_reliees(graphe2)
        #                 print(tab)
        #                 if len(tab) > 0 :
        #                     compteur +=1
                            
            
                        
                         
                        graphe_commun_max, sim_max = comparaison(graphe1, graphe2, "petit rat")   
                         
                        #dico_sim.update({(liste_fichiers[i][8:len(liste_fichiers[i])-7], liste_fichiers[j][8:len(liste_fichiers[j])-7]) : sim_max})
                        print("sim max")
                        print(sim_max)
                        print("graphe commun max")
                        print(graphe_commun_max.nodes.data())
                        print(graphe_commun_max.edges.data())
                        
                        dico_graphe.update({(liste_fichiers[i][8:len(liste_fichiers[i])-7], liste_fichiers[j][8:len(liste_fichiers[j])-7]) : {"sim" : sim_max, "graphe" : graphe_commun_max}})
                        
                        compteur += 1
            #break
    
    with open("Nouvelles_donnees/Resultats/dico_graphe_sim_new_algo.pickle", 'wb') as fichier_graphe :
        mon_pickler = pickle.Pickler(fichier_graphe)
        mon_pickler.dump(dico_graphe)
            
#     with open("Nouvelles_donnees/Resultats/dico_sim_new_algo.pickle", 'wb') as fichier_ecriture :
#         mon_pickler = pickle.Pickler(fichier_ecriture)
#         mon_pickler.dump(dico_sim)
 
'''09/09/2019
effectuer les comparaisons au sein des groupes d'homologues '''  
def comparaison_homologues(num_ARN1, num_ARN2, num1, num2):     
    
    #compteur = 0
    compter_partie = 1    
    with open("groupes_%s_homologues.pickle"%num_ARN1, 'rb') as fichier1_homologues :
        mon_depickler_1 = pickle.Unpickler(fichier1_homologues)
        groupes1_homologues = mon_depickler_1.load()
        
    with open("groupes_%s_homologues.pickle"%num_ARN2, 'rb') as fichier2_homologues :
        mon_depickler_2 = pickle.Unpickler(fichier2_homologues)
        groupes2_homologues = mon_depickler_2.load()
        
        
        compteur1 = 0
        compteur2 = 0
        dico_graphe = {}
        for groupe1 in groupes1_homologues :
            if compteur1 == num1 :
                
                for i in range(len(groupe1)) :
                    #print(liste_fichiers[i])
                    for groupe2 in groupes2_homologues :
                        if compteur2 == num2 :
                            for j in range(len(groupe2)) :
                        #if compteur < 15 :
                        #print(liste_fichiers[j])
                        
                        #if liste_fichiers[i] == "fichier_4V9F_0_48_2_3.pickle" and liste_fichiers[j] == "fichier_3JCS_1_25_16_3.pickle" :
                        #if compteur < 10000 and (liste_fichiers[i][8:len(liste_fichiers[i])-7], liste_fichiers[j][8:len(liste_fichiers[j])-7]) not in dico_graphe.keys() :  
                            #if "fichier_"+str(groupe[i][0])+"_"+str(groupe[i][1])+".pickle" in os.listdir(NEW_EXTENSION_PATH_TAILLE) and "fichier_"+str(groupe[j][0])+"_"+str(groupe[j][1])+".pickle" in os.listdir(NEW_EXTENSION_PATH_TAILLE) :
                                with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+str(groupe1[i][0])+"_"+str(groupe1[i][1])+".pickle", 'rb') as fichier1 :
                                    mon_depickler = pickle.Unpickler(fichier1)
                                    graphe1 = mon_depickler.load()
                                
                                
                                
                                with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+str(groupe2[j][0])+"_"+str(groupe2[j][1])+".pickle", 'rb') as fichier2 :
                                    mon_depickler = pickle.Unpickler(fichier2)
                                    graphe2 = mon_depickler.load()
                    #                 tab = chaines_reliees(graphe2)
                    #                 print(tab)
                    #                 if len(tab) > 0 :
                    #                     compteur +=1
                                        
                        
                                    
                                     
                                    graphe_commun_max, sim_max = comparaison(graphe1, graphe2, "petit rat")   
                                     
                                    #dico_sim.update({(liste_fichiers[i][8:len(liste_fichiers[i])-7], liste_fichiers[j][8:len(liste_fichiers[j])-7]) : sim_max})
                                    print("sim max")
                                    print(sim_max)
                                    print("graphe commun max")
                                    print(graphe_commun_max.nodes.data())
                                    print(graphe_commun_max.edges.data())
                                    
                                    dico_graphe.update({(str(groupe1[i][0])+"_"+str(groupe1[i][1]), str(groupe2[j][0])+"_"+str(groupe2[j][1])) : {"sim" : sim_max, "graphe" : graphe_commun_max}})
                                    
                                    if len(dico_graphe) > 8000 :
                                            with open("Nouvelles_donnees/Resultats/dico_sim_new_algo_%s_%s_homologues_groupe_%s_groupe_%s_part_%s.pickle"%(num_ARN1, num_ARN2, num1, num2, compter_partie), 'wb') as fichier_sim :
                                                    mon_pickler = pickle.Pickler(fichier_sim)
                                                    mon_pickler.dump(dico_graphe)
                                            dico_graphe.clear()
                                            compter_partie += 1 
                        compteur2 += 1
            compteur1 += 1
                                
                    
                                    
            
#                 with open("Nouvelles_donnees/Resultats/dico_graphe_sim_new_algo_%s_homologues_groupe_%s.pickle"%(num_ARN, compteur), 'wb') as fichier_graphe :
#                     mon_pickler = pickle.Pickler(fichier_graphe)
#                     mon_pickler.dump(dico_graphe)
        with open("Nouvelles_donnees/Resultats/dico_sim_new_algo_%s_%s_homologues_groupe_%s_groupe_%s_part_%s.pickle"%(num_ARN1,num_ARN2, num1, num2, compter_partie), 'wb') as fichier_sim :
            mon_pickler = pickle.Pickler(fichier_sim)
            mon_pickler.dump(dico_graphe)
                #compteur += 1
            
#     with open("Nouvelles_donnees/Resultats/dico_sim_new_algo.pickle", 'wb') as fichier_ecriture :
#         mon_pickler = pickle.Pickler(fichier_ecriture)
#         mon_pickler.dump(dico_si

'''09/09/19
chercher les vrais identiques : d'un point de vue extension aussi '''
def obs_resultats(num_ARN):
#     with open("Nouvelles_donnees/Resultats/dico_graphe_sim_new_algo_23s_homologues_groupe_%s.pickle"%5, 'rb') as fichier_graphe :
#         mon_depickler = pickle.Unpickler(fichier_graphe)
#         dico_graphe_sim = mon_depickler.load()
#         
#         print(len(dico_graphe_sim))
#         #print(dico_graphe_sim.keys())
#         
#         compter_idem = 0
#         tab_sim = []
#         for cle in dico_graphe_sim.keys() :
#             tab_sim.append(dico_graphe_sim[cle]["sim"]) 
#             if dico_graphe_sim[cle]["sim"] ==  1.0 :
# #                 print(cle)
# #                 print(dico_graphe_sim[cle])
#                 #print(dico_graphe_sim[cle]["graphe"].edges.data())
#                 compter_idem += 1
#         plt.plot(tab_sim)
#         plt.show()
        #print(compter_idem)
        dico_graphe_sim_complet = {}
        compteur = 0
        while "dico_graphe_sim_new_algo_%s_homologues_groupe_%s.pickle"%(num_ARN, compteur) in os.listdir("/media/coline/Maxtor/Resultats/%s"%num_ARN) :
            with open("/media/coline/Maxtor/Resultats/%s/dico_graphe_sim_new_algo_%s_homologues_groupe_%s.pickle"%(num_ARN, num_ARN, compteur), 'rb') as fichier_graphe :
                mon_depickler = pickle.Unpickler(fichier_graphe)
                dico_graphe_sim = mon_depickler.load()
                
                dico_graphe_sim_complet.update(dico_graphe_sim)
                print(compteur)
                compteur += 1
            
        
        compteur_nb_groupes_vraiment_identiques = 0
        with open("groupes_%s_homologues_sequences.pickle"%num_ARN, 'rb') as fichier_homologues :
            mon_depickler_1 = pickle.Unpickler(fichier_homologues)
            groupes_homologues = mon_depickler_1.load()
            print(len(groupes_homologues))
            for groupe in groupes_homologues :
                print(len(groupe))
                print(groupe)
            #print(groupes_homologues)
             
            with open("groupes_%s_identiques_sequences.pickle"%num_ARN, 'rb') as fichier_identiques :
                mon_depickler_2 = pickle.Unpickler(fichier_identiques)
                groupes_identiques = mon_depickler_2.load()
                
                print(len(groupes_identiques))
                print(groupes_identiques)
                 
                with open("all_aminor.pickle", 'rb') as fichier_all_aminor :
                    mon_depickler = pickle.Unpickler(fichier_all_aminor)
                    all_aminor = mon_depickler.load()

                    compter_vraiment_idem = 0
                    compteur = 0
                    groupes_a_rassembler = []
                    for groupes in groupes_identiques :
                        
                        
                        #if compteur > 5 :
#                             if "dico_graphe_sim_new_algo_%s_homologues_groupe_%s.pickle"%(num_ARN, compteur) in os.listdir("/media/coline/Maxtor/Resultats/%s"%num_ARN) :
#                                 with open("/media/coline/Maxtor/Resultats/%s/dico_graphe_sim_new_algo_%s_homologues_groupe_%s.pickle"%(num_ARN, num_ARN, compteur), 'rb') as fichier_graphe :
#                                     mon_depickler = pickle.Unpickler(fichier_graphe)
#                                     dico_graphe_sim = mon_depickler.load()
                                    
                                
                                    print(len(groupes))
                                    for groupe in groupes :
                                        groupes_vraiment_identiques = []
                                        print(len(groupe))
                                        print(groupe)
            #                             for i in range(len(groupe)) :
            #                                 for j in range(i+1, len(groupe)) :
            #                                     if (str(groupe[i][0])+"_"+str(groupe[i][1]), str(groupe[j][0])+"_"+str(groupe[j][1])) in dico_graphe_sim.keys() :
            #                                         print(dico_graphe_sim[(str(groupe[i][0])+"_"+str(groupe[i][1]), str(groupe[j][0])+"_"+str(groupe[j][1]))])
            #                                     else  :
            #                                         print(dico_graphe_sim[(str(groupe[j][0])+"_"+str(groupe[j][1]), str(groupe[i][0])+"_"+str(groupe[i][1]))])
                                        for i in range(len(groupe)) :
                                            #if groupe[i] not in deux_chaines :
                                                for j in range(i+1, len(groupe)) :
                                                    #print((groupe[i], groupe[j]))
                                                    if (str(groupe[i][0])+"_"+str(groupe[i][1]), str(groupe[j][0])+"_"+str(groupe[j][1])) in dico_graphe_sim_complet.keys() :
                                                        if dico_graphe_sim_complet[(str(groupe[i][0])+"_"+str(groupe[i][1]), str(groupe[j][0])+"_"+str(groupe[j][1]))]["sim"] == 1.0 :
                                                            #print((groupe[i], groupe[j]))
                                                        #else :
                                                            ok = False
                                                            compter_vraiment_idem += 1
                                                            for elt in groupes_vraiment_identiques :
                                                                if groupe[j] in elt and groupe[i] not in elt :
                                                                    elt.append(groupe[i])
                                                                    ok = True
                                                                elif groupe[i] in elt and groupe[j] not in elt :
                                                                    elt.append(groupe[j])
                                                                    ok = True
                                                                elif groupe[i] in elt and groupe[j] in elt :
                                                                    ok = True
                                                                    
                                                            if not ok :
                                                                groupes_vraiment_identiques.append([groupe[i], groupe[j]])
                                                                
                                                        else : 
                                                            print(groupe[i])
                                                            print(groupe[j])
                                                            print(dico_graphe_sim_complet[(str(groupe[i][0])+"_"+str(groupe[i][1]), str(groupe[j][0])+"_"+str(groupe[j][1]))])
                                                            #print(dico_graphe_sim[(str(groupe[i][0])+"_"+str(groupe[i][1]), str(groupe[j][0])+"_"+str(groupe[j][1]))]["graphe"].edges.data())
                                                            
                                                            with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+str(groupe[i][0])+"_"+str(groupe[i][1])+".pickle", 'rb') as fichier1 :
                                                                mon_depickler_1 = pickle.Unpickler(fichier1)
                                                                graphe1 = mon_depickler_1.load()
                                                                
                                                            with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+str(groupe[j][0])+"_"+str(groupe[j][1])+".pickle", 'rb') as fichier2 :
                                                                mon_depickler_2 = pickle.Unpickler(fichier2)
                                                                graphe2 = mon_depickler_2.load()
                                                                
                                                            print(graphe1.nodes.data())
                                                            print(graphe2.nodes.data())
                                                            
                                                            for cle in all_aminor.keys() :
                                                                if cle == str(groupe[i][0]) :
                                                                    c = 1
                                                                    for elt in all_aminor[cle] :
                                                                        if c == groupe[i][1] :
                                                                            print(elt.nodes.data())
                                                                        c += 1
                                                                        
                                                                if cle == str(groupe[j][0]) :
                                                                    c = 1
                                                                    for elt in all_aminor[cle] :
                                                                        if c == groupe[j][1] :
                                                                            print(elt.nodes.data())
                                                                        c += 1
                                                                        
                                                                        
                                                            
                                                    elif  (str(groupe[j][0])+"_"+str(groupe[j][1]), str(groupe[i][0])+"_"+str(groupe[i][1])) in dico_graphe_sim_complet.keys() :
                                                        if dico_graphe_sim_complet[(str(groupe[j][0])+"_"+str(groupe[j][1]), str(groupe[i][0])+"_"+str(groupe[i][1]))]["sim"] == 1.0 :
                #                                             print((groupe[j], groupe[i]))
                #                                         else :
                                                            ok = False
                                                            compter_vraiment_idem += 1
                                                            for elt in groupes_vraiment_identiques :
                                                                if groupe[j] in elt and groupe[i] not in elt :
                                                                    elt.append(groupe[i])
                                                                    ok = True
                                                                elif groupe[i] in elt and groupe[j] not in elt :
                                                                    elt.append(groupe[j])
                                                                    ok = True
                                                                elif groupe[i] in elt and groupe[j] in elt :
                                                                    ok = True
                                                                    
                                                            if not ok :
                                                                groupes_vraiment_identiques.append([groupe[i], groupe[j]])
                                                                
                                                        
                                                        else : 
                                                            print(groupe[i])
                                                            print(groupe[j])
                                                            print(dico_graphe_sim_complet[(str(groupe[j][0])+"_"+str(groupe[j][1]), str(groupe[i][0])+"_"+str(groupe[i][1]))])
                                                            #print(dico_graphe_sim[(str(groupe[j][0])+"_"+str(groupe[j][1]), str(groupe[i][0])+"_"+str(groupe[i][1]))]["graphe"].edges.data())
                                                    else :
                                                        print("bizarre")
                                                            
                                        print(groupes_vraiment_identiques)
    #                                     print("nb de groupes vraiment identiques")
    #                                     print(len(groupes_vraiment_identiques))
    #                                     with open(NEW_EXTENSION_PATH_TAILLE+"Resultats/groupes_vraiment_identiques_%s"%num_ARN, 'wb') as fichier_vraiment_id :
    #                                         mon_pickler = pickle.Pickler(fichier_vraiment_id)
    #                                         mon_pickler.dump(groupes_vraiment_identiques)
                                        compteur_nb_groupes_vraiment_identiques += len(groupes_vraiment_identiques)
    #                                     print("compteur etape 1 ")
    #                                     print(compteur_nb_groupes_vraiment_identiques)
                                        
    #                                     with open(NEW_EXTENSION_PATH_TAILLE+"groupe_%s.pickle"%num_ARN, 'rb') as fichier_sortie :
    #                                         mon_depickler = pickle.Unpickler(fichier_sortie)
    #                                         liste = mon_depickler.load()
                                        #compteur_seuls = 0
    #                                         print("liste")
    #                                         print(len(liste))
    #                                     for elt in groupe :
    #                                             y_est = False
    #                                             for groupes in groupes_vraiment_identiques :
    #                                                 for elt_id in groupes :
    #                                                     if elt == elt_id :
    #                                                         y_est = True
    #                                             if not y_est :
    #                                                 compteur_seuls += 1
                                        
                                        for g in groupes_vraiment_identiques :
                                            groupes_a_rassembler.append(list(g))
                                        
                                        
    #                                     print("tout seul")
    #                                     print(compteur_seuls)
                                        
                                        #compteur_nb_groupes_vraiment_identiques += compteur_seuls
                                        print("compteur etape 2 ")
                                        print(compteur_nb_groupes_vraiment_identiques)
                                        for i in range(len(groupes_vraiment_identiques)) :
                                            for j in range(i+1, len(groupes_vraiment_identiques)) :
                                                for elt1 in groupes_vraiment_identiques[i] :
                                                    for elt2 in groupes_vraiment_identiques[j] : 
                                                        if elt1 == elt2 :
                                                            print("gros rat")
                                                            print(elt1)
                                    
                                 
                                    compteur += 1
        with open(NEW_EXTENSION_PATH_TAILLE+"groupe_%s.pickle"%num_ARN, 'rb') as fichier_sortie :
            mon_depickler = pickle.Unpickler(fichier_sortie)
            liste = mon_depickler.load()
            print("liste")
            print(len(liste))
            print(liste)
            #print(groupes_vraiment_identiques)
            
            for groupes in groupes_a_rassembler : 
                print(len(groupes))
        
            compteur_seuls = 0
            for elt in liste :
                y_est = False
                for groupes in groupes_a_rassembler :
                    for elt_id in groupes :
                        if elt == elt_id :
                            y_est = True
                if not y_est :
                    compteur_seuls += 1
        print("nb de groupes d'identiques")
        print(compteur_nb_groupes_vraiment_identiques)
        print("compteur_seuls")
        print(compteur_seuls)
        compteur_nb_groupes_vraiment_identiques += compteur_seuls
        print(compter_vraiment_idem)
        print("nb vraiment identiques")
        print(compteur_nb_groupes_vraiment_identiques)
        with open(NEW_EXTENSION_PATH_TAILLE+"Resultats/groupes_vraiment_identiques_sequences_%s.pickle"%num_ARN, 'wb') as fichier_vraiment_id :
            mon_pickler = pickle.Pickler(fichier_vraiment_id)
            mon_pickler.dump(groupes_a_rassembler)
            

#         plt.plot(tab_sim)
#         plt.show()

'''11/09/19'''
def ecarter_doublons(num_ARN):
    with open("Nouvelles_donnees/script_deplacer_doublons_%s.sh"%num_ARN, 'w') as fichier :
        
        with open(NEW_EXTENSION_PATH_TAILLE+"Resultats/groupes_vraiment_identiques_%s.pickle"%num_ARN, 'rb') as fichier_vraiment_id :
            mon_depickler = pickle.Unpickler(fichier_vraiment_id)
            groupes_a_rassembler = mon_depickler.load()
            
            with open(NEW_EXTENSION_PATH_TAILLE+"groupe_%s.pickle"%num_ARN, 'rb') as fichier_sortie :
                mon_depickler = pickle.Unpickler(fichier_sortie)
                liste = mon_depickler.load()
                
            
                compteur_seuls = 0
                for elt in liste :
                    y_est = False
                    for groupes in groupes_a_rassembler :
                        for elt_id in groupes :
                            if elt == elt_id :
                                y_est = True
                    if not y_est :
                        compteur_seuls += 1
            
            #print(groupes_a_rassembler)
            print(len(groupes_a_rassembler))
            print(compteur_seuls)
            print(compteur_seuls+len(groupes_a_rassembler))
            
            for groupe in groupes_a_rassembler :
                for i in range(1, len(groupe)) :
                    fichier.write("mv fichier_%s_%s.pickle Doublons/fichier_%s_%s.pickle\n"%(groupe[i][0],groupe[i][1],groupe[i][0],groupe[i][1]))
                    fichier.write("mv fichier_%s_%s.png Doublons/fichier_%s_%s.png\n"%(groupe[i][0],groupe[i][1],groupe[i][0],groupe[i][1]))

'''13/09/19 '''
def distribution_homologues(num_ARN, c):
    with open("groupes_%s_homologues.pickle"%num_ARN, 'rb') as fichier_homologues :
            mon_depickler_1 = pickle.Unpickler(fichier_homologues)
            groupes_homologues = mon_depickler_1.load()
            
            compteur = 0
            for groupe in groupes_homologues :
                if compteur == c : 
                    with open("Nouvelles_donnees/Resultats/dico_sim_new_algo_%s_homologues_groupe_%s.pickle"%(num_ARN, compteur), 'rb') as fichier_graphe :
                        mon_depickler = pickle.Unpickler(fichier_graphe)
                        dico_graphe_sim = mon_depickler.load()
                        
                        plt.plot(dico_graphe_sim.values())
                compteur += 1
            plt.title("Distribution des similarites au sein d'un groupe d'homologues d'ARNr 23S")
            ax = plt.gca()
            ax.set_ylabel("Valeur de similarite")
            plt.show()
       
'''11/09/19'''
def groupe_entre_deux_chaines(num_ARN):
    with open("all_aminor.pickle", 'rb') as fichier_all_aminor :
        mon_depickler = pickle.Unpickler(fichier_all_aminor)
        all_aminor = mon_depickler.load()
        with open("groupes_%s_identiques.pickle"%num_ARN, 'rb') as fichier_identiques :
            mon_depickler_1 = pickle.Unpickler(fichier_identiques)
            groupes_identiques = mon_depickler_1.load() 
            
            for groupes in groupes_identiques :
                for groupe in groupes :
                    nb_entre_deux_chaines = 0
                    for elt in groupe :
                        for cle in all_aminor.keys() :
                            if cle == elt[0] :
                                c = 1
                                for graphe in all_aminor[cle] :
                                    if c == elt[1] :
                                        liste_nums_chaines_par_graphe = []
                                        for noeud in graphe.nodes() :
                                            if noeud[0] not in liste_nums_chaines_par_graphe :
                                                liste_nums_chaines_par_graphe.append(noeud[0])
                                        if len(liste_nums_chaines_par_graphe) == 3 :
                                            nb_entre_deux_chaines +=1
                                    c += 1
                    if nb_entre_deux_chaines != 0 : 
                        print(groupe)
#                         for e in groupe :
#                          
#                             for cle in all_aminor.keys() :
#                                 if cle == e[0] :
#                                     c = 1
#                                     for graphe in all_aminor[cle] :
#                                         if c == e[1] :
#                                             print(cle)
#                                             print(graphe.nodes.data())
#                                         c += 1
                    if nb_entre_deux_chaines != len(groupe) and nb_entre_deux_chaines != 0 :
                        print(nb_entre_deux_chaines)
                        print(len(groupe))
                        print(groupe)
                        print("bizarre")
                        
                        for e in groupe :
                        
                            for cle in all_aminor.keys() :
                                if cle == e[0] :
                                    c = 1
                                    for graphe in all_aminor[cle] :
                                        if c == e[1] :
                                            print(cle)
                                            print(graphe.nodes.data())
                                        c += 1

