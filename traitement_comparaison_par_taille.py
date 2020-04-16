'''
Created on 10 avr. 2019

@author: coline
'''
import csv
import networkx as nx
import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
from recup_data.regroupement_graphes import creation_dico_comp
from recup_data.calcul_sim import stockage_sim, distribution_similarite
from recup_data.clustering import construction_graphe_complet_pondere,\
    premier_clustering
from recup_data.constantes import EXTENSION_PATH, EXTENSION_TOUTES_ARETES,\
    GROUPES_TOUTES_ARETES_MAX_4_10_07, GROUPE_ARICH_DE_07, GROUPE_DENSE_DE_07
from recup_data.etude_composantes_connexes import draw_composantes,\
    stats_composantes_connexes, compter_composantes_differentes,\
    comparaison_cluster, compter_composantes_differentes_par_groupe,\
    draw_groupes_par_type, non_homologues_dans_cluster,\
    comparaison_cluster_identiques_plus_de_2, afficher_groupes,\
    draw_grandes_composantes, draw_groupes_par_type, clustering_composante,\
    retrouver_types, retrouver_lien_chaines_1_3, calcul_groupe_etendu_seuil_sup,\
    draw_groupe_gephi_unique, draw_une_composante, draw_une_grande_composante
from recup_data.draw_isomorphism import draw_isomorphisme

def obtention(taille_ext):
    creation_dico_comp(taille_ext) 
    stockage_sim("toutes_aretes_coeff_all1", "toutes_aretes_coeff_all1", "extensions", taille_ext)
      #stockage_sim("toutes_aretes_coeff_all1", "toutes_aretes_coeff_all1_par_k", "extensions", taille_ext)
    graphe_complet = construction_graphe_complet_pondere(taille_ext, "toutes_aretes_coeff_all1", "toutes_aretes_coeff_all1")
    print(graphe_complet.number_of_nodes())
    print(graphe_complet.number_of_edges())
     
#     with open(EXTENSION_PATH%taille_ext+"sim_extensions_toutes_aretes_coeff_all1_taille_%s.pickle"%taille_ext, 'rb') as fichier_sim_new :
#         mon_depickler_new = pickle.Unpickler(fichier_sim_new)
#         dico_sim_new = mon_depickler_new.load()
#      
#         with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/result_graphes_extension_%s/sim_extensions_toutes_aretes_coeff_all1_taille_%s.pickle"%(taille_ext, taille_ext), 'rb') as fichier_sim :
#             mon_depickler = pickle.Unpickler(fichier_sim)
#             dico_sim = mon_depickler.load()
#              
#             print(len(dico_sim))
#             print(dico_sim.keys())
#              
#             for noeud1, data1 in graphe_complet.nodes(data=True) :
#                 for noeud2, data2 in graphe_complet.nodes(data=True) :
#                     if noeud1 != noeud2 :
#                         if (noeud1,noeud2) not in graphe_complet.edges() :
#                             if (data1["nom"], data2["nom"]) in dico_sim.keys() :
#                                 graphe_complet.add_edge(noeud1, noeud2, poids = dico_sim[(data1["nom"], data2["nom"])])
#                                 dico_sim_new.update({(data1["nom"], data2["nom"]) : dico_sim[(data1["nom"], data2["nom"])] })
#                             else :
#                                 graphe_complet.add_edge(noeud1, noeud2, poids = dico_sim[(data2["nom"], data1["nom"])])
#                                 dico_sim_new.update({(data2["nom"], data1["nom"]) : dico_sim[(data2["nom"], data1["nom"])] })
#  
#     print(len(dico_sim_new))
#     print(graphe_complet.number_of_edges())
#     
#         
#     with open(EXTENSION_PATH%taille_ext+"graphe_complet_pondere_sim_%s_taille_%s.pickle"%("toutes_aretes_coeff_all1",taille_ext), 'wb') as fichier_graphe_complet :
#         mon_pickler_complet = pickle.Pickler(fichier_graphe_complet)
#         mon_pickler_complet.dump(graphe_complet)      
#          
#     with open(EXTENSION_PATH%taille_ext+"sim_extensions_toutes_aretes_coeff_all1_taille_%s.pickle"%taille_ext, 'wb') as fichier_sim_new_new :
#         mon_pickler_new = pickle.Pickler(fichier_sim_new_new)
#         mon_pickler_new.dump(dico_sim_new)     
       
    with open(EXTENSION_PATH%taille_ext+"fichier_petits_groupes_toutes_aretes_taille_%s.txt"%taille_ext, 'w') as fichier_ecriture :
        i = 0.1
        while i <= 1.0 :
            print(i)
            fichier_ecriture.write("Valeur: "+str(i)+"\n")
            premier_clustering(graphe_complet, i, fichier_ecriture,"extensions", "toutes_aretes_coeff_all1", taille_ext, taille_ext)
            i = i+0.1
    draw_composantes("toutes_aretes_coeff_all1", "extensions", taille_ext, taille_ext, 1, type_comp="extensions_toutes_aretes")
    distribution_similarite("extensions", "toutes_aretes_coeff_all1", taille_ext)
    stats_composantes_connexes(taille_ext, 'extensions', 'toutes_aretes_coeff_all1', 1)
    
    
    #distribution_similarite("extensions", "toutes_aretes_coeff_all1_par_k", taille_ext)
    #draw_isomorphisme(taille_ext)
    
#     stockage_sim("toutes_aretes_coeff_all1", "toutes_aretes_coeff_all1_chaines_1_3", "extensions", taille_ext)
#     distribution_similarite("extensions", "toutes_aretes_coeff_all1_chaines_1_3", taille_ext)
#     graphe_complet = construction_graphe_complet_pondere(taille_ext, "toutes_aretes_coeff_all1", "toutes_aretes_coeff_all1_chaines_1_3")
#     with open(EXTENSION_PATH%taille_ext+"fichier_petits_groupes_toutes_aretes_chaines_1_3_taille_%s.txt"%taille_ext, 'w') as fichier_ecriture :
#         i = 0.1
#         while i <= 1.0 :
#             print(i)
#             fichier_ecriture.write("Valeur: "+str(i)+"\n")
#             premier_clustering(graphe_complet, i, fichier_ecriture,"extensions", "toutes_aretes_coeff_all1_chaines_1_3", taille_ext)
#             i = i+0.1
#     stats_composantes_connexes(taille_ext, 'extensions', 'toutes_aretes_coeff_all1_chaines_1_3', 1)
#     draw_composantes("toutes_aretes_coeff_all1_chaines_1_3", "extensions", taille_ext, 1, type_comp="extensions_toutes_aretes")
    
def observation(taille_ext):
    stats_composantes_connexes(taille_ext, 'extensions', 'toutes_aretes_coeff_all1', 1)
    compter_composantes_differentes(taille_ext)
    #draw_groupes_par_type_motif(taille_ext, ['1FJG_A_271_1', '4V9F_0_30_23', '4V9F_0_48_16', '1FJG_A_109_6'], "dans_10")
    #non_homologues_dans_cluster(taille_ext, 'extensions', 'toutes_aretes_coeff_all1')

def comparaison(taille_ext_1, taille_ext_2):
    groupe_1, groupe_2 = comparaison_cluster(taille_ext_1, taille_ext_2)
    return groupe_1, groupe_2

def generation_distrib_csv(val_min, val_max, depart, typ_graphe, typ_sim):
    with open(EXTENSION_TOUTES_ARETES+"csv_distrib_sim_%s_%s_%s.csv"%(typ_sim, val_min, val_max), 'w', newline='') as fichier_csv :
        csvwriter = csv.writer(fichier_csv)
        for i in range(val_min,val_max+1) :
            dico_sim = stockage_sim(typ_graphe, typ_sim, depart, i)
    
            nombre_intervalles = [0]*20
            for cle in dico_sim.keys() :
                if dico_sim[cle] > 4.25 :
                    nombre_intervalles[19] += 1
                else :
                    nombre_intervalles[int(dico_sim[cle]/0.225)] += 1
            
            proportion_intervalles = [0]*20
            for i in range(20) :
                proportion_intervalles[i] = nombre_intervalles[i]/4005

            
            csvwriter.writerow(proportion_intervalles) 

            
def traitement_recherche_kmax(val_min, val_max, depart, typ_graphe, typ_sim):
    dico_crible = {}
    maxi_tot = -1.0
    for i in range(val_min,val_max+1) :
        with open(EXTENSION_PATH%i+"sim_%s_%s_taille_%s.pickle"%(depart, typ_sim, i), 'rb') as fichier_sim_new_new :
            mon_depickler_new = pickle.Unpickler(fichier_sim_new_new)
            dico_sim = mon_depickler_new.load()   
            #dico_sim = stockage_sim(typ_graphe, typ_sim, depart, i)
            with open(EXTENSION_PATH%i+"graphe_complet_pondere_sim_%s_taille_%s.pickle"%(typ_sim,i), 'rb') as fichier_graphe_complet :
                mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
                graphe_complet = mon_depickler_complet.load() 
        #graphe_complet = construction_graphe_complet_pondere(i, typ_graphe, typ_sim)
            
                for cle in dico_sim.keys() :
                    if dico_sim[cle] > 1 :
                        print(i)
                        print(dico_sim[cle])
#                     if i == 10 and cle[0] in ['5FDU_1A_272_1', '3JCS_1_25_46'] and cle[1] in ['5FDU_1A_272_1', '3JCS_1_25_46'] :
#                         print(dico_sim[cle])
#                         exit()
                    if cle not in dico_crible.keys() and  (cle[1], cle[0]) not in dico_crible.keys() :
                        dico_crible.update({cle : [dico_sim[cle]]})
                    elif (cle[1], cle[0]) in dico_crible.keys() :
                        dico_crible[(cle[1], cle[0])].append(dico_sim[cle])
                    else :
                        dico_crible[cle].append(dico_sim[cle])
                    if dico_sim[cle] > maxi_tot :
                        maxi_tot = dico_sim[cle]
    print("ramou")
    print(maxi_tot)
    #print(dico_crible)

    dico_crible_max_avec_k = {}
    dico_crible_max = {}
    for cle in dico_crible.keys() :
        maxi = -1.0
        pos = -1
        compteur = val_min
        for elt in dico_crible[cle] :
            if elt > maxi :
                maxi = elt
                pos = compteur
            compteur += 1 
        if cle not in dico_crible_max.keys() :
            dico_crible_max_avec_k.update({cle : (maxi, pos)})
            dico_crible_max.update({cle : maxi})
        else :
            print("bizarre")
    #print(dico_crible_max)
    #print(dico_crible)
     
    nombre_intervalles = [0]*20
    nombre_k = [0]*val_max
    for cle in dico_crible_max_avec_k.keys() :
        #print(dico_crible[cle])
        #print(dico_crible_max_avec_k[cle])
        nombre_k[dico_crible_max_avec_k[cle][1]-val_min] += 1    
        if dico_crible_max_avec_k[cle][0] > maxi_tot - (maxi_tot/20) :
            nombre_intervalles[19] += 1
        else :
            nombre_intervalles[int(dico_crible_max_avec_k[cle][0]/(maxi_tot/20))] += 1
         
    proportion_intervalles = [0]*20
    for i in range(20) :
        proportion_intervalles[i] = nombre_intervalles[i]/4005
         
    with open(EXTENSION_TOUTES_ARETES+"csv_distrib_sim_%s_%d_%d.csv"%(typ_sim, val_min, val_max), 'w', newline='') as fichier_csv :
        csvwriter = csv.writer(fichier_csv)
        csvwriter.writerow(proportion_intervalles)
    print(len(dico_crible))
    print(len(dico_crible_max_avec_k))
    print(nombre_k)   
     
    graphe_complet_max = graphe_complet.copy()
    for u,v in graphe_complet.edges() :
        if (graphe_complet.nodes[u]["nom"], graphe_complet.nodes[v]["nom"]) in dico_crible_max.keys() :
            graphe_complet_max.edges[u,v]["poids"] = dico_crible_max[(graphe_complet.nodes[u]["nom"], graphe_complet.nodes[v]["nom"])]
        else :
            graphe_complet_max.edges[u,v]["poids"] = dico_crible_max[(graphe_complet.nodes[v]["nom"], graphe_complet.nodes[u]["nom"])]
                 
    with open(EXTENSION_PATH%"taille_max"+"sim_%s_%s_max_taille_taille_max.pickle"%(depart, typ_sim), 'wb') as fichier_sim :
        mon_depickler = pickle.Pickler(fichier_sim)
        mon_depickler.dump(dico_crible_max)
        
    with open(EXTENSION_PATH%"taille_max"+"sim_%s_%s_max_taille_taille_max_avec_val_k.pickle"%(depart, typ_sim), 'wb') as fichier_sim :
        mon_depickler = pickle.Pickler(fichier_sim)
        mon_depickler.dump(dico_crible_max_avec_k)
         
        with open(EXTENSION_TOUTES_ARETES+"fichier_petits_groupes_%s_taille_max.txt"%(typ_sim), 'w') as fichier_ecriture :
            i = 0.1*maxi_tot
            while i <= 1.0*maxi_tot+0.000000005 :
                print(i)
                fichier_ecriture.write("Valeur: "+str(i)+"\n")
                premier_clustering(graphe_complet, i, fichier_ecriture, depart, typ_sim+"_max", "taille_max", "taille_max")
                i = i+(0.1*maxi_tot)
                     
    with open(EXTENSION_PATH%"taille_max"+"graphe_complet_pondere_sim_%s_max_taille_taille_max.pickle"%(typ_sim), 'wb') as fichier_graphe :
        mon_depickler = pickle.Pickler(fichier_graphe)
        mon_depickler.dump(graphe_complet_max)            
     
        draw_composantes(typ_sim+"_max", depart, "taille_max", "taille_max", maxi_tot, type_comp="extensions_toutes_aretes")
        stats_composantes_connexes("taille_max", depart, typ_sim+"_max", maxi_tot)
    print(len(dico_crible))
    print(len(dico_crible_max_avec_k))
    print(nombre_k) 
    print(maxi_tot)

def comparaison_sein_taille(taille_ext, depart, typ1, typ2, chaine):
    with open(EXTENSION_PATH%taille_ext+"sim_%s_%s_taille_%s.pickle"%(depart,typ1,taille_ext), 'rb') as fichier_sim_1 :
        mon_depickler_1 = pickle.Unpickler(fichier_sim_1)
        dico_sim_1 = mon_depickler_1.load()   
        
        with open(EXTENSION_PATH%taille_ext+"sim_%s_%s_taille_%s.pickle"%(depart,typ2,taille_ext), 'rb') as fichier_sim_2 :
            mon_depickler_2 = pickle.Unpickler(fichier_sim_2)
            dico_sim_2 = mon_depickler_2.load() 
            
            tab_diff_ecart_type = [[],[],[],[],[],[],[],[],[],[],[]]
            print(tab_diff_ecart_type)
            tab_difference = [0]*11
            tab_taille = [0]*11
            for cle in dico_sim_1.keys() :
#                 print(cle)
#                 print(dico_sim_1[cle])
                tab_diff_ecart_type[int(dico_sim_1[cle]*10)].append(dico_sim_1[cle] - dico_sim_2[cle])
                #print(tab_diff_ecart_type)
                tab_difference[int(dico_sim_1[cle]*10)] += dico_sim_1[cle] - dico_sim_2[cle] 
                tab_taille[int(dico_sim_1[cle]*10)] += 1
            
            tab_moyenne_difference = [0]*11
            compteur = 0
            for intervalle in tab_difference : 
                print(intervalle)
                if intervalle != 0 :
                    tab_moyenne_difference[compteur] = intervalle/tab_taille[compteur]
                compteur += 1
            print(tab_moyenne_difference)
            x = [0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05]
            y = [0.0]*11
            plt.plot(x,y, color="black")
            plt.scatter(x, tab_moyenne_difference)
            axes = plt.gca()
            axes.xaxis.set_ticks(np.arange(0,1.1,0.1))
            plt.xlabel("Intervalle de valeurs de similarite pour les chaines %d et %d"%(chaine[0], chaine[1]))
            plt.ylabel("Moyenne des différences")
            
            plt.title("Répartition de la moyenne des différences par tranche entre les valeurs \n de similarité pour les chaînes %d et %d par rapport à la similarité globale"%(chaine[0], chaine[1]))
            plt.savefig(EXTENSION_PATH%taille_ext+"moyenne_diff_chaines_%d_%d_avec_global.png"%(chaine[0], chaine[1]))
              
            
            tab_ecart_type = [0]*11
            compteur = 0
            for elt in tab_diff_ecart_type :
                print(elt)
                if len(elt) > 0 :
                    tab_ecart_type[compteur] = np.std(np.array(elt))
                compteur += 1
            print(tab_ecart_type)
#              
#             for j in range(len(tab_moyenne_difference)) :
#                 plt.annotate(str(round(tab_ecart_type[j],2)), xy = (j/10+0.05, tab_moyenne_difference[j]), xytext = (j/10, tab_moyenne_difference[j]-0.02))
# #                
        
            plt.errorbar([0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05], tab_moyenne_difference, yerr = tab_ecart_type,fmt = 'none', capsize = 5, ecolor = 'red', elinewidth = 1, capthick = 4) 
            axes.yaxis.set_ticks(np.arange(-0.4,0.5,0.1))
            plt.savefig(EXTENSION_PATH%taille_ext+"moyenne_diff_chaines_%d_%d_avec_global_avec_ecart_type_barre_erreur.png"%(chaine[0], chaine[1]))
            plt.show()
#             new_val_densite = list(densite_par_taille[tab_taille[j]-1])[0] + tab_densite[j]
#                         new_val_nombre = list(densite_par_taille[tab_taille[j]-1])[1] + 1
#                         
#                         densite_par_taille[tab_taille[j]-1] = (new_val_densite, new_val_nombre)
#             
        #print(tab_difference)
        
def afficher_sim_k_max(cle_1, cle_2, depart, typ_sim, rep):
    with open(EXTENSION_PATH%rep+"sim_%s_%s_taille_taille_max_avec_val_k.pickle"%(depart, typ_sim), 'rb') as fichier_sim :
        mon_depickler = pickle.Unpickler(fichier_sim)
        dico_crible_max_avec_k = mon_depickler.load()
        
        print(cle_1)
        print(cle_2)
        if (cle_1, cle_2) in dico_crible_max_avec_k.keys() :
            print(dico_crible_max_avec_k[(cle_1, cle_2)])
        else :
            print(dico_crible_max_avec_k[(cle_2, cle_1)])

def composantes_seuil_069(rep, depart, typ_sim, taille_ext):            
    with open(EXTENSION_PATH%rep+"graphe_complet_pondere_sim_%s_taille_%s.pickle"%(typ_sim, taille_ext), 'rb') as fichier_graphe :
        mon_depickler = pickle.Unpickler(fichier_graphe)
        graphe_complet = mon_depickler.load() 
        
        with open(EXTENSION_PATH%rep+"fichier_petits_groupes_%s_taille_%s_069.txt"%(typ_sim, taille_ext), 'w') as fichier_ecriture :
            premier_clustering(graphe_complet, 0.69, fichier_ecriture, depart, typ_sim, taille_ext, rep)
          
    
if __name__ == '__main__':

    
    print("ramousnif")
    #for i in range(7, 10) :
    #draw_isomorphisme(4)
    #draw_isomorphisme(8)
    for i in range(1,11) :
        creation_dico_comp(i) 
        #obtention(i)
#
    #draw_composantes("toutes_aretes_coeff_all1_max", "extensions","taille_max/result_k_max_4_10_toutes_aretes", "taille_max", 1, type_comp="extensions_toutes_aretes")
    #traitement_recherche_kmax(4, 10, "extensions", "toutes_aretes_coeff_all1", "toutes_aretes_coeff_all1")
    
    #afficher_sim_k_max("5J7L_DA_30_15", "5J7L_DA_25_10", "extensions", "toutes_aretes_coeff_all1_max", "taille_max/result_k_max_4_10_toutes_aretes")
    
#     composantes_seuil_069("taille_max/result_k_max_4_10_toutes_aretes", "extensions", "toutes_aretes_coeff_all1_max", "taille_max")
#     draw_une_composante("toutes_aretes_coeff_all1_max", "extensions", "taille_max", 1.0, 0.69, "taille_max/result_k_max_4_10_toutes_aretes" )
    #draw_une_grande_composante("toutes_aretes_coeff_all1", "extensions", 4, 4, 0.7)
    
    
    #traitement_recherche_kmax(4, 10, "extensions", "toutes_aretes_coeff_all1", "toutes_aretes_coeff_all1")
#     for i in range(1,11) :
#         draw_isomorphisme(i)
    #draw_isomorphisme(10)
    #obtention(5)
#     for i in range(2,11) :
#         comparaison_sein_taille(i, "extensions", "toutes_aretes_coeff_all1_chaines_2_4", "toutes_aretes_coeff_all1",[2,4] )
#     
    #observation(4)
    #groupe_1, groupe_2 = comparaison(9,10)
#     compter_composantes_differentes_par_groupe(groupe_1)
#     compter_composantes_differentes_par_groupe(groupe_2)
#     comparaison(5,10)
#     groupe_1, groupe_2 = comparaison(6,10)
#     compter_composantes_differentes_par_groupe(groupe_1)
#     comparaison(8,10)
#     tab_idem, tab_inclus_dans = comparaison_cluster_identiques_plus_de_2([4,5,6,7,8,-1,10])
#     tab_idem_567, tab_inclus_dans_567 = comparaison_cluster_identiques_plus_de_2([7,8,-1,10])
#     print("petit rat")
#     for groupe in tab_idem_567 :
#         existe = False 
#         for elt in groupe :
#             for groupe_tout in tab_idem :
#                 if elt in groupe_tout :
#                     existe = True
#         if existe == False :
#             print(groupe)
#             
#     print("petit rat")
#     for groupe in tab_inclus_dans_567 :
#         existe = False 
#         for elt in groupe :
#             for groupe_tout in tab_inclus_dans :
#                 if elt in groupe_tout :
#                     existe = True
#         if existe == False :
#             print(groupe)
#     comparaison_cluster_identiques_plus_de_2([6,7])
    #comparaison_cluster_identiques_plus_de_2([7,8,-1,10])
    #comparaison_cluster_identiques_plus_de_2([4,5])
    
#     for i in range(2,11) :
#         print("i : %d"%i)
#         dico_sim = stockage_sim("toutes_aretes_coeff_all1", "toutes_aretes_coeff_all1_par_k", "extensions", i)
#         groupe1, groupe2, comp_1, comp_2 = comparaison_cluster(i, i, "extensions", "toutes_aretes_coeff_all1", "extensions", "toutes_aretes_coeff_all1_par_k", 1, max(dico_sim.values()))
#         compter_composantes_differentes_par_groupe(groupe1)
#         compter_composantes_differentes_par_groupe(groupe2)
#         compter_composantes_differentes_par_groupe(comp_1)
#         compter_composantes_differentes_par_groupe(comp_2)
    
#     dico_sim = stockage_sim("toutes_aretes_coeff_all1", "toutes_aretes_coeff_all1_par_k", "extensions", 2)   
#     print(max(dico_sim.values()))
#     nombre = 0
#     for elt in dico_sim.keys() :
#         if dico_sim[elt] == max(dico_sim.values()) :
#             nombre += 1
#             print(elt)
#     print(nombre)
#     print(dico_sim[("5J7L_DA_25_10", '5DM6_X_25_34')])
#     
#     dico_sim = stockage_sim("toutes_aretes_coeff_all1", "toutes_aretes_coeff_all1", "extensions", 2)   
#     print(max(dico_sim.values()))
#     nombre = 0
#     for elt in dico_sim.keys() :
#         if dico_sim[elt] == max(dico_sim.values()) :
#             nombre += 1
#             print(elt)
#     print(nombre)

#
#     for i in range(2,7) :
#         groupe1, groupe2 = comparaison_cluster(i,"taille_max", "extensions", "toutes_aretes_coeff_all1", "extensions", "toutes_aretes_coeff_all1_par_k_max", 1, 4.5)   
#         compter_composantes_differentes_par_groupe(groupe1)
#         compter_composantes_differentes_par_groupe(groupe2)
    
#     comparaison_cluster_identiques_plus_de_2([4,5,6,7], "extensions", "toutes_aretes_coeff_all1_chaines_1_3")
    #draw_isomorphisme(2)
    
#     groupe1, groupe2, comp_1, comp_2 = comparaison_cluster("taille_max","taille_max", "taille_max/result_k_max_3_10_toutes_aretes", "taille_max/result_k_max_3_10",  "extensions", "toutes_aretes_coeff_all1_max", "extensions", "toutes_aretes_coeff_all1_par_k_max", 1, 4.5)   
#     #groupe1, groupe2, comp_1, comp_2 = comparaison_cluster("taille_max","taille_max", "taille_max/result_k_max_4_10_toutes_aretes", "taille_max/result_k_max_4_10",  "extensions", "toutes_aretes_coeff_all1_max", "extensions", "toutes_aretes_coeff_all1_par_k_max", 1, 4.5)   
#     
#     compter_composantes_differentes_par_groupe(groupe1)
#     compter_composantes_differentes_par_groupe(groupe2)
#     compter_composantes_differentes_par_groupe(comp_1)
#     compter_composantes_differentes_par_groupe(comp_2)
    
    #compter_composantes_differentes("extensions", "toutes_aretes_coeff_all1_par_k_max", "taille_max/result_k_max_4_10", "taille_max", 4.5)
    #afficher_groupes("extensions", "toutes_aretes_coeff_all1_par_k_max", "taille_max/result_k_max_2_10", "taille_max", 4.5)
    
    #draw_grandes_composantes("toutes_aretes_coeff_all1_max", "extensions", "taille_max", "taille_max/result_k_max_4_10_toutes_aretes", 1.0)
    
#     for i in range(len(GROUPES_TOUTES_ARETES_MAX_4_10_07)) :
#         draw_groupes_par_type("taille_max/result_k_max_4_10_toutes_aretes", "taille_max", "extensions", "toutes_aretes_coeff_all1_max", True, GROUPES_TOUTES_ARETES_MAX_4_10_07[i], "toutes_aretes_max_4_10_0.7_%d"%i, 0.7)
    
#     for i in range(len(GROUPES_TOUTES_ARETES_MAX_4_10_07[0])) :
#         afficher_sim_k_max(GROUPES_TOUTES_ARETES_MAX_4_10_07[0][i], "4V9F_0_25_56", "extensions", "toutes_aretes_coeff_all1_max", "taille_max/result_k_max_4_10_toutes_aretes" )
    
    #clustering_composante("taille_max/result_k_max_4_10_toutes_aretes", "extensions", "toutes_aretes_coeff_all1_max", "taille_max", 0.7)
    
    #retrouver_types("extensions", "toutes_aretes_coeff_all1_max", "taille_max/result_k_max_4_10_toutes_aretes", "taille_max", 1, "GNRA")
    #retrouver_lien_chaines_1_3("extensions", "toutes_aretes_coeff_all1_max", "taille_max/result_k_max_4_10_toutes_aretes", "taille_max", 1)
    #draw_groupes_par_type("taille_max/result_k_max_4_10_toutes_aretes", "taille_max", "extensions", "toutes_aretes_coeff_all1_max", True, GROUPE_ARICH_DE_07, "Groupe A-rich loop toutes aretes max seuil=0.7", 0.7)
    
#     groupe_etendu = calcul_groupe_etendu_seuil_sup("toutes_aretes_coeff_all1_max", "taille_max/result_k_max_4_10_toutes_aretes", "taille_max", GROUPES_TOUTES_ARETES_MAX_4_10_07[2], 0.6)
#     print(groupe_etendu)
#     print(len(groupe_etendu))
#     
#     
#     draw_groupe_gephi_unique(groupe_etendu, GROUPES_TOUTES_ARETES_MAX_4_10_07[2], "taille_max/result_k_max_4_10_toutes_aretes", "toutes_aretes_coeff_all1_max", "taille_max", 0.6, "groupe_4_10_toutes_aretes_0.7_55_2_etendu_0.6")
#     
#     with open(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes/"+ "sim_extensions_toutes_aretes_coeff_all1_max_taille_taille_max_avec_val_k.pickle", 'rb') as fichier_sim_1 :
#         mon_depickler_1 = pickle.Unpickler(fichier_sim_1)
#         dico_sim_1 = mon_depickler_1.load()  
#         
#         print(dico_sim_1)
#         if ('5FDU_1A_25_68', '4V88_A5_25_47') in dico_sim_1.keys() :
#             print(dico_sim_1[('5FDU_1A_25_68', '4V88_A5_25_47')])
#         else :
#             print(dico_sim_1[('4V88_A5_25_47', '5FDU_1A_25_68')])
#             
#         
#         if ('5J7L_DA_25_12', '4V88_A5_25_47') in dico_sim_1.keys() :
#             print(dico_sim_1[('5J7L_DA_25_12', '4V88_A5_25_47')])
#         else :
#             print(dico_sim_1[('4V88_A5_25_47', '5J7L_DA_25_12')])

    