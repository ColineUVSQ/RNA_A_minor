'''
Created on 9 janv. 2019

@author: coline

Etudes des composantes connexes et autres groupes qu'on peut obtenir avec les donnees de CaRNAval 
(affichage, statistiques de taille et de nombre, clustering...)
'''
import pickle
import networkx as nx
import matplotlib.pyplot as plt 
import numpy as np  
import os
import csv
from recup_data.calcul_sim import calcul_sim_non_cov_sans_motif, \
    calcul_sim_avec_poids, calcul_sim_non_cov_sans_motif_par_chaine, \
    calcul_sim_avec_poids_par_chaine, calcul_sim_aretes_avec_coeff_par_chaine
from recup_data.sous_graphe_commun_max_version_grands_graphes import ajout_attribut_chaine, \
    recup_chaines
from recup_data.constantes import EXTENSION_PATH, GROUPE_GBULGE, GROUPE_GNRA, \
    GROUPE_ARICH, EXTENSION_PATH_TAILLE, HOMOLOGUES
from sklearn.cluster import KMeans
from sklearn.cluster.dbscan_ import DBSCAN
from sklearn.cluster import AgglomerativeClustering
from recup_data.retrouver_GNRA import lien_chaines_1_3

''' dans le graphe global passe en parametre on ajoute les numeros de chaine en attributs '''
def ajout_chaines_grands_graphes(graphe, occ_a_minor): 
    nx.set_node_attributes(graphe, -1, "chaine")
    chaine = recup_chaines(graphe, occ_a_minor)
    
    for i in range(4) :
        for e in chaine[i] :
#             print(e)
            if graphe.nodes[e]["chaine"] == -1 :
                graphe.nodes[e]["chaine"] = []
            if i + 1 not in graphe.nodes[e]["chaine"] :

                graphe.nodes[e]["chaine"].append(i + 1)

#             print("voisins")
            for voisin in graphe[e] :
#                 print(voisin)
                if graphe.nodes[voisin]["chaine"] == -1 :
                    graphe.nodes[voisin]["chaine"] = []
                if i + 1 not in graphe.nodes[voisin]["chaine"] and voisin not in occ_a_minor :
                    graphe.nodes[voisin]["chaine"].append(i + 1)
    
    for noeud in graphe.nodes() :
        if graphe.nodes[noeud]["chaine"] == -1 :
            print(noeud)
            
    return graphe
        
''' renvoie le nombre d'arcs et d'aretes du graphe passe en parametre '''
def nombre_arc_arete_graphe(graphe):
    compteur_arc = 0
    compteur_arete = 0
    for (u, v, t) in graphe.edges(data="label") :
        if t == 'B53' :
            compteur_arc += 1
        else :
            compteur_arete += 1
    return compteur_arc, compteur_arete


''' renvoie les noeuds en commun entre les graphe1 et graphe2 chaine par chaine (dans un tableau à deux dimensions de 4 lignes, une par chaine)
graphe_comp : le graphe complet contenant toutes les occurrences d'Aminor (donnees CaRNAval) et avec aretes ponderees par les similarites
u : numero dans graphe_comp du premier graphe considere 
v : numero dans graphe_comp du deuxieme graphe considere
graphe1 : premier graphe d'extension considere
graphe2 : deuxieme graphe d'extension considere
fichier_graphe_commun : nom du fichier dans lequel est stocke le dictionnaire avec tous les graphes communs '''
def recherche_chaines_par_morceau(u, v, depart, graphe_comp, graphe1, graphe2, fichier_graphe_commun):
    if depart == "extensions" :
#         with open("graphes_extension/fichier_"+graphe_comp.nodes[u]["nom"]+".pickle", 'rb') as fichier_graphe1 :
#             mon_depickler_graphe1 = pickle.Unpickler(fichier_graphe1)
#             graphe1 = mon_depickler_graphe1.load()
#             with open("graphes_extension/fichier_"+graphe_comp.nodes[v]["nom"]+".pickle", 'rb') as fichier_graphe2 :
#                 mon_depickler_graphe2 = pickle.Unpickler(fichier_graphe2)
#                 graphe2 = mon_depickler_graphe2.load()
                
                chaines_1 = [[], [], [], []]
                
                for noeud, ch in graphe1.nodes(data="chaine") :

                    if 1 in ch :
                        chaines_1[0].append(noeud)
                    if 2 in ch :
                        chaines_1[1].append(noeud)
                    if 3 in ch :
                        chaines_1[2].append(noeud)
                    if 4 in ch :
                        chaines_1[3].append(noeud)
         
                chaines_2 = [[], [], [], []]
                
                for noeud, ch in graphe2.nodes(data="chaine") :
                    
                    if 1 in ch :
                        chaines_2[0].append(noeud)
                    if 2 in ch :
                        chaines_2[1].append(noeud)
                    if 3 in ch :
                        chaines_2[2].append(noeud)
                    if 4 in ch :
                        chaines_2[3].append(noeud)
                        
                with open(fichier_graphe_commun, 'rb') as fichier_commun :
                    mon_depickler_commun = pickle.Unpickler(fichier_commun)
                    dico_graphe = mon_depickler_commun.load()
                    
                    chaines_commun = [[], [], [], []]
                    if ("fichier_" + graphe_comp.nodes[u]["nom"], "fichier_" + graphe_comp.nodes[v]["nom"]) in dico_graphe.keys() :
                        graphe_commun = dico_graphe[("fichier_" + graphe_comp.nodes[u]["nom"], "fichier_" + graphe_comp.nodes[v]["nom"])]
                        
                        for noeud in graphe_commun.nodes() :
                        
                            for j in range(4) :
                                if noeud[0] in chaines_1[j] and noeud[1] in chaines_2[j] :
                                    chaines_commun[j].append(noeud)

                        return chaines_1, chaines_2, chaines_commun
                    else :
                        graphe_commun = dico_graphe[("fichier_" + graphe_comp.nodes[v]["nom"], "fichier_" + graphe_comp.nodes[u]["nom"])]
                        for noeud in graphe_commun.nodes() :
                        
                            for j in range(4) :
                                if noeud[0] in chaines_2[j] and noeud[1] in chaines_1[j] :
                                    chaines_commun[j].append(noeud)

                        return chaines_1, chaines_2, chaines_commun
                    
    else :
#         with open("grands_graphes.pickle", 'rb') as fichier_graphes :
#             mon_depickler_graphes = pickle.Unpickler(fichier_graphes)
#             graphes = mon_depickler_graphes.load()
#             
            elt1 = (graphe_comp.nodes[u]["nom"].split("_")[0], graphe_comp.nodes[u]["nom"].split("_")[1], int(graphe_comp.nodes[u]["nom"].split("_")[2]), int(graphe_comp.nodes[u]["nom"].split("_")[3]))
            elt2 = (graphe_comp.nodes[v]["nom"].split("_")[0], graphe_comp.nodes[v]["nom"].split("_")[1], int(graphe_comp.nodes[v]["nom"].split("_")[2]), int(graphe_comp.nodes[v]["nom"].split("_")[3]))
   
            chaines_1 = [[], [], [], []]
                
            for noeud, ch in graphe1.nodes(data="chaine") :
                if 1 in ch :
                    chaines_1[0].append(noeud)
                if 2 in ch :
                    chaines_1[1].append(noeud)
                if 3 in ch :
                    chaines_1[2].append(noeud)
                if 4 in ch :
                    chaines_1[3].append(noeud)
     
            chaines_2 = [[], [], [], []]
            
            for noeud, ch in graphe2.nodes(data="chaine") :
                if 1 in ch :
                    chaines_2[0].append(noeud)
                if 2 in ch :
                    chaines_2[1].append(noeud)
                if 3 in ch :
                    chaines_2[2].append(noeud)
                if 4 in ch :
                    chaines_2[3].append(noeud)

            with open(fichier_graphe_commun, 'rb') as fichier_commun :
                mon_depickler_commun = pickle.Unpickler(fichier_commun)
                dico_graphe = mon_depickler_commun.load()
                
                chaines_commun = [[], [], [], []]
#                 print(dico_graphe.keys())
                if (elt1, elt2) in dico_graphe.keys() :
                    graphe_commun = dico_graphe[(elt1, elt2)]
                    
                    for noeud in graphe_commun.nodes() :
                        for j in range(4) :
                            if noeud[0] in chaines_1[j] and noeud[1] in chaines_2[j] :
                                chaines_commun[j].append(noeud)
                    
                    return chaines_1, chaines_2, chaines_commun
                else :
                    graphe_commun = dico_graphe[(elt2, elt1)]
                
                    for noeud in graphe_commun.nodes() :
                        for j in range(4) :
                            if noeud[0] in chaines_2[j] and noeud[1] in chaines_1[j] :
                                chaines_commun[j].append(noeud)
                    
                    return chaines_2, chaines_1, chaines_commun
   

''' calcule les similarites associees a chaque chaine pour toutes les paires de graphe
et les stocke dans un dictionnaire '''   
def calcul_sim_morceau(typ, depart, fichier_graphe_commun) :
    with open("Extensions/Metrique_toutes_aretes/graphe_complet_pondere_sim.pickle", 'rb') as fichier_graphe_complet :
        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
        graphe_comp = mon_depickler_complet.load()
        dico_sim_par_chaine = {}
        if depart == "extensions" :
            
#             i = 0.1
#             #i = 0.4
#             while i <= 1.0 :
#                 with open("composantes_connexes_"+type_comp+"_"++str(i)+".pickle", 'rb') as fichier_comp:
#                     mon_depickler_comp = pickle.Unpickler(fichier_comp)
#                     composantes_connexes = mon_depickler_comp.load()
#                     
#                     for composante in composantes_connexes :
#                         if len(composante) < 10 and len(composante) > 2 :
            # graphe_comp = graphe_complet.subgraph(groupe)
            compteur = 0
            for u, v in graphe_comp.edges() :

                if graphe_comp.nodes[u]["nom"] in ['1FJG_A_48_11', '3UCZ_R_62_15', '4V9F_0_48_26', '4V9F_0_44_3', '5J7L_DA_50_21', '1FJG_A_58_23', '1U9S_A_58_11'] and graphe_comp.nodes[v]["nom"] in ['1FJG_A_48_11', '3UCZ_R_62_15', '4V9F_0_48_26', '4V9F_0_44_3', '5J7L_DA_50_21', '1FJG_A_58_23', '1U9S_A_58_11'] :
                    print(graphe_comp.nodes[u]["nom"])
                    print(graphe_comp.nodes[v]["nom"])
                    print(compteur)
                    compteur = compteur + 1
                    with open("Extensions/Metrique_toutes_aretes/graphes_extension/fichier_" + graphe_comp.nodes[u]["nom"] + ".pickle", 'rb') as fichier_graphe1 :
                        mon_depickler_graphe1 = pickle.Unpickler(fichier_graphe1)
                        graphe1 = mon_depickler_graphe1.load()
                        with open("Extensions/Metrique_toutes_aretes/graphes_extension/fichier_" + graphe_comp.nodes[v]["nom"] + ".pickle", 'rb') as fichier_graphe2 :
                            mon_depickler_graphe2 = pickle.Unpickler(fichier_graphe2)
                            graphe2 = mon_depickler_graphe2.load()
                            chaines_1, chaines_2, chaines_commun = recherche_chaines_par_morceau(u, v, depart, graphe_comp, graphe1, graphe2, fichier_graphe_commun)
                            
                            with open(fichier_graphe_commun, 'rb') as fichier_commun :
                                mon_depickler_commun = pickle.Unpickler(fichier_commun)
                                dico_graphe = mon_depickler_commun.load()
                                
                                tab_sim = []
                                if ("fichier_" + graphe_comp.nodes[u]["nom"], "fichier_" + graphe_comp.nodes[v]["nom"]) in dico_graphe.keys() :
                                    cle = ("fichier_" + graphe_comp.nodes[u]["nom"], "fichier_" + graphe_comp.nodes[v]["nom"])
                                    graphe_commun = dico_graphe[cle]
#                                     for j in range(4) :
#                                         if typ == "longue_distance" :
#                                             sim = calcul_sim_non_cov_sans_motif(graphe1, graphe2, graphe_commun, chaines_1, chaines_2, chaines_commun, j)
#                                         if typ == "non_cov" :
#                                             sim = calcul_sim_avec_poids_par_chaine(graphe1, graphe2, graphe_commun, chaines_1, chaines_2, chaines_commun, j)
#                                         if typ == "toutes_aretes" :
#                                             sim = calcul_sim_aretes_avec_coeff_par_chaine(graphe1, graphe2, graphe_commun, chaines_1, chaines_2, chaines_commun, cle, 1, 1, 1, j)
#                                         #print(sim)                    
#                                         
#                                         tab_sim.append(sim)
                                    if typ == "toutes_aretes" :
                                        sim = calcul_sim_aretes_avec_coeff_par_chaine(graphe1, graphe2, graphe_commun, chaines_1, chaines_2, chaines_commun, cle, 1, 1, 1, [0, 2])
                                        tab_sim.append(sim)
                                    
                                    dico_sim_par_chaine.update({(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]): tab_sim})
                                    
                                else :
                                    cle = ("fichier_" + graphe_comp.nodes[v]["nom"], "fichier_" + graphe_comp.nodes[u]["nom"])
                                    graphe_commun = dico_graphe[cle]
#                                     for j in range(4) :
#                                         if typ == "longue_distance" :
#                                             sim = calcul_sim_non_cov_sans_motif(graphe2, graphe1, graphe_commun, chaines_2, chaines_1, chaines_commun, j)
#                                         if typ == "non_cov" :
#                                             sim = calcul_sim_avec_poids_par_chaine(graphe2, graphe1, graphe_commun, chaines_2, chaines_1, chaines_commun, j)
#                                         #print(sim)   
#                                         if typ == "toutes_aretes" :
#                                             sim = calcul_sim_aretes_avec_coeff_par_chaine(graphe2, graphe1, graphe_commun, chaines_2, chaines_1, chaines_commun, cle, 1, 1, 1, j)
#                      
#                                         tab_sim.append(sim)
                                    if typ == "toutes_aretes" :
                                        sim = calcul_sim_aretes_avec_coeff_par_chaine(graphe2, graphe1, graphe_commun, chaines_2, chaines_1, chaines_commun, cle, 1, 1, 1, [0, 2])
                                        tab_sim.append(sim)
                                    dico_sim_par_chaine.update({(graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"]): tab_sim})
                                            
#                 i = i+0.1
        else :
#             i = 0.1
#             #i = 0.4
#             while i <= 1.0 :
#                 with open("composantes_connexes_"+type_comp+"_"+str(i)+".pickle", 'rb') as fichier_comp:
#                     mon_depickler_comp = pickle.Unpickler(fichier_comp)
#                     composantes_connexes = mon_depickler_comp.load()
#                     
#                     for composante in composantes_connexes :
#                         if len(composante) < 10 and len(composante) > 2 :
            # graphe_comp = graphe_complet.subgraph(composante)
            for u, v in graphe_comp.edges() :
#                                 if (graphe_comp.nodes[u]["nom"] == '4V9F_0_48_26' and graphe_comp.nodes[v]["nom"] == '5J7L_DA_48_30') or  (graphe_comp.nodes[v]["nom"] == '4V9F_0_48_26' and graphe_comp.nodes[u]["nom"] == '5J7L_DA_48_30') :
#                                 
                    chaines_1, chaines_2, chaines_commun = recherche_chaines_par_morceau(u, v, graphe_comp, depart)
                    print(u)
                    print(v)
                    print(chaines_commun)
                    with open("grands_graphes.pickle", 'rb') as fichier_graphes :
                        mon_depickler_graphes = pickle.Unpickler(fichier_graphes)
                        graphes = mon_depickler_graphes.load()
                        
                        with open("fichier_comp_grands_graphes_V2.pickle", 'rb') as fichier_commun :
                            mon_depickler_commun = pickle.Unpickler(fichier_commun)
                            dico_graphe = mon_depickler_commun.load()
                            
                            elt1 = (graphe_comp.nodes[u]["nom"].split("_")[0], graphe_comp.nodes[u]["nom"].split("_")[1], int(graphe_comp.nodes[u]["nom"].split("_")[2]), int(graphe_comp.nodes[u]["nom"].split("_")[3]))
                            elt2 = (graphe_comp.nodes[v]["nom"].split("_")[0], graphe_comp.nodes[v]["nom"].split("_")[1], int(graphe_comp.nodes[v]["nom"].split("_")[2]), int(graphe_comp.nodes[v]["nom"].split("_")[3]))
                            
                            if (elt1, elt2) in dico_graphe.keys() :
                                cle = (elt1, elt2)
                                cle_sim = (graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])
                            else :
                                cle = (elt2, elt1)
                                cle_sim = (graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"])
                                
                            if cle_sim not in dico_sim_par_chaine.keys() :
                                
                                chaines_1, chaines_2, chaines_commun = recherche_chaines_par_morceau(u, v, depart, graphe_comp, graphes[cle[0]], graphes[cle[1]], fichier_graphe_commun)
                                
                                graphe_commun = dico_graphe[cle]
                                tab_sim = []
                                for j in range(4) :
                                    if typ == "non_cov" :
                                        sim = calcul_sim_non_cov_sans_motif_par_chaine(graphes[cle[0]], graphes[cle[1]], graphe_commun, chaines_1, chaines_2, chaines_commun, j)
                                        if (elt1 == ('4V9F', '0', 48, 26) and elt2 == ('5J7L', 'DA', 48, 30)) or (elt2 == ('4V9F', '0', 48, 26) and elt1 == ('5J7L', 'DA', 48, 30))  :
                                            print(cle[0])
                                            print(cle[1])
                                            print(j)
                                            print(sim)
                                    # print(sim)                    
                                    tab_sim.append(sim)
                                    dico_sim_par_chaine.update({cle_sim : tab_sim})
                                            
#                 i = i+0.1
        
    with open("Extensions/Metrique_toutes_aretes/sim_" + depart + "_" + typ + "_par_chaine.pickle", 'wb') as fichier_graphe_complet_plus :
        mon_pickler_complet = pickle.Pickler(fichier_graphe_complet_plus)
        mon_pickler_complet.dump(dico_sim_par_chaine)


''' affiche chaque composante connexes separement dans un fichier .png avec les aretes ponderees par la sim
en rouge les aretes ponderees par une valeur au-dessus du seuil (composantes de taille < 10 et > 2)
type_sim, depart_sim : infos pour choisir le fichier a considerer pour la similarite
taille_ext : taille de l'extensions consideree
max_val : valeur maximale de la similarite (pour savoir comment on a decoupe les seuils)
rep : repertoire dans lequel on stocke le fichier .png'''
def draw_composantes(typ_sim, depart_sim, rep, taille_ext, max_val, *args, **kwargs) :
    par_chaine = kwargs.get('par_chaine', False)
    type_comp = kwargs.get('type_comp', "")
    
    if par_chaine == False :
#         with open("Extensions/Metrique_toutes_aretes/sim_"+depart_sim+"_"+typ_sim+".pickle", 'rb') as fichier_sim :
#             mon_depickler = pickle.Unpickler(fichier_sim)
#             dico_sim = mon_depickler.load()
#             print(dico_sim.keys())
            
            with open(EXTENSION_PATH % rep + "sim_%s_%s_taille_%s_avec_val_k.pickle" % (depart_sim, typ_sim, taille_ext), 'rb') as fichier_sim :
                mon_depickler = pickle.Unpickler(fichier_sim)
                dico_sim_1 = mon_depickler.load()
                print(dico_sim_1.keys())
            
#                 with open("Extensions/Metrique_toutes_aretes/sim_extensions_toutes_aretes_coeffn1_a0_c0_superposition_all1.pickle", 'rb') as fichier_sim :
#                     mon_depickler = pickle.Unpickler(fichier_sim)
#                     dico_sim_2 = mon_depickler.load()
#                     print(dico_sim_2.keys())
            
                with open(EXTENSION_PATH % rep + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (typ_sim, taille_ext), 'rb') as fichier_graphe_complet :
                        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
                        graphe_complet = mon_depickler_complet.load()
                        
                        i = 0.1 * max_val
                        while i <= 1.0 * max_val + 0.000005 :
                            # i = 0.6
                            with open(EXTENSION_PATH % rep + "composantes_connexes_" + depart_sim + "_" + typ_sim + "_taille_" + str(taille_ext) + "_" + str(round(i, 2)) + ".pickle", 'rb') as fichier_comp:
                            # with open("Extensions/Metrique_toutes_aretes/composantes_connexes_extensions_toutes_aretes_coeffn1_a1_c1_0.6_0.6.pickle", 'rb') as fichier_comp:
                            
                                mon_depickler_comp = pickle.Unpickler(fichier_comp)
                                composantes_connexes = mon_depickler_comp.load()
                                
                                compteur = 1
                                for composante in composantes_connexes :
                                    
                                    if len(composante) < 10 and len(composante) > 2 :
                                        graphe_comp = nx.Graph()
                                        graphe_comp = graphe_complet.subgraph(composante)
                                        nx.set_node_attributes(graphe_comp, (33, 33), "coordonnees")
                                        pos = nx.circular_layout(graphe_comp)
                                        
                                        node_labels = dict([(u, (d["nom"]))for u, d in graphe_comp.nodes(data=True)])
                
                                        for u, v, d in graphe_comp.edges(data=True) :
                                            
                                            if (graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]) in dico_sim_1.keys() :
                                                graphe_comp.edges[u, v]["poids"] = (round(dico_sim_1[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])][0], 2), dico_sim_1[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])][1])
                                            else :
                                                graphe_comp.edges[u, v]["poids"] = (round(dico_sim_1[(graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"])][0], 2), dico_sim_1[(graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"])][1])
                                        
#                                         for u,v,d in graphe_comp.edges(data=True) :
#                                             
#                                             if (graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]) in dico_sim_2.keys() :
#                                                 graphe_comp.edges[u,v]["poids"].append(round(dico_sim_2[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])],2))
#                                             else :
#                                                 graphe_comp.edges[u,v]["poids"].append(round(dico_sim_2[(graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"])],2))
#         
                                            
                                        edge_labels = dict([((u, v), d["poids"])for u, v, d in graphe_comp.edges(data=True)])
                                        
                                        red_edges = []
                                        
                                        for u, v, d in graphe_comp.edges(data=True) :
                                            if d["poids"][0] >= i :  # and d["poids"][1] >= i:
                                                red_edges.append((u, v))
                                        
                                        edge_colors = ['red' if edge in red_edges else 'black' for edge in graphe_comp.edges()]
                                        densite = len(red_edges) / ((len(composante) * (len(composante) - 1)) / 2)
                                        
                                        plt.figure(figsize=(9, 7))
                                        # plt.title(depart_sim)
                                        print(i / max_val)
                                        plt.title("poids >= " + str(round(i, 2)) + " " + str(int((i / max_val) * 100)) + "% densite : " + str(round(densite, 2)))
                                        plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.9)
                                        
                                        nx.draw_networkx_nodes(graphe_comp, pos)
                                        nx.draw_networkx_labels(graphe_comp, pos, labels=node_labels, font_size=8)
                                        nx.draw_networkx_edges(graphe_comp, pos, edge_color=edge_colors)
                                        nx.draw_networkx_edge_labels(graphe_comp, pos, edge_labels=edge_labels, label_pos=0.3)
                                        plt.axis('off')
                                        
                                        # plt.show()
                                        # plt.savefig("Extensions/Metrique_toutes_aretes/composantes_connexes_"+depart_sim+"_"+typ_sim+"/comp_"+type_comp+"_"+str(i)+"_"+str(i)+"_"+str(compteur)+".png") # save as png
                                        os.makedirs(EXTENSION_PATH % taille_ext + "composantes_connexes_%s_%s_taille_%s/" % (depart_sim, typ_sim, taille_ext), exist_ok=True)
                                        plt.savefig(EXTENSION_PATH % taille_ext + "composantes_connexes_%s_%s_taille_%s/comp_" % (depart_sim, typ_sim, taille_ext) + type_comp + "_" + str(round(i, 2)) + "_" + str(compteur) + ".png")  # save as png

                                        plt.close()
                                        compteur += 1
                                    
                            i = i + (0.1 * max_val)
    else :
        with open("sim_" + depart_sim + "_" + typ_sim + "par_chaine.pickle", 'rb') as fichier_sim :
            mon_depickler = pickle.Unpickler(fichier_sim)
            dico_sim = mon_depickler.load()
            print(dico_sim.keys())
            
            with open("graphe_complet_pondere_sim.pickle", 'rb') as fichier_graphe_complet :
                mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
                graphe_complet = mon_depickler_complet.load()
                
                i = 0.1
                while i <= 1.0 :
                    with open("composantes_connexes_" + type_comp + "_" + str(i) + ".pickle", 'rb') as fichier_comp:
                        mon_depickler_comp = pickle.Unpickler(fichier_comp)
                        composantes_connexes = mon_depickler_comp.load()
                        
                        compteur = 1
                        for composante in composantes_connexes :
                            
                            if len(composante) < 10 and len(composante) > 2 :
                                graphe_comp = nx.Graph()
                                graphe_comp = graphe_complet.subgraph(composante)
                                nx.set_node_attributes(graphe_comp, (33, 33), "coordonnees")
                                pos = nx.circular_layout(graphe_comp)
                                
                                node_labels = dict([(u, (d["nom"]))for u, d in graphe_comp.nodes(data=True)])
        
                                for u, v, d in graphe_comp.edges(data=True) :
                                    
                                    if (graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]) in dico_sim.keys() :
                                        for j in range(4) :
                                            print(dico_sim[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])][j])
                                        graphe_comp.edges[u, v]["poids"] = (dico_sim[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])][0][0] + dico_sim[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])][1][0] + dico_sim[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])][2][0] + dico_sim[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])][3][0]) / 4
                                    else :
                                        graphe_comp.edges[u, v]["poids"] = (dico_sim[(graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"])][0][0] + dico_sim[(graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"])][1][0] + dico_sim[(graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"])][2][0] + dico_sim[(graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"])][3][0]) / 4
                                    
                                edge_labels = dict([((u, v), (round(d["poids"], 2)))for u, v, d in graphe_comp.edges(data=True)])
                                
                                red_edges = []
                                
                                for u, v, d in graphe_comp.edges(data=True) :
                                    if d["poids"] >= i :
                                        red_edges.append((u, v))
                                
                                edge_colors = ['red' if edge in red_edges else 'black' for edge in graphe_comp.edges()]
                                densite = len(red_edges) / ((len(composante) * (len(composante) - 1)) / 2)
                                
                                plt.figure(figsize=(9, 7))
                                plt.title(depart_sim)
                                # plt.title("poids > "+str(round(i,1)) + " densite : "+ str(round(densite,2)))
                                plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.9)
                                
                                nx.draw_networkx_nodes(graphe_comp, pos)
                                nx.draw_networkx_labels(graphe_comp, pos, labels=node_labels, font_size=8)
                                nx.draw_networkx_edges(graphe_comp, pos, edge_color=edge_colors)
                                nx.draw_networkx_edge_labels(graphe_comp, pos, edge_labels=edge_labels, label_pos=0.3)
                                plt.axis('off')
                                
                                # plt.show()
                                plt.savefig("composantes_connexes_" + depart_sim + "_" + typ_sim + "par_chaine/comp_" + type_comp + "_" + str(round(i, 1)) + "_" + str(compteur) + ".png")  # save as png
                                plt.close()
                                compteur += 1
                            
                    i = i + 0.1


''' affiche la ou les composantes connexes associees a un seuil (dans des fichier .png, avec les aretes ponderees par la sim
en rouge les aretes ponderees par une valeur au-dessus du seuil (composantes de taille < 10 et > 2))
type_sim, depart_sim : infos pour choisir le fichier a considerer pour la similarite 
taille_ext : taille de l'extensions consideree
max_val : valeur maximale de la similarite (pour savoir comment on a decoupe les seuils)
rep : repertoire dans lequel on stocke le fichier .png
i : seuil a considerer'''
def draw_une_composante(typ_sim, depart_sim, taille_ext, max_val, i, rep) :
    with open(EXTENSION_PATH % rep + "sim_%s_%s_taille_%s.pickle" % (depart_sim, typ_sim, taille_ext), 'rb') as fichier_sim :
                mon_depickler = pickle.Unpickler(fichier_sim)
                dico_sim_1 = mon_depickler.load()
                print(dico_sim_1.keys())
            
#                 with open("Extensions/Metrique_toutes_aretes/sim_extensions_toutes_aretes_coeffn1_a0_c0_superposition_all1.pickle", 'rb') as fichier_sim :
#                     mon_depickler = pickle.Unpickler(fichier_sim)
#                     dico_sim_2 = mon_depickler.load()
#                     print(dico_sim_2.keys())
            
                with open(EXTENSION_PATH % rep + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (typ_sim, taille_ext), 'rb') as fichier_graphe_complet :
                        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
                        graphe_complet = mon_depickler_complet.load()
                        
                        with open(EXTENSION_PATH % rep + "composantes_connexes_" + depart_sim + "_" + typ_sim + "_taille_" + str(taille_ext) + "_" + str(round(i, 2)) + ".pickle", 'rb') as fichier_comp:
                            # with open("Extensions/Metrique_toutes_aretes/composantes_connexes_extensions_toutes_aretes_coeffn1_a1_c1_0.6_0.6.pickle", 'rb') as fichier_comp:
                            
                                mon_depickler_comp = pickle.Unpickler(fichier_comp)
                                composantes_connexes = mon_depickler_comp.load()
                                
                                compteur = 1
                                for composante in composantes_connexes :
                                    
                                    if len(composante) < 10 and len(composante) > 2 :
                                        graphe_comp = nx.Graph()
                                        graphe_comp = graphe_complet.subgraph(composante)
                                        nx.set_node_attributes(graphe_comp, (33, 33), "coordonnees")
                                        pos = nx.circular_layout(graphe_comp)
                                        
                                        node_labels = dict([(u, (d["nom"]))for u, d in graphe_comp.nodes(data=True)])
                
                                        for u, v, d in graphe_comp.edges(data=True) :
                                            
                                            if (graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]) in dico_sim_1.keys() :
                                                graphe_comp.edges[u, v]["poids"] = round(dico_sim_1[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])], 2)
                                            else :
                                                graphe_comp.edges[u, v]["poids"] = round(dico_sim_1[(graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"])], 2)
                                        
#                                         for u,v,d in graphe_comp.edges(data=True) :
#                                             
#                                             if (graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]) in dico_sim_2.keys() :
#                                                 graphe_comp.edges[u,v]["poids"].append(round(dico_sim_2[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])],2))
#                                             else :
#                                                 graphe_comp.edges[u,v]["poids"].append(round(dico_sim_2[(graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"])],2))
#         
                                            
                                        edge_labels = dict([((u, v), d["poids"])for u, v, d in graphe_comp.edges(data=True)])
                                        
                                        red_edges = []
                                        
                                        for u, v, d in graphe_comp.edges(data=True) :
                                            if d["poids"] >= i :  # and d["poids"][1] >= i:
                                                red_edges.append((u, v))
                                        
                                        edge_colors = ['red' if edge in red_edges else 'black' for edge in graphe_comp.edges()]
                                        densite = len(red_edges) / ((len(composante) * (len(composante) - 1)) / 2)
                                        
                                        plt.figure(figsize=(9, 7))
                                        # plt.title(depart_sim)
                                        print(i / max_val)
                                        plt.title("poids >= " + str(round(i, 2)) + " " + str(int((i / max_val) * 100)) + "% densite : " + str(round(densite, 2)))
                                        plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.9)
                                        
                                        nx.draw_networkx_nodes(graphe_comp, pos)
                                        nx.draw_networkx_labels(graphe_comp, pos, labels=node_labels, font_size=8)
                                        nx.draw_networkx_edges(graphe_comp, pos, edge_color=edge_colors)
                                        nx.draw_networkx_edge_labels(graphe_comp, pos, edge_labels=edge_labels, label_pos=0.3)
                                        plt.axis('off')
                                        
                                        # plt.show()
                                        # plt.savefig("Extensions/Metrique_toutes_aretes/composantes_connexes_"+depart_sim+"_"+typ_sim+"/comp_"+type_comp+"_"+str(i)+"_"+str(i)+"_"+str(compteur)+".png") # save as png
                                        # os.makedirs(EXTENSION_PATH%taille_ext+"composantes_connexes_%s_%s_taille_%s/"%(depart_sim, typ_sim, taille_ext), exist_ok=True)
                                        plt.savefig(EXTENSION_PATH % rep + "composantes_connexes_%s_%s_taille_%s/comp_" % (depart_sim, typ_sim, taille_ext) + "_" + str(round(i, 3)) + "_" + str(compteur) + ".png")  # save as png

                                        plt.close()
                                        compteur += 1

''' affiche les grandes composantes connexes (de taille > 10)
stocke les informations dans des fichiers csv utilisables par Gephi
calcule le clustering DBSCAN, kmeans et agglomerative clustering pour chaque composante et stocke l'information dans le fichier csv des noeuds '''                                    
def draw_grandes_composantes(typ_sim, depart_sim, taille_ext, rep, max_val, *args, **kwargs):
    
    homologues = [['4V9F_0_62_12', '5J7L_DA_62_5', '5DM6_X_328_2', '5FDU_1A_62_14'],
                  ['4V9F_0_30_4', '5J7L_DA_30_15', '5FDU_1A_30_17'],
                  ['4V9F_0_48_13', '5J7L_DA_48_20', '5DM6_X_48_28', '5FDU_1A_48_25'],
                  ['4V9F_0_207_3', '5J7L_DA_272_2', '5FDU_1A_272_1'],
                  ['4V9F_0_25_56', '5J7L_DA_25_12', '5FDU_1A_25_68', '4V88_A5_25_30'],
                  ['4V9F_0_137_5', '5J7L_DA_48_1', '5FDU_1A_137_6', '5DM6_X_48_10', '4V88_A5_48_3'],
                  ['4V9F_0_134_5', '5FDU_1A_74_7', '5DM6_X_127_7'],
                  ['4V9F_0_48_21', '5J7L_DA_197_4', '5FDU_1A_197_3', '5DM6_X_48_9'],
                  ['4V9F_0_127_6', '5J7L_DA_134_1', '5FDU_1A_134_3', '5DM6_X_134_2'],
                  ['4V9F_0_287_2', '5DM6_X_25_15'],
                  ['5J7L_DA_25_10', '5DM6_X_25_34', '4V88_A5_25_47', '5FDU_1A_25_78'],
                  ['1FJG_A_48_17', '5J5B_BA_48_14'],
                  ['1FJG_A_48_8', '5J5B_BA_48_23'],
                  ['1FJG_A_294_1', '5J5B_BA_294_2'],
                  ['1FJG_A_58_23', '5J5B_BA_58_3'],
                  ['1FJG_A_138_3', '5J5B_BA_138_2'],
                  ['5J5B_BA_48_7', '4V88_A6_48_12']]
    
    par_chaine = kwargs.get('par_chaine', False)
    type_comp = kwargs.get('type_comp', "")
    
    if par_chaine == False :
#         with open("Extensions/Metrique_toutes_aretes/sim_"+depart_sim+"_"+typ_sim+".pickle", 'rb') as fichier_sim :
#             mon_depickler = pickle.Unpickler(fichier_sim)
#             dico_sim = mon_depickler.load()
#             print(dico_sim.keys())
            
#             with open(EXTENSION_PATH%taille_ext+"sim_%s_%s_taille_%s.pickle"%(depart_sim, typ_sim, taille_ext), 'rb') as fichier_sim :
#                 mon_depickler = pickle.Unpickler(fichier_sim)
#                 dico_sim_1 = mon_depickler.load()
#                 print(dico_sim_1.keys())
            
#                 with open("Extensions/Metrique_toutes_aretes/sim_extensions_toutes_aretes_coeffn1_a0_c0_superposition_all1.pickle", 'rb') as fichier_sim :
#                     mon_depickler = pickle.Unpickler(fichier_sim)
#                     dico_sim_2 = mon_depickler.load()
#                     print(dico_sim_2.keys())
            
                with open(EXTENSION_PATH % rep + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (typ_sim, taille_ext), 'rb') as fichier_graphe_complet :
                        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
                        graphe_complet = mon_depickler_complet.load()
                        
                        i = 0.1 * max_val
                        while i <= 1.0 * max_val + 0.000005 :
                            # i = 0.6
                            with open(EXTENSION_PATH % rep + "composantes_connexes_" + depart_sim + "_" + typ_sim + "_taille_" + str(taille_ext) + "_" + str(round(i, 2)) + ".pickle", 'rb') as fichier_comp:
                            # with open("Extensions/Metrique_toutes_aretes/composantes_connexes_extensions_toutes_aretes_coeffn1_a1_c1_0.6_0.6.pickle", 'rb') as fichier_comp:
                            
                                mon_depickler_comp = pickle.Unpickler(fichier_comp)
                                composantes_connexes = mon_depickler_comp.load()
                                
                                compteur = 1
                                for composante in composantes_connexes :
                                    
                                    if len(composante) > 8  :
                                        graphe_comp = nx.Graph()
                                        graphe_comp = graphe_complet.subgraph(composante)
                                        nx.set_node_attributes(graphe_comp, (33, 33), "coordonnees")
                                        
                                        graphe_comp_aretes_sim = nx.Graph()
                                        graphe_comp_aretes_sim.add_nodes_from(graphe_comp.nodes(data=True))
                                        
                                        print(len(composante))
                                        
                                        with open(EXTENSION_PATH % rep + "fichier_csv_grandes_composantes_%s_%s.csv" % (round(i, 2), compteur), 'w') as fichier_csv:
                                            csvwriter = csv.writer(fichier_csv)
                                            csvwriter.writerow(["source", "target", "label"])
                                        
                                            for u, v, data in graphe_comp.edges(data=True) :
                                                if data["poids"] >= i :
                                                    graphe_comp_aretes_sim.add_edge(u, v, **data)
                                                    csvwriter.writerow([u, v, round(data["poids"], 2)])
                                        
                                        with open(EXTENSION_PATH % rep + "fichier_csv_grandes_composantes_noeuds_%s_%s.csv" % (round(i, 2), compteur), 'w') as fichier_csv:
                                            csvwriter = csv.writer(fichier_csv)
                                            csvwriter.writerow(["id", "label", "homologues", "clustering_kmeans", "clustering_dbscan", "clustering_agglo"])
                                            
                                            if i == 0.7 or i == 0.6 :
                                                tab_clustering_kmeans = clustering_composante(rep, depart_sim, typ_sim, taille_ext, i, 'kmeans')
                                                tab_clustering_dbscan = clustering_composante(rep, depart_sim, typ_sim, taille_ext, i, 'dbscan')
                                                tab_clustering_agglo = clustering_composante(rep, depart_sim, typ_sim, taille_ext, i, 'agglo')
                                                print("nombre de noeuds")
                                                print(graphe_comp_aretes_sim.number_of_nodes())
                                                
                                            for noeud, data in graphe_comp_aretes_sim.nodes(data=True) :
                                                
                                                numero_homologues = 0
                                                for j in range(len(homologues)) :
                                                    
                                                    if data["nom"] in homologues[j] :
                                                        numero_homologues = j + 1
                                                        if data["nom"] == '5DM6_X_328_2' and i == 0.7 :
                                                            print("ramou")
                                                
                                                numero_clustering_kmeans = -1   
                                                numero_clustering_dbscan = -1         
                                                numero_clustering_agglo = -1     
                                                if i == 0.7 or i == 0.6:  
                                                    for j in range(len(tab_clustering_kmeans)) :
                                                        if data["nom"] in tab_clustering_kmeans[j] :
                                                            numero_clustering_kmeans = j
                                                    
                                                    for j in range(len(tab_clustering_dbscan)) :
                                                        if data["nom"] in tab_clustering_dbscan[j] :
                                                            numero_clustering_dbscan = j
                                                    
                                                    for j in range(len(tab_clustering_agglo)) :
                                                        if data["nom"] in tab_clustering_agglo[j] :
                                                            numero_clustering_agglo = j
                                                
                                                csvwriter.writerow([noeud, data["nom"], numero_homologues, numero_clustering_kmeans, numero_clustering_dbscan, numero_clustering_agglo])
                                        
                                        compteur += 1
                                    
                            i = i + (0.1 * max_val)
    else :
        with open("sim_" + depart_sim + "_" + typ_sim + "par_chaine.pickle", 'rb') as fichier_sim :
            mon_depickler = pickle.Unpickler(fichier_sim)
            dico_sim = mon_depickler.load()
            print(dico_sim.keys())
            
            with open("graphe_complet_pondere_sim.pickle", 'rb') as fichier_graphe_complet :
                mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
                graphe_complet = mon_depickler_complet.load()
                
                i = 0.1
                while i <= 1.0 :
                    with open("composantes_connexes_" + type_comp + "_" + str(i) + ".pickle", 'rb') as fichier_comp:
                        mon_depickler_comp = pickle.Unpickler(fichier_comp)
                        composantes_connexes = mon_depickler_comp.load()
                        
                        compteur = 1
                        for composante in composantes_connexes :
                            
                            if len(composante) < 10 and len(composante) > 2 :
                                graphe_comp = nx.Graph()
                                graphe_comp = graphe_complet.subgraph(composante)
                                nx.set_node_attributes(graphe_comp, (33, 33), "coordonnees")
                                pos = nx.circular_layout(graphe_comp)
                                
                                node_labels = dict([(u, (d["nom"]))for u, d in graphe_comp.nodes(data=True)])
        
                                for u, v, d in graphe_comp.edges(data=True) :
                                    
                                    if (graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]) in dico_sim.keys() :
                                        for j in range(4) :
                                            print(dico_sim[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])][j])
                                        graphe_comp.edges[u, v]["poids"] = (dico_sim[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])][0][0] + dico_sim[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])][1][0] + dico_sim[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])][2][0] + dico_sim[(graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"])][3][0]) / 4
                                    else :
                                        graphe_comp.edges[u, v]["poids"] = (dico_sim[(graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"])][0][0] + dico_sim[(graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"])][1][0] + dico_sim[(graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"])][2][0] + dico_sim[(graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"])][3][0]) / 4
                                    
                                edge_labels = dict([((u, v), (round(d["poids"], 2)))for u, v, d in graphe_comp.edges(data=True)])
                                
                                red_edges = []
                                
                                for u, v, d in graphe_comp.edges(data=True) :
                                    if d["poids"] >= i :
                                        red_edges.append((u, v))
                                
                                edge_colors = ['red' if edge in red_edges else 'black' for edge in graphe_comp.edges()]
                                densite = len(red_edges) / ((len(composante) * (len(composante) - 1)) / 2)
                                
                                plt.figure(figsize=(9, 7))
                                plt.title(depart_sim)
                                # plt.title("poids > "+str(round(i,1)) + " densite : "+ str(round(densite,2)))
                                plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.9)
                                
                                nx.draw_networkx_nodes(graphe_comp, pos)
                                nx.draw_networkx_labels(graphe_comp, pos, labels=node_labels, font_size=8)
                                nx.draw_networkx_edges(graphe_comp, pos, edge_color=edge_colors)
                                nx.draw_networkx_edge_labels(graphe_comp, pos, edge_labels=edge_labels, label_pos=0.3)
                                plt.axis('off')
                                
                                # plt.show()
                                plt.savefig("composantes_connexes_" + depart_sim + "_" + typ_sim + "par_chaine/comp_" + type_comp + "_" + str(round(i, 1)) + "_" + str(compteur) + ".png")  # save as png
                                plt.close()
                                compteur += 1
                            
                    i = i + 0.1


''' affiche les grandes composantes connexes (de taille > 10) associees a un seuil i
stocke les informations dans des fichiers csv utilisables par Gephi
calcule le clustering DBSCAN, kmeans et agglomerative clustering pour chaque composante et stocke l'information dans le fichier csv des noeuds ''' 
def draw_une_grande_composante(typ_sim, depart_sim, taille_ext, rep, i):

    with open(EXTENSION_PATH % rep + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (typ_sim, taille_ext), 'rb') as fichier_graphe_complet :
                        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
                        graphe_complet = mon_depickler_complet.load()
                      
                        with open(EXTENSION_PATH % rep + "composantes_connexes_" + depart_sim + "_" + typ_sim + "_taille_" + str(taille_ext) + "_" + str(round(i, 3)) + ".pickle", 'rb') as fichier_comp:
                            # with open("Extensions/Metrique_toutes_aretes/composantes_connexes_extensions_toutes_aretes_coeffn1_a1_c1_0.6_0.6.pickle", 'rb') as fichier_comp:
                            
                                mon_depickler_comp = pickle.Unpickler(fichier_comp)
                                composantes_connexes = mon_depickler_comp.load()
                                
                                compteur = 1
                                for composante in composantes_connexes :
                                    
                                    if len(composante) > 8  :
                                        graphe_comp = nx.Graph()
                                        graphe_comp = graphe_complet.subgraph(composante)
                                        nx.set_node_attributes(graphe_comp, (33, 33), "coordonnees")
                                        
                                        graphe_comp_aretes_sim = nx.Graph()
                                        graphe_comp_aretes_sim.add_nodes_from(graphe_comp.nodes(data=True))
                                        
                                        with open(EXTENSION_PATH % rep + "fichier_csv_grandes_composantes_%s_%s.csv" % (round(i, 3), compteur), 'w') as fichier_csv:
                                            csvwriter = csv.writer(fichier_csv)
                                            csvwriter.writerow(["source", "target", "label"])
                                        
                                            for u, v, data in graphe_comp.edges(data=True) :
                                                if data["poids"] >= i :
                                                    graphe_comp_aretes_sim.add_edge(u, v, **data)
                                                    csvwriter.writerow([u, v, round(data["poids"], 2)])
                                        
                                        with open(EXTENSION_PATH % rep + "fichier_csv_grandes_composantes_noeuds_%s_%s.csv" % (round(i, 3), compteur), 'w') as fichier_csv:
                                            csvwriter = csv.writer(fichier_csv)
                                            csvwriter.writerow(["id", "label", "homologues", "clustering_kmeans", "clustering_dbscan", "clustering_agglo"])
                                            
                                            if i == 0.7 or i == 0.6 :
                                                tab_clustering_kmeans = clustering_composante(rep, depart_sim, typ_sim, taille_ext, i, 'kmeans')
                                                tab_clustering_dbscan = clustering_composante(rep, depart_sim, typ_sim, taille_ext, i, 'dbscan')
                                                tab_clustering_agglo = clustering_composante(rep, depart_sim, typ_sim, taille_ext, i, 'agglo')
                                                print("nombre de noeuds")
                                                print(graphe_comp_aretes_sim.number_of_nodes())
                                                
                                            for noeud, data in graphe_comp_aretes_sim.nodes(data=True) :
                                                
                                                numero_homologues = 0
                                                for j in range(len(HOMOLOGUES)) :
                                                    
                                                    if data["nom"] in HOMOLOGUES[j] :
                                                        numero_homologues = j + 1
                                                        if data["nom"] == '5DM6_X_328_2' and i == 0.7 :
                                                            print("ramou")
                                                
                                                numero_clustering_kmeans = -1   
                                                numero_clustering_dbscan = -1         
                                                numero_clustering_agglo = -1     
                                                if i == 0.7 or i == 0.6:  
                                                    for j in range(len(tab_clustering_kmeans)) :
                                                        if data["nom"] in tab_clustering_kmeans[j] :
                                                            numero_clustering_kmeans = j
                                                    
                                                    for j in range(len(tab_clustering_dbscan)) :
                                                        if data["nom"] in tab_clustering_dbscan[j] :
                                                            numero_clustering_dbscan = j
                                                    
                                                    for j in range(len(tab_clustering_agglo)) :
                                                        if data["nom"] in tab_clustering_agglo[j] :
                                                            numero_clustering_agglo = j
                                                
                                                csvwriter.writerow([noeud, data["nom"], numero_homologues, numero_clustering_kmeans, numero_clustering_dbscan, numero_clustering_agglo])
                                        
                                        compteur += 1

''' renvoie une matrice de distance associee aux elements de la composante passee en parametre
pour les paires de graphe ayant une arete, on associe la valeur de distance 1 - sim
pour les paires de graphe n'ayant aucune arete, on associe la valeur de distance 50
rep, type_sim, taille_ext : infos pour retrouver le nom du fichier dans lequel est stocke le graphe complet pour avoir acces aux valeurs de sim'''
def construction_matrice_distance(composante, rep, typ_sim, taille_ext, seuil):
    with open(EXTENSION_PATH % rep + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (typ_sim, taille_ext), 'rb') as fichier_graphe_complet :
        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
        graphe_complet = mon_depickler_complet.load()
        
        compteur = 0
        graphe_comp_aretes_sim = nx.Graph()
        print(len(composante))
        for elt in composante :
            graphe_comp_aretes_sim.add_node(compteur, **graphe_complet.nodes[elt])
            compteur += 1
        
        print(graphe_comp_aretes_sim.nodes.data())
        print(graphe_complet.nodes.data())
        for u, v, data in graphe_complet.edges(data=True) :
            # if data["poids"] >= seuil :
                u_new = -1
                v_new = -1
                for noeud, data_noeud in graphe_comp_aretes_sim.nodes(data=True) :
#                     print(data_noeud)
#                     print(graphe_complet.nodes[u])
                    if graphe_complet.nodes[u]["nom"] == data_noeud["nom"] :
                        u_new = noeud
                    if graphe_complet.nodes[v]["nom"] == data_noeud["nom"] :
                        v_new = noeud
                if u_new != -1  and v_new != -1 :
                    graphe_comp_aretes_sim.add_edge(u_new, v_new, **data)
                
        matrice = [[0] * graphe_comp_aretes_sim.number_of_nodes() for _ in range(graphe_comp_aretes_sim.number_of_nodes())]
        print(matrice)
        
        compteur = 0
        for i in range(graphe_comp_aretes_sim.number_of_nodes()) :
            for j in range(graphe_comp_aretes_sim.number_of_nodes()) :
                if i != j :
                    if (i, j) in graphe_comp_aretes_sim.edges() or (j, i) in graphe_comp_aretes_sim.edges() :
                        matrice[i][j] = 1 - graphe_comp_aretes_sim.edges[i, j]["poids"]
                        compteur += 1
                    else :
                        matrice[i][j] = 50.0
                
        print(matrice)
        print(graphe_comp_aretes_sim.number_of_edges())
        print(compteur)
        
        # print(matrice[0][32])
        return matrice, graphe_comp_aretes_sim

''' calcule une methode de clustering (passee en parametre : dbscan, kmeans, ou agglo) sur les grandes composantes connexes (taille > 10) associees a un seuil
pour cela, cree la matrice de distance, appelle les fonctions de clustering du package sklearn et renvoie un tableau a deux dimensions dans lequel chaque ligne correspond à un cluster
rep, depart_sim, type_sim, taille_ext : pour retrouver le fichier ou sont stockes les composantes 
'''
def clustering_composante(rep, depart_sim, typ_sim, taille_ext, seuil, methode_clustering):

    with open(EXTENSION_PATH % rep + "composantes_connexes_" + depart_sim + "_" + typ_sim + "_taille_" + str(taille_ext) + "_" + str(round(seuil, 2)) + ".pickle", 'rb') as fichier_comp:
        # with open("Extensions/Metrique_toutes_aretes/composantes_connexes_extensions_toutes_aretes_coeffn1_a1_c1_0.6_0.6.pickle", 'rb') as fichier_comp:
        
            mon_depickler_comp = pickle.Unpickler(fichier_comp)
            composantes_connexes = mon_depickler_comp.load()
            
            for composante in composantes_connexes :
                
                if len(composante) > 10  :
                    matrice, graphe_comp_sim_aretes = construction_matrice_distance(composante, rep, typ_sim, taille_ext, seuil)
                    
                    # init = np.array([matrice[16], matrice[28], matrice[39], matrice[44], matrice[53]], np.float64)
                    # init = init.reshape(-1,1)
                    if methode_clustering == 'kmeans' :
                        print("kmeans")
                        clustering = KMeans(n_clusters=7, n_init=50).fit(matrice)  # init=init, n_init=1).fit(matrice)
                    # print(kmeans.labels_)
                    
                    if methode_clustering == "dbscan" :
                        print("dbscan")
                        clustering = DBSCAN(min_samples=2, eps=0.277, metric='precomputed').fit(matrice)
                    
                    if methode_clustering == "agglo" :
                        print("agglo")
                        clustering = AgglomerativeClustering(n_clusters=7).fit(matrice)
                    
                    tab_clustering = [['']]
                    print(len(tab_clustering))
                    print(len(clustering.labels_))
                    print(clustering.labels_)
                    compteur = 0
                    for elt in clustering.labels_ :
                        if elt != -1 :
                            for _ in range(len(tab_clustering), elt + 1) :
                                    tab_clustering.append([''])
                            tab_clustering[elt].append(graphe_comp_sim_aretes.nodes[compteur]["nom"])
                        compteur += 1
                        
                    for i in range(len(tab_clustering)) :
                        print(len(tab_clustering[i]) - 1)
                        print(tab_clustering[i])
                    # print(dbscan.cluster_centers_)
    return tab_clustering


''' dans le graphe complet stocke dans le fichier graphe_complet_pondere_sim.pickle, stocke les valeurs de similarites de la methode passee en parametre
(non_cov_sans_motif ou avec_poids)
pas terrible comme fonction '''
def choix_sim_stockee(choix):
        with open("graphe_complet_pondere_sim.pickle", 'rb') as fichier_graphe_complet :
            mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
            graphe_complet = mon_depickler_complet.load()
            
            for u, v in graphe_complet.edges() :
                elt1 = graphe_complet.nodes[u]["nom"]
                elt2 = graphe_complet.nodes[v]["nom"]
                if "4V9F_0_25_4" not in elt1 and "4V9F_0_25_4" not in elt2 :
                    print(elt1)
                    print(elt2)
                    with open("graphes_extension/fichier_" + elt1 + ".pickle", 'rb') as fichier_graphe1 :
                        mon_depickler_graphe1 = pickle.Unpickler(fichier_graphe1)
                        graphe1 = mon_depickler_graphe1.load()
                        with open("graphes_extension/fichier_" + elt2 + ".pickle", 'rb') as fichier_graphe2 :
                            mon_depickler_graphe2 = pickle.Unpickler(fichier_graphe2)
                            graphe2 = mon_depickler_graphe2.load()
                            
                            with open("dico_graphe_epure_en_tout_test.pickle", 'rb') as fichier_commun :
                                mon_depickler_commun = pickle.Unpickler(fichier_commun)
                                dico_graphe = mon_depickler_commun.load()
                                
                                if choix == "non_cov_sans_motif" :
                                    if ("fichier_" + elt1, "fichier_" + elt2) in dico_graphe.keys() :
                                        sim = calcul_sim_non_cov_sans_motif(graphe1, graphe2, dico_graphe[("fichier_" + elt1, "fichier_" + elt2)])
                                    else :
                                        sim = calcul_sim_non_cov_sans_motif(graphe2, graphe1, dico_graphe[("fichier_" + elt2, "fichier_" + elt1)])
                                else :
                                    if ("fichier_" + elt1, "fichier_" + elt2) in dico_graphe.keys() :
                                        sim = calcul_sim_avec_poids(graphe1, graphe2, dico_graphe[("fichier_" + elt1, "fichier_" + elt2)])
                                    else :
                                        sim = calcul_sim_avec_poids(graphe2, graphe1, dico_graphe[("fichier_" + elt2, "fichier_" + elt1)])
                                        
                                graphe_complet.edges[u, v]["poids"] = sim
                            
        with open("graphe_complet_pondere_sim.pickle", 'wb') as fichier_graphe_complet_refait :
            mon_pickler_complet = pickle.Pickler(fichier_graphe_complet_refait)
            mon_pickler_complet.dump(graphe_complet)                    
    
''' ne sais pas trop a quoi ca sert'''
def calcul_sim_globale_par_groupe(composantes_connexes, graphe_complet_pondere):
    with open("dico_graphe_epure_en_tout.pickle", 'rb') as fichier_graphe :
        mon_depickler = pickle.Unpickler(fichier_graphe)
        dico_graphe = mon_depickler.load()
    
        for composante in composantes_connexes :
            if len(composante) < 10 and len(composante) > 2 :
                graphe_comp = graphe_complet_pondere.subgraph(composante)
                
                max_nombre_arc = -1
                max_nombre_arete = -1
                max_nombre_sommet = -1
                for noeud, data in graphe_comp.nodes(data=True) :
                    with open("graphes_extension/" + data["nom"] + ".pickle", 'rb') as fichier :
                        mon_depickler = pickle.Unpickler(fichier)
                        graphe = mon_depickler.load()
                        
                        nombre_arc, nombre_arete = nombre_arc_arete_graphe(graphe)
                        if nombre_arc < max_nombre_arc and nombre_arete < max_nombre_arete and max_nombre_sommet < graphe.number_of_nodes() :
                            max_nombre_arc = nombre_arc
                            max_nombre_arete = nombre_arete
                            max_nombre_sommet = graphe.number_of_nodes()
                
                aretes_idem = []
                print(graphe_comp.edges())            
                for u, v in graphe_comp.edges() :
                    nom_1 = graphe.nodes[u]["nom"]
                    nom_2 = graphe.nodes[v]["nom"]
                    if (nom_1, nom_2) in dico_graphe.keys() :
                        cle_1 = (nom_1, nom_2)
                    else :
                        cle_1 = (nom_2, nom_1)
                        
                    for w, z in graphe_comp.edges() :
                        nom_3 = graphe.nodes[w]["nom"]
                        nom_4 = graphe.nodes[z]["nom"]
                        if (nom_3, nom_4) in dico_graphe.keys() :
                            cle_2 = (nom_3, nom_4)
                        else :
                            cle_2 = (nom_4, nom_3) 
                            
                        if cle_1 != cle_2 :
                            for e1, e2, data_1 in dico_graphe[cle_1].edges(data=True) :
                                print("ramousnif")
                            
''' affiche et stocke dans un fichier .png un graphe complet pondere par les valeurs de sim correspondant a un groupe defini dans les parametres
rep, taille_ext, typ : pour retrouver le graphe complet total
groupe : liste des noms des elements du groupe
groupe_nom : nom du groupe
seuil : pour choisir quelles aretes seront dessinees en rouge
avec_val_k : pour les valeurs de sim, version en divisant le nombre d'aretes en commun par la taille de l'extension (mais c'est pas terrible en fait) ''' 
def draw_groupes_par_type(rep, taille_ext, depart, typ, avec_val_k, groupe, groupe_nom, seuil):
    with open(EXTENSION_PATH % rep + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (typ, taille_ext), 'rb') as fichier_graphe_complet :
        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
        graphe_complet = mon_depickler_complet.load()          
        
#         with open("graphe_complet_pondere_sim_grands_graphes_metrique_sommets_aretes.pickle", 'rb') as fichier_graphe_complet_grands_graphes :
#             mon_depickler_complet_2 = pickle.Unpickler(fichier_graphe_complet_grands_graphes)
#             graphe_complet_grands_graphes = mon_depickler_complet_2.load() 
        if avec_val_k :    
            fichier_sim_a_ouvrir = EXTENSION_PATH % rep + "sim_%s_%s_taille_%s_avec_val_k.pickle" % (depart, typ, taille_ext)
        else :
            fichier_sim_a_ouvrir = EXTENSION_PATH % rep + "sim_%s_%s_taille_%s.pickle" % (depart, typ, taille_ext)
        
        with open(fichier_sim_a_ouvrir, 'rb') as fichier_sim:

                mon_depickler_sim = pickle.Unpickler(fichier_sim)
                dico_sim = mon_depickler_sim.load()
                print(dico_sim)
                
                # print(graphe_complet.number_of_nodes())
                groupe_noeuds = []
                for noeud, nom in graphe_complet.nodes(data="nom") :
                    if nom in groupe :
                        groupe_noeuds.append(noeud)
                
                print(groupe_noeuds)        
                
                graphe_groupe = graphe_complet.subgraph(groupe_noeuds)
                pos = nx.circular_layout(graphe_groupe)
                
                node_labels = dict([(u, (d["nom"]))for u, d in graphe_groupe.nodes(data=True)])
            
#             for u,v,d in graphe_gnra.edges(data=True) :
#                 for u2,v2,d2 in graphe_complet_grands_graphes.edges(data=True) :
#                     nom_1 = (graphe_complet.nodes[u]["nom"].split("_")[0], graphe_complet.nodes[u]["nom"].split("_")[1], int(graphe_complet.nodes[u]["nom"].split("_")[2]), int(graphe_complet.nodes[u]["nom"].split("_")[3]))
#                     nom_2 = (graphe_complet.nodes[v]["nom"].split("_")[0], graphe_complet.nodes[v]["nom"].split("_")[1], int(graphe_complet.nodes[v]["nom"].split("_")[2]), int(graphe_complet.nodes[v]["nom"].split("_")[3]))
#     
#                     if (graphe_complet_grands_graphes.nodes[u2]["nom"] == nom_1 and graphe_complet_grands_graphes.nodes[v2]["nom"] == nom_2) or (graphe_complet_grands_graphes.nodes[v2]["nom"] == nom_1 and graphe_complet_grands_graphes.nodes[u2]["nom"] == nom_2) :
#                         graphe_gnra.edges[u,v]["poids"] = d2["poids"]
            
                for u, v, d in graphe_groupe.edges(data=True) :
                    for cle in dico_sim.keys() :
    #                     for i in range(4) : 
    #                         dico_sim[cle][i] = round(dico_sim[cle][i], 2)
                        if (graphe_groupe.nodes[u]["nom"] == cle[0] and graphe_groupe.nodes[v]["nom"] == cle[1]) or (graphe_groupe.nodes[u]["nom"] == cle[1] and graphe_groupe.nodes[v]["nom"] == cle[0]) :
                            if len(dico_sim[cle]) == 1 : 
                                graphe_groupe.edges[u, v]["poids"] = dico_sim[cle]
                            else :
                                graphe_groupe.edges[u, v]["poids"] = dico_sim[cle][0]
                                graphe_groupe.edges[u, v]["k"] = dico_sim[cle][1]
                print(dico_sim)  
                print(graphe_groupe.edges.data())  
                
                red_colors = []
                
    #             for noeud, data in graphe_gnra.nodes(data=True) :
    #                 if data["nom"] in ['5J5B_BA_48_23', '5J7L_DA_48_15', '1FJG_A_48_8', '5FDU_1A_48_19', '1U9S_A_58_11'] :
    #                     red_colors.append(noeud)
    #             node_colors = ['red' if noeud in red_colors else 'blue' for noeud in graphe_gnra.nodes()]
                 
                # edge_labels = dict([((u,v), (round(d["poids"])))for u,v,d in graphe_gnra.edges(data=True)])
    
                edge_labels = dict([((u, v), (d["k"], round(d["poids"], 2)))for u, v, d in graphe_groupe.edges(data=True)])
                red_edges = []
                for u, v, d in graphe_groupe.edges(data=True) :
                    if d["poids"] >= seuil :
                        red_edges.append((u, v))
                
                densite = len(red_edges) / graphe_groupe.number_of_edges()
                edge_colors = ['red' if edge in red_edges else 'black' for edge in graphe_groupe.edges()]
                
                plt.figure(figsize=(9, 7))
                # plt.title("GNRA extensions αN=1 αCI=1 αA=1 taille 5")
                plt.title(groupe_nom + " densite : %f" % (round(densite, 2)))
                plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.9)
                
                nx.draw_networkx_nodes(graphe_groupe, pos, node_color="red")
                nx.draw_networkx_labels(graphe_groupe, pos, labels=node_labels, font_size=8)
                nx.draw_networkx_edges(graphe_groupe, pos, edge_color=edge_colors)
                nx.draw_networkx_edge_labels(graphe_groupe, pos, edge_labels=edge_labels, label_pos=0.4)
                plt.axis('off')
                
                plt.savefig(EXTENSION_PATH % rep + "%s.png" % groupe_nom)  # save as png
        
                # plt.show()
                plt.close()                       


''' affiche dans la console des infos sur les composantes connexes contenant au moins un element du type cherche (GNRA, A-rich ou G-bulge) '''
def retrouver_types(depart, type_sim, rep, taille_ext, val_max, type_cherche):
    i = 0.1 * val_max
    
    if type_cherche == "GBULGE" :
        groupe = GROUPE_GBULGE
    elif type_cherche == 'GNRA' :
        groupe = GROUPE_GNRA
    elif type_cherche == 'ARICH' :
        groupe = GROUPE_ARICH
    
    while i < 1.0 + 0.00000005 :
            with open(EXTENSION_PATH % rep + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (type_sim, str(taille_ext)), 'rb') as fichier_graphe_complet_1 :
                mon_depickler_complet_1 = pickle.Unpickler(fichier_graphe_complet_1)
                graphe_complet = mon_depickler_complet_1.load()
                
                with open(EXTENSION_PATH % rep + "composantes_connexes_%s_%s_taille_%s_%s.pickle" % (depart, type_sim, str(taille_ext), str(round(i, 2))), 'rb') as fichier_comp:
                    mon_depickler_comp = pickle.Unpickler(fichier_comp)
                    composantes_connexes = mon_depickler_comp.load()
                    
                    num = 0
                    for composante in composantes_connexes :
                        compteur = 0
                        for elt in composante : 
                            if graphe_complet.nodes[elt]["nom"] in groupe :
                                compteur += 1
                                print(graphe_complet.nodes[elt]["nom"])
                        if compteur > 0 :
                            print("i= " + str(i))
                            print("nombre du type = " + str(compteur))
                            print("numero composante = " + str(num))
                            print("taille composante = " + str(len(composante)))
                            
                        num += 1
            i = i + 0.1
                        
''' affiche dans la console des infos sur les composantes connexes contenant au moins un element dans lequel il y a un lien entre la chaine 1 et la chaine 3 '''
def retrouver_lien_chaines_1_3(depart, type_sim, rep, taille_ext, val_max):
    i = 0.1 * val_max
    
    while i < 1.0 + 0.00000005 :
            with open(EXTENSION_PATH % rep + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (type_sim, str(taille_ext)), 'rb') as fichier_graphe_complet_1 :
                mon_depickler_complet_1 = pickle.Unpickler(fichier_graphe_complet_1)
                graphe_complet = mon_depickler_complet_1.load()
                
                with open(EXTENSION_PATH % rep + "composantes_connexes_%s_%s_taille_%s_%s.pickle" % (depart, type_sim, str(taille_ext), str(round(i, 2))), 'rb') as fichier_comp:
                    mon_depickler_comp = pickle.Unpickler(fichier_comp)
                    composantes_connexes = mon_depickler_comp.load()
                    
                    num = 0
                    for composante in composantes_connexes :
                        lien_chaine_1_3 = -1
                        print("i= " + str(i))
                        print("numero composante = " + str(num))
                        print("taille composante = " + str(len(composante)))
                        for elt in composante : 
                            with open(EXTENSION_PATH_TAILLE % 10 + "fichier_%s.pickle" % (graphe_complet.nodes[elt]["nom"]), 'rb') as fichier_graphe:
                                mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
                                graphe = mon_depickler_graphe.load()
                                lien_chaine_1_3 = lien_chaines_1_3(graphe)
                                
                                if lien_chaine_1_3 :
                                    print("nom = " + str(graphe_complet.nodes[elt]["nom"]))
                            
                        num += 1
            i = i + 0.1

''' affiche dans la console (et stocke dans un fichier csv) les noms des graphes de chaque composante connexe et le nombre de composantes differentes '''
def afficher_groupes(depart, typ, rep, taille_ext, val_max):                       
    i = 1.0 * val_max
    compteur = 0
    
    meme_groupe_deja_vu = []
    
    with open(EXTENSION_PATH % rep + "fichier_csv_groupes_differents_par_seuil.csv", 'w', newline="") as fichier_csv :
        csvwriter = csv.writer(fichier_csv)
    
        while i > 0.0 + 0.00000005 :
            with open(EXTENSION_PATH % rep + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (typ, str(taille_ext)), 'rb') as fichier_graphe_complet_1 :
                mon_depickler_complet_1 = pickle.Unpickler(fichier_graphe_complet_1)
                graphe_complet = mon_depickler_complet_1.load()
                
                with open(EXTENSION_PATH % rep + "composantes_connexes_%s_%s_taille_%s_%s.pickle" % (depart, typ, str(taille_ext), str(round(i, 2))), 'rb') as fichier_comp:
                    mon_depickler_comp = pickle.Unpickler(fichier_comp)
                    composantes_connexes = mon_depickler_comp.load()
                    
                    for composante in composantes_connexes :
                        
                        if len(composante) < 10 and len(composante) > 2 :
                            deja_vu = False
                            for groupe in meme_groupe_deja_vu :
                                nb_deja_vu = 0
                                for elt2 in composante :
        #                             print(groupe)
        #                             print(elt2)
                                    if elt2 in groupe :
                                        nb_deja_vu += 1
                                if nb_deja_vu == len(composante) and nb_deja_vu == len(groupe) :  
                                    deja_vu = True
        #                         print(nb_deja_vu)
        #                     print(deja_vu)
                            if deja_vu == False :
                                
                                print("composante numero : %d i : %f" % (compteur + 1, i / val_max))
                                str_comp = ""
                                for elt in composante :
                                    print(graphe_complet.nodes[elt]["nom"])
                                    str_comp += "\n" + graphe_complet.nodes[elt]["nom"] 
                                csvwriter.writerow([round(i / val_max, 2), str_comp])
                                meme_groupe_deja_vu.append(composante)
                                compteur += 1 
                            
            i = i - (0.1 * val_max)        
    print("Nombre de composantes differentes : " + str(compteur))

''' renvoie une liste dans laquelle sont stockees les elements du groupe passe en parametre + leurs voisins si les aretes sont au-dessus du seuil passe en parametre '''
def calcul_groupe_etendu_seuil_sup(typ, rep, taille_ext, groupe, seuil):
    new_groupe = []
    with open(EXTENSION_PATH % rep + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (typ, str(taille_ext)), 'rb') as fichier_graphe_complet :
        mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
        graphe_complet = mon_depickler_complet.load()
        
        for noeud, data in graphe_complet.nodes(data=True) :
            if data["nom"] in groupe :
                if data["nom"] not in new_groupe :
                    new_groupe.append(data["nom"])
                    
                for voisin in graphe_complet[noeud] :
                    if graphe_complet.edges[noeud, voisin]["poids"] > seuil :
                        if graphe_complet.nodes[voisin]["nom"] not in new_groupe :
                            new_groupe.append(graphe_complet.nodes[voisin]["nom"])
                            
    return new_groupe
                    
''' cree les fichiers csv pour Gephi pour une grande composante connexe dans le nom est passe en parametre
groupe : liste des elements du groupe
groupe_base : liste des elements du groupe appartenant au groupe directement (sans passer par le voisinage a un seuil plus bas, voir fonction precedente) 
seuil : ne stocke que les aretes au-dessus de ce seuil
'''
def draw_groupe_gephi_unique(groupe, groupe_base, rep, typ, taille_ext, seuil, nom_groupe):
    with open(EXTENSION_PATH % rep + "fichier_csv_grandes_composantes_%s.csv" % (nom_groupe), 'w') as fichier_csv:
        csvwriter = csv.writer(fichier_csv)
        csvwriter.writerow(["source", "target", "label"])
        
        with open(EXTENSION_PATH % rep + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (typ, str(taille_ext)), 'rb') as fichier_graphe_complet :
            mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
            graphe_complet = mon_depickler_complet.load()
            
            for u, v, data in graphe_complet.edges(data=True) :
                if graphe_complet.nodes[u]["nom"] in groupe and graphe_complet.nodes[v]["nom"] in groupe :
                    if data["poids"] > seuil :
                        csvwriter.writerow([u, v, round(data["poids"], 2)])
            
            with open(EXTENSION_PATH % rep + "fichier_csv_grandes_composantes_noeuds_%s.csv" % (nom_groupe), 'w') as fichier_csv:
                csvwriter = csv.writer(fichier_csv)
                csvwriter.writerow(["id", "label", "homologues", "dans_groupe"])
                    
                for noeud, data in graphe_complet.nodes(data=True) :
                    if data["nom"] in groupe :
                        numero_homologues = 0
                        for j in range(len(HOMOLOGUES)) :
                            
                            if data["nom"] in HOMOLOGUES[j] :
                                numero_homologues = j + 1
                                
                        if data["nom"] in groupe_base :
                            dans_groupe = 1
                        else :
                            dans_groupe = 0
                        csvwriter.writerow([noeud, data["nom"], numero_homologues, dans_groupe])

''' affiche dans la console le nombre de composantes avec des elements differents '''   
def compter_composantes_differentes(depart, typ, rep, taille_ext, val_max):
    i = 0.1 * val_max
    compteur = 0
    
    meme_groupe_deja_vu = []
    while i <= 1.0 * val_max :
        with open(EXTENSION_PATH % rep + "composantes_connexes_%s_%s_taille_%s_%s.pickle" % (depart, typ, str(taille_ext), str(round(i, 2))), 'rb') as fichier_comp:
            mon_depickler_comp = pickle.Unpickler(fichier_comp)
            composantes_connexes = mon_depickler_comp.load()
            
            for composante in composantes_connexes :
                
                if len(composante) < 10 and len(composante) > 2 :
                    deja_vu = False
                    for groupe in meme_groupe_deja_vu :
                        nb_deja_vu = 0
                        for elt2 in composante :
#                             print(groupe)
#                             print(elt2)
                            if elt2 in groupe :
                                nb_deja_vu += 1
                        if nb_deja_vu == len(composante) :  
                            deja_vu = True
#                         print(nb_deja_vu)
#                     print(deja_vu)
                    if deja_vu == False :
                        meme_groupe_deja_vu.append(composante)
                        compteur += 1 
                        
            i = i + (0.1 * val_max)        
    print("Nombre de composantes differentes : " + str(compteur))

''' affiche dans la console le nombre de composantes differentes dans un groupe de composantes passees en parametre  '''
def compter_composantes_differentes_par_groupe(groupe):

    compteur = 0
    
    meme_groupe_deja_vu = []

    for composante in groupe :
                
                    deja_vu = False
                    for groupe in meme_groupe_deja_vu :
                        nb_deja_vu = 0
                        for elt2 in composante :
#                             print(groupe)
#                             print(elt2)
                            if elt2 in groupe :
                                nb_deja_vu += 1
                        if nb_deja_vu == len(composante) :  
                            deja_vu = True
#                         print(nb_deja_vu)
#                     print(deja_vu)
                    if deja_vu == False :
                        meme_groupe_deja_vu.append(composante)
                        compteur += 1 
      
    print("Nombre de composantes differentes : " + str(compteur))

                                     
'''comparaison des clusters obtenus avec les tailles d'extension 4 et 10 (petites composantes connexes de taille comprise entre 3 et 9) 
fait des tas de trucs compliques pour les comparer'''
def comparaison_cluster(taille_ext_1, taille_ext_2, rep1, rep2, depart1, typ1, depart2, typ2, max_val_1, max_val_2):
    tab_idem = []
    nombre_comp_1 = [0] * 10
    nombre_comp_2 = [0] * 10
    tab_inclus_dans_1 = []
    tab_inclus_dans_2 = []
    comp_1 = []
    comp_2 = []
    with open(EXTENSION_PATH % rep1 + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (typ1, str(taille_ext_1)), 'rb') as fichier_graphe_complet_1 :
            mon_depickler_complet_1 = pickle.Unpickler(fichier_graphe_complet_1)
            graphe_complet_1 = mon_depickler_complet_1.load()
                        
            with open(EXTENSION_PATH % rep2 + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (typ2, str(taille_ext_2)), 'rb') as fichier_graphe_complet_2 :
                mon_depickler_complet_2 = pickle.Unpickler(fichier_graphe_complet_2)
                graphe_complet_2 = mon_depickler_complet_2.load()
                        
                i = 0.1
                while i <= 1.0 :
                            
                            with open(EXTENSION_PATH % rep1 + "composantes_connexes_%s_%s_taille_%s_%s.pickle" % (depart1, typ1, str(taille_ext_1), str(round(i * max_val_1, 2))), 'rb') as fichier_comp_1:
                                mon_depickler_comp_1 = pickle.Unpickler(fichier_comp_1)
                                composantes_connexes_1 = mon_depickler_comp_1.load()
                                
                                nombre_1 = 1
                                for composante_1 in composantes_connexes_1 :
                                            
                                        if len(composante_1) < 10 and len(composante_1) > 2 :
                                            graphe_comp_1 = nx.Graph()
                                            graphe_comp_1 = graphe_complet_1.subgraph(composante_1)
                                            
                                            comp_1.append(composante_1)
                                            j = 0.1
                                            while j <= 1.0 :
                                                    
                                                with open(EXTENSION_PATH % rep2 + "composantes_connexes_%s_%s_taille_%s_%s.pickle" % (depart2, typ2, str(taille_ext_2), str(round(j * max_val_2, 2))), 'rb') as fichier_comp_2:
                                                    mon_depickler_comp_2 = pickle.Unpickler(fichier_comp_2)
                                                    composantes_connexes_2 = mon_depickler_comp_2.load()
                                                    
                                                    nombre_2 = 1
                                                    for composante_2 in composantes_connexes_2 :
                                                            if len(composante_2) < 10 and len(composante_2) > 2 :
                                                                graphe_comp_2 = nx.Graph()
                                                                graphe_comp_2 = graphe_complet_2.subgraph(composante_2)
                                                                if composante_2 not in comp_2 :
                                                                    comp_2.append(composante_2)
                                                                compteur = 0
                                                                for noeud_1, data_1 in graphe_comp_1.nodes(data=True) :
                                                                    for noeud_2, data_2 in graphe_comp_2.nodes(data=True) :
                                                                        if data_1["nom"] == data_2["nom"] : 
                                                                            compteur += 1 
                                                                if compteur == len(composante_1) and compteur == len(composante_2) :
                                                                    # print("%f,%f"%(i,j))
                                                                    # print("%d, %d\n"%(nombre_1, nombre_2))
                                                                    tab_idem.append(((round(i, 1), nombre_1), (round(j, 1), nombre_2)))
                                                                else :
                                                                    if compteur == len(composante_1) :
                                                                        tab_inclus_dans_2.append(((round(i, 1), nombre_1), (round(j, 1), nombre_2)))
                                                                        
                                                                    if compteur == len(composante_2) :
                                                                        tab_inclus_dans_1.append(((round(i, 1), nombre_1), (round(j, 1), nombre_2)))
                                                                nombre_2 += 1
                                                    nombre_comp_2[round(j * 10) - 1] = nombre_2 - 1 
                                                j = j + 0.1
                                            nombre_1 += 1
                                        
                                nombre_comp_1[round(i * 10) - 1] = nombre_1 - 1       
                            i = i + 0.1
                                    
#                 print(nombre_comp_1)
#                 print(nombre_comp_2)
                tab_vu_1 = []
                tab_vu_2 = []
                for elt in tab_idem :
                    if elt[0] not in tab_vu_1 :
                        tab_vu_1.append(elt[0])
                    if elt[1] not in tab_vu_2 :
                        tab_vu_2.append(elt[1])
                
                tab_non_idem_1 = []
                tab_non_idem_2 = []        
                for i in range(10) :
                    for j in range(1, nombre_comp_1[i] + 1) :
                        est_la = False
                        for elt in tab_vu_1 :
                            if elt[0] == (i + 1) / 10 and elt[1] == j :
                                est_la = True
                        if est_la == False :
                            tab_non_idem_1.append(((i + 1) / 10, j))
                
                for i in range(10) :
                    for j in range(1, nombre_comp_2[i] + 1) :
                        est_la = False
                        for elt in tab_vu_2 :
                            if elt[0] == (i + 1) / 10 and elt[1] == j :
                                est_la = True
                        if est_la == False :
                            tab_non_idem_2.append(((i + 1) / 10, j))   
                print(tab_idem)
#                 print(tab_non_idem_1)
#                 print(tab_non_idem_2)     
                
                tab_vu_1 = []
                tab_vu_2 = []
                for elt in tab_inclus_dans_1 :
                    if elt[0] not in tab_vu_1 :
                        tab_vu_1.append(elt[0])
                    if elt[1] not in tab_vu_2 :
                        tab_vu_2.append(elt[1])

                for elt in tab_inclus_dans_2 :
                    if elt[0] not in tab_vu_1 :
                        tab_vu_1.append(elt[0])
                    if elt[1] not in tab_vu_2 :
                        tab_vu_2.append(elt[1])
                
                # print(tab_inclus_dans_1)
#                 print(tab_vu_1)
#                 print(tab_vu_2)        
                for i in range(10) :
                    for j in range(1, nombre_comp_1[i] + 1) :
                        est_la = False
                        for elt in tab_vu_1 :
                            if elt[0] == (i + 1) / 10 and elt[1] == j :
                                est_la = ((i + 1) / 10, j)
                        if est_la != False and est_la in tab_non_idem_1 :
                            # print(tab_non_idem_1)
                            # print(est_la)
                            tab_non_idem_1.remove(((i + 1) / 10, j))
                            # print(tab_non_idem_1)
                
                for i in range(10) :
                    for j in range(1, nombre_comp_2[i] + 1) :
                        est_la = False
                        for elt in tab_vu_2 :
                            if elt[0] == (i + 1) / 10 and elt[1] == j :
                                est_la = ((i + 1) / 10, j)
                        if est_la != False and est_la in tab_non_idem_2 :
                            tab_non_idem_2.remove(((i + 1) / 10, j))   
                            
#                 print(tab_inclus_dans_1)
#                 print(tab_inclus_dans_2)
                print(tab_non_idem_1)
                print(tab_non_idem_2)
                groupe_1 = []
                for elt in tab_non_idem_1 : 
                    with open(EXTENSION_PATH % rep1 + "composantes_connexes_%s_%s_taille_%s_%s.pickle" % (depart1, typ1, str(taille_ext_1), str(round(elt[0] * max_val_1, 2))), 'rb') as fichier_comp_1:
                        mon_depickler_comp_1 = pickle.Unpickler(fichier_comp_1)
                        composantes_connexes_1 = mon_depickler_comp_1.load()
                        
                        nombre_1 = 1
                        for composante_1 in composantes_connexes_1 :
                                            
                            if len(composante_1) < 10 and len(composante_1) > 2 :
                                if nombre_1 == elt[1]: 
                                    groupe_1.append(composante_1)
                                nombre_1 += 1 
                groupe_2 = []
                for elt in tab_non_idem_2 : 
                    with open(EXTENSION_PATH % rep2 + "composantes_connexes_%s_%s_taille_%s_%s.pickle" % (depart2, typ2, str(taille_ext_2), str(round(elt[0] * max_val_2, 2))), 'rb') as fichier_comp_2:
                        mon_depickler_comp_2 = pickle.Unpickler(fichier_comp_2)
                        composantes_connexes_2 = mon_depickler_comp_2.load()
                        
                        nombre_2 = 1
                        for composante_2 in composantes_connexes_2 :
                                            
                            if len(composante_2) < 10 and len(composante_2) > 2 :
                                if nombre_2 == elt[1]: 
                                    groupe_2.append(composante_2)        
                                nombre_2 += 1
                        
                return groupe_1, groupe_2, comp_1, comp_2


'''comparaison des clusters obtenus avec les tailles d'extension 4 et 10 (petites composantes connexes de taille comprise entre 3 et 9) 
fait des tas de trucs compliques pour les comparer'''
def comparaison_cluster_identiques_plus_de_2(tab_taille_ext, depart, typ):                
    
    tab_idem = []
    nombre_comp_1 = [0] * 10
    nombre_comp_2 = [0] * 10
    tab_inclus_dans = []
    # tab_inclus_dans_2 = []
    
    with open(EXTENSION_PATH % tab_taille_ext[0] + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (typ, str(tab_taille_ext[0])), 'rb') as fichier_graphe_complet_1 :
        mon_depickler_complet_1 = pickle.Unpickler(fichier_graphe_complet_1)
        graphe_complet_1 = mon_depickler_complet_1.load()
        
        # recherche d'identiques ou d'inclus au sein d'une taille 
        i = 0.1
        while i <= 1.0 :

            with open(EXTENSION_PATH % tab_taille_ext[0] + "composantes_connexes_%s_%s_taille_%s_%s.pickle" % (depart, typ, str(tab_taille_ext[0]), str(round(i, 1))), 'rb') as fichier_comp_1:
                mon_depickler_comp_1 = pickle.Unpickler(fichier_comp_1)
                composantes_connexes_1 = mon_depickler_comp_1.load()

                nombre_1 = 1
                for composante_1 in composantes_connexes_1 :   
                    if len(composante_1) < 10 and len(composante_1) > 2 :
                        graphe_comp_1 = nx.Graph()
                        graphe_comp_1 = graphe_complet_1.subgraph(composante_1)
                        
                        trouve = False
                        for elt in tab_idem :
                            if (tab_taille_ext[0], round(i, 1), nombre_1) in elt :
                                trouve = True
                        if trouve == False :
                            tab_idem.append([(tab_taille_ext[0], round(i, 1), nombre_1)])
                        
                        trouve = False
                        for elt in tab_inclus_dans :
                            if (tab_taille_ext[0], round(i, 1), nombre_1) in elt :
                                trouve = True
                        if trouve == False :
                            tab_inclus_dans.append([(tab_taille_ext[0], round(i, 1), nombre_1)])   

                        j = 0.1
                        while j <= 1.0 :
                            with open(EXTENSION_PATH % tab_taille_ext[0] + "composantes_connexes_%s_%s_taille_%s_%s.pickle" % (depart, typ, str(tab_taille_ext[0]), str(round(j, 1))), 'rb') as fichier_comp_2:
                                mon_depickler_comp_2 = pickle.Unpickler(fichier_comp_2)
                                composantes_connexes_2 = mon_depickler_comp_2.load()

                                nombre_2 = 1
                                for composante_2 in composantes_connexes_2 :
                                        if i != j or composante_1 != composante_2 :
                                            if len(composante_2) < 10 and len(composante_2) > 2 :
                                                nombre_idem = 0
                                                for elt in tab_idem :
                                                    if (tab_taille_ext[0], round(i, 1), nombre_1) in elt :
                                                        graphe_comp_2 = nx.Graph()
                                                        graphe_comp_2 = graphe_complet_1.subgraph(composante_2)
                                                         
                                                        compteur = 0
                                                        for noeud_1, data_1 in graphe_comp_1.nodes(data=True) :
                                                            for noeud_2, data_2 in graphe_comp_2.nodes(data=True) :
                                                                if data_1["nom"] == data_2["nom"] : 
                                                                    compteur += 1 
                                                        if compteur == len(composante_1) and compteur == len(composante_2) :
                                                            # print("%f,%f"%(i,j))
                                                            # print("%d, %d\n"%(nombre_1, nombre_2))
                                                            # print(tab_idem[nombre_idem])
                                                            # print(nombre_2)
                                                            if (tab_taille_ext[0], round(j, 1), nombre_2) not in tab_idem[nombre_idem] :
                                                                tab_idem[nombre_idem].append((tab_taille_ext[0], round(j, 1), nombre_2))     
    #                                                         if compteur == len(composante_2) :
    #                                                             tab_inclus_dans.append(((round(i,1),nombre_1),(round(j,1), nombre_2)))
    #                                              
                                                    nombre_idem += 1
                                                
                                                nombre_inclus = 0    
                                                for elt in tab_inclus_dans :
                                                    if (tab_taille_ext[0], round(i, 1), nombre_1) in elt :
                                                        graphe_comp_2 = nx.Graph()
                                                        graphe_comp_2 = graphe_complet_1.subgraph(composante_2)
                                                         
                                                        compteur = 0
                                                        for noeud_1, data_1 in graphe_comp_1.nodes(data=True) :
                                                            for noeud_2, data_2 in graphe_comp_2.nodes(data=True) :
                                                                if data_1["nom"] == data_2["nom"] : 
                                                                    compteur += 1 
                                                        if compteur == len(composante_1) or compteur == len(composante_2) :
                                                            if (tab_taille_ext[0], round(j, 1), nombre_2) not in tab_inclus_dans[nombre_inclus] :
                                                                tab_inclus_dans[nombre_inclus].append((tab_taille_ext[0], round(j, 1), nombre_2))
                                                    nombre_inclus += 1             
    #                                                         if compteur == len(composante_2) :
    #                                                             tab_inclus_dans.append(((round(i,1),nombre_1),(round(j,1), nombre_2)))
    #                                              
                                                    nombre_idem += 1
                                                nombre_2 += 1
                                            
                                nombre_comp_2[round(j * 10) - 1] = nombre_2 - 1 
                            j = j + 0.1

                        nombre_1 += 1
            i += 0.1

        for k in range(1, len(tab_taille_ext)) :    
            if tab_taille_ext[k] != -1 :
                    with open(EXTENSION_PATH % tab_taille_ext[k] + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (typ, str(tab_taille_ext[k])), 'rb') as fichier_graphe_complet_2 :
                        mon_depickler_complet_2 = pickle.Unpickler(fichier_graphe_complet_2)
                        graphe_complet_2 = mon_depickler_complet_2.load()
                     
                        i = 0.1
                        while i <= 1.0 :
 
                                with open(EXTENSION_PATH % tab_taille_ext[0] + "composantes_connexes_%s_%s_taille_%s_%s.pickle" % (depart, typ, str(tab_taille_ext[0]), str(round(i, 1))), 'rb') as fichier_comp_1:
                                    mon_depickler_comp_1 = pickle.Unpickler(fichier_comp_1)
                                    composantes_connexes_1 = mon_depickler_comp_1.load()
                    
                                    nombre_1 = 1
                                    for composante_1 in composantes_connexes_1 :   
                                        if len(composante_1) < 10 and len(composante_1) > 2 :
                                            graphe_comp_1 = nx.Graph()
                                            graphe_comp_1 = graphe_complet_1.subgraph(composante_1)
    
                                            j = 0.1
                                            while j <= 1.0 :
                                                with open(EXTENSION_PATH % tab_taille_ext[k] + "composantes_connexes_%s_%s_taille_%s_%s.pickle" % (depart, typ, str(tab_taille_ext[k]), str(round(j, 1))), 'rb') as fichier_comp_2:
                                                    mon_depickler_comp_2 = pickle.Unpickler(fichier_comp_2)
                                                    composantes_connexes_2 = mon_depickler_comp_2.load()
                    
                                                    nombre_2 = 1
                                                    for composante_2 in composantes_connexes_2 :
                                                            if i != j or composante_1 != composante_2 :
                                                                if len(composante_2) < 10 and len(composante_2) > 2 :
                                                                    nombre_idem = 0
                                                                    for elt in tab_idem :
                                                                        if (tab_taille_ext[0], round(i, 1), nombre_1) in elt :
                                                                            graphe_comp_2 = nx.Graph()
                                                                            graphe_comp_2 = graphe_complet_2.subgraph(composante_2)
                                                                             
                                                                            compteur = 0
                                                                            for noeud_1, data_1 in graphe_comp_1.nodes(data=True) :
                                                                                for noeud_2, data_2 in graphe_comp_2.nodes(data=True) :
                                                                                    if data_1["nom"] == data_2["nom"] : 
                                                                                        compteur += 1 
                                                                            if compteur == len(composante_1) and compteur == len(composante_2) :
                                                                                # print("%f,%f"%(i,j))
                                                                                # print("%d, %d\n"%(nombre_1, nombre_2))
                                                                                # print(tab_idem[nombre_idem])
                                                                                # print(nombre_2)
                                                                                if (tab_taille_ext[k], round(j, 1), nombre_2) not in tab_idem[nombre_idem] :
                                                                                    tab_idem[nombre_idem].append((tab_taille_ext[k], round(j, 1), nombre_2))     
                        #                                                         if compteur == len(composante_2) :
                        #                                                             tab_inclus_dans.append(((round(i,1),nombre_1),(round(j,1), nombre_2)))
                        #                                              
                                                                        nombre_idem += 1
                                                                    
                                                                    nombre_inclus = 0    
                                                                    for elt in tab_inclus_dans :
                                                                        if (tab_taille_ext[0], round(i, 1), nombre_1) in elt :
                                                                            graphe_comp_2 = nx.Graph()
                                                                            graphe_comp_2 = graphe_complet_1.subgraph(composante_2)
                                                                             
                                                                            compteur = 0
                                                                            for noeud_1, data_1 in graphe_comp_1.nodes(data=True) :
                                                                                for noeud_2, data_2 in graphe_comp_2.nodes(data=True) :
                                                                                    if data_1["nom"] == data_2["nom"] : 
                                                                                        compteur += 1 
                                                                            if compteur == len(composante_1) or compteur == len(composante_2) :
                                                                                if (tab_taille_ext[k], round(j, 1), nombre_2) not in tab_inclus_dans[nombre_inclus] :
                                                                                    tab_inclus_dans[nombre_inclus].append((tab_taille_ext[k], round(j, 1), nombre_2))
                                                                        nombre_inclus += 1             
                        #                                                         if compteur == len(composante_2) :
                        #                                                             tab_inclus_dans.append(((round(i,1),nombre_1),(round(j,1), nombre_2)))
                        #                                              
                                                                        nombre_idem += 1
                                                                    nombre_2 += 1
                                                                
                                                    nombre_comp_2[round(j * 10) - 1] = nombre_2 - 1 
                                                j = j + 0.1
                    
                                            nombre_1 += 1
                                i += 0.1
    
    print(tab_idem) 
    print(len(tab_idem))
    print(tab_inclus_dans) 
    print(len(tab_inclus_dans))  
    
    a_enlever = []
    for groupe_idem in tab_idem :
        nombre_par_taille = [0] * len(tab_taille_ext)
        for elt in groupe_idem :
#             print(elt[0])
            nombre_par_taille[elt[0] - tab_taille_ext[0]] += 1
        # print(nombre_par_taille)
        compteur = 0
        for elt in nombre_par_taille :
            # print(tab_taille_ext[elt])
            if elt == 0 and tab_taille_ext[compteur] != -1:
                if groupe_idem not in a_enlever :
                    a_enlever.append(groupe_idem)
            compteur += 1
    for elt in a_enlever :
        tab_idem.remove(elt)
    
#     print(a_enlever)    
    
    a_enlever = []
    for groupe_inclus in tab_inclus_dans :
        nombre_par_taille = [0] * len(tab_taille_ext)
        for elt in groupe_inclus :
            nombre_par_taille[elt[0] - tab_taille_ext[0]] += 1
        
        compteur = 0
        for elt in nombre_par_taille :
            if elt == 0 and tab_taille_ext[compteur] != -1:
                if groupe_inclus not in a_enlever :
                    a_enlever.append(groupe_inclus)
            compteur += 1
    for elt in a_enlever :
        tab_inclus_dans.remove(elt)
    
#     print(a_enlever)
    
    for elt in tab_idem :
        print(elt)                                   
    print(len(tab_idem))
    for elt in tab_inclus_dans :
        print(elt)
    print(len(tab_inclus_dans)) 
    return tab_idem, tab_inclus_dans
                   

''' affiche des graphiques du nombre de composantes connexes en fonction de la taille de la composante pour chaque seuil
et les stocke dans des fichiers .png '''       
def stats_composantes_connexes(taille_ext, depart, typ, max_val):
    plt.figure(figsize=(20, 20))
    plt.suptitle("Distribution des composantes connexes")
    plt.gcf().subplots_adjust(left=0.25, bottom=0.05,
                       right=0.75, top=0.90, wspace=0, hspace=1.0)
    i = 0.1 * max_val
    while i <= 1.0 * max_val + 0.0000005 :
        with open(EXTENSION_PATH % taille_ext + "composantes_connexes_" + depart + "_" + typ + "_taille_" + str(taille_ext) + "_" + str(round(i, 2)) + ".pickle", 'rb') as fichier_comp:
            mon_depickler = pickle.Unpickler(fichier_comp)
            composantes_connexes = mon_depickler.load()
            print(fichier_comp)
            
            with open(EXTENSION_PATH % taille_ext + "graphe_complet_pondere_sim_%s_taille_%s.pickle" % (typ, taille_ext), 'rb') as fichier_graphe_complet :
                mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
                graphe_complet = mon_depickler_complet.load()
                
                tab_densite = []
                tab_taille = []
                print(composantes_connexes)
                for composante in composantes_connexes :
                    if len(composante) < 2 :
                        densite_comp = None
                    else :
                        nombre_aretes = 0
                        for u, v, data in graphe_complet.subgraph(composante).edges(data=True) :
                            if data["poids"] >= i :
                                nombre_aretes += 1
                        densite_comp = nombre_aretes / ((len(composante) * (len(composante) - 1)) / 2)
                    tab_densite.append(densite_comp)
                    tab_taille.append(len(composante))
                
                # print(tab_taille)
                # print(tab_densite)
                nb_par_taille = [0] * 90
                densite_par_taille = [(0, 0)] * 90
                # print(densite_par_taille)
                
                for j in range(len(tab_taille)) :
                    nb_par_taille[tab_taille[j] - 1] += 1
                    if tab_densite[j] != None :
                        # print(densite_par_taille[tab_taille[j]-1])
                        new_val_densite = list(densite_par_taille[tab_taille[j] - 1])[0] + tab_densite[j]
                        new_val_nombre = list(densite_par_taille[tab_taille[j] - 1])[1] + 1
                        
                        densite_par_taille[tab_taille[j] - 1] = (new_val_densite, new_val_nombre)
                    else :
                        densite_par_taille[tab_taille[j] - 1] = (None, None)
                print(nb_par_taille)
                print(densite_par_taille)   
                
                densite_moy_par_taille = [0] * 90
                for j in range(len(densite_par_taille)) :
                    if list(densite_par_taille[j])[1] != 0 :
                        if list(densite_par_taille[j])[1] != None :
                            densite_moy_par_taille[j] = list(densite_par_taille[j])[0] / list(densite_par_taille[j])[1]
                        else :
                            densite_moy_par_taille[j] = None
                    
                print(densite_moy_par_taille)
                
                ax1 = plt.subplot(10, 1, round(i / max_val * 10))
                plt.title("pour les arêtes de poids >= %1.1f (%d %%)" % (round(i, 2), i / max_val * 100))
                print(i / max_val * 10)
                ax1.xaxis.set_ticks(np.arange(0, 90, 5))
                ax1.set_xlabel("Taille")
                ax1.set_ylabel("Nombre")
                ax1.xaxis.set_label_coords(1.0, -0.3)
                for j in range(len(densite_moy_par_taille)) :
                    if densite_moy_par_taille[j] != 0 and densite_moy_par_taille[j] != None :
                        if j < 10 :
                            dec = 1
                        else :
                            dec = 0
                        if i / max_val < 0.6 :
                            plt.annotate(str(round(densite_moy_par_taille[j], 2)), xy=(j, 0), xytext=(dec + j, 1))
                        elif i / max_val == 0.6 :
                            if j < 3 : 
                                plt.annotate(str(round(densite_moy_par_taille[j], 2)), xy=(j, 0), xytext=(dec + j, 2))
                            else :
                                plt.annotate(str(round(densite_moy_par_taille[j], 2)), xy=(j, 0), xytext=(dec + j, 1))
                        elif i / max_val > 0.6 :
                            plt.annotate(str(round(densite_moy_par_taille[j], 2)), xy=(j, 0), xytext=(dec + j, (j % 3) * 10 + 5))       
                plt.hist(tab_taille, range=(1, 90), bins=91, edgecolor='black')
        i = i + (0.1 * max_val)
    # plt.show()
    plt.savefig(EXTENSION_PATH % taille_ext + "distribution_comp_connexes_%s_%s.png" % (depart, typ))

''' affiche dans la console les composantes qui ne contiennent pas que des homologues '''
def non_homologues_dans_cluster(taille_ext, depart, typ):
    homologues = [['4V9F_0_62_12', '5J7L_DA_62_5', '5DM6_X_328_2', '5FDU_1A_62_14'],
                  ['4V9F_0_30_4', '5J7L_DA_30_15', '5FDU_1A_30_17'],
                  ['4V9F_0_48_13', '5J7L_DA_48_20', '5DM6_X_48_28', '5FDU_1A_48_25'],
                  ['4V9F_0_207_3', '5J7L_DA_272_2', '5FDU_1A_272_1'],
                  ['4V9F_0_25_56', '5J7L_DA_25_12', '5FDU_1A_25_68', '4V88_A5_25_30'],
                  ['4V9F_0_137_5', '5J7L_DA_48_1', '5FDU_1A_137_6', '5DM6_X_48_10', '4V88_A5_48_3'],
                  ['4V9F_0_134_5', '5FDU_1A_74_7', '5DM6_X_127_7'],
                  ['4V9F_0_48_21', '5J7L_DA_197_4', '5FDU_1A_197_3', '5DM6_X_48_9'],
                  ['4V9F_0_127_6', '5J7L_DA_134_1', '5FDU_1A_134_3', '5DM6_X_134_2'],
                  ['4V9F_0_287_2', '5DM6_X_25_15'],
                  ['5J7L_DA_25_10', '5DM6_X_25_34', '4V88_A5_25_47', '5FDU_1A_25_78'],
                  ['1FJG_A_48_17', '5J5B_BA_48_14'],
                  ['1FJG_A_48_8', '5J5B_BA_48_23'],
                  ['1FJG_A_294_1', '5J5B_BA_294_2'],
                  ['1FJG_A_58_23', '5J5B_BA_58_3'],
                  ['1FJG_A_138_3', '5J5B_BA_138_2'],
                  ['5J5B_BA_48_7', '4V88_A6_48_12'],
                  ['4YAZ_R_36_25', '3UCZ_R_62_15']]
    
    i = 0.1
    while i <= 1.0 :
        with open(EXTENSION_PATH % taille_ext + "composantes_connexes_" + depart + "_" + typ + "_taille_" + str(taille_ext) + "_" + str(round(i, 1)) + ".pickle", 'rb') as fichier_comp:
            mon_depickler = pickle.Unpickler(fichier_comp)
            composantes_connexes = mon_depickler.load()
            
            with open(EXTENSION_PATH % taille_ext + "graphe_complet_pondere_sim_coeff_all1_taille_%s.pickle" % taille_ext, 'rb') as fichier_graphe_complet :
                mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
                graphe_complet = mon_depickler_complet.load()
                
                compteur = 1
                for composante in composantes_connexes :
                    tous_homologues = False
                    if len(composante) > 2 and len(composante) < 10 :
                        for groupe in homologues :
                            nombre_dans_groupe = 0
                            for elt in composante :
                                if graphe_complet.nodes[elt]["nom"] in groupe :
                                    nombre_dans_groupe += 1
                            if nombre_dans_groupe == len(composante) :
                                tous_homologues = True
                        if tous_homologues == False :
                            print(i)
                            print(compteur)
                        compteur += 1
            
        i = i + 0.1
        
                    
liste = ['5J7L_DA_191_3', '5J7L_DA_191_4', '5FDU_1A_301_1', '5J7L_DA_301_2', '5DM6_X_334_1', '5FDU_1A_334_2', '4V9F_0_335_1', '5J7L_DA_335_2', '3JCS_1_137_4', '4V88_A5_290_1', '4V88_A6_314_2', '5J7L_DA_218_3', '4V9F_0_251_2', '1FJG_A_62_8', '5J7L_DA_137_1', '4V9F_0_118_1', '4V9F_0_62_2', '5J7L_DA_271_2', '4V9F_0_224_1', '5DM6_X_197_1', '3GX5_A_138_6', '1FJG_A_317_2', '5J5B_BA_317_1', '1FJG_A_326_1', '5DM6_X_137_3', '5J5B_BA_314_1', '4V9F_0_134_6', '4V9F_0_328_1', '4V9F_0_197_2', '4V9F_0_62_16', '5J7L_DA_282_2', '4V88_A5_137_2', '5FDU_1A_224_3', '5J7L_DA_326_2']

if __name__ == '__main__':
    comparaison_cluster()
    # compter_composantes_differentes()
    
    # draw_composantes("toutes_aretes", "extensions", type_comp="extensions_toutes_aretes")
    
    # calcul_sim_morceau("toutes_aretes", "extensions", "Extensions/Metrique_toutes_aretes/dico_comp_complet_metrique_toutes_aretes_coeffn1_a1_c1.pickle")
    
    # draw_groupes_par_type_motif()
#     with open("graphe_complet_pondere_sim.pickle", 'rb') as fichier_graphe_complet :
#                 mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
#                 graphe_complet = mon_depickler_complet.load()
#                  
#                 i = 0.1
#                 while i <= 1.0 :
#                     with open("composantes_connexes_extensions_non_cov_nouvelle_metrique_"+str(i)+".pickle", 'rb') as fichier_comp:
#                         mon_depickler_comp = pickle.Unpickler(fichier_comp)
#                         composantes_connexes = mon_depickler_comp.load()
#                         compteur = 1
#                         for composante in composantes_connexes :
#                              
#                             if len(composante) < 10 and len(composante) > 2 :
#                                 print(i)
#                                 est_la = False
#                                 for elt in composante :
#                                     if graphe_complet.nodes[elt]["nom"] in ['5FDU_1A_272_1', '5J7L_DA_272_2', '4V9F_0_207_3', '4FAU_A_207_4'] :
#                                         print(graphe_complet.nodes[elt]["nom"] )
#                                         est_la = True
#                                 if est_la == True :
#                                     print(i)
#                                     print(compteur)
#                                     print(composante)
#                                 compteur += 1
#                     i = i+0.1
    
#             print(graphe_complet_grands_graphes.nodes.data())
#             print(graphe_complet_grands_graphes.edges.data())

#     with open("sim_extensions_non_covpar_chaine.pickle", 'rb') as fichier_graphe_complet_plus :
#         mon_pickler_complet = pickle.Unpickler(fichier_graphe_complet_plus)
#         dico_sim_par_chaine = mon_pickler_complet.load()
#         print(len(dico_sim_par_chaine.keys()))
# #         print(dico_sim_par_chaine.keys())
#         print(dico_sim_par_chaine[('fichier_1FJG_A_109_6', 'fichier_4V9F_0_30_23')])
#         for cle in dico_sim_par_chaine.keys() :
#             for i in range(4) :
#                 if dico_sim_par_chaine[cle][i][0] < 0.0 or dico_sim_par_chaine[cle][i][0] > 1.0 :
#                         print(cle)
# #                         print(i)
#                         print(dico_sim_par_chaine[cle])
#          
#         with open("grands_graphes.pickle", 'rb') as fichier_graphes :
#             mon_depickler_graphes = pickle.Unpickler(fichier_graphes)
#             graphes = mon_depickler_graphes.load()       
#             print(graphes[('4V9F','0',48,26)].nodes.data())         
#                         
#     with open("graphe_complet_pondere_sim.pickle", 'rb') as fichier_graphe_complet :
#         mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
#         graphe_complet = mon_depickler_complet.load()
#         for x,data in graphe_complet.nodes(data=True) :
#             if data["nom"] == '4V9F_0_48_26' :
#                 u = x
#             if data["nom"] == '5J7L_DA_48_30' :
#                 v = x
#          
#         chaines_1, chaines_2, chaines_commun = recherche_chaines_par_morceau(u,v, graphe_complet, "graphe_global")
#         print(u)
#         print(v)
#         print(chaines_commun) 
#         with open("grands_graphes.pickle", 'rb') as fichier_graphes :
#             mon_depickler_graphes = pickle.Unpickler(fichier_graphes)
#             graphes = mon_depickler_graphes.load()
#             with open("fichier_comp_grands_graphes_V2.pickle", 'rb') as fichier_commun :
#                 mon_depickler_commun = pickle.Unpickler(fichier_commun)
#                 dico_graphe = mon_depickler_commun.load()
#                  
#                 elt1 = ('4V9F', '0', 48, 26)
#                 elt2 = ('5J7L', 'DA', 48, 30)
#                  
#                 if (elt1, elt2) in dico_graphe.keys() :
#                     cle = (elt1, elt2)
#                 else :
#                     cle = (elt2, elt1)
#                 graphe_commun = dico_graphe[cle]
#                 tab_sim = []
#                 for j in range(4) :
#                     print(j)
#                     sim = calcul_sim_non_cov_sans_motif_par_chaine(graphes[cle[0]], graphes[cle[1]], graphe_commun, chaines_1, chaines_2, chaines_commun, j)
#                     print(sim)
          
#         with open("graphe_complet_pondere_sim.pickle", 'rb') as fichier_graphe_complet :
#             mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
#             graphe_complet = mon_depickler_complet.load()
#             print(graphe_complet.nodes.data())
#             for u,v in graphe_complet.edges() :
#                 if (graphe_complet.nodes[u]["nom"] == '5FJC_A_138_1' and graphe_complet.nodes[v]["nom"] == '5DM6_X_227_2') or (graphe_complet.nodes[v]["nom"] == '5FJC_A_138_1' and graphe_complet.nodes[u]["nom"] == '5DM6_X_227_2'):
#                     chaines_1, chaines_2, chaines_commun = recherche_chaines_par_morceau(u,v, graphe_complet, "graphe_global")
#                     print(chaines_1)
#                     print(chaines_2)
#                     print(chaines_commun)

#     with open("grands_graphes.pickle", 'rb') as fichier :
#         mon_depickler = pickle.Unpickler(fichier)
#         dico_graphes = mon_depickler.load()
#         print(dico_graphes[('5FJC', 'A',138,1)].edges.data())
        
#           
#         with open("fichiers_pickle/a-minor_test2.pickle", 'rb') as fichier_pickle :
#             mon_depickler_aminor = pickle.Unpickler(fichier_pickle)
#             tab_aminor = mon_depickler_aminor.load()
#             
#             for occ in tab_aminor :
#                 est_la = False
#                 for elt in liste :
#                     if occ["num_PDB"]+"_"+ occ["num_ch"]+"_"+ str(occ["num_motif"])+"_"+ str(occ["num_occ"]) == elt :
#                         est_la = True
#                         
#                 if est_la == False :
#                 
#                     cle = (occ["num_PDB"], occ["num_ch"], occ["num_motif"], occ['num_occ'])
#                     print(cle)
#                     ajout_chaines_grands_graphes(dico_graphes[cle], occ["a_minor"])
#                     print(dico_graphes[cle].nodes.data())
#                     
#     with open("grands_graphes.pickle", 'wb') as fichier_ec :
#         mon_pickler = pickle.Pickler(fichier_ec)
#         mon_pickler.dump(dico_graphes)
            
    # choix_sim_stockee("par_poids_sans_motif")
    # calcul_sim_morceau("non_cov", "extensions")
    
    # draw_composantes("non_cov", "graphe_global", type_comp="graphe_global_non_cov")
#     
#     
#     with open("graphe_complet_pondere_sim.pickle", 'rb') as fichier_graphe_complet :
#         mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
#         graphe_complet = mon_depickler_complet.load()
#         
#         tab_proportion = [0,0,0,0]
#         nb_aretes = 0
#         
#         i = 0.1
#         #i = 0.4
#         graphes_deja_faits = []
#         while i <= 1.0 :
#             with open("composantes_connexes_"+str(i)+".pickle", 'rb') as fichier_comp:
#                 mon_depickler_comp = pickle.Unpickler(fichier_comp)
#                 composantes_connexes = mon_depickler_comp.load()
#                 
#                 for composante in composantes_connexes :
#                     if len(composante) < 10 and len(composante) > 2 :
#                         graphe_comp = graphe_complet.subgraph(composante)
#                         for u,v,t in graphe_comp.edges(data="poids") :
#                             if (graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]) not in graphes_deja_faits and (graphe_comp.nodes[v]["nom"], graphe_comp.nodes[u]["nom"]) not in graphes_deja_faits :
#                                 if t < i :
#                                     for j in range(4) :
#                                         if graphe_complet.edges[u,v]["sim_par_chaine"][j] > i :
#                                             tab_proportion[j] += 1
#                                     nb_aretes += 1
#                                 graphes_deja_faits.append((graphe_comp.nodes[u]["nom"], graphe_comp.nodes[v]["nom"]))
#         
#             i = i + 0.1
#         for j in range(4) :
#             tab_proportion[j] = tab_proportion[j] / nb_aretes   
#         print(tab_proportion)
#         for i in range(4) :
#             tab_proportion[i] = tab_proportion[i] / nb_aretes
                       
