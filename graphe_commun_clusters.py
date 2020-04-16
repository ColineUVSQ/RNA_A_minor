'''
Created on 28 mars 2019

@author: coline

Fonctions de calculs du sous-graphe commun moyen a un cluster
(version CaRNAval et new data)
'''

import networkx as nx
import pickle
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Circle
from recup_data.constantes import GROUPES_TOUTES_ARETES_MAX_4_10_07,\
    EXTENSION_PATH, GROUPE_DENSE_DE_07, GROUPE_ARICH_DE_07,\
    EXTENSION_PATH_TAILLE, GROUPE_ETOILE_GNRA,\
    CLUSTERING_PEREZ_VERSION_NON_CAN_2, NEW_EXTENSION_PATH_TAILLE
from recup_data.draw_like_carnaval import motif_commun_to_2D

# def noeud_deja_dans_graphe(graphe, noeud, num_elt):
#     for n, data in graphe.nodes(data=True) :
#         if data["num_elt"] == num_elt and data["num_noeud"] == noeud :
#             return True
#     return False

''' renvoie le numero du sommet suivant au noeud passe en parametre si celui-ci 
est aussi dans l'ensemble des cliques correspondant au sous-graphe commun 
noeud : le noeud courant
graphe : le graphe courant
liste_sommets : l'ensemble des cliques correspondant au sous-graphe commun
num : dans les cliques, le numero correspondant au graphe courant '''
def suivant_existe(noeud, graphe, liste_sommets, num):
#     print(liste_sommets)
#     print(num)
    for clique in liste_sommets : 
        for elt in clique :  
            #print("gros gros rat")
            #print(elt)
            if elt[1] == num :
                #print("gros gros rat")
                #print(graphe.nodes[elt[0]]["position"][0])
                #print(graphe.nodes[noeud]["position"][1]+1)
#                 print(elt[0])
#                 print(noeud)
#                 print(num)
#                 print(clique)
                if graphe.nodes[elt[0]]["position"][0] == graphe.nodes[noeud]["position"][1]+1 :
                    return elt[0]
    return -1

''' recherche le plus grand sous-graphe commun dans un cluster a l'aide de cliques
(version CaRNAval)
rep_comparaison : chemin du fichier ou se trouve le dictionnaire des graphes communs
rep_extension : chemin du repertoire ou se trouvent les fichiers de graphes d'extension '''
def commun_cluster_clique(cluster, rep_comparaison, rep_extension):

    graphe_commun = nx.MultiDiGraph()
    print(len(cluster))

    if len(cluster) > 2 :
        
        ## Construction du graphe dans lequel on recherchera des cliques ##
        with open(rep_comparaison, 'rb') as fichier_comp :
            mon_depickler = pickle.Unpickler(fichier_comp)
            dico_comp = mon_depickler.load()
            graphe_cluster_clique = nx.Graph()
            
            for i in range(len(cluster)) :
                for j in range(i+1, len(cluster)) :
                    if ("fichier_"+cluster[i], "fichier_"+cluster[j]) in dico_comp.keys() :
                        cle_0 = i
                        cle_1 = j 
                    else :
                        cle_0 = j 
                        cle_1 = i
                        
                    for noeud in dico_comp[("fichier_"+cluster[cle_0], "fichier_"+cluster[cle_1])].nodes() :
                        #print(noeud)
                        if noeud not in [(1,1),(2,2),(3,3),(4,4),(5,5)] :
                            if (noeud[0],cle_0) not in graphe_cluster_clique.nodes() :
                                graphe_cluster_clique.add_node((noeud[0],cle_0))
                            if (noeud[1],cle_1) not in graphe_cluster_clique.nodes() :
                                graphe_cluster_clique.add_node((noeud[1],cle_1))
                            
                            graphe_cluster_clique.add_edge((noeud[0],cle_0), (noeud[1],cle_1))
        #print(graphe_cluster_clique.nodes.data())
        #print(graphe_cluster_clique.edges.data())
        
        ## Recherche des cliques de toutes tailles ##
        cliques = list(nx.enumerate_all_cliques(graphe_cluster_clique))
        #print(cliques)
        
        ## Recuperation des cliques de taille len(cluster) s'il y en a (ce sont les cliques qui nous interessent) ##
        liste_sommets_motifs = []
        liste_cliques = []
        for elt in cliques :
            if len(elt) == len(cluster) :
                #print(elt)
                if elt[0][1] != 0 :
                    print("probleme")
                for indice in elt :
                    if indice[1] == 0 :
                        liste_sommets_motifs.append(indice[0])
                liste_cliques.append(elt)
        for i in range(1,6) :
            elt = []
            for j in range(len(cluster)) :
                elt.append((i,j))
            liste_cliques.append(elt)
        liste_sommets_motifs.extend([1,2,3,4,5])
        #print(liste_sommets_motifs)
        
        ### Construction du graphe commun moyen associe (avec les attributs necessaires : espacement_motif, poids, chaine, position (pas les vraies positions, 
        ### juste pour l'affichage)  ###
        with open(rep_extension+"/"+"fichier_"+cluster[0]+".pickle", 'rb') as fichier_ext :
            mon_depickler_ext = pickle.Unpickler(fichier_ext)
            graphe = mon_depickler_ext.load()
            
            for noeud in liste_sommets_motifs :
                graphe.nodes[noeud].update({'num_seq' : [(cluster[0], graphe.nodes[noeud]["position"])]})
                graphe_commun.add_node(noeud, **graphe.nodes[noeud])
            
            for u,v,data in graphe.edges(data=True) :
                if u in graphe_commun.nodes() and v in graphe_commun.nodes() and data["label"] != 'B53' :
                    graphe_commun.add_edge(u,v,**data)
    #     
            nx.set_node_attributes(graphe_commun, [], "espacement_motif")
    #         noeuds_isoles_a_enlever = []
    
    
            chaines = [[1]]
            for i in range(1,5) :
                    compteur = i
                    if i != 1 : chaines.append([i])
                    liaison_B53 = True
                    while liaison_B53 :
                        liaison_B53 = False
                        temp = compteur
                        if i == 1 or i == 4 :
                            for voisin in graphe.successors(compteur) :
                                for arc in graphe[compteur][voisin] :
                                    if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and graphe[compteur][voisin][arc]["label"] == 'B53' :
                                        liaison_B53 = True
                                        temp = voisin
                                        chaines[len(chaines)-1].append(voisin)
                        else :            
                            for voisin in graphe.predecessors(compteur) :
                                for arc in graphe[voisin][compteur] :
                                    if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and graphe[voisin][compteur][arc]["label"] == 'B53' :
                                        liaison_B53 = True
                                        temp = voisin
                                        chaines[len(chaines)-1].append(voisin)
                        compteur = temp
            noeuds_isoles_a_enlever = []            
            for noeud in graphe_commun.nodes() :
                    
                    liaison_b53_a_mettre = []
                    
                    liaison_b53_a_mettre.append(suivant_existe(noeud, graphe, liste_cliques, 0)) 
                    #print("gros rat")
                    #print(liaison_b53_a_mettre)
                    
                    a_enlever = []
                    espacement_motif = []
                     
                    mini_poids = graphe.nodes[noeud]["poids"]
                    chaine = graphe.nodes[noeud]["chaine"]
                    liste_poids = []
                     
                    voisins = graphe_commun[noeud]
            #                 print(voisins)
    #                 liaison_autre_que_b53 = False
    #                 if len(voisins) == 1 :
    #                     for edge in graphe_commun[noeud][next(iter(voisins))] :
    #                         if graphe_commun[noeud][next(iter(voisins))][edge]["label"] != 'B53' :
    #                             liaison_autre_que_b53 = True
    #                 if noeud == (5,5) :
    #                     print("petit rat")
    #                     print(len(voisins))
    #                     print(liaison_autre_que_b53)           
                    if len(voisins) >= 1  : ## si pas de voisin, on l enleve
                        compteur = 1
                        for ch in chaines :
                            if noeud in ch :
                                espacement_motif.append(abs(graphe.nodes[noeud]["position"][0]-graphe.nodes[compteur]["position"][0]))
                            compteur += 1
             
                        for i in range(1,len(cluster)) :
                            with open(rep_extension+"/"+"fichier_"+cluster[i]+".pickle", 'rb') as fichier_ext_2 :
                                mon_depickler_ext_2 = pickle.Unpickler(fichier_ext_2)
                                graphe2 = mon_depickler_ext_2.load()
                                if ("fichier_"+cluster[0], "fichier_"+cluster[i]) in dico_comp.keys() :
                                    place_i = 1
                                else :
                                    place_i = 0
                                   
                                for noeud2 in dico_comp[("fichier_"+cluster[abs(place_i-1)*i], "fichier_"+cluster[place_i*i])].nodes() :
                                    if noeud == noeud2[abs(place_i-1)] :
                                        graphe_commun.nodes[noeud]["num_seq"].append((cluster[i], graphe2.nodes[noeud2[place_i]]["position"]))
                                        
                                        for elt in chaine :
                                            if elt not in graphe2.nodes[noeud2[place_i]]["chaine"] :
                                                a_enlever.append(elt)
                                        if mini_poids > graphe2.nodes[noeud2[place_i]]["poids"] :
                                            mini_poids = graphe2.nodes[noeud2[place_i]]["poids"]
                                            #print(mini_poids)
                                            #print(cluster[0])
                                            #print(cluster[i])
                                        if graphe2.nodes[noeud2[place_i]]["poids"] not in liste_poids :
                                            liste_poids.append(graphe2.nodes[noeud2[place_i]]["poids"])
                                        #print("tout petit rat")  
                                        compteur = 1
            #                                     print("rat")
            #                                     print(noeud[place_0])
            #                                     print(chaines)
                                        for ch in chaines :
                                            if noeud in ch :
                                                espacement_motif.append(abs(graphe2.nodes[noeud2[place_i]]["position"][0]-graphe2.nodes[compteur]["position"][0]))
                                                 
                                            compteur += 1
                                        #print(noeud)
                                        liaison_b53_a_mettre.append(suivant_existe(noeud2[place_i], graphe2, liste_cliques, i)) 
            #                     print("ramous")
            #                     print(espacement_motif)
                        graphe_commun.nodes[noeud]["espacement_motif"] = list(espacement_motif)
            #                     print(a_enlever)
                        for elt in a_enlever :
                            if elt in chaine :
                                chaine.remove(elt)
                        #print("gros rat")
                        #print(liste_poids)
            #                 a_enlever = []
            #                 if len(liste_poids) > 1 :
            #                     for u,v,key,data in digraphe_commun.edges(data=True, keys=True) :
            #                         if u == noeud :
            #                             if data["label"] == 'B53' :
            #                                 a_enlever.append((u,v,key))
            #                 print("gros rat")
            #                 print(a_enlever)
            #                 for elt in a_enlever :
            #                     digraphe_commun.remove_edge(elt[0], elt[1], key=elt[2])
                         
             
                        graphe_commun.nodes[noeud]["poids"] = mini_poids          
                        graphe_commun.nodes[noeud]["chaine"] = chaine
                        
                        liaison_b53_a_ajouter = True
                        for elt in liaison_b53_a_mettre :
                            if elt == -1  :
                                liaison_b53_a_ajouter = False
                        #print("ramou")
                        #print(liaison_b53_a_mettre)
                        #print(liaison_b53_a_ajouter)
                        if liaison_b53_a_ajouter : 
                            graphe_commun.add_edge(noeud, liaison_b53_a_mettre[0], label='B53', long_range=False)
                        
                        chaine_position = []
                        for i in range(4) :
                            for j in range(len(chaines[i])) :
                                if noeud == chaines[i][j] :
            #                             print(num_noeud)
            #                             print(chaines[i][j])
                                    if i == 0 or i == 3 :
                                        chaine_position.append(j)
                                    else :
                                        chaine_position.append(10-j)
                        graphe_commun.nodes[noeud]["position"] = chaine_position
                            
                    else : 
                        noeuds_isoles_a_enlever.append(noeud)
            for elt in noeuds_isoles_a_enlever :
                graphe_commun.remove_node(elt)
            for i in range(1,6) :
                graphe_commun.nodes[noeud]["positon"] = [i]
    #         
            #print(graphe_commun.nodes.data())
    #         
            return graphe_commun, liste_cliques  
        
''' recherche le plus grand sous-graphe commun dans un cluster a l'aide de cliques
(version toutes donnees PDB)
rep_comparaison : chemin du fichier ou se trouve le dictionnaire des graphes communs
rep_extension : chemin du repertoire ou se trouvent les fichiers de graphes d'extension '''
def commun_cluster_clique_new_data(cluster, rep_comparaison, rep_extension):
    graphe_commun = nx.MultiDiGraph()
    print(len(cluster))
    with open(rep_comparaison, 'rb') as fichier_comp :
        mon_depickler = pickle.Unpickler(fichier_comp)
        dico_comp = mon_depickler.load()
    
        if len(cluster) > 1 :
        ## construction graphe des cliques
            graphe_cluster_clique = nx.Graph()
            
            for i in range(len(cluster)) :
                for j in range(i+1, len(cluster)) :
                    if (cluster[i], cluster[j]) in dico_comp.keys() :
                        cle_0 = i
                        cle_1 = j 
                    else :
                        cle_0 = j 
                        cle_1 = i
                        
                    for noeud in dico_comp[(cluster[cle_0], cluster[cle_1])]["graphe"].nodes() :
                        #print(noeud)
                        #print(cle_0)
                        if cluster[cle_0] == ('1mms', 1) or cluster[cle_1] == ('1mms', 1) :
                            print("rap")
                            print(cluster[cle_0])
                            print(cluster[cle_1])
                            print(noeud)
                        if noeud not in [(1,1),(2,2),(3,3),(4,4),(5,5)] :
                            if (noeud[0],cle_0) not in graphe_cluster_clique.nodes() :
                                graphe_cluster_clique.add_node((noeud[0],cle_0))
                            if (noeud[1],cle_1) not in graphe_cluster_clique.nodes() :
                                graphe_cluster_clique.add_node((noeud[1],cle_1))
                            
                            graphe_cluster_clique.add_edge((noeud[0],cle_0), (noeud[1],cle_1))
        #print(graphe_cluster_clique.nodes.data())
        #print(graphe_cluster_clique.edges.data())
        
            ## recherche des cliques de taille len(cluster)
            cliques = list(nx.enumerate_all_cliques(graphe_cluster_clique))
            #print(cliques)
            print(len(cliques))
            liste_sommets_motifs = []
            liste_cliques = []
            for elt in cliques :
                if len(elt) == len(cluster) :
                    print(elt)
                    if elt[0][1] != 0 :
                        print("probleme")
                    for indice in elt :
                        if indice[1] == 0 :
                            liste_sommets_motifs.append(indice[0])
                    liste_cliques.append(elt)
            for i in range(1,6) :
                elt = []
                for j in range(len(cluster)) :
                    elt.append((i,j))
                liste_cliques.append(elt)
            liste_sommets_motifs.extend([1,2,3,4,5])
            print(liste_sommets_motifs)
            print(len(liste_cliques))
            #return liste_cliques
            #exit(0)
            
#             for i in range(1,len(cluster)) :
#                 with open(rep_extension+"/"+"fichier_"+cluster[i][0]+"_"+str(cluster[i][1])+"_5.pickle", 'rb') as fichier_ext_2 :
#                     mon_depickler_ext_2 = pickle.Unpickler(fichier_ext_2)
#                     graphe2 = mon_depickler_ext_2.load()
                    

            
            ## remplissage du graphe commun
            with open(rep_extension+"/"+"fichier_"+cluster[0][0]+"_"+str(cluster[0][1])+"_2.pickle", 'rb') as fichier_ext :
                mon_depickler_ext = pickle.Unpickler(fichier_ext)
                graphe = mon_depickler_ext.load()
    
                for noeud in liste_sommets_motifs :
                    graphe.nodes[noeud].update({'num_seq' : [(cluster[0], graphe.nodes[noeud]["position"], graphe.nodes[noeud]["num_ch"])]})
                    graphe_commun.add_node(noeud, **graphe.nodes[noeud])
                
                print(graphe_commun.nodes.data())
               
                ## ajout des aretes apparaissant dans tous les graphes
                liste_aretes = []
                for u,v,data in graphe.edges(data=True) :
                    if u in graphe_commun.nodes() and v in graphe_commun.nodes() and data["label"] != 'B53' :
                        if (u,v, data) not in liste_aretes : 
                            liste_aretes.append((u,v,data))
                 
                for i in range(1,len(cluster)) :
                    with open(rep_extension+"/"+"fichier_"+cluster[i][0]+"_"+str(cluster[i][1])+"_2.pickle", 'rb') as fichier_ext_2 :
                        mon_depickler_ext_2 = pickle.Unpickler(fichier_ext_2)
                        graphe2 = mon_depickler_ext_2.load()
                        print(cluster[i][0])
                        if (cluster[0], cluster[i]) in dico_comp.keys() :
                            place_i = 1
                        else :
                            place_i = 0
                        
                        graphe_comp = dico_comp[(cluster[abs(place_i-1)*i], cluster[place_i*i])]["graphe"]
                        compteur = 0
                        a_enlever = []
                        for u,v,data in liste_aretes :
                            ok = False
                            u2 = -1
                            v2 = -1
                            for noeud in graphe_comp.nodes() :
                                if noeud[abs(place_i-1)] == u :
                                    u2 = noeud[place_i]
                                if noeud[abs(place_i-1)] == v :
                                    v2 = noeud[place_i]
                            if u2 != -1 and v2 != -1 :
                                if (u2,v2) in graphe2.edges() :
                                    
                                    for edge in graphe2[u2][v2] :
                                        if graphe2[u2][v2][edge]["label"] == data["label"] :
                                            ok = True
                            else :
                                print("probleme !!!")
                            if not ok :
                                a_enlever.append((u,v,data))
                        for elt in a_enlever :
                            del(liste_aretes[liste_aretes.index(elt)])
                
                for u,v,data in liste_aretes :
                    graphe_commun.add_edge(u,v,**data)
                
                print(graphe_commun.nodes.data())
                print(graphe_commun.number_of_nodes())
                
        #     
                nx.set_node_attributes(graphe_commun, [], "espacement_motif")
        #         noeuds_isoles_a_enlever = []
                print("gros ramou")
                print(graphe_commun.edges.data())
                #exit(0)
                
                ## gestion des attributs des noeuds
    
                chaines = [[1]]
                for i in range(1,5) :
                        compteur = i
                        if i != 1 : chaines.append([i])
                        liaison_B53 = True
                        while liaison_B53 :
                            liaison_B53 = False
                            temp = compteur
                            if i == 1 or i == 4 :
                                for voisin in graphe.successors(compteur) :
                                    for arc in graphe[compteur][voisin] :
                                        if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and graphe[compteur][voisin][arc]["label"] == 'B53' :
                                            liaison_B53 = True
                                            temp = voisin
                                            chaines[len(chaines)-1].append(voisin)
                            else :            
                                for voisin in graphe.predecessors(compteur) :
                                    for arc in graphe[voisin][compteur] :
                                        if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and graphe[voisin][compteur][arc]["label"] == 'B53' :
                                            liaison_B53 = True
                                            temp = voisin
                                            chaines[len(chaines)-1].append(voisin)
                            compteur = temp
                noeuds_isoles_a_enlever = []            
                for noeud in graphe_commun.nodes() :
                        print("noeuds")
                        print(noeud)
                        liaison_b53_a_mettre = []
                        
                        liaison_b53_a_mettre.append(suivant_existe(noeud, graphe, liste_cliques, 0)) 
                        #print("gros rat")
                        #print(liaison_b53_a_mettre)
                        
                        a_enlever = []
                        espacement_motif = []
                         
                        mini_poids = graphe.nodes[noeud]["poids"]
                        chaine = list(graphe.nodes[noeud]["chaine"])
                        liste_poids = []
#                         print("chaine")
#                         print(chaine)
                        
                         
                        voisins = graphe_commun[noeud]
                #                 print(voisins)
        #                 liaison_autre_que_b53 = False
        #                 if len(voisins) == 1 :
        #                     for edge in graphe_commun[noeud][next(iter(voisins))] :
        #                         if graphe_commun[noeud][next(iter(voisins))][edge]["label"] != 'B53' :
        #                             liaison_autre_que_b53 = True
        #                 if noeud == (5,5) :
        #                     print("petit rat")
        #                     print(len(voisins))
        #                     print(liaison_autre_que_b53)           
                        if len(voisins) >= 1  : ## si pas de voisin, on l enleve
                            compteur = 1
                            for ch in chaines :
                                if noeud in ch :
                                    espacement_motif.append(abs(graphe.nodes[noeud]["position"][0]-graphe.nodes[compteur]["position"][0]))
                                compteur += 1
                 
                            for i in range(1,len(cluster)) :
                            
                                with open(rep_extension+"/"+"fichier_"+cluster[i][0]+"_"+str(cluster[i][1])+"_2.pickle", 'rb') as fichier_ext_2 :
                                    mon_depickler_ext_2 = pickle.Unpickler(fichier_ext_2)
                                    graphe2 = mon_depickler_ext_2.load()
                                    
                                    
                                    if (cluster[0], cluster[i]) in dico_comp.keys() :
                                        place_i = 1
                                    else :
                                        place_i = 0
                                       
                                    for noeud2 in dico_comp[(cluster[abs(place_i-1)*i], cluster[place_i*i])]["graphe"].nodes() :
                                        if noeud == noeud2[abs(place_i-1)] :
                                            graphe_commun.nodes[noeud]["num_seq"].append((cluster[i], graphe2.nodes[noeud2[place_i]]["position"], graphe2.nodes[noeud2[place_i]]["num_ch"]))
                                            
                                            for elt in chaine :
                                                if elt not in graphe2.nodes[noeud2[place_i]]["chaine"] :
    
                                                    a_enlever.append(elt)
#                                                     print("gros tas")
#                                                     print(cluster[i])
#                                                     print(elt)
#                                                     print(noeud2[place_i])
#                                                     print(graphe2.nodes[noeud2[place_i]]["chaine"])
                                            if mini_poids > graphe2.nodes[noeud2[place_i]]["poids"] :
                                                mini_poids = graphe2.nodes[noeud2[place_i]]["poids"]
                                                #print(mini_poids)
                                                #print(cluster[0])
                                                #print(cluster[i])
                                            if graphe2.nodes[noeud2[place_i]]["poids"] not in liste_poids :
                                                liste_poids.append(graphe2.nodes[noeud2[place_i]]["poids"])
                                            #print("tout petit rat")  
                                            compteur = 1
                #                                     print("rat")
                #                                     print(noeud[place_0])
                #                                     print(chaines)
                                            for ch in chaines :
                                                if noeud in ch :
                                                    espacement_motif.append(abs(graphe2.nodes[noeud2[place_i]]["position"][0]-graphe2.nodes[compteur]["position"][0]))
                                                     
                                                compteur += 1
                                            #print(noeud)
                                            print(cluster[i])
                                            print(noeud2[place_i])
                                            print(cluster)
                                            liaison_b53_a_mettre.append(suivant_existe(noeud2[place_i], graphe2, liste_cliques, i)) 
                #                     print("ramous")
                #                     print(espacement_motif)
                            graphe_commun.nodes[noeud]["espacement_motif"] = list(espacement_motif)
                #                     print(a_enlever)
                            
                            for elt in a_enlever :
                                if elt in chaine :
                                    chaine.remove(elt)
#                             print("chaines2")
#                             print(chaine)
                            #if noeud == 20  :
#                             print(graphe_commun.nodes[noeud])
                            #print("gros rat")
                            #print(liste_poids)
                #                 a_enlever = []
                #                 if len(liste_poids) > 1 :
                #                     for u,v,key,data in digraphe_commun.edges(data=True, keys=True) :
                #                         if u == noeud :
                #                             if data["label"] == 'B53' :
                #                                 a_enlever.append((u,v,key))
                #                 print("gros rat")
                #                 print(a_enlever)
                #                 for elt in a_enlever :
                #                     digraphe_commun.remove_edge(elt[0], elt[1], key=elt[2])
                             
                 
                            graphe_commun.nodes[noeud]["poids"] = mini_poids          
                            graphe_commun.nodes[noeud]["chaine"] = list(chaine)
                            
                            liaison_b53_a_ajouter = True
                            for elt in liaison_b53_a_mettre :
                                if elt == -1  :
                                    liaison_b53_a_ajouter = False
                            #print("ramou")
                            #print(liaison_b53_a_mettre)
                            #print(liaison_b53_a_ajouter)
                            if liaison_b53_a_ajouter : 
                                graphe_commun.add_edge(noeud, liaison_b53_a_mettre[0], label='B53', long_range=False)
                            
                            chaine_position = []
                            for i in range(4) :
                                for j in range(len(chaines[i])) :
                                    if noeud == chaines[i][j] :
                #                             print(num_noeud)
                #                             print(chaines[i][j])
                                        if i == 0 or i == 3 :
                                            chaine_position.append(j)
                                        else :
                                            chaine_position.append(10-j)
                            graphe_commun.nodes[noeud]["position"] = chaine_position
                                
                        else : 
                            noeuds_isoles_a_enlever.append(noeud)
                for elt in noeuds_isoles_a_enlever :
                    graphe_commun.remove_node(elt)
                for i in range(1,6) :
                    graphe_commun.nodes[noeud]["positon"] = [i]
        #       
                
          
                print(graphe_commun.nodes.data())
                print(graphe_commun.number_of_nodes())
                for noeud, data in graphe_commun.nodes(data=True) :
                    print(noeud, data)
                #exit()
                return graphe_commun, liste_cliques  
#         else :
#             if (cluster[0], cluster[1]) in dico_comp.keys() :
#                 return dico_comp[(cluster[0], cluster[1])]["graphe"], []
#             else :
#                 return dico_comp[(cluster[1], cluster[0])]["graphe"], []
        

''' recherche le plus grand sous-graphe commun dans un cluster V1 (en prenant seulement les sommets superposes entre deux graphes, puis en regardant quels sommets gardes avec un troisieme graphe etc...)
(mais le resultat change avec l'ordre de traitement des graphes du cluster)
(version CaRNAval)
rep_comparaison : chemin du fichier ou se trouve le dictionnaire des graphes communs
rep_extension : chemin du repertoire ou se trouvent les fichiers de graphes d'extension
rep_sauve : chemin du repertoire ou stocker le graphe commun final 
num : numero du cluster 
seuil : seuil de sim considere
taille_comp : taille de l'extension
commut : si vrai, les noeuds du graphe commun sont convertis en entiers de 1 a la taille du graphe '''
def commun_cluster(cluster, rep_extension, rep_comparaison, num, taille_comp, seuil, rep_sauve, commut): ## cluster : liste de nom de graphes
    graphe_commun = nx.MultiDiGraph()
    with open(rep_comparaison, 'rb') as fichier_comp :
        mon_depickler = pickle.Unpickler(fichier_comp)
        dico_comp = mon_depickler.load()
        
#         print(dico_comp.keys())
        place_0 = -1
        if ("fichier_"+cluster[0], "fichier_"+cluster[1]) in dico_comp.keys() :
            #graphe_commun.add_nodes_from(dico_comp[(cluster[0], cluster[1])].nodes(data=True))
            graphe_commun = dico_comp[("fichier_"+cluster[0], "fichier_"+cluster[1])].copy()
#             print("fichier_"+cluster[0])
#             print("fichier_"+cluster[1])
            place_0 = 0
        else :
            #graphe_commun.add_nodes_from(dico_comp[(cluster[1], cluster[0])].nodes(data=True))
            graphe_commun = dico_comp[("fichier_"+cluster[1], "fichier_"+cluster[0])].copy()
#             print("fichier_"+cluster[1])
#             print("fichier_"+cluster[0])
            place_0 = 1
        
#         print(graphe_commun.nodes.data())
#         print(graphe_commun.edges.data())
            
        for i in range(2,len(cluster)) :
            
            a_enlever_noeuds = []
            a_enlever_aretes = []
            if ("fichier_"+cluster[0], "fichier_"+cluster[i]) in dico_comp.keys() :
                place_i = 1
            else :
                place_i = 0
#             print("fichier_"+cluster[abs(place_i-1)*i])
#             print("fichier_"+cluster[place_i*i])
#             print(dico_comp[("fichier_"+cluster[abs(place_i-1)*i], "fichier_"+cluster[place_i*i])].nodes.data())
#             print(dico_comp[("fichier_"+cluster[abs(place_i-1)*i], "fichier_"+cluster[place_i*i])].edges.data())
#             print(a_enlever_aretes)
            for noeud1 in graphe_commun.nodes() :
                noeud1_existe = False
                for noeud2 in  dico_comp[("fichier_"+cluster[abs(place_i-1)*i], "fichier_"+cluster[place_i*i])].nodes() :
                    if noeud1[place_0] == noeud2[abs(place_i-1)] :
                        noeud1_existe = True
                        
                if noeud1_existe == False :
                    a_enlever_noeuds.append(noeud1)
#                     for edge in graphe_commun.in_edges(noeud1) :
#                         if edge not in a_enlever_aretes :
#                             a_enlever_aretes.append(edge)
#                     for edge in graphe_commun.out_edges(noeud1) :
#                         if edge not in a_enlever_aretes :
#                             a_enlever_aretes.append(edge)
                    for edge in graphe_commun.edges(noeud1) :
                        a_enlever_aretes.append((edge))
                else : 
                    for voisin in graphe_commun[noeud1] :
                        voisin_existe = False
                        for noeud2 in  dico_comp[("fichier_"+cluster[abs(place_i-1)*i], "fichier_"+cluster[place_i*i])].nodes() :
                            if voisin[place_0] == noeud2[abs(place_i-1)] :
                                voisin_existe = True
                        if voisin_existe == False :
                            a_enlever_noeuds.append(voisin)
                            
                            for edge in graphe_commun.edges(voisin) :
                                a_enlever_aretes.append((edge))
#                                 print(edge)
                        else :
                            liaison_existe = False
                            for u,v,e in dico_comp[("fichier_"+cluster[abs(place_i-1)*i], "fichier_"+cluster[place_i*i])].edges(keys=True) :
                                if (u[abs(place_i-1)] == voisin[place_0] and v[abs(place_i-1)] == noeud1[place_0]) or (v[abs(place_i-1)] == voisin[place_0] and u[abs(place_i-1)] == noeud1[place_0]) : 
#                                     print(noeud1[place_0])
#                                     print(voisin[place_0])
                                    
#                                     if dico_comp[("fichier_"+cluster[abs(place_i-1)*i], "fichier_"+cluster[place_i*i])][u][v][e]["type"] != 'B53' :
#                                     
#                                         if (noeud1[place_0],voisin[place_0]) in graphe1.edges() :
#                                             for edge in graphe1[noeud1[place_0]][voisin[place_0]] :
#                                                 if graphe1[noeud1[place_0]][voisin[place_0]][edge]["label"] != 'B53' :
#                                                     label = graphe1[noeud1[place_0]][voisin[place_0]][edge]["label"]
#                                             print(graphe_commun[noeud1][voisin])
#                                             for edge in graphe_commun[noeud1][voisin] :
#                                                 if label == graphe_commun[noeud1][voisin][edge]["type"] :
#                                                     liaison_existe = True
#                                     else:
                                    liaison_existe = True
                                        
                                    
                            if liaison_existe == False :
                                a_enlever_aretes.append((noeud1, voisin))
#                                 print((noeud1, voisin))
                            
#             print(a_enlever_aretes)
#             print("ramou")
#             print(graphe_commun.edges())
            for elt in a_enlever_aretes :
                if elt in graphe_commun.edges() :
                    graphe_commun.remove_edge(elt[0], elt[1])
            for elt in a_enlever_noeuds :
                if elt in graphe_commun.nodes() :
                    graphe_commun.remove_node(elt)
            
#         print(graphe_commun.nodes.data())
#         print(graphe_commun.edges.data())
        
        with open(EXTENSION_PATH%rep_sauve+"/graphe_commun_%s_%s_taille_%s.pickle"%(seuil, num, taille_comp), 'wb') as fichier : 
            mon_pickler = pickle.Pickler(fichier)
            mon_pickler.dump(graphe_commun)

        with open(rep_extension+"/"+"fichier_"+cluster[0]+".pickle", 'rb') as fichier_ext :
            mon_depickler_ext = pickle.Unpickler(fichier_ext)
            graphe = mon_depickler_ext.load()
            
            
            digraphe_commun = nx.MultiDiGraph()
            digraphe_commun.add_nodes_from(graphe_commun.nodes(data=True))
            
            for u,v,e,data in graphe_commun.edges(keys=True, data=True) :
                label = ""
                if data["type"] != 'B53' : 
## plus besoin, comme maintenant les labels des graphes communs sont les memes que ceux des graphes d extensions
#                     if data["type"] == 'CAN' :
#                         label = 'CWW'
#                     elif data["type"] == 'ART' :
#                         label = '0'
#                     else :
#                         label = 'NON_CAN'
                    label = data["type"]
                    del(data["type"])
                    data.update({"label" : label})
                    if data["label"] == '0' :
                        digraphe_commun.add_edge(u, v, **data)
                        
                        ## si graphe commun pas multigraph
                        #digraphe_commun.add_edge(v, u, **data)
                    else :   
#                         print(u[place_0])
#                         print(v[place_0])
#                         print(e)
                        
                        ## si graphe commun pas multigraph :
#                         if len(graphe[u[place_0]][v[place_0]]) > e :
#                             if graphe[u[place_0]][v[place_0]][e]["label"] == label : 
#                                 digraphe_commun.add_edge(u, v, **data)
#                                 label = label[0] + label[2] + label[1] + label[3:]
#                                 data["label"] = label
#                                 digraphe_commun.add_edge(v, u, **data)
#                             elif graphe[v[place_0]][u[place_0]][e]["label"] == label :
#                                 digraphe_commun.add_edge(v, u, **data)
#                                 label = label[0] + label[2] + label[1] + label[3:]
#                                 data["label"] = label
#                                 digraphe_commun.add_edge(u, v, **data)
#                             else :
#                                 print("bizarre")

                        ## si graphe commun deja multigraph :
                            digraphe_commun.add_edge(u, v, **data)
                            
                    
                else :
                    label = 'B53'
                    del(data["type"])
                    data.update({"label" : label})
#                     print("petit rat")
#                     print((u[place_0], v[place_0]))
                    if (u[place_0], v[place_0]) in graphe.edges() :
                        est_dans_le_bon_sens = False
                        for edge in graphe[u[place_0]][v[place_0]] :
                            if graphe[u[place_0]][v[place_0]][edge]["label"] == 'B53' :
                                est_dans_le_bon_sens = True
                        if est_dans_le_bon_sens :
                            digraphe_commun.add_edge(u, v, **data)
                        else :
                            digraphe_commun.add_edge(v, u, **data)
                    else :
                        digraphe_commun.add_edge(v, u, **data)
#                     print(digraphe_commun.edges())
#             print(graphe[27])
            chaines = [[1]]
            for i in range(1,5) :
                    compteur = i
                    if i != 1 : chaines.append([i])
                    liaison_B53 = True
                    while liaison_B53 :
                        liaison_B53 = False
                        temp = compteur
                        if i == 1 or i == 4 :
                            for voisin in graphe.successors(compteur) :
                                for arc in graphe[compteur][voisin] :
                                    if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and graphe[compteur][voisin][arc]["label"] == 'B53' :
                                        liaison_B53 = True
                                        temp = voisin
                                        chaines[len(chaines)-1].append(voisin)
                        else :            
                            for voisin in graphe.predecessors(compteur) :
                                for arc in graphe[voisin][compteur] :
                                    if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and graphe[voisin][compteur][arc]["label"] == 'B53' :
                                        liaison_B53 = True
                                        temp = voisin
                                        chaines[len(chaines)-1].append(voisin)
                        compteur = temp

            
            nx.set_node_attributes(digraphe_commun, -2, "type")
            nx.set_node_attributes(digraphe_commun, -1, "poids")
            nx.set_node_attributes(digraphe_commun, [], "chaine")
            nx.set_node_attributes(digraphe_commun, [], "position")
            nx.set_node_attributes(digraphe_commun, [], "espacement_motif")
            noeuds_isoles_a_enlever = []
            for noeud in digraphe_commun.nodes() :
                a_enlever = []
                espacement_motif = []
                num_noeud = noeud[place_0]
                digraphe_commun.nodes[noeud]["type"] = graphe.nodes[num_noeud]["type"]
                
                mini_poids = graphe.nodes[num_noeud]["poids"]
                chaine = graphe.nodes[num_noeud]["chaine"]
                liste_poids = []
                
                voisins = digraphe_commun[noeud]
#                 print(voisins)
                liaison_autre_que_b53 = False
                if len(voisins) == 1 :
                    for edge in digraphe_commun[noeud][next(iter(voisins))] :
                        if digraphe_commun[noeud][next(iter(voisins))][edge]["label"] != 'B53' :
                            liaison_autre_que_b53 = True
                if noeud == (5,5) :
                    print("petit rat")
                    print(len(voisins))
                    print(liaison_autre_que_b53)           
                if len(voisins) > 1 or (len(voisins) == 1 and liaison_autre_que_b53) : ## si pas de voisin, on l enleve
                    compteur = 1
                    for ch in chaines :
                        if noeud[place_0] in ch :
                            espacement_motif.append(abs(graphe.nodes[noeud[place_0]]["position"][0]-graphe.nodes[compteur]["position"][0]))
                        compteur += 1
    
                    for i in range(1,len(cluster)) :
                        with open(rep_extension+"/"+"fichier_"+cluster[i]+".pickle", 'rb') as fichier_ext_2 :
                            mon_depickler_ext_2 = pickle.Unpickler(fichier_ext_2)
                            graphe2 = mon_depickler_ext_2.load()
                            if ("fichier_"+cluster[0], "fichier_"+cluster[i]) in dico_comp.keys() :
                                place_i = 1
                            else :
                                place_i = 0
                                
                            for noeud2 in dico_comp[("fichier_"+cluster[abs(place_i-1)*i], "fichier_"+cluster[place_i*i])].nodes() :
                                if noeud[place_0] == noeud2[abs(place_i-1)] :
                                    for elt in chaine :
                                        if elt not in graphe2.nodes[noeud2[place_i]]["chaine"] :
                                            a_enlever.append(elt)
                                    if mini_poids > graphe2.nodes[noeud2[place_i]]["poids"] :
                                        mini_poids = graphe2.nodes[noeud2[place_i]]["poids"]
                                        print(mini_poids)
                                        print(cluster[0])
                                        print(cluster[i])
                                    if graphe2.nodes[noeud2[place_i]]["poids"] not in liste_poids :
                                        liste_poids.append(graphe2.nodes[noeud2[place_i]]["poids"])
                                    
                                    compteur = 1
#                                     print("rat")
#                                     print(noeud[place_0])
#                                     print(chaines)
                                    for ch in chaines :
                                        if noeud[place_0] in ch :
                                            espacement_motif.append(abs(graphe2.nodes[noeud2[place_i]]["position"][0]-graphe2.nodes[compteur]["position"][0]))
                                            
                                        compteur += 1
#                     print("ramous")
#                     print(espacement_motif)
                    digraphe_commun.nodes[noeud]["espacement_motif"] = list(espacement_motif)
#                     print(a_enlever)
                    for elt in a_enlever :
                        if elt in chaine :
                            chaine.remove(elt)
                    #print("gros rat")
                    #print(liste_poids)
    #                 a_enlever = []
    #                 if len(liste_poids) > 1 :
    #                     for u,v,key,data in digraphe_commun.edges(data=True, keys=True) :
    #                         if u == noeud :
    #                             if data["label"] == 'B53' :
    #                                 a_enlever.append((u,v,key))
    #                 print("gros rat")
    #                 print(a_enlever)
    #                 for elt in a_enlever :
    #                     digraphe_commun.remove_edge(elt[0], elt[1], key=elt[2])
                    

                    digraphe_commun.nodes[noeud]["poids"] = mini_poids          
                    digraphe_commun.nodes[noeud]["chaine"] = chaine
                    
                    
                    chaine_position = []
                    for i in range(4) :
                        for j in range(len(chaines[i])) :
                            if num_noeud == chaines[i][j] :
    #                             print(num_noeud)
    #                             print(chaines[i][j])
                                if i == 0 or i == 3 :
                                    chaine_position.append(j)
                                else :
                                    chaine_position.append(10-j)
                    digraphe_commun.nodes[noeud]["position"] = chaine_position
                else : 
                    noeuds_isoles_a_enlever.append(noeud)
            
#             print(chaines)    
#             print(digraphe_commun.nodes.data())   
#             print(digraphe_commun.edges.data())
            print(noeuds_isoles_a_enlever)
            for elt in noeuds_isoles_a_enlever :
                digraphe_commun.remove_node(elt)
            
            
            
#             print(digraphe_commun.nodes.data())
#             print(digraphe_commun.edges.data()) 
            
            if commut :
                digraphe_commun_integer = nx.convert_node_labels_to_integers(digraphe_commun, first_label=1, ordering='default')
            else :
                digraphe_commun_integer = digraphe_commun.copy()
            with open(EXTENSION_PATH%rep_sauve+"digraphe_commun_%s_%s_taille_%s.pickle"%(seuil, num, taille_comp), 'wb') as fichier : 
                mon_pickler = pickle.Pickler(fichier)
                mon_pickler.dump(digraphe_commun_integer)      
            
#             print(espacement_motif)
                
            return digraphe_commun_integer     


''' recherche le plus grand sous-graphe commun dans un cluster V1 (en prenant seulement les sommets superposes entre deux graphes, puis en regardant quels sommets gardes avec un troisieme graphe etc...)
(mais le resultat change avec l'ordre de traitement des graphes du cluster)
(version CaRNAval) 
difference avec la fonction du dessus : on recherche un sous-graphe commun connexe autour du motif
rep_comparaison : chemin du fichier ou se trouve le dictionnaire des graphes communs
rep_extension : chemin du repertoire ou se trouvent les fichiers de graphes d'extension
num : numero du cluster 
taille_comp : taille de l'extension'''
def commun_cluster_raccourci(cluster, rep_extension, rep_comparaison, num, taille_comp): ## cluster : liste de nom de graphes
    graphe_commun = nx.MultiDiGraph()
    with open(rep_comparaison, 'rb') as fichier_comp :
        mon_depickler = pickle.Unpickler(fichier_comp)
        dico_comp = mon_depickler.load()
        
        print(dico_comp.keys())
        place_0 = -1
        if ("fichier_"+cluster[0], "fichier_"+cluster[1]) in dico_comp.keys() :
            #graphe_commun.add_nodes_from(dico_comp[(cluster[0], cluster[1])].nodes(data=True))
            graphe_commun = dico_comp[("fichier_"+cluster[0], "fichier_"+cluster[1])].copy()
            print(graphe_commun.nodes.data())
            place_0 = 0
        else :
            #graphe_commun.add_nodes_from(dico_comp[(cluster[1], cluster[0])].nodes(data=True))
            graphe_commun = dico_comp[("fichier_"+cluster[1], "fichier_"+cluster[0])].copy()
            place_0 = 1
            print(graphe_commun.nodes.data())
         
            
        for i in range(2,len(cluster)) :
            a_enlever_noeuds = []
            a_enlever_aretes = []
            if ("fichier_"+cluster[0], "fichier_"+cluster[i]) in dico_comp.keys() :
                place_i = 1
            else :
                place_i = 0
            print("fichier_"+cluster[abs(place_i-1)*i])
            print("fichier_"+cluster[place_i*i])
            print(dico_comp[("fichier_"+cluster[abs(place_i-1)*i], "fichier_"+cluster[place_i*i])].nodes.data())
            print(dico_comp[("fichier_"+cluster[abs(place_i-1)*i], "fichier_"+cluster[place_i*i])].edges.data())
            print(a_enlever_aretes)
            for noeud1 in graphe_commun.nodes() :
                noeud1_existe = False
                for noeud2 in  dico_comp[("fichier_"+cluster[abs(place_i-1)*i], "fichier_"+cluster[place_i*i])].nodes() :
                    if noeud1[place_0] == noeud2[abs(place_i-1)] :
                        noeud1_existe = True
                        
                if noeud1_existe == False :
                    a_enlever_noeuds.append(noeud1)
#                     for edge in graphe_commun.in_edges(noeud1) :
#                         if edge not in a_enlever_aretes :
#                             a_enlever_aretes.append(edge)
#                     for edge in graphe_commun.out_edges(noeud1) :
#                         if edge not in a_enlever_aretes :
#                             a_enlever_aretes.append(edge)
                    for edge in graphe_commun.edges(noeud1) :
                        a_enlever_aretes.append((edge))
                else : 
                    for voisin in graphe_commun[noeud1] :
                        voisin_existe = False
                        for noeud2 in  dico_comp[("fichier_"+cluster[abs(place_i-1)*i], "fichier_"+cluster[place_i*i])].nodes() :
                            if voisin[place_0] == noeud2[abs(place_i-1)] :
                                voisin_existe = True
                        if voisin_existe == False :
                            a_enlever_noeuds.append(voisin)
                            
                            for edge in graphe_commun.edges(voisin) :
                                a_enlever_aretes.append((edge))
                                print(edge)
                        else :
                            liaison_existe = False
                            for u,v in dico_comp[("fichier_"+cluster[abs(place_i-1)*i], "fichier_"+cluster[place_i*i])].edges() :
                                if (u[abs(place_i-1)] == voisin[place_0] and v[abs(place_i-1)] == noeud1[place_0]) or (v[abs(place_i-1)] == voisin[place_0] and u[abs(place_i-1)] == noeud1[place_0]) : 
                                    liaison_existe = True
                                    
                            if liaison_existe == False :
                                a_enlever_aretes.append((noeud1, voisin))
                                print((noeud1, voisin))
                            
            print(a_enlever_aretes)
            print(graphe_commun.edges())
            for elt in a_enlever_aretes :
                if elt in graphe_commun.edges() :
                    graphe_commun.remove_edge(elt[0], elt[1])
            for elt in a_enlever_noeuds :
                if elt in graphe_commun.nodes() :
                    graphe_commun.remove_node(elt)
                    
        print(graphe_commun.nodes.data())
        print(graphe_commun.edges.data())
        
#         with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/result_graphes_extension_taille_max/graphe_commun_0.7_55_%s_taille_%s.pickle"%(num, taille_comp), 'wb') as fichier : 
#             mon_pickler = pickle.Pickler(fichier)
#             mon_pickler.dump(graphe_commun)
        
        
            
       
        with open(rep_extension+"/"+"fichier_"+cluster[0]+".pickle", 'rb') as fichier_ext :
            mon_depickler_ext = pickle.Unpickler(fichier_ext)
            graphe = mon_depickler_ext.load()
            
            
            digraphe_commun = nx.MultiDiGraph()
            digraphe_commun.add_nodes_from(graphe_commun.nodes(data=True))
            
            for u,v,data in graphe_commun.edges(data=True) :
                label = ""
                if data["type"] != 'COV' : 
                    if data["type"] == 'CAN' :
                        label = 'CWW'
                    elif data["type"] == 'ART' :
                        label = '0'
                    else :
                        label = 'NON_CAN'
                        
                    del(data["type"])
                    data.update({"label" : label})    
                    digraphe_commun.add_edge(u, v, **data)
                    digraphe_commun.add_edge(v, u, **data)
                else :
                    label = 'B53'
                    del(data["type"])
                    data.update({"label" : label})
                    if (u[place_0], v[place_0]) in graphe.edges() :
                        est_dans_le_bon_sens = False
                        for edge in graphe[u[place_0]][v[place_0]] :
                            if graphe[u[place_0]][v[place_0]][edge]["label"] == 'B53' :
                                est_dans_le_bon_sens = True
                        if est_dans_le_bon_sens :
                            digraphe_commun.add_edge(u, v, **data)
                        else :
                            digraphe_commun.add_edge(v, u, **data)
                    else :
                        digraphe_commun.add_edge(v, u, **data)
            
#             print(graphe[27])
            chaines = [[1]]
            for i in range(1,5) :
                    compteur = i
                    if i != 1 : chaines.append([i])
                    liaison_B53 = True
                    while liaison_B53 :
                        liaison_B53 = False
                        temp = compteur
                        if i == 1 or i == 4 :
                            for voisin in graphe.successors(compteur) :
                                for arc in graphe[compteur][voisin] :
                                    if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and graphe[compteur][voisin][arc]["label"] == 'B53' :
                                        liaison_B53 = True
                                        temp = voisin
                                        chaines[len(chaines)-1].append(voisin)
                        else :            
                            for voisin in graphe.predecessors(compteur) :
                                for arc in graphe[voisin][compteur] :
                                    if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and graphe[voisin][compteur][arc]["label"] == 'B53' :
                                        liaison_B53 = True
                                        temp = voisin
                                        chaines[len(chaines)-1].append(voisin)
                        compteur = temp

           
            nx.set_node_attributes(digraphe_commun, -2, "type")
            nx.set_node_attributes(digraphe_commun, -1, "poids")
            nx.set_node_attributes(digraphe_commun, [], "chaine")
            nx.set_node_attributes(digraphe_commun, [], "position")
            nx.set_node_attributes(digraphe_commun, [], "espacement_motif")
            a_enlever_n = []
            for noeud in digraphe_commun.nodes() :
                a_enlever = []
                
                espacement_motif = []
                num_noeud = noeud[place_0]
                digraphe_commun.nodes[noeud]["type"] = graphe.nodes[num_noeud]["type"]
                
                mini_poids = graphe.nodes[num_noeud]["poids"]
                chaine = graphe.nodes[num_noeud]["chaine"]
                liste_poids = []
                
                compteur = 1
                for ch in chaines :
                    if noeud[place_0] in ch :
                        espacement_motif.append(abs(graphe.nodes[noeud[place_0]]["position"][0]-graphe.nodes[compteur]["position"][0]))
                    compteur += 1

                for i in range(1,len(cluster)) :
                    with open(rep_extension+"/"+"fichier_"+cluster[i]+".pickle", 'rb') as fichier_ext_2 :
                        mon_depickler_ext_2 = pickle.Unpickler(fichier_ext_2)
                        graphe2 = mon_depickler_ext_2.load()
                        if ("fichier_"+cluster[0], "fichier_"+cluster[i]) in dico_comp.keys() :
                            place_i = 1
                        else :
                            place_i = 0
                            
                        for noeud2 in dico_comp[("fichier_"+cluster[abs(place_i-1)*i], "fichier_"+cluster[place_i*i])].nodes() :
                            if noeud[place_0] == noeud2[abs(place_i-1)] :
                                for elt in chaine :
                                    if elt not in graphe2.nodes[noeud2[place_i]]["chaine"] :
                                        a_enlever.append(elt)
                                if mini_poids > graphe2.nodes[noeud2[place_i]]["poids"] :
                                    mini_poids = graphe2.nodes[noeud2[place_i]]["poids"]
                                if graphe2.nodes[noeud2[place_i]]["poids"] not in liste_poids :
                                    liste_poids.append(graphe2.nodes[noeud2[place_i]]["poids"])
                                
                                compteur = 1
                                print("rat")
                                print(noeud[place_0])
                                print(chaines)
                                for ch in chaines :
                                    if noeud[place_0] in ch :
                                        print()
                                        espacement_motif.append(abs(graphe2.nodes[noeud2[place_i]]["position"][0]-graphe2.nodes[compteur]["position"][0]))
                                        
                                    compteur += 1
                print("ramous")
                print(espacement_motif)
                
                if len(espacement_motif) > 0 :
                    val = espacement_motif[0]
                    for elt in espacement_motif :
                        if elt != val :
                            if noeud not in a_enlever_n : 
                                a_enlever_n.append(noeud)
                            
                digraphe_commun.nodes[noeud]["espacement_motif"] = list(espacement_motif)
                for elt in a_enlever :
                    chaine.remove(elt)
                #print("gros rat")
                #print(liste_poids)
                a_enlever = []
                if len(liste_poids) > 1 :
                    for u,v,key,data in digraphe_commun.edges(data=True, keys=True) :
                        if u == noeud :
                            if data["label"] == 'B53' :
                                a_enlever.append((u,v,key))
                #print(a_enlever)
                for elt in a_enlever :
                    digraphe_commun.remove_edge(elt[0], elt[1], key=elt[2])
                
                digraphe_commun.nodes[noeud]["poids"] = mini_poids          
                digraphe_commun.nodes[noeud]["chaine"] = chaine
                
                
                chaine_position = []
                for i in range(4) :
                    for j in range(len(chaines[i])) :
                        if num_noeud == chaines[i][j] :
#                             print(num_noeud)
#                             print(chaines[i][j])
                            if i == 0 or i == 3 :
                                chaine_position.append(j)
                            else :
                                chaine_position.append(10-j)
                digraphe_commun.nodes[noeud]["position"] = chaine_position
#            
#             print(chaines)    
#             print(digraphe_commun.nodes.data())   
#             print(digraphe_commun.edges.data())
            print(a_enlever_n)
            for elt in a_enlever_n :
                digraphe_commun.remove_node(elt)
            print(digraphe_commun.nodes.data())
            print(digraphe_commun.edges.data()) 
            with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/result_graphes_extension_taille_max/digraphe_commun_0.7_55_%s_taille_%s_raccourci.pickle"%(num, taille_comp), 'wb') as fichier : 
                mon_pickler = pickle.Pickler(fichier)
                mon_pickler.dump(digraphe_commun)      
            
            print(espacement_motif)
                
            return digraphe_commun      

''' petite fonction utilitaire pour afficher des aretes courbes
G : graphe d'extension
edges_sans_labels : aretes a mettre en courbe
pos : coordonnees des noeuds '''
def draw_network(G,edges_sans_label, pos,ax,sg=None):
    
    #print(edges_sans_label)
    nodes = []
    for elt in edges_sans_label :
        if elt[0] not in nodes :
            nodes.append(elt[0])
        if elt[1] not in nodes :
            nodes.append(elt[1])
    
    for n in nodes:
        c=Circle(pos[n],radius=0.02,alpha=0.5)
        ax.add_patch(c)
        G.node[n]['patch']=c
        x,y=pos[n]
    seen={}
    lab = 0
    for u,v,data in G.edges(data=True) :
        #print(u)
        #print(v)
        #lab = edges_sans_label.index((u,v))
        try :
            lab = edges_sans_label.index((u,v))
        except ValueError :
            lab = None
        #print(lab)
        #print(data["label"])
        if lab != None :

            n1=G.node[u]['patch']
            n2=G.node[v]['patch']
            rad=0.1
            if (u,v) in seen:
                rad=seen.get((u,v))
                rad=(rad+np.sign(rad)*0.1)*-1
            alpha=1.0
            #color='k'
            if data["label"] == 'CWW' :
                color = 'blue'
            elif data["label"] == '0':
                color = 'blue'
            elif data["label"] == 'B53' :
                color = 'green'
            else :
                color = 'black'
            
            e = FancyArrowPatch(n1.center,n2.center,patchA=n1,patchB=n2,
                                arrowstyle='<|-|>',
                                connectionstyle='arc3,rad=%s'%rad,
                                mutation_scale=10.0,
                                lw=2,
                                alpha=alpha,
                                edgecolor=color, 
                                label=data["label"])
#             if G.nodes[u]["coordonnees"][0] == G.nodes[v]["coordonnees"][0] :
#                 if G.nodes[u]["coordonnees"][1] > G.nodes[v]["coordonnees"][1] :
#                     plt.text(G.nodes[u]["coordonnees"][0]-0.15, G.nodes[v]["coordonnees"][1] + 0.25, data["label"], fontsize=8)
#                 else :
#                     plt.text(G.nodes[u]["coordonnees"][0]-0.15, G.nodes[u]["coordonnees"][1] + 0.25, data["label"], fontsize=8)
#             elif G.nodes[u]["coordonnees"][1] == G.nodes[v]["coordonnees"][1] :
#                 if G.nodes[u]["coordonnees"][0] > G.nodes[v]["coordonnees"][0] :
#                     plt.text(G.nodes[v]["coordonnees"][0]+0.25, G.nodes[u]["coordonnees"][1]-0.15, data["label"], fontsize=8)
#                 else :    
#                     plt.text(G.nodes[u]["coordonnees"][0]+0.25, G.nodes[u]["coordonnees"][1]-0.15, data["label"], fontsize=8)
            seen[(u,v)]=rad
            ax.add_patch(e)
    return e

''' affiche le sous-graphe commun G et le stocke en format png dans rep_sauve
(version CaRNAval) '''
def draw(G, rep_sauve, num, taille_comp, commut):

    print(G.nodes.data())
    
    nx.set_node_attributes(G, (33,33), "coordonnees")
    if commut :
        G.nodes[1]["coordonnees"] = (0.0,0.5)
        G.nodes[2]["coordonnees"] = (2.0,0.5)
        G.nodes[3]["coordonnees"] = (0.0,0.0)
        G.nodes[4]["coordonnees"] = (2.0,0.0)
        G.nodes[5]["coordonnees"] = (3.0,0.5)
    else :
        G.nodes[(1,1)]["coordonnees"] = (0.0,0.5)
        G.nodes[(2,2)]["coordonnees"] = (2.0,0.5)
        G.nodes[(3,3)]["coordonnees"] = (0.0,0.0)
        G.nodes[(4,4)]["coordonnees"] = (2.0,0.0)
        G.nodes[(5,5)]["coordonnees"] = (3.0,0.5)
     
#                 fichier.write(str(element)+"\n") 
#                 fichier.write(str(G.number_of_nodes())+"\n") 
    print(G.edges.data())

    nodes_list = [u for u,d in G.nodes(data=True)]#and len(G[u]) > 0] 
    print(nodes_list)

    coordonnees = []
    for noeud in nodes_list :
        #voisins = G[noeud]
        print(noeud)
        if noeud not in [1,2,3,4,5] :
            if G.nodes[noeud]["type"] != None and G.nodes[noeud]["type"] != -1 :
                    chaine = G.nodes[noeud]["chaine"][0]
                    #print(chaine)
                    if commut : 
                        if chaine == 1 or chaine == 3 :
                            G.nodes[noeud]["coordonnees"] =  (G.nodes[chaine]["coordonnees"][0], G.nodes[chaine]["coordonnees"][1] + (G.nodes[noeud]["position"][0] - G.nodes[chaine]["position"][0])/2) 
                        else :
                            G.nodes[noeud]["coordonnees"] =  (G.nodes[chaine]["coordonnees"][0], G.nodes[chaine]["coordonnees"][1] + (G.nodes[chaine]["position"][0] - G.nodes[noeud]["position"][0])/2) 
                    else :
                        if chaine == 1 or chaine == 3 :
                            G.nodes[noeud]["coordonnees"] =  (G.nodes[(chaine,chaine)]["coordonnees"][0], G.nodes[(chaine,chaine)]["coordonnees"][1] + (G.nodes[noeud]["position"][0] - G.nodes[(chaine,chaine)]["position"][0])/2) 
                        else :
                            G.nodes[noeud]["coordonnees"] =  (G.nodes[(chaine,chaine)]["coordonnees"][0], G.nodes[(chaine,chaine)]["coordonnees"][1] + (G.nodes[(chaine,chaine)]["position"][0] - G.nodes[noeud]["position"][0])/2) 
        
                    coordonnees.append(G.nodes[noeud]["coordonnees"])
            else :
                    voisin = list(G[noeud])
                    print(len(voisin))
                    print(voisin)
                    if (G.nodes[voisin[0]]["coordonnees"][0] + 0.5, G.nodes[voisin[0]]["coordonnees"][1]) not in coordonnees :
                        G.nodes[noeud]["coordonnees"] = (G.nodes[voisin[0]]["coordonnees"][0] + 0.5, G.nodes[voisin[0]]["coordonnees"][1])
                        coordonnees.append(G.nodes[noeud]["coordonnees"])
                    elif  (G.nodes[voisin[0]]["coordonnees"][0] - 0.5, G.nodes[voisin[0]]["coordonnees"][1]) not in coordonnees :
                        G.nodes[noeud]["coordonnees"] = (G.nodes[voisin[0]]["coordonnees"][0] - 0.5, G.nodes[voisin[0]]["coordonnees"][1])  
                        coordonnees.append(G.nodes[noeud]["coordonnees"])  
                    else :
                        print("probleme")
                        

    print("pos")            
    #pos = nx.get_node_attributes(G, 'coordonnees')
    pos = dict([(u,(d["coordonnees"]))for u,d in G.nodes(data=True) if u in nodes_list])

    print(pos)
    for elt in pos :
        if pos[elt][0] > 30 or pos[elt][1] > 30 :
            print(elt)
     
    fig = plt.figure(figsize =(5,12))
    courbes = []      
    for noeud in nodes_list :
        courbes_x = []
        courbes_y = []  
        coordonnees_x = []
        coordonnees_y = []
        for voisin in G[noeud] :
            if noeud == (46,40) :
                print(voisin)
                print(G.nodes[voisin]["coordonnees"][0])
                print(G.nodes[noeud]["coordonnees"][0])
            #if (noeud, voisin) in G.edges() and (voisin, noeud) in G.edges() :
#                 print(voisin)
#                 print(G.nodes[voisin]["coordonnees"][1])
#                 print(noeud)
#                 print(G.nodes[noeud]["coordonnees"][1])
            if G.nodes[voisin]["coordonnees"][0] == G.nodes[noeud]["coordonnees"][0] :
                    for elt in coordonnees_x :
                        if G.nodes[voisin]["coordonnees"][0] == elt[1] :
                            courbes_x.append((elt[0], voisin))
            if G.nodes[voisin]["coordonnees"][1] == G.nodes[noeud]["coordonnees"][1] :
                    print("ramou")
                    for elt in coordonnees_y :
                        if G.nodes[voisin]["coordonnees"][1] == elt[1] :
                            courbes_y.append((elt[0], voisin))
            coordonnees_x.append((voisin, G.nodes[voisin]["coordonnees"][0]))
            coordonnees_y.append((voisin, G.nodes[voisin]["coordonnees"][1]))
        #if noeud == (33,28) :
        print(noeud)
        print(courbes_x)
        print(courbes_y)
        if len(courbes_x) > 0 :
            for elt in courbes_x :
                if (G.nodes[noeud]["coordonnees"][1] < G.nodes[elt[0]]["coordonnees"][1] and G.nodes[noeud]["coordonnees"][1] < G.nodes[elt[1]]["coordonnees"][1]) or (G.nodes[noeud]["coordonnees"][1] > G.nodes[elt[0]]["coordonnees"][1] and G.nodes[noeud]["coordonnees"][1] > G.nodes[elt[1]]["coordonnees"][1]) :
                    if abs(G.nodes[elt[0]]["coordonnees"][1] - G.nodes[noeud]["coordonnees"][1]) > abs(G.nodes[elt[1]]["coordonnees"][1] - G.nodes[noeud]["coordonnees"][1]) :
                        if (noeud, elt[0]) not in courbes and (elt[0], noeud) not in courbes  :
                            courbes.append((noeud, elt[0]))
                        
                    else :
                        if (noeud, elt[1]) not in courbes and (elt[1], noeud) not in courbes :
                            courbes.append((noeud, elt[1]))
                        
        if len(courbes_y) > 0 :
            for elt in courbes_y :
                if (G.nodes[noeud]["coordonnees"][0] < G.nodes[elt[0]]["coordonnees"][0] and G.nodes[noeud]["coordonnees"][0] < G.nodes[elt[1]]["coordonnees"][0]) or (G.nodes[noeud]["coordonnees"][0] > G.nodes[elt[0]]["coordonnees"][0] and G.nodes[noeud]["coordonnees"][0] > G.nodes[elt[1]]["coordonnees"][0]) :
                    if abs(G.nodes[elt[0]]["coordonnees"][0] - G.nodes[noeud]["coordonnees"][0]) > abs(G.nodes[elt[1]]["coordonnees"][0] - G.nodes[noeud]["coordonnees"][0]) :
                        if (noeud, elt[0]) not in courbes and (elt[0], noeud) not in courbes  :
                            courbes.append((noeud, elt[0]))
                            
                    else :
                        if (noeud, elt[1]) not in courbes and (elt[1], noeud) not in courbes :
                            courbes.append((noeud, elt[1]))
    
    for u,v,key,data in G.edges(data=True, keys=True) :
        if u in nodes_list and v in nodes_list :
            if G.nodes[u]["type"] != None and G.nodes[v]["type"] != None and data["label"] != 'B53' :
                if (u,v) not in courbes and (v,u) not in courbes :
                    courbes.append((u,v))
                            
    print(courbes)
    green_edges = []
    blue_edges = []
    black_edges = []
    grey_edges = []
    #black_edges = [edge for edge in G.edges() if edge not in red_edges]
    for u,v,edata, in G.edges(data=True) :
            if edata["label"] == "B53" :
                green_edges.append((u,v))
            elif edata["label"] == "CWW" :
                blue_edges.append((u,v)) 
            elif edata["label"] == '0' :
                grey_edges.append((u,v))
            else :
                black_edges.append((u,v))

    #edge_labels=dict([((u,v,),d["label"])for u,v,d in G.edges(data=True) if d["label"] != 'B53' and d["label"] != 'CWW' and ((u,v) not in courbes and (v,u) not in courbes)])
    #print(edge_labels)
    node_labels=dict([(u,(d["type"], d["poids"]))for u,d in G.nodes(data=True) if u in nodes_list and d["type"] != -1 and d["type"] != None])## if d["type"] != None])
    #node_labels=dict([(u, (u,d["type"], d["poids"])) if d["type"] != None else (u, (u)) for u,d in G.nodes(data=True) ])
    #node_labels=dict([(u, (u)) for u,d in G.nodes(data=True) ])
    #print(node_labels)
    nx.draw_networkx_nodes(G, pos, node_size=150, nodelist=nodes_list, node_color="orange")
                #nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), 
                #           node_color = values, node_size = 500)
    nx.draw_networkx_labels(G, pos, labels = node_labels, font_size = 8)
    #nx.draw_networkx_edges(G, pos, edgelist=red_edges, edge_color='r', edge_labels = edge_labels)
    #nx.draw_networkx_edges(G, pos, edgelist=blue_edges, edge_color='b', edge_labels = edge_labels)
    #nx.draw_networkx_edges(G, pos, edgelist=green_edges, edge_color='g', edge_labels = edge_labels)
    #nx.draw_networkx_edges(G, pos, edgelist=black_edges, edge_labels = edge_labels)
    
    edge_colors = ['black' if (u,v) in black_edges else 'blue' if (u,v) in blue_edges else 'green' if (u,v) in green_edges else 'grey' for u,v,d in G.edges(data=True) if (u,v) not in courbes and (v,u) not in courbes and d["label"] != '0' and u in nodes_list and v in nodes_list]
#         print(black_edges)
#         print(red_edges)
#         print(edge_colors)
    #nx.draw_networkx_edge_labels(G,pos)
    edges_list = [(u,v,)for u,v,d in G.edges(data=True) if ((u,v) not in courbes and (v,u) not in courbes) and d["label"] != '0' and u in nodes_list and v in nodes_list] 
    edgelabels = dict([((u,v,), (d["label"])) for u,v,d in G.edges(data=True) if (u,v) in edges_list])
    #print(edges_list)
    #nx.draw_networkx_edge_labels(G,pos, edge_labels=edgelabels, font_size=8)
    nx.draw_networkx_edges(G,pos, edgelist = edges_list, edge_color=edge_colors)
    
    ax=plt.gca()
    draw_network(G, courbes, pos,ax)
    print("petit rat")
    plt.axis('off')
    plt.savefig(EXTENSION_PATH%rep_sauve+"graphe_commun_toutes_aretes_max_0.7_%s_taille_%d.png"%(num, taille_comp)) # save as png
    #plt.savefig("graphes_extension/fichier_1FJG_A_48_8.png") # save as png
    print(courbes)
    print(nodes_list)
    for noeud, data, in G.nodes(data=True) :
        if noeud in nodes_list :
            print(noeud)
            print(data)
    #plt.show()
    
    with open(EXTENSION_PATH%rep_sauve+"graphe_commun_toutes_aretes_max_0.7_%s_taille_%d_avec_coord.pickle"%(num, taille_comp), 'wb') as fichier_ecriture :
        mon_pickler = pickle.Pickler(fichier_ecriture)
        mon_pickler.dump(G)
    plt.clf()
    plt.close()
    ''' affiche le sous-graphe commun G et le stocke en format png dans rep_sauve
(version CaRNAval) '''    

''' affiche le sous-graphe commun G et le stocke en format png dans rep_sauve
(version toutes donnees PDB) 
num : numero du cluster
taille_comp : taille d'extension
commut : depend de G'''
def draw_new_data(G, rep_sauve, num, taille_comp, commut, vrai_graphe, **kwargs):
    motif = kwargs.get("motif", [])
    num_graphe = kwargs.get("num_graphe", [])
    print(G.nodes.data())
    print(G.edges.data())
    
    nx.set_node_attributes(G, (33,33), "coordonnees")
    if not vrai_graphe :
        if commut :
            G.nodes[1]["coordonnees"] = (0.0,0.5)
            G.nodes[2]["coordonnees"] = (2.0,0.5)
            G.nodes[3]["coordonnees"] = (0.0,0.0)
            G.nodes[4]["coordonnees"] = (2.0,0.0)
            G.nodes[5]["coordonnees"] = (3.0,0.5)
        else :
            G.nodes[(1,1)]["coordonnees"] = (0.0,0.5)
            G.nodes[(2,2)]["coordonnees"] = (2.0,0.5)
            G.nodes[(3,3)]["coordonnees"] = (0.0,0.0)
            G.nodes[(4,4)]["coordonnees"] = (2.0,0.0)
            G.nodes[(5,5)]["coordonnees"] = (3.0,0.5)
    else :
        G.nodes[motif[0]]["coordonnees"] = (0.0,0.5)
        G.nodes[motif[1]]["coordonnees"] = (2.0,0.5)
        G.nodes[motif[2]]["coordonnees"] = (0.0,0.0)
        G.nodes[motif[3]]["coordonnees"] = (2.0,0.0)
        G.nodes[motif[4]]["coordonnees"] = (3.0,0.5)
        
     
#                 fichier.write(str(element)+"\n") 
#                 fichier.write(str(G.number_of_nodes())+"\n") 
    print(G.edges.data())
    if not vrai_graphe :
        nodes_list = [1,2,3,4,5]
    else :
        nodes_list = [motif[0], motif[1], motif[2], motif[3], motif[4]]
    nodes_list.extend([u for u,d in G.nodes(data=True) if u not in [1,2,3,4,5]])#and len(G[u]) > 0] 
    print(nodes_list)

    coordonnees = []
    for noeud in nodes_list :
        #voisins = G[noeud]
            print(noeud)
        #if noeud not in [1,2,3,4,5] :
#             if vrai_graphe :
#                 for elt in G.nodes[noeud]["num_seq"] :
#                     if elt[0] == num_graphe :
#                         for noeud, data in dico_graphe[0].nodes(data=True) :
#                             if data["fr3d"] == 
        
            
            if G.nodes[noeud]["type"] != None and G.nodes[noeud]["type"] != -1  :
                if G.nodes[noeud]["coordonnees"] == (33,33) :
                    if not vrai_graphe :
                        chaine = G.nodes[noeud]["chaine"][0]
                    else :
                        chaine = motif[G.nodes[noeud]["chaine"][0]-1]
                    #print(chaine)
                    if commut : 
                        if G.nodes[noeud]["chaine"][0] == 1 or G.nodes[noeud]["chaine"][0] == 3 :
                            print(noeud)
                            print(chaine)
                            G.nodes[noeud]["coordonnees"] =  (G.nodes[chaine]["coordonnees"][0], G.nodes[chaine]["coordonnees"][1] + (G.nodes[noeud]["position"][0] - G.nodes[chaine]["position"][0])/2) 
                        else :
                            G.nodes[noeud]["coordonnees"] =  (G.nodes[chaine]["coordonnees"][0], G.nodes[chaine]["coordonnees"][1] + (G.nodes[chaine]["position"][0] - G.nodes[noeud]["position"][0])/2) 
                    else :
                        if G.nodes[noeud]["chaine"][0] == 1 or G.nodes[noeud]["chaine"][0] == 3 :
                            G.nodes[noeud]["coordonnees"] =  (G.nodes[(chaine,chaine)]["coordonnees"][0], G.nodes[(chaine,chaine)]["coordonnees"][1] + (G.nodes[noeud]["position"][0] - G.nodes[(chaine,chaine)]["position"][0])/2) 
                        else :
                            G.nodes[noeud]["coordonnees"] =  (G.nodes[(chaine,chaine)]["coordonnees"][0], G.nodes[(chaine,chaine)]["coordonnees"][1] + (G.nodes[(chaine,chaine)]["position"][0] - G.nodes[noeud]["position"][0])/2) 
        
                    coordonnees.append(G.nodes[noeud]["coordonnees"])
                print(noeud)
                print("voisins")
                for voisin in G[noeud] :
                        print(voisin)
                        if G.nodes[voisin]["coordonnees"] == (33,33) :
                            if G.nodes[voisin]["type"] == None or G.nodes[voisin]["type"] == -1 :
                                if (G.nodes[noeud]["coordonnees"][0] + 0.5, G.nodes[noeud]["coordonnees"][1]) not in coordonnees :
                                    G.nodes[voisin]["coordonnees"] = (G.nodes[noeud]["coordonnees"][0] + 0.5, G.nodes[noeud]["coordonnees"][1])
                                    coordonnees.append(G.nodes[voisin]["coordonnees"])
                                elif  (G.nodes[noeud]["coordonnees"][0] - 0.5, G.nodes[noeud]["coordonnees"][1]) not in coordonnees :
                                    G.nodes[voisin]["coordonnees"] = (G.nodes[noeud]["coordonnees"][0] - 0.5, G.nodes[noeud]["coordonnees"][1])  
                                    coordonnees.append(G.nodes[voisin]["coordonnees"])  
                                elif (G.nodes[noeud]["coordonnees"][0] + 0.25, G.nodes[noeud]["coordonnees"][1]+0.25) not in coordonnees :
                                    G.nodes[voisin]["coordonnees"] = (G.nodes[noeud]["coordonnees"][0] + 0.25, G.nodes[noeud]["coordonnees"][1]+0.25)  
                                    coordonnees.append(G.nodes[voisin]["coordonnees"])  
                                elif (G.nodes[noeud]["coordonnees"][0] - 0.25, G.nodes[noeud]["coordonnees"][1]-0.25) not in coordonnees :
                                    G.nodes[voisin]["coordonnees"] = (G.nodes[noeud]["coordonnees"][0] - 0.25, G.nodes[noeud]["coordonnees"][1]-0.25)  
                                    coordonnees.append(G.nodes[voisin]["coordonnees"])  
                                else :
                                    print("probleme")
                    

#             else :
#                     voisin = list(G[noeud])
#                     print(len(voisin))
#                     print(voisin)
#                     if (G.nodes[voisin[0]]["coordonnees"][0] + 0.5, G.nodes[voisin[0]]["coordonnees"][1]) not in coordonnees :
#                         G.nodes[noeud]["coordonnees"] = (G.nodes[voisin[0]]["coordonnees"][0] + 0.5, G.nodes[voisin[0]]["coordonnees"][1])
#                         coordonnees.append(G.nodes[noeud]["coordonnees"])
#                     elif  (G.nodes[voisin[0]]["coordonnees"][0] - 0.5, G.nodes[voisin[0]]["coordonnees"][1]) not in coordonnees :
#                         G.nodes[noeud]["coordonnees"] = (G.nodes[voisin[0]]["coordonnees"][0] - 0.5, G.nodes[voisin[0]]["coordonnees"][1])  
#                         coordonnees.append(G.nodes[noeud]["coordonnees"])  
#                     elif (G.nodes[voisin[0]]["coordonnees"][0] + 0.25, G.nodes[voisin[0]]["coordonnees"][1]+0.25) not in coordonnees :
#                         G.nodes[noeud]["coordonnees"] = (G.nodes[voisin[0]]["coordonnees"][0] + 0.25, G.nodes[voisin[0]]["coordonnees"][1]+0.25)  
#                         coordonnees.append(G.nodes[noeud]["coordonnees"])  
#                     elif (G.nodes[voisin[0]]["coordonnees"][0] - 0.25, G.nodes[voisin[0]]["coordonnees"][1]-0.25) not in coordonnees :
#                         G.nodes[noeud]["coordonnees"] = (G.nodes[voisin[0]]["coordonnees"][0] - 0.25, G.nodes[voisin[0]]["coordonnees"][1]-0.25)  
#                         coordonnees.append(G.nodes[noeud]["coordonnees"])  
#                     else :
#                         print("probleme")
                        

    print("pos")            
    #pos = nx.get_node_attributes(G, 'coordonnees')
    pos = dict([(u,(d["coordonnees"]))for u,d in G.nodes(data=True) if u in nodes_list])
    
    print(pos)
    for elt in pos :
        if pos[elt][0] > 30 or pos[elt][1] > 30 :
            print(elt)
     
    fig = plt.figure(figsize =(5,12))
    courbes = []      
    for noeud in nodes_list :
        courbes_x = []
        courbes_y = []  
        coordonnees_x = []
        coordonnees_y = []
        for voisin in G[noeud] :
            if noeud == (46,40) :
                print(voisin)
                print(G.nodes[voisin]["coordonnees"][0])
                print(G.nodes[noeud]["coordonnees"][0])
            #if (noeud, voisin) in G.edges() and (voisin, noeud) in G.edges() :
#                 print(voisin)
#                 print(G.nodes[voisin]["coordonnees"][1])
#                 print(noeud)
#                 print(G.nodes[noeud]["coordonnees"][1])
            if G.nodes[voisin]["coordonnees"][0] == G.nodes[noeud]["coordonnees"][0] :
                    for elt in coordonnees_x :
                        if G.nodes[voisin]["coordonnees"][0] == elt[1] :
                            courbes_x.append((elt[0], voisin))
            if G.nodes[voisin]["coordonnees"][1] == G.nodes[noeud]["coordonnees"][1] :
                    print("ramou")
                    for elt in coordonnees_y :
                        if G.nodes[voisin]["coordonnees"][1] == elt[1] :
                            courbes_y.append((elt[0], voisin))
            coordonnees_x.append((voisin, G.nodes[voisin]["coordonnees"][0]))
            coordonnees_y.append((voisin, G.nodes[voisin]["coordonnees"][1]))
        #if noeud == (33,28) :
        print(noeud)
        print(courbes_x)
        print(courbes_y)
        if len(courbes_x) > 0 :
            for elt in courbes_x :
                if (G.nodes[noeud]["coordonnees"][1] < G.nodes[elt[0]]["coordonnees"][1] and G.nodes[noeud]["coordonnees"][1] < G.nodes[elt[1]]["coordonnees"][1]) or (G.nodes[noeud]["coordonnees"][1] > G.nodes[elt[0]]["coordonnees"][1] and G.nodes[noeud]["coordonnees"][1] > G.nodes[elt[1]]["coordonnees"][1]) :
                    if abs(G.nodes[elt[0]]["coordonnees"][1] - G.nodes[noeud]["coordonnees"][1]) > abs(G.nodes[elt[1]]["coordonnees"][1] - G.nodes[noeud]["coordonnees"][1]) :
                        if (noeud, elt[0]) not in courbes and (elt[0], noeud) not in courbes  :
                            courbes.append((noeud, elt[0]))
                        
                    else :
                        if (noeud, elt[1]) not in courbes and (elt[1], noeud) not in courbes :
                            courbes.append((noeud, elt[1]))
                        
        if len(courbes_y) > 0 :
            for elt in courbes_y :
                if (G.nodes[noeud]["coordonnees"][0] < G.nodes[elt[0]]["coordonnees"][0] and G.nodes[noeud]["coordonnees"][0] < G.nodes[elt[1]]["coordonnees"][0]) or (G.nodes[noeud]["coordonnees"][0] > G.nodes[elt[0]]["coordonnees"][0] and G.nodes[noeud]["coordonnees"][0] > G.nodes[elt[1]]["coordonnees"][0]) :
                    if abs(G.nodes[elt[0]]["coordonnees"][0] - G.nodes[noeud]["coordonnees"][0]) > abs(G.nodes[elt[1]]["coordonnees"][0] - G.nodes[noeud]["coordonnees"][0]) :
                        if (noeud, elt[0]) not in courbes and (elt[0], noeud) not in courbes  :
                            courbes.append((noeud, elt[0]))
                            
                    else :
                        if (noeud, elt[1]) not in courbes and (elt[1], noeud) not in courbes :
                            courbes.append((noeud, elt[1]))
    
    for u,v,key,data in G.edges(data=True, keys=True) :
        if u in nodes_list and v in nodes_list :
            if G.nodes[u]["type"] != None and G.nodes[v]["type"] != None and data["label"] != 'B53' :
                if (u,v) not in courbes and (v,u) not in courbes :
                    courbes.append((u,v))
                            
    print(courbes)
    green_edges = []
    blue_edges = []
    black_edges = []
    grey_edges = []
    #black_edges = [edge for edge in G.edges() if edge not in red_edges]
    for u,v,edata, in G.edges(data=True) :
            if edata["label"] == "B53" :
                green_edges.append((u,v))
            elif edata["label"] == "CWW" :
                blue_edges.append((u,v)) 
            elif edata["label"] == '0' :
                grey_edges.append((u,v))
            else :
                black_edges.append((u,v))

    #edge_labels=dict([((u,v,),d["label"])for u,v,d in G.edges(data=True) if d["label"] != 'B53' and d["label"] != 'CWW' and ((u,v) not in courbes and (v,u) not in courbes)])
    #print(edge_labels)
    node_labels=dict([(u,(u))for u,d in G.nodes(data=True) ])#if u in nodes_list and d["type"] != -1 and d["type"] != None])## if d["type"] != None])
    #node_labels=dict([(u, (u,d["type"], d["poids"])) if d["type"] != None else (u, (u)) for u,d in G.nodes(data=True) ])
    #node_labels=dict([(u, (u)) for u,d in G.nodes(data=True) ])
    print(node_labels)
    nx.draw_networkx_nodes(G, pos, node_size=150, nodelist=nodes_list, node_color="orange")
                #nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), 
                #           node_color = values, node_size = 500)
    nx.draw_networkx_labels(G, pos, labels = node_labels, font_size = 12)
    #nx.draw_networkx_edges(G, pos, edgelist=red_edges, edge_color='r', edge_labels = edge_labels)
    #nx.draw_networkx_edges(G, pos, edgelist=blue_edges, edge_color='b', edge_labels = edge_labels)
    #nx.draw_networkx_edges(G, pos, edgelist=green_edges, edge_color='g', edge_labels = edge_labels)
    #nx.draw_networkx_edges(G, pos, edgelist=black_edges, edge_labels = edge_labels)
    
    edge_colors = ['black' if (u,v) in black_edges else 'blue' if (u,v) in blue_edges else 'green' if (u,v) in green_edges else 'blue' for u,v,d in G.edges(data=True) if (u,v) not in courbes and (v,u) not in courbes and d["label"] != '0' and u in nodes_list and v in nodes_list]
#         print(black_edges)
#         print(red_edges)
#         print(edge_colors)
    #nx.draw_networkx_edge_labels(G,pos)
    edges_list = [(u,v,)for u,v,d in G.edges(data=True) if ((u,v) not in courbes and (v,u) not in courbes) and d["label"] != '0' and u in nodes_list and v in nodes_list] 
    edgelabels = dict([((u,v,), (d["label"])) for u,v,d in G.edges(data=True) if (u,v) in edges_list])
    #print(edges_list)
    #nx.draw_networkx_edge_labels(G,pos, edge_labels=edgelabels, font_size=8)
    nx.draw_networkx_edges(G,pos, edgelist = edges_list, edge_color=edge_colors)
    
    ax=plt.gca()
    draw_network(G, courbes, pos,ax)
    print("petit rat")
    plt.axis('off')
    #plt.savefig(NEW_EXTENSION_PATH_TAILLE+rep_sauve+"/graphe_commun_global_cluster_%s_%s_%s.png"%(num, taille_comp, num_graphe)) # save as png
    #plt.savefig("graphes_extension/fichier_1FJG_A_48_8.png") # save as png
    print(courbes)
    print(nodes_list)
    for noeud, data, in G.nodes(data=True) :
        if noeud in nodes_list :
            print(noeud)
            print(data)
    #plt.show()
    
#     with open(NEW_EXTENSION_PATH_TAILLE+rep_sauve+"/graphe_commun_global_cluster_%s_%s_avec_coord.pickle"%(num, taille_comp), 'wb') as fichier_ecriture :
#         mon_pickler = pickle.Pickler(fichier_ecriture)
#         mon_pickler.dump(G)
    plt.show()
    plt.clf()
    plt.close()
    

''' calcule le nombre d'aretes de tous les sous-graphes communs à un cluster des differentes tailles d'extension et renvoie le maximum avec la taille associee
(version CaRNAval mais peut etre adaptee pour les nouvelles donnees) '''
def plus_grand_nombre_aretes(rep, num):
    maximum = -1
    taille_maximum = -1
    for i in range(4,11) :
        with open(EXTENSION_PATH%rep+"digraphe_commun_0.7_%s_taille_%s.pickle"%(num, i), 'rb') as fichier : 
            mon_depickler = pickle.Unpickler(fichier)
            digraphe_commun = mon_depickler.load()  
            
            print(digraphe_commun.nodes.data())
            print(digraphe_commun.edges.data())
            compter_aretes = 0
            for u,v,key,data in digraphe_commun.edges(data=True, keys=True) :
                if data["label"] != 'B53' :
                    compter_aretes += digraphe_commun.nodes[u]["poids"]
            compter_aretes = compter_aretes/2 - 4
            
            
            if compter_aretes > maximum :
                maximum = compter_aretes
                taille_maximum = i
    
    print(maximum)
    print(taille_maximum)


''' renvoie l'espacement maximal dans un sous-graphe commun à cluster pour une certaine taille d'extension 
cad la difference de distance au motif la plus grande entre les noeuds associes des differents graphes du cluster
(version CaRNAval mais peut etre adaptee pour les nouvelles donnees) '''
def espacement_max(rep, taille_ext, num):
    with open(EXTENSION_PATH%rep+"digraphe_commun_0.7_55_%s_taille_%s.pickle"%(num, taille_ext), 'rb') as fichier : 
        mon_depickler = pickle.Unpickler(fichier)
        digraphe_commun = mon_depickler.load()  
        
        diff_max = -1
        for noeud, data in digraphe_commun.nodes(data=True) :
            print(data)
            for i in range(len(data["espacement_motif"])) :
                for j in range(i+1, len(data["espacement_motif"])) :
                    if abs(data["espacement_motif"][i] - data["espacement_motif"][j]) > diff_max :
                        diff_max = abs(data["espacement_motif"][i] - data["espacement_motif"][j])
        print(taille_ext)
        print(diff_max)
            

if __name__ == '__main__':
    with open("/media/coline/Maxtor/clustering_perez_tot_new_data_100320_sim_par_branche_0.65.pickle", 'rb') as fichier_sortie :
        mon_depickler = pickle.Unpickler(fichier_sortie)
        clusters = mon_depickler.load()
        
        digraphe_commun, liste_cliques = commun_cluster_clique_new_data(clusters[55], "/media/coline/Maxtor/dico_new_100320_sim_par_branche_0.65.pickle", NEW_EXTENSION_PATH_TAILLE)

        draw_new_data(digraphe_commun, "Graphes_communs_avril_2020", 56, 4, True)
        
        exit()
        

    
    for j in range(len(CLUSTERING_PEREZ_VERSION_NON_CAN_2)) :
        for i in range(4,11) :
            graphe_commun = commun_cluster_clique(CLUSTERING_PEREZ_VERSION_NON_CAN_2[j], EXTENSION_PATH%i+"dico_comp_complet_metrique_toutes_aretes_coeff_all1_taille_%s.pickle"%(i), EXTENSION_PATH_TAILLE%i)
            if graphe_commun != None :
                draw(graphe_commun,"taille_max/result_k_max_4_10_toutes_aretes", "groupe_clustering_perez_%d"%j, i, True)
    #     for i in range(4,11) :
#         digraphe_commun = commun_cluster(GROUPE_ETOILE_GNRA, EXTENSION_PATH_TAILLE%i, EXTENSION_PATH%i+"dico_comp_complet_metrique_toutes_aretes_coeff_all1_taille_%s.pickle"%(i), "etoile_gnra", i, 0.7)                
#         draw(digraphe_commun,"taille_max/result_k_max_4_10_toutes_aretes", "etoile_gnra", i) 
#     for j in range(4,11) :
#         for i in range(len(GROUPES_TOUTES_ARETES_MAX_4_10_07)) :
#             digraphe_commun = commun_cluster_raccourci(GROUPES_TOUTES_ARETES_MAX_4_10_07[i], "Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/graphes_extension_taille_%s"%j, "Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/result_graphes_extension_%s/dico_comp_complet_metrique_toutes_aretes_coeff_all1_taille_%s.pickle"%(j,j), i, j)                
# 
# #             print(digraphe_commun.nodes.data())
#             draw(digraphe_commun,"taille_max/result_k_max_4_10_toutes_aretes", i, j)   

    #plus_grand_nombre_aretes("taille_max/result_k_max_4_10_toutes_aretes/sous_graphe_commun_etoile_gnra_5_elts", "etoile_gnra")
    
    #for i in range(4,11) :
    #espacement_max("taille_max/result_k_max_4_10_toutes_aretes", 10, 5.1)
            
#     with open(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes"+"digraphe_commun_0.7_55_%s_taille_%s.pickle"%(0, 4), 'rb') as fichier : 
#         mon_depickler = pickle.Unpickler(fichier)
#         digraphe_commun = mon_depickler.load() 
#         
#         print(digraphe_commun.nodes.data())
#         
#         with open(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes"+"digraphe_commun_0.7_55_%s_taille_%s_raccourci.pickle"%(0, 4), 'rb') as fichier_rac : 
#             mon_depickler = pickle.Unpickler(fichier_rac)
#             digraphe_commun_rac = mon_depickler.load() 
#             
#             print(digraphe_commun_rac.nodes.data())
            
        
             
    