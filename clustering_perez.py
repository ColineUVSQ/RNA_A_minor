'''
Created on 12 juin 2019

@author: coline

Implementation de l'algorithme de clustering Perez
(utilisé dans les version CaRNAval et données toutes PDB)

'''

import networkx as nx
from networkx import algorithms
from recup_data.constantes import EXTENSION_PATH
import pickle
import matplotlib.pyplot as plt
import csv
import seaborn as sns
import numpy as np

dico_intra_sim = {}

'''Calcul de la similarite moyenne intra du wsgraphe
cad la similarite moyenne au sein du wsgraphe '''
def approx_intra_sim(wsgraphe, centre) :
    
    if wsgraphe.number_of_nodes() < 2 :
        return 0
    else : 
        somme_sim = 0
#         for noeud in wsgraphe.nodes() :
#             if noeud != centre :
#                 somme_sim += wsgraphe.edges[centre,noeud]["poids"]
        for u,v,data in wsgraphe.edges(data=True) :
            somme_sim += data["sim"]
                
        return somme_sim/((wsgraphe.number_of_nodes() * (wsgraphe.number_of_nodes() - 1))/2)

'''stocke toutes les intra_sim de tous les sommets du graphe depart dans une variable 
(pour eviter d avoir a les recalculer a chaque fois) '''
def stockage_intra_sim(graphe_depart):
    for noeud in graphe_depart.nodes() :
        wsgraphe = construction_wsgraphe(graphe_depart, noeud)
        dico_intra_sim.update({noeud : approx_intra_sim(wsgraphe, noeud)})

''' Calcul de la densite relative du sommet cad
la proportion de ses voisins qui ont un degre plus eleve que lui '''
def calcul_relative_density(graphe_depart, sommet):
    
    if len(graphe_depart[sommet]) == 0 :
        return 0
    else : 
        density = 0
        for voisin in graphe_depart[sommet] :
            if graphe_depart.degree(voisin) <= graphe_depart.degree(sommet) :
                density += 1
               
        return density/len(list(graphe_depart[sommet]))

''' Calcul de la compaction relative du sommet
cad la proportion de ses voisins qui ont une similarite moyenne intra (la similarite moyenne intra du wsgraphe de centre ce sommet) superieure a la sienne '''
def calcul_relative_compactness(graphe_depart, sommet):
    if len(graphe_depart[sommet]) == 0 :
        return 0
    else : 
        compactness = 0
        for voisin in graphe_depart[sommet] :
            if dico_intra_sim[sommet] >= dico_intra_sim[voisin] :
                compactness += 1
        return compactness/len(list(graphe_depart[sommet]))

''' Calcul de la pertinence du sommet
cad la moyenne entre la densite relative et la compaction relative de ce sommet '''
def calcul_relevance(graphe_depart, sommet):
    relevance = (calcul_relative_density(graphe_depart, sommet) + calcul_relative_compactness(graphe_depart, sommet))/2
    print(sommet)
    print(calcul_relative_compactness(graphe_depart, sommet))
    print(relevance)
    return relevance

''' Renvoie le wsgraphe associe a un sommet 
cad le sous-graphe du graphe de depart comprenant le sommet et tous ses voisins (et donc toutes les aretes qui les relient entre eux) '''
def construction_wsgraphe(graphe_depart, sommet):
    noeuds = list(graphe_depart[sommet])
    noeuds.append(sommet)
    
    wsgraphe = graphe_depart.subgraph(noeuds)

    return wsgraphe

''' Indique si le sommet appartient deja a un cluster (est centre et/ou adjacent a un centre) '''
def sommet_couvert(graphe_depart, centres, sommet):
    for elt in centres :
        if sommet in graphe_depart[elt] or sommet == elt :
            return True
    return False

''' Calcul de la densite d'un cluster '''
def calcul_densite(cluster, graphe_depart):
    composante = graphe_depart.subgraph(cluster)
    if composante.number_of_nodes() <= 2 :
        return 1
    else :
        return composante.number_of_edges()/((composante.number_of_nodes()*(composante.number_of_nodes()-1))/2)

''' Algo principal
avec comme modif par rapport a l'original 
de ne selectionner que des centres dont la densite du wsgraphe associee est superieure a un seuil'''
def algo_principal(graphe_depart):
    
    stockage_intra_sim(graphe_depart)
    print(dico_intra_sim)
    
    ## Couverture initiale
    nx.set_node_attributes(graphe_depart, "satellite", "type")
    dico_relevance = {}
    
    centres = []
    
    for noeud in graphe_depart.nodes() :
        dico_relevance.update({noeud : calcul_relevance(graphe_depart, noeud)})
    
    dico_relevance_ordonne = sorted(dico_relevance.items(), key=lambda t: t[1], reverse=True)
    print(dico_relevance)
    
    for elt in dico_relevance.keys() :
        #print(graphe_depart.nodes[elt]["nom"])
        print(dico_relevance[elt])
        print(elt)
    
    for elt in dico_relevance_ordonne :
        if elt[1] > 0 :
            centre_ok = False
            if sommet_couvert(graphe_depart, centres, elt[0]) :
                voisins_non_couverts = 0
                for voisin in graphe_depart[elt[0]] :
                    if not sommet_couvert(graphe_depart, centres, voisin) :
                        voisins_non_couverts += 1
                if voisins_non_couverts > 0 :
                    centre_ok = True
                    #centres.append(elt[0])
            else :
                centre_ok = True
                #centres.append(elt[0])   
            ## version verif de la densite au moment de la selection des centres
            if centre_ok :
                cluster = [elt[0]]
                cluster.extend(graphe_depart[elt[0]])
                if calcul_densite(cluster, graphe_depart) > 0.7 :
                    centres.append(elt[0])
            
    print(centres)
    print(len(centres))
    
    clusters = [[]]
    for elt in centres :
        cluster = [elt]
        cluster.extend(graphe_depart[elt]) ## ajout au cluster de tous les sommets adjacents au sommet seed
        clusters.append(cluster)
    print(clusters)
    print(len(clusters))
    
    a_enlever = []
    for cluster in clusters :
        print("cluster")
        ## version supprimer les clusters de densite trop faible
#         if calcul_densite(cluster, graphe_depart) < 0.77 :
#             a_enlever.append(cluster)
    
        #if cluster == [86, 9, 16, 34, 50]:
        for c in cluster :
                print(c)
    
#     for elt in a_enlever :
#         clusters.remove(elt)
#         centres.remove(elt[0])
        
    ## Elimination des centres inutiles
    dico_degres_centres = {}
    dico_analyses_centres = {}
    dico_linked_centres = {}
    dico_seed_centres = {}
    for elt in centres :
        dico_degres_centres.update({elt : graphe_depart.degree(elt)})
        dico_analyses_centres.update({elt : "non_analyse"})
        dico_linked_centres.update({elt : []})
        dico_seed_centres.update({elt : "non_seed"})
       
    dico_degres_ordonnes = sorted(dico_degres_centres.items(), key=lambda t: t[1], reverse=True)
    for elt in dico_degres_ordonnes :
        if elt[0] in centres : ## pour verifier que l element courant est encore dans les centres 
            print("centre")
            print(elt[0])
            print("voisin")
            for voisin in graphe_depart[elt[0]] :
                if voisin in centres and dico_analyses_centres[voisin] == "non_analyse" :
                    print(voisin)
                    nb_adjacents_partages = 0
                    nb_adjacents_non_partages = 0
                    adjacents_non_partages = []
                    for adj in graphe_depart[voisin] :
                        print(adj)
                        autres_centres = list(centres)
                        autres_centres.remove(voisin)
                        
                        if sommet_couvert(graphe_depart, autres_centres, adj) :
                            nb_adjacents_partages += 1
                            print("petit rat")
                        else :
                            nb_adjacents_non_partages += 1
                            adjacents_non_partages.append(adj)
                            print("gros rat")
                    
                    print(nb_adjacents_non_partages)
                    print(nb_adjacents_partages)
                           
                    if nb_adjacents_partages > nb_adjacents_non_partages :
                        print("ramousnif")
                        centres.remove(voisin)
                        dico_analyses_centres.pop(voisin, None) ## pour supprimer une cle
                        dico_linked_centres.pop(voisin, None)
                        
                        dico_linked_centres[elt[0]].extend(adjacents_non_partages)
                    else :
                        dico_analyses_centres[voisin] = "analyse"
                        
            dico_seed_centres[elt[0]] = "seed"
    
    print(dico_seed_centres)   
    
    ## Creation clusters
    
    clusters = [[]]
    for cle in dico_seed_centres.keys() :
        if dico_seed_centres[cle] == "seed" :
            cluster = [cle]
            cluster.extend(graphe_depart[cle]) ## ajout au cluster de tous les sommets adjacents au sommet seed
            cluster.extend(dico_linked_centres[cle]) ## ajout au cluster de tous les sommets linked potentiels
            
            clusters.append(cluster)
    print(clusters)
    print(len(clusters))
    
    tab_clusters = []
    for cluster in clusters :
        print("cluster")
        print(calcul_densite(cluster, graphe_depart))
        tab_clusters.append([])
        for c in cluster :
            tab_clusters[len(tab_clusters)-1].append(c)
#             if graphe_depart.nodes[c]["nom"] == "1U9S_A_58_11" :
#                 print(c)
            if c == 34 :
                print("petit ramou")
            #print(graphe_depart.nodes[c]["nom"])
            
    print(tab_clusters)   
    print(len(tab_clusters))
    
    for noeud, data in graphe_depart.nodes(data=True) :
        est_dans_un_cluster = False
        for cluster in tab_clusters :
            if noeud in cluster :
                est_dans_un_cluster = True  
        #print(est_dans_un_cluster)
        if not est_dans_un_cluster :
            tab_clusters.append([noeud])
    
    
    print(tab_clusters)   
    print(len(tab_clusters))
    return tab_clusters, dico_relevance


def noeuds_a_garder_centrality(num_ARN):
    with open("/media/coline/Maxtor/graphe_dico_sim_sans_doublons_avec_manque_%s_new_noms_res_3a.pickle"%num_ARN, 'rb') as fichier_graphe :
        mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
        graphe_complet = mon_depickler_graphe.load() 
        liste = [(u,v) for u,v in graphe_complet.edges() if u == v]
        #print(liste)
        graphe_complet.remove_edges_from([(u,v) for u,v in graphe_complet.edges() if u == v])
        print(graphe_complet.edges())
        print(len(graphe_complet.edges()))
        print(graphe_complet.number_of_edges())
        #print(graphe_complet.edges.data())
        a_enlever = []
        for u,v,data in graphe_complet.edges(data=True) : 
            if data["sim"] < 0.6 : 
                a_enlever.append((u,v))
          
        for elt in a_enlever :
            graphe_complet.remove_edge(elt[0], elt[1])
             
        print(graphe_complet.number_of_edges())
        
        clusters = algo_principal(graphe_complet)
        print(len(clusters))
        print(clusters)
        
        with open("/media/coline/Maxtor/clustering_perez_%s_new_data_res_inf_3.pickle"%num_ARN, 'wb') as fichier_sortie :
            mon_pickler = pickle.Pickler(fichier_sortie)
            mon_pickler.dump(clusters)
        
        nx.set_node_attributes(graphe_complet, [] , "clustering_perez")
        compteur = 1
        for cluster in clusters :
            for elt in cluster :
                #print(elt)
#                 for noeud, data in graphe_complet.nodes(data=True) :
#                     if data["nom"] == elt :
#                         print(data["nom"])
                        
                    liste = list(graphe_complet.nodes[elt]["clustering_perez"])
                    liste.append(compteur)
                    graphe_complet.nodes[elt]["clustering_perez"] = list(liste)
                        
                        #exit(0)
                    
            compteur += 1
        print(graphe_complet.nodes.data())
        #node_labels=dict([(u, (d["clustering_perez"])) for u,d in graphe_complet.nodes(data=True)]) #else (u, (u)) for u,d in G.nodes(data=True) ])
        
        liste_noeuds_a_garder = []
        
        nx.set_node_attributes(graphe_complet, -1 , "centrality_cluster")
        for cluster in clusters :
            
            ok_2_elts =True
            if len(cluster) == 2 :
                for elt in cluster :
                    #for noeud, data in graphe_complet.nodes(data=True) :
                        #if data["nom"] == elt :
                            if len(graphe_complet[elt]) > 1 :
                                ok_2_elts = False
                                
            if len(cluster) > 2 or ok_2_elts :
                liste_elt = []
                for elt in cluster :
                    #for noeud, data in graphe_complet.nodes(data=True) :
                        #if data["nom"] == elt :
                            if len(graphe_complet.nodes[elt]["clustering_perez"]) == 1 :
                                liste_elt.append(elt)
                            else :
                                liste_noeuds_a_garder.append(elt)
                if len(liste_elt) > 0 :
                    subgraph = nx.subgraph(graphe_complet, liste_elt)
                    centrality = nx.eigenvector_centrality(subgraph)
                    print(centrality)
                    
                    noeud_max = -1
                    maxi = 0
                    
                    for noeud in centrality :
                        if maxi < centrality[noeud] :
                            maxi = centrality[noeud]
                            noeud_max = noeud
                        graphe_complet.nodes[noeud]["centrality_cluster"] = round(centrality[noeud], 2)
                        
                    liste_noeuds_a_garder.append(noeud_max)
            else :
                for elt in cluster :
                    #for noeud, data in graphe_complet.nodes(data=True) :
                        #if data["nom"] == elt :
                            liste_noeuds_a_garder.append(elt)
        
        print(len(liste_noeuds_a_garder))
                       
        node_labels=dict([(u, (d["centrality_cluster"])) for u,d in graphe_complet.nodes(data=True)]) #else (u, (u)) for u,d in G.nodes(data=True) ])
        
  
        nx.draw(graphe_complet, node_size = 3, labels = node_labels, font_size=20)
        plt.plot()
        plt.show()
        
        
        print(graphe_complet.nodes.data())
        with open("/media/coline/Maxtor/fichier_csv_new_data_seuil_0.6_avec_clustering_centrality_%s_avec_res_3.csv"%num_ARN,'w') as fichier_csv:
                    csvwriter = csv.writer(fichier_csv)
                    csvwriter.writerow(["source", "target", "label"])
                     
                    #seuil = 0.6
                    for u,v,data in graphe_complet.edges(data=True) :
                        #if data["sim"] > seuil :
                        csvwriter.writerow([u,v,round(data["sim"],2)])
                         
                    with open("/media/coline/Maxtor/fichier_csv_new_data_seuil_0.6_avec_clustering_centrality_noeud_%s_avec_res_3.csv"%num_ARN,'w') as fichier_csv:
                            csvwriter2 = csv.writer(fichier_csv)
                            csvwriter2.writerow(["id", "label", "nom_extension", "clusters", "garde"])
                                 
                            for noeud,data in graphe_complet.nodes(data=True) :
                                num_clusters = []
                                compteur = 1
                                for cluster in clusters:
                                    if noeud in cluster :
                                        num_clusters.append(compteur)
                                    compteur += 1
                                
                                a_garder = 0
                                if noeud in liste_noeuds_a_garder :
                                    a_garder = 1
                                
                                csvwriter2.writerow([noeud, graphe_complet.nodes[noeud]["centrality_cluster"], noeud, num_clusters, a_garder])
#         liste_a_garder_noms = []
#         for elt in liste_noeuds_a_garder :
#             liste_a_garder_noms.append(graphe_complet.nodes[elt]["nom"])
#         print(liste_a_garder_noms)
        
        with open("/media/coline/Maxtor/noeuds_a_garder_selon_centrality_%s_avec_res_3.pickle"%num_ARN, 'wb') as fichier_sortie :
            mon_pickler = pickle.Pickler(fichier_sortie)
            mon_pickler.dump(liste_noeuds_a_garder)
            
            
def clustering_pour_complet(liste_num_ARN, seuil):
    liste_tout = []
    for elt in liste_num_ARN :
        with open("resolutions.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            resolutions = mon_depickler.load()
             
            with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
                    mon_depickler = pickle.Unpickler(fichier_sortie)
                    liste_a_garder_noms = mon_depickler.load()    
                       
                    for element in liste_a_garder_noms :
                        if resolutions[element[0]] <= 3 :
                            liste_tout.append((elt, element))
    graphe_complet = nx.Graph()
    for elt in liste_tout :
        graphe_complet.add_node(elt[1], type=elt[0], nom=elt[1])
    
    with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_612.pickle", 'rb') as fichier_graphe :
        mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
        dico_complet = mon_depickler_graphe.load() 
        #print(liste)
        
        for cle in dico_complet.keys() :
                if isinstance(cle[0], tuple) :
                    if isinstance(dico_complet[cle], dict) :
                        graphe_complet.add_edge(cle[0], cle[1], sim=dico_complet[cle]["sim"])
                    else :
                        graphe_complet.add_edge(cle[0], cle[1], sim=dico_complet[cle])
                else :
                    cle0 = (cle[0].split("_")[0], int(cle[0].split("_")[1]))
                    cle1 = (cle[1].split("_")[0], int(cle[1].split("_")[1]))
#                     if (cle0, cle1) in dico_manque.keys() :
#                         print(dico_manque[(cle0, cle1)])
#                         print(dico_manque[cle])
#                     elif (cle1, cle0) in dico_manque.keys() :
#                         print(dico_manque[(cle1, cle0)])
#                         print(dico_manque[cle])
                    
                    if isinstance(dico_complet[cle], dict) :
                        graphe_complet.add_edge(cle0, cle1, sim=dico_complet[cle]["sim"])
                    else :
                        graphe_complet.add_edge(cle0, cle1, sim=dico_complet[cle])
        
        graphe_complet.remove_edges_from([(u,v) for u,v in graphe_complet.edges() if u == v])
        #print(graphe_complet.edges())
        print(len(graphe_complet.edges()))
        #print(graphe_complet.number_of_edges())
        #print(graphe_complet.edges.data())
        
        a_enlever = []
        for u,v,data in graphe_complet.edges(data=True) : 
            print(data["sim"])
            if isinstance(data["sim"], dict) :
                if data["sim"]["sim"] < seuil : 
                    a_enlever.append((u,v))
                else :
                    data["sim"] = data["sim"]["sim"]
            else :
                if data["sim"] < seuil : 
                    a_enlever.append((u,v)) 
          
        for elt in a_enlever :
            graphe_complet.remove_edge(elt[0], elt[1])
             
        print(graphe_complet.number_of_edges())
        
        
        
        clusters, dico_relevance = algo_principal(graphe_complet)
        print(len(clusters))
        
        tailles_clusters = []
        somme_tailles = 0
        for cluster in clusters :
            tailles_clusters.append(len(cluster))
            print(len(cluster))
            somme_tailles += len(cluster)
        print(somme_tailles)
        
        axs = sns.distplot(tailles_clusters, kde=False, color="blue", bins=35)
        axs.set_xticks([x for x in np.arange(0,max(tailles_clusters)+1, 2)])
        axs.set_xlabel("Cluster size")
        axs.set_ylabel("Number of clusters")
        #plt.savefig("/home/coline/Documents/Extensions/poster_macim/distrib_taille_cluster.svg", format='svg', transparent=True)
        plt.show()
        #print(algorithms.smallworld.sigma(graphe_complet, niter=1, nrand=1))
        print("petit rat")
       
        #print(clusters)
        
        with open("/media/coline/Maxtor/clustering_perez_%s_new_data_res_inf_3_avec_modif_612.pickle"%liste_num_ARN, 'wb') as fichier_sortie :
            mon_pickler = pickle.Pickler(fichier_sortie)
            mon_pickler.dump(clusters)
         
        with open("/media/coline/Maxtor/clustering_perez_%s_new_data_res_inf_3_avec_modif_relevance_612.pickle"%liste_num_ARN, 'wb') as fichier_relevance :
            mon_pickler = pickle.Pickler(fichier_relevance)
            mon_pickler.dump(dico_relevance)
        
        return len(clusters)
        
if __name__ == '__main__':
    types = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16S_arnm", "arnt_16S"]
#     distrib_nb_clusters = []
#     for seuil in np.arange(0.0, 1.1, 0.1) :
    clustering_pour_complet(types, 0.6)
#         print(seuil)
#     print(distrib_nb_clusters)
#     sns.distplot(distrib_nb_clusters)
    #plt.show()
    
#     for typ in types :
#         with open("/media/coline/Maxtor/clustering_perez_%s_new_data.pickle"%typ, 'rb') as fichier_sortie :
#                 mon_depickler = pickle.Unpickler(fichier_sortie)
#                 clusters = mon_depickler.load()
#                 print(typ)
#                 print(len(clusters))
        #noeuds_a_garder_centrality(typ)
#         graphe_test = nx.Graph()
#         graphe_test.add_edge(1,2, poids = 0.15)
#         graphe_test.add_edge(2,3, poids = 0.20)
#         graphe_test.add_edge(2,4, poids = 0.30)
#         graphe_test.add_edge(3,4, poids = 0.20)
#         graphe_test.add_edge(4,5, poids = 0.40)
#         graphe_test.add_edge(4,6, poids = 0.30)
#         graphe_test.add_edge(4,7, poids = 0.10)
#         graphe_test.add_edge(6,7, poids = 0.10)
#         graphe_test.add_edge(7,8, poids = 0.10)
#         graphe_test.add_edge(7,9, poids = 0.10)
#         graphe_test.add_edge(7,10, poids = 0.10)
#         graphe_test.add_edge(8,11, poids = 0.70)
#         graphe_test.add_edge(8,12, poids = 0.50)
#         graphe_test.add_edge(8,13, poids = 0.30)
#         graphe_test.add_edge(8,9, poids = 0.20)
#         graphe_test.add_edge(8,14, poids = 0.16)
#         graphe_test.add_edge(14,15, poids = 0.32)
#         graphe_test.add_edge(14,16, poids = 0.40)
#         graphe_test.add_edge(14,17, poids = 0.20)
#         graphe_test.add_edge(14,18, poids = 0.25)
#         graphe_test.add_edge(9,10, poids = 0.20)
#         graphe_test.add_edge(10,19, poids = 0.20)
#         graphe_test.add_edge(10,20, poids = 0.20)
#         graphe_test.add_edge(10,21, poids = 0.20)
#         graphe_test.add_edge(10,22, poids = 0.50)
#         graphe_test.add_edge(21,22, poids = 0.20)
#         graphe_test.add_edge(22,23, poids = 0.30)
#         
#         algo_principal(graphe_test)
        
                            
                
     


        
        
         
    