'''
Created on 13 mai 2019

@author: coline

Implementation du clustering MCODE
(utilisé pour la version CaRNAval mai 2019)

'''
import networkx as nx
import pickle
import csv
import numpy as np
from recup_data.constantes import EXTENSION_PATH,\
    GROUPES_TOUTES_ARETES_MAX_4_10_07

''' renvoie le plus grand k-core du graphe
cad le sous-graphe central le plus densement connecte
et la valeur de k (dans un k-core, tous les noeuds sont de degre au moins k)'''
def trouve_highest_kcore(graphe):
#     print(nx.is_frozen(graphe))
    k = 1
    degree_sequence = sorted([d for n, d in graphe.degree()], reverse=True)
    dmax = max(degree_sequence)
    
#     print(graphe.nodes.data())
#     print(dmax)
    while graphe.number_of_nodes() > 0 and k <= dmax :
#         print(k)
        graphe_k_avant = graphe.copy()
        nombre_noeuds = graphe.number_of_nodes()
        nombre_noeuds_new = -1
        while nombre_noeuds != nombre_noeuds_new :
            nombre_noeuds = graphe.number_of_nodes()
            for noeud in graphe.nodes() :
#                 print("noeud")
#                 print(noeud)
                if graphe.degree(noeud) < k :
                    #print("ramousnif")
#                     print(graphe.nodes.data())
                    graphe.remove_node(noeud)
                    break
            nombre_noeuds_new = graphe.number_of_nodes()
        k = k+1
    
    if graphe.number_of_nodes() == 0 :
#         print(k-2)
#         print(graphe_k_avant.nodes.data())
#         print(graphe_k_avant.edges.data())
        return k-2, graphe_k_avant
    else :
#         print('ramousnif')
#         print(k-1)
#         print(graphe.nodes.data())
#         print(graphe.edges.data()) 
        return k-1, graphe  
        
''' renvoie la densite du graphe passe en entree'''
def densite(graphe):
#     print(graphe.number_of_edges())
    return 2*graphe.number_of_edges()/(graphe.number_of_nodes()*(graphe.number_of_nodes()-1))     

''' renvoie un dictionnaire contenant les noeuds du graphe ordonnes par densite '''
def ordonner_par_densite(graphe):
    dico_noeuds = {}
    for noeud, data in graphe.nodes(data=True) :
        dico_noeuds[noeud] = data["density"]
    
    dico_noeuds_ordonnes = sorted(dico_noeuds.items(), key=lambda t: t[1], reverse=True) 
    return dico_noeuds_ordonnes

''' Fonction recursive 
Ajoute le noeud au cluster courant si  la densite du noeud est superieure a un seuil defini'''
def ajout_cluster(graphe, pourcentage_d, noeud, deja_vu, cluster_c) :
    if noeud in deja_vu : 
        return
    else :
#         print("ramousnif")
        deja_vu.append(noeud)
        for voisin in graphe[noeud] :
#             print(graphe.nodes[voisin]["density"])
            if voisin not in deja_vu and graphe.nodes[voisin]["density"] > graphe.nodes[noeud]["density"]*(1-pourcentage_d) :
#                 print(voisin)
                cluster_c.append(voisin)
                ajout_cluster(graphe, pourcentage_d, voisin, deja_vu, cluster_c)
                #if voisin not in deja_vu :
                #    deja_vu.append(voisin)
        #deja_vu.append(noeud)
        

''' Algo principal '''
def clustering_mcode(graphe, pourcentage_d, fluff, haircut, pourcentage_f) :
    ## Phase 1 : ponderation des noeuds
    nx.set_node_attributes(graphe, 0, "density")
    for noeud in graphe.nodes() :
        voisins = graphe[noeud]
        noeud_plus_voisins = list(voisins)
        noeud_plus_voisins.append(noeud)
#         print(noeud_plus_voisins)
        sous_graphe = graphe.subgraph(noeud_plus_voisins).copy()
#         print(nx.is_frozen(sous_graphe))
        val_k_core, graphe_k_core = trouve_highest_kcore(sous_graphe)
        densite_k_core = densite(graphe_k_core)
#         print(densite_k_core)
        poids_noeud = densite_k_core * val_k_core
#         print(poids_noeud)
        graphe.nodes[noeud]["density"] = poids_noeud
        
    ## Phase 2 : trouver les clusters
#     print("petit rat")
    deja_vu = []
    dico_noeuds_ordonnes = ordonner_par_densite(graphe)
#     print(dico_noeuds_ordonnes)
    clusters = []
    for elt in dico_noeuds_ordonnes :
        if elt[0] not in deja_vu :
            cluster_c = [elt[0]]
#             print(deja_vu)
            vieux_deja_vu = list(deja_vu)
            ajout_cluster(graphe, pourcentage_d, elt[0], deja_vu, cluster_c)
            clusters.append(cluster_c)
            deja_vu = list(vieux_deja_vu)
            deja_vu.append(elt)
#     print(deja_vu)        
#     print(clusters)
#     print(clusters)
    ## Phase 3 : post-processing
    a_enlever = []
    
    for c in clusters :
        k, g = trouve_highest_kcore(graphe.subgraph(c).copy())
        if k < 2 :
            a_enlever.append(c)
    
    for elt in a_enlever :
        clusters.remove(elt)
    
    a_enlever = []   
    a_ajouter = []     
    for c in clusters :
        if fluff :
                for elt in c :
                    voisins = list(graphe[elt])
                    voisins.append(elt)
                    densite_c = densite(graphe.subgraph(voisins))
#                     print("gros rat")
#                     print(deja_vu)
#                     print(len(deja_vu))
#                     print(densite_c)    
                    if densite_c > pourcentage_f :
                        for voisin in graphe[elt] : 
                            est_dans_un_cluster = False
                            for c2 in clusters :
                                if voisin in c2 :
                                    est_dans_un_cluster = True
                            if not est_dans_un_cluster :
                                c.append(voisin)
        if haircut : ## enlever les sommets de degre 1
                groupe_noeuds = []
                for elt in c :
                    if graphe.degree(elt) > 1 :
                        groupe_noeuds.append(elt)
                if len(groupe_noeuds) < graphe.number_of_nodes() :
                    a_enlever.append(c)
                    a_ajouter.append(groupe_noeuds)
                                            
    
    for elt in a_enlever :
        clusters.remove(elt)
    for elt in a_ajouter :
        clusters.append(elt)

#     print(len(clusters)) 
# #     print(len(clusters[0])) 
#     print(clusters)
    return clusters


''' Applique l'algo de clustering mcode a toutes les composantes connexes 
obtenues a partir d'un seuil defini '''
def clustering_mcode_composante(rep, depart_sim, typ_sim, taille_ext, seuil):
    with open(EXTENSION_PATH%rep+"composantes_connexes_"+depart_sim + "_" + typ_sim+ "_taille_"+str(taille_ext)+"_"+str(round(seuil,2))+".pickle", 'rb') as fichier_comp:
        #with open("Extensions/Metrique_toutes_aretes/composantes_connexes_extensions_toutes_aretes_coeffn1_a1_c1_0.6_0.6.pickle", 'rb') as fichier_comp:
    
            mon_depickler_comp = pickle.Unpickler(fichier_comp)
            composantes_connexes = mon_depickler_comp.load()
            
            with open(EXTENSION_PATH%rep+"graphe_complet_pondere_sim_%s_taille_%s.pickle"%(typ_sim, taille_ext), 'rb') as fichier_graphe_complet :
                mon_depickler_complet = pickle.Unpickler(fichier_graphe_complet)
                graphe_complet = mon_depickler_complet.load()
            
                nb_comp = 1
                for composante in composantes_connexes :
                    graphe_comp = nx.Graph()
                    for elt in composante :
                        graphe_comp.add_node(elt, **graphe_complet.nodes[elt])
                    
                    for u,v,data in graphe_complet.edges(data=True) :
                        if data["poids"] >= seuil :
                            u_new = -1
                            v_new = -1
                            for noeud, data_noeud in graphe_comp.nodes(data=True) :
            #                     print(data_noeud)
            #                     print(graphe_complet.nodes[u])
                                if graphe_complet.nodes[u]["nom"] == data_noeud["nom"] :
                                    u_new = noeud
                                if graphe_complet.nodes[v]["nom"] == data_noeud["nom"] :
                                    v_new = noeud
                            if u_new != -1  and v_new != -1 :
                                graphe_comp.add_edge(u_new,v_new,**data)
                    if len(composante) > 10  :
                        for elt in composante : 
                            print(graphe_comp.nodes[elt]["nom"])
                        max_taille_clusters = -1
                        clusters_max = []
                        val_d = -1
                        val_f = -1
                        clusters = [[]]
                        for i1 in np.arange(0.1, 0.15, 0.1) :
                                print(i1)
                            #for i2 in np.arange(0.9, 1.05, 0.1) :
#                                 print("i2")
#                                 print(i2)
                                clusters_prec = list(clusters)
                                #print(graphe_comp.nodes.data())
                                clusters = clustering_mcode(graphe_comp, i1, False, True, -1)
    #                             print(clusters)
                                for c in clusters :
                                    print("cluster")
                                    print(densite(graphe_comp.subgraph(c)))
#                                     for elt in c :
#                                         print(graphe_comp.nodes[elt]["nom"])
          
                                with open(EXTENSION_PATH%rep+"clustering_mcode/fichier_csv_grandes_composantes_clustering_mcode_%s_%s_%s.csv"%(round(seuil, 3),i1, nb_comp),'w') as fichier_csv:
                                    csvwriter = csv.writer(fichier_csv)
                                    csvwriter.writerow(["source", "target", "label"])
                                     
                                    for u,v,data in graphe_comp.edges(data=True) :
                                        csvwriter.writerow([u,v,round(data["poids"],2)])
                                     
                                with open(EXTENSION_PATH%rep+"clustering_mcode/fichier_csv_grandes_composantes_noeuds_clustering_mcode_%s_%s_%s.csv"%(round(seuil, 3),i1, nb_comp),'w') as fichier_csv:
                                    csvwriter = csv.writer(fichier_csv)
                                    csvwriter.writerow(["id", "label", "clustering_mcode"])
                            
                                    for noeud,data in graphe_comp.nodes(data=True) :
                                             
                                        num_cluster = []
                                        compteur = 1
                                        for c in clusters :
                                            if noeud in c :
                                                num_cluster.append(compteur)
                                            compteur += 1
                                             
                                        csvwriter.writerow([noeud, data["nom"], num_cluster])
                                
                                
                                different = False
                                for j in range(len(clusters)) :
                                    for l in range(len(clusters[j])) :
                                        if clusters[j][l] not in clusters_prec[j] :
                                            different = True
                                            break
                                    if different  :
                                        break
                                #if different :
                                print(clusters)        
                                
                                nombre_par_groupe = [0,0,0,0]
                                tot = 0
                                for c in clusters :
                                    compteur = 0
                                    for groupe in GROUPES_TOUTES_ARETES_MAX_4_10_07[:4] :
                                        compter = 0
                                        for elt in c :
                                            if graphe_comp.nodes[elt]["nom"] in groupe :
                                                compter += 1
                                        if compter > nombre_par_groupe[compteur] :
                                            nombre_par_groupe[compteur] = compter
     
                                        
                                        compteur += 1
                                    #print(nombre_par_groupe)
                                    
                                for elt in nombre_par_groupe :
                                    tot += elt
                                
                                #print(nombre_par_groupe)                
                                if max_taille_clusters < tot :
                                    max_taille_clusters = tot
                                    val_d  = i1
                                    #val_f = i2
                                    clusters_max = clusters 
                                
                            
                        print(val_d)
                        print(val_f)
                        print(clusters_max)
                        print(max_taille_clusters)
                        
                        for cluster in clusters : 
                            print("cluster :")
                            for c in cluster : 
                                print(graphe_comp.nodes[c]["nom"])
#                         print(graphe_comp.nodes[40]["nom"])
#                         print(graphe_comp.nodes[31]["nom"])
                        nb_comp += 1
                            
                        
                        
if __name__ == '__main__':
#     graphe = nx.dense_gnm_random_graph(5,5)
# #     graphe = nx.Graph()
# #     graphe.add_nodes_from([0,1,2,3,4])
# #     graphe.add_edges_from([(0,1), (0,2), (0,3), (0,4), (1,2), (1,4), (2,3), (3,4)])
#     print(graphe.nodes.data())
#     print(graphe.edges.data())
#     for u in graphe.nodes() :
#         graphe.nodes[u]["nom"] = "ramou"
#     
#     
#     print(graphe.nodes.data())
#     dico_noeuds = sorted(graphe.nodes.data()[1].items(), key=lambda t: t[0])
#     clustering_mcode(graphe,0.1, True, True, 0.7)

    clustering_mcode_composante("taille_max/result_k_max_4_8_toutes_aretes", "extensions", "toutes_aretes_coeff_all1_max", "taille_max", 0.7)
    
    
      