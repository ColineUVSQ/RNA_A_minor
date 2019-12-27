'''
Created on 16 juil. 2019

@author: coline

Comparaison des clusterings obtenus avec notre metrique et avec la RMSD a l'aide de l'index de Jaccard
(deux methodes de clustering testees : kmeans et recouvrant (Perez))

(version CaRNAval)
peut-etre a tester sur toutes donnees PDB
'''
from recup_data.constantes import EXTENSION_PATH, PATH_MMCIF
import pickle
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, IndexLocator,\
    FixedFormatter, FixedLocator
import seaborn as sns
from sklearn.cluster import KMeans
import numpy as np
from recup_data.interface_web import alter_base_extension_clustering,\
    alter_base_rmsd, ajout_clustering_perez
from recup_data import clustering_perez
from formatter import NullFormatter

''' calcul de l index de Jaccard entre le clustering1 et le clustering2 (clusterings recouvrants)
nb de paires d'elements clusterises dans le meme cluster dans les deux clustering / (nb de paires d'elements clusterises dans le meme cluster dans les deux clustering + nb de paires d'elements clusterises dans le meme cluster dans un des clusterings mais pas dans l'autre)'''
def calcul_sim_jaccard_pour_clustering_recouvrant(clustering1, clustering2, liste_elts) :
    nb_meme_cluster_ccprime = 0
    nb_meme_cluster_c = 0
    nb_meme_cluster_cprime = 0
    
    for i in range(len(liste_elts)) :
        for j in range(i+1, len(liste_elts)) :
            groupe_i1 = -1
            groupe_i2 = -1
            groupe_j1 = -1
            groupe_j2 = -1
            compteur = 0
            meme_groupe_1 = False
            meme_groupe_2 = False
            
            continuer = True
            while continuer and compteur < len(clustering1) : 
                groupe1 = clustering1[compteur]
                i_ok = False
                j_ok = False
                for i1 in range(len(groupe1)) :
                    if groupe1[i1] == liste_elts[i] :
                        i_ok = True
                    if groupe1[i1] == liste_elts[j] :
                        j_ok = True
                if i_ok and j_ok  :
                    
                    continuer = False
                    meme_groupe_1 = True
                compteur += 1
            
            compteur = 0
            continuer = True
            while continuer and compteur < len(clustering2) : 
                groupe2 = clustering2[compteur]
                i_ok = False
                j_ok = False
                for i2 in range(len(groupe2)) :
                    if groupe2[i2] == liste_elts[i] :
                        i_ok = True
                    if groupe2[i2] == liste_elts[j] :
                        j_ok = True
                        
#                 print(i+1, j+1)
#                 print(i_ok)
#                 print(j_ok)
                if i_ok and j_ok  :
#                     print(i+1, j+1)
                    continuer = False
                    meme_groupe_2 = True
                compteur += 1
                        
            if meme_groupe_1 and meme_groupe_2 :
                nb_meme_cluster_ccprime += 1
                #print(i+1,j+1)
            elif meme_groupe_1 :
                nb_meme_cluster_c += 1
            elif meme_groupe_2 :
                nb_meme_cluster_cprime += 1
                    
#     print(nb_meme_cluster_c)
#     print(nb_meme_cluster_cprime)
#     print(nb_meme_cluster_ccprime)
    if (nb_meme_cluster_ccprime+nb_meme_cluster_c+nb_meme_cluster_cprime) != 0 :
        tab_res = nb_meme_cluster_ccprime/(nb_meme_cluster_ccprime+nb_meme_cluster_c+nb_meme_cluster_cprime)
    else : 
        tab_res = 0
    return tab_res, (nb_meme_cluster_c, nb_meme_cluster_cprime, nb_meme_cluster_ccprime)


''' calcul de l index de Jaccard entre le clustering1 et le clustering2 (clusterings non recouvrants)
nb de paires d'elements clusterises dans le meme cluster dans les deux clustering / (nb de paires d'elements clusterises dans le meme cluster dans les deux clustering + nb de paires d'elements clusterises dans le meme cluster dans un des clusterings mais pas dans l'autre)'''
def calcul_sim_jaccard(clustering1, clustering2):
    nb_meme_cluster_ccprime = 0
    nb_meme_cluster_c = 0
    nb_meme_cluster_cprime = 0
    
    for groupe1 in clustering1 :
        for i1 in range(len(groupe1)) :
            for j1 in range(i1+1, len(groupe1)) :
                groupe_i2 = -1
                groupe_j2 = -1
                compteur_2 = 1
                for groupe2 in clustering2 :
                    for i2 in range(len(groupe2)) :
                        if groupe1[i1] == groupe2[i2] :
                            groupe_i2 = compteur_2
                        if groupe1[j1] == groupe2[i2] : 
                            groupe_j2 = compteur_2
                    compteur_2 += 1
                
                if groupe_i2 != -1 and groupe_j2 != -1 :
                    if groupe_i2 == groupe_j2 :
                        nb_meme_cluster_ccprime += 1
                    else :
                        nb_meme_cluster_c += 1
                else :
                    print("bizarre")
                    print(i1, j1)
                    #exit()
                    
    for groupe2 in clustering2 :
        for i2 in range(len(groupe2)) :
            for j2 in range(i2+1, len(groupe2)) :
                groupe_i1 = -1
                groupe_j1 = -1
                compteur_1 = 1
                for groupe1 in clustering1 :
                    for i1 in range(len(groupe1)) :
                        if groupe2[i2] == groupe1[i1] :
                            groupe_i1 = compteur_1
#                             print(groupe2[i2])
                        if groupe2[j2] == groupe1[i1] : 
                            groupe_j1 = compteur_1
#                             print(groupe2[j2])
                    compteur_1 += 1
                
                if groupe_i1 != -1 and groupe_j1 != -1 :
                    if groupe_i1 != groupe_j1 :
                        nb_meme_cluster_cprime += 1
                else :
                    print("bizarre")
                    #exit()
                    
#     print(nb_meme_cluster_c)
#     print(nb_meme_cluster_cprime)
#     print(nb_meme_cluster_ccprime)
    if (nb_meme_cluster_ccprime+nb_meme_cluster_c+nb_meme_cluster_cprime) != 0 :
        tab_res = nb_meme_cluster_ccprime/(nb_meme_cluster_ccprime+nb_meme_cluster_c+nb_meme_cluster_cprime)
    else : 
        tab_res = 0
    return tab_res, (nb_meme_cluster_c, nb_meme_cluster_cprime, nb_meme_cluster_ccprime)
                        
    
'''Renvoie le clustering kmeans de la matrice en entree sous forme d'une liste de listes (version sim)'''                     
def clustering_kmeans_sim(graphe_complet_sans_deux_elts, matrice, nb_clusters):
        
        clustering = KMeans(n_clusters=nb_clusters).fit(matrice)
        #print(clustering.labels_)
          
        compteur = 0
        tab_clustering_sim = [[]]
        for elt in clustering.labels_ :
            if elt != -1 :
                for _ in range(len(tab_clustering_sim), elt+1) :
                    tab_clustering_sim.append([])
                tab_clustering_sim[elt].append(graphe_complet_sans_deux_elts.nodes[compteur]["nom"])
            compteur += 1
#         print(len(tab_clustering_sim))    
#         print(tab_clustering_sim)
        
        return tab_clustering_sim

'''Renvoie le clustering kmeans de la matrice en entree sous forme d'une liste de listes (version rmsd mais ca a l'air d'etre les memes)'''
def clustering_kmeans_rmsd(graphe_complet_sans_deux_elts, matrice, nb_clusters):
    
#                             exit()
#             print(len(matrice))
#             print(matrice)                 
            clustering = KMeans(n_clusters=nb_clusters).fit(matrice)
            #print(clustering.labels_)
             
            compteur = 0
            tab_clustering_rmsd = [[]]
            for elt in clustering.labels_ :
                if elt != -1 :
                    for _ in range(len(tab_clustering_rmsd), elt+1) :
                        tab_clustering_rmsd.append([])
                    tab_clustering_rmsd[elt].append(graphe_complet_sans_deux_elts.nodes[compteur]["nom"])
                compteur += 1
#             print(len(tab_clustering_rmsd))    
#             print(tab_clustering_rmsd)
        
            return tab_clustering_rmsd
             
#             compteur = 1
#             for groupe in tab_clustering :
#                 for elt in groupe :
#                     if elt != '' :
#                         alter_base_extension_clustering(elt, compteur,"clustering_rmsd")
#                 compteur += 1
        
'''Effectue les clusterings kmeans (pour un nb de clusters variant de 2 à 14) avec les valeurs de sim et les valeurs de RMSD, puis calcule l'index de Jaccard entre les deux clusterings
stocke tous les resultats (nb de clusters, valeur de Jaccard, clusterings) dans un fichier .pickle
on choisit la taille d'extension dans la fonction'''        
def comp_kmeans() :
    taille_ext = 5
    with open(EXTENSION_PATH%taille_ext+"graphe_complet_pondere_sim_toutes_aretes_coeff_all1_taille_%d.pickle"%taille_ext, 'rb') as fichier_graphes :
        mon_depickler = pickle.Unpickler(fichier_graphes)
        graphe_complet = mon_depickler.load()
#         
        print(graphe_complet.number_of_nodes())
        print(graphe_complet.number_of_edges())
#        
        a_enlever = []
        for noeud, data in graphe_complet.nodes(data=True) :
            if data["nom"] in ['2XD0_V_36_21', '3JCS_1_282_1', '4FAU_A_207_4'] :
                a_enlever.append(noeud)
                     
        for elt in a_enlever :
            graphe_complet.remove_node(elt)
             
        graphe_complet_sans_deux_elts = nx.convert_node_labels_to_integers(graphe_complet)
        print(graphe_complet_sans_deux_elts.number_of_nodes())
        print(graphe_complet_sans_deux_elts.edges.data())
             
            
        matrice_sim = [[0] *graphe_complet_sans_deux_elts.number_of_nodes() for _ in range(graphe_complet_sans_deux_elts.number_of_nodes())]
#         print(matrice)
#         print(graphe_complet.edges.data())
#         print(graphe_complet.edges[1,7])
              
        compteur = 0
        for i in range(graphe_complet_sans_deux_elts.number_of_nodes()) :
            for j in range(graphe_complet_sans_deux_elts.number_of_nodes()) :
                if i != j :
                    if (i,j) in graphe_complet_sans_deux_elts.edges() or (j,i) in graphe_complet_sans_deux_elts.edges() :
                        matrice_sim[i][j] = 1 - graphe_complet_sans_deux_elts.edges[i,j]["poids"]
                        compteur += 1
                    else :
#                         print(i,j)
                        matrice_sim[i][j] = 50.0
#         print(compteur)
#         print(matrice)
            
             
        with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_1_normalise.pickle"%taille_ext, 'rb') as fichier_rmsd :
            mon_depickler = pickle.Unpickler(fichier_rmsd)
            dico_rmsd = mon_depickler.load()
                
#             mini = min(list(dico_rmsd.values()))
#             maxi = max(list(dico_rmsd.values()))
#             for cle in dico_rmsd.keys() :
#                 cle_1 = cle[0][8:len(cle[0])-13]
#                 cle_2 = cle[1][8:len(cle[1])-13]
#                 if dico_rmsd[cle] != None :
#                     alter_base_rmsd(cle_1, cle_2, dico_rmsd[cle])
                 
            matrice_rmsd = [[0] *graphe_complet_sans_deux_elts.number_of_nodes() for _ in range(graphe_complet_sans_deux_elts.number_of_nodes())]
#             print(matrice)
    #         print(graphe_complet.edges.data())
    #         print(graphe_complet.edges[1,7])
                 
            compteur = 0
            for i in range(graphe_complet_sans_deux_elts.number_of_nodes()) :
                for j in range(graphe_complet_sans_deux_elts.number_of_nodes()) :
                    if i != j :
                        if ("fichier_"+graphe_complet_sans_deux_elts.nodes[i]["nom"]+"_taille_%d.pdb"%taille_ext, "fichier_"+graphe_complet_sans_deux_elts.nodes[j]["nom"]+"_taille_%d.pdb"%taille_ext) in dico_rmsd.keys() and dico_rmsd[("fichier_"+graphe_complet_sans_deux_elts.nodes[i]["nom"]+"_taille_%d.pdb"%taille_ext, "fichier_"+graphe_complet_sans_deux_elts.nodes[j]["nom"]+"_taille_%d.pdb"%taille_ext)] != None :
                            matrice_rmsd[i][j] = dico_rmsd[("fichier_"+graphe_complet_sans_deux_elts.nodes[i]["nom"]+"_taille_%d.pdb"%taille_ext, "fichier_"+graphe_complet_sans_deux_elts.nodes[j]["nom"]+"_taille_%d.pdb"%taille_ext)]
                            compteur += 1
                        elif ("fichier_"+graphe_complet_sans_deux_elts.nodes[j]["nom"]+"_taille_%d.pdb"%taille_ext, "fichier_"+graphe_complet_sans_deux_elts.nodes[i]["nom"]+"_taille_%d.pdb"%taille_ext) in dico_rmsd.keys() and dico_rmsd[("fichier_"+graphe_complet_sans_deux_elts.nodes[j]["nom"]+"_taille_%d.pdb"%taille_ext, "fichier_"+graphe_complet_sans_deux_elts.nodes[i]["nom"]+"_taille_%d.pdb"%taille_ext)] != None :
                            matrice_rmsd[i][j] = dico_rmsd[("fichier_"+graphe_complet_sans_deux_elts.nodes[j]["nom"]+"_taille_%d.pdb"%taille_ext, "fichier_"+graphe_complet_sans_deux_elts.nodes[i]["nom"]+"_taille_%d.pdb"%taille_ext)]
                            compteur += 1
                        else :
                            matrice_rmsd[i][j] = 100.0 
            
        tab_res = []
        maxi_jaccard = [0, 0, 0]
        tab_maxi_tot = []
             
        for i in range(2, 15) :
            for j in range(2, 15) :
                print(i)
                print(j)
                maxi_jaccard = [0, 0, 0]
                for _ in range(200) :
                    tab_clustering_sim = clustering_kmeans_sim(graphe_complet_sans_deux_elts, matrice_sim, i)
                    tab_clustering_rmsd = clustering_kmeans_rmsd(graphe_complet_sans_deux_elts, matrice_rmsd, j)
                         
                    res = calcul_sim_jaccard(tab_clustering_rmsd, tab_clustering_sim)
                    if res[0] > maxi_jaccard[2] :
                        maxi_jaccard = list([i,j, res[0], tab_clustering_sim, tab_clustering_rmsd, res[1]])
                tab_maxi_tot.append(maxi_jaccard)   
                print(tab_maxi_tot)
                         
        with open("fichier_comp_sim_rmsd_par_nb_clusters_tot_normalise_taille_5.pickle", 'wb') as fichier :
            mon_pickler = pickle.Pickler(fichier)
            mon_pickler.dump(tab_maxi_tot)    
 
'''Effectue les clusterings recouvrants (pour un seuil variant de 0 a 1 avec un pas de 0.1) avec les valeurs de sim et les valeurs de RMSD, puis calcule l'index de Jaccard entre les deux clusterings
stocke tous les resultats (nb de clusters, valeur de Jaccard, clusterings) dans un fichier .pickle
on choisit la taille d'extension dans la fonction'''   
def comp_perez():
    taille_ext = 4
    with open(EXTENSION_PATH%taille_ext+"graphe_complet_pondere_sim_toutes_aretes_coeff_all1_taille_%d.pickle"%taille_ext, 'rb') as fichier_graphes :
        mon_depickler = pickle.Unpickler(fichier_graphes)
        graphe_complet = mon_depickler.load()
#         
        print(graphe_complet.number_of_nodes())
        print(graphe_complet.number_of_edges())
#        
        liste_elts = []
        a_enlever = []
        for noeud, data in graphe_complet.nodes(data=True) :
            if data["nom"] in ['2XD0_V_36_21', '3JCS_1_282_1', '4FAU_A_207_4'] :
                a_enlever.append(noeud)
            else :
                liste_elts.append(data["nom"])
                     
        for elt in a_enlever :
            graphe_complet.remove_node(elt)
             
        graphe_complet_sans_deux_elts = nx.convert_node_labels_to_integers(graphe_complet)
        print(graphe_complet_sans_deux_elts.number_of_nodes())
        print(graphe_complet_sans_deux_elts.edges.data())
             
#         print(compteur)
#         print(matrice)
        
        nx.set_edge_attributes(graphe_complet_sans_deux_elts, -1, "rmsd")    
             
        with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_1_normalise.pickle"%taille_ext, 'rb') as fichier_rmsd :
            mon_depickler = pickle.Unpickler(fichier_rmsd)
            dico_rmsd = mon_depickler.load()
                
                 
            compteur = 0
            for i in range(graphe_complet_sans_deux_elts.number_of_nodes()) :
                for j in range(graphe_complet_sans_deux_elts.number_of_nodes()) :
                    if i != j :
                        if ("fichier_"+graphe_complet_sans_deux_elts.nodes[i]["nom"]+"_taille_%d.pdb"%taille_ext, "fichier_"+graphe_complet_sans_deux_elts.nodes[j]["nom"]+"_taille_%d.pdb"%taille_ext) in dico_rmsd.keys() and dico_rmsd[("fichier_"+graphe_complet_sans_deux_elts.nodes[i]["nom"]+"_taille_%d.pdb"%taille_ext, "fichier_"+graphe_complet_sans_deux_elts.nodes[j]["nom"]+"_taille_%d.pdb"%taille_ext)] != None :
                            graphe_complet_sans_deux_elts.edges[i,j]["rmsd"] = dico_rmsd[("fichier_"+graphe_complet_sans_deux_elts.nodes[i]["nom"]+"_taille_%d.pdb"%taille_ext, "fichier_"+graphe_complet_sans_deux_elts.nodes[j]["nom"]+"_taille_%d.pdb"%taille_ext)]
                            compteur += 1
                        elif ("fichier_"+graphe_complet_sans_deux_elts.nodes[j]["nom"]+"_taille_%d.pdb"%taille_ext, "fichier_"+graphe_complet_sans_deux_elts.nodes[i]["nom"]+"_taille_%d.pdb"%taille_ext) in dico_rmsd.keys() and dico_rmsd[("fichier_"+graphe_complet_sans_deux_elts.nodes[j]["nom"]+"_taille_%d.pdb"%taille_ext, "fichier_"+graphe_complet_sans_deux_elts.nodes[i]["nom"]+"_taille_%d.pdb"%taille_ext)] != None :
                            graphe_complet_sans_deux_elts.edges[i,j]["rmsd"] = dico_rmsd[("fichier_"+graphe_complet_sans_deux_elts.nodes[j]["nom"]+"_taille_%d.pdb"%taille_ext, "fichier_"+graphe_complet_sans_deux_elts.nodes[i]["nom"]+"_taille_%d.pdb"%taille_ext)]
                            compteur += 1
                        else :
                            graphe_complet_sans_deux_elts.edges[i,j]["rmsd"] = 100.0 
            
        tab_maxi_tot = []
             
        for i in np.arange(0,1.1,0.1) :
            graphe_seuil_i = graphe_complet_sans_deux_elts.copy()
            a_enlever = []
            for u,v,data in graphe_complet_sans_deux_elts.edges(data=True) : 
                if data["poids"] < i : 
                    a_enlever.append((u,v))
          
            for elt in a_enlever :
                graphe_seuil_i.remove_edge(elt[0], elt[1])
            
            for j in np.arange(1.0,-0.01,-0.1) :
                
                graphe_seuil_j = graphe_complet_sans_deux_elts.copy()
                a_enlever = []
                for u,v,data in graphe_complet_sans_deux_elts.edges(data=True) : 
                    if data["rmsd"] > j or j < 0.1: 
                        a_enlever.append((u,v))
              
                for elt in a_enlever :
                    graphe_seuil_j.remove_edge(elt[0], elt[1])
                print(i)
                print(j)
                print(graphe_seuil_i.number_of_edges())
                print(graphe_seuil_j.number_of_edges())
                
                clustering_sim = clustering_perez.algo_principal(graphe_seuil_i)
                clustering_rmsd = clustering_perez.algo_principal(graphe_seuil_j)
                
#                 for noeud, data in graphe_seuil_j.nodes(data=True) :
#                     print(len(graphe_seuil_j[noeud]))
#                     if len(graphe_seuil_j[noeud]) < 4  :
#                         print(data["nom"])
                         
                res = calcul_sim_jaccard_pour_clustering_recouvrant(clustering_sim, clustering_rmsd, liste_elts)
                tab_maxi_tot.append([i,j, res[0], clustering_sim, clustering_rmsd, res[1]])   
                print(tab_maxi_tot)
                print(liste_elts)
                         
        with open("fichier_comp_sim_rmsd_par_seuil_tot_normalise_taille_4_singleton.pickle", 'wb') as fichier :
            mon_pickler = pickle.Pickler(fichier)
            mon_pickler.dump(tab_maxi_tot)  
 
''' Affiche la distribution des index de Jaccard pour les comparaisons des clusterings recouvrants pour toutes les paires de seuils testees '''
def distrib_perez():
    with open("fichier_comp_sim_rmsd_par_seuil_tot_normalise_taille_4_singleton.pickle", 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        tab_max_par_nb_clusters = mon_depickler.load()
        print(tab_max_par_nb_clusters)
        liste_sim_jaccard = []
        compteur = 1
        for elt in tab_max_par_nb_clusters :
            if compteur < 15 :
                print(elt)
                print(len(elt[3]))
                #print(len(elt[3][1]))
                print(len(elt[4]))
                #print(len(elt[4][1]))
            compteur += 1
                 
            liste_sim_jaccard.append(elt[2])
             
        liste_x = []
        liste_locator_x = []
        
        compteur = 0
        for i in range(11) :
            k = 0
            for j in np.arange(1.0,-0.1,-0.1) :
                if k%2 == 0 :
                    liste_x.append(round(j,1))
                    #if k != 0 or compteur == 0 :
                    liste_locator_x.append(k+11*compteur)
                k = k+1
            #liste_locator_x.append(11+10*compteur)
            compteur += 1
        print(liste_x)
        
        matplotlib.rc('xtick', labelsize=6)
        ax = plt.gca()
        #fig, ax = plt.subplots()
        ax.stem(liste_sim_jaccard, markerfmt=' ')
        
        #ax.set_xticklabels(list(liste_x))
        
        ax.xaxis.set_major_locator(FixedLocator(liste_locator_x))
        ax.xaxis.set_major_formatter(FixedFormatter(list(liste_x)))
        
        ax.xaxis.set_minor_locator(FixedLocator(range(0, 121, 1)))
        #ax.xaxis.set_minor_formatter(NullFormatter())
        #ax.xaxis.set_minor_locator(MultipleLocator(0.25))
        
        #ax.set_xticks(range(121))
        #ax.set_xticklabels(list(liste_x))
        
        matplotlib.rc('xtick', labelsize=10)
        ax2 = ax.twiny()
  
        # Decide the ticklabel position in the new x-axis,
        # then convert them to the position in the old x-axis # labels of the xticklabels: the position in the new x-axis
          
        #newpos   = [k2degc(t) for t in newlabel]   # position of the xticklabels in the old x-axis
        
        
        ax2.set_xticks(np.arange(0,121,11))
        ax2.set_xticklabels([round(x,1) for x in np.arange(0,1.1, 0.1)])
          
#         ax2.xaxis.set_major_locator(np.arange(0,169,13))
#         ax2.xaxis.set_major_formatter(newlabel)
          
        # For the minor ticks, use no labels; default NullFormatter.
        #ax2.xaxis.set_minor_locator(range(169))
          
          
        ax2.xaxis.set_ticks_position('bottom') # set the position of the second x-axis to bottom
        ax2.xaxis.set_label_position('bottom') # set the position of the second x-axis to bottom
        ax2.spines['bottom'].set_position(('outward', 36))
        ax2.set_xlabel('Seuil de similarité')
        ax2.set_xlim(ax.get_xlim())
        
        ax.set_ylim((0,1.0))
             
        ax.set_xlabel("Seuil de RMSD")
        ax.set_ylabel("Indice de Jaccard")
        plt.title("Comparaison clustering sim/RMSD avec seuils et clustering recouvrant")
        plt.show()

''' Affiche la distribution des index de Jaccard pour les comparaisons des clusterings kmeans pour toutes les paires de nb de clusters testees '''
def distrib_kmeans():
    with open("fichier_comp_sim_rmsd_par_nb_clusters_tot_normalise_taille_5.pickle", 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        tab_max_par_nb_clusters = mon_depickler.load()
        print(tab_max_par_nb_clusters)
        liste_sim_jaccard = []
        compteur = 1
        for elt in tab_max_par_nb_clusters :
            if compteur < 15 :
                print(elt)
                print(len(elt[3]))
                #print(len(elt[3][1]))
                print(len(elt[4]))
                #print(len(elt[4][1]))
            compteur += 1
                 
            liste_sim_jaccard.append(elt[2])
             
        liste_x = []
        liste_locator_x = []
        
        compteur = 0
        for i in range(2,15) :
            k = 0
            for j in np.arange(2,15) :
                if k%2 == 0 :
                    liste_x.append(j)
                    #if k != 0 or compteur == 0 :
                    liste_locator_x.append(k+13*compteur)
                k = k+1
            #liste_locator_x.append(11+10*compteur)
            compteur += 1
        print(liste_x)
        
        matplotlib.rc('xtick', labelsize=6)
        ax = plt.gca()
        #fig, ax = plt.subplots()
        ax.stem(liste_sim_jaccard, markerfmt=' ')
        
        #ax.set_xticklabels(list(liste_x))
        
        ax.xaxis.set_major_locator(FixedLocator(liste_locator_x))
        ax.xaxis.set_major_formatter(FixedFormatter(list(liste_x)))
        
        ax.xaxis.set_minor_locator(FixedLocator(range(169)))
        #ax.xaxis.set_minor_formatter(NullFormatter())
        #ax.xaxis.set_minor_locator(MultipleLocator(0.25))
        
        #ax.set_xticks(range(121))
        #ax.set_xticklabels(list(liste_x))
        
        matplotlib.rc('xtick', labelsize=10)
        ax2 = ax.twiny()
  
        # Decide the ticklabel position in the new x-axis,
        # then convert them to the position in the old x-axis # labels of the xticklabels: the position in the new x-axis
          
        #newpos   = [k2degc(t) for t in newlabel]   # position of the xticklabels in the old x-axis
        
        
        ax2.set_xticks(np.arange(0,169,13))
        ax2.set_xticklabels(range(2,15))
          
#         ax2.xaxis.set_major_locator(np.arange(0,169,13))
#         ax2.xaxis.set_major_formatter(newlabel)
          
        # For the minor ticks, use no labels; default NullFormatter.
        #ax2.xaxis.set_minor_locator(range(169))
          
          
        ax2.xaxis.set_ticks_position('bottom') # set the position of the second x-axis to bottom
        ax2.xaxis.set_label_position('bottom') # set the position of the second x-axis to bottom
        ax2.spines['bottom'].set_position(('outward', 36))
        ax2.set_xlabel('Nombre de clusters clustering similarité')
        ax2.set_xlim(ax.get_xlim())
        
        ax.set_ylim((0,1.0))
             
        ax.set_xlabel("Nombre de clusters clustering RMSD")
        ax.set_ylabel("Indice de Jaccard")
        plt.title("Comparaison clustering sim/RMSD avec kmeans")
        plt.show()
        
''' 22/08/19
cherche la distribution en valeurs de sim des deux groupes trouves avec kmeans sur les valeurs de RMSD '''
def distrib_kmeans_sim_rmsd_2paquets():
    with open("fichier_comp_sim_rmsd_par_nb_clusters_tot_normalise_taille_4.pickle", 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        tab_max_par_nb_clusters = mon_depickler.load()
        
        for elt in tab_max_par_nb_clusters : 
            if elt[1] == 2 : 
                clustering_rmsd = elt[4]
                print(len(clustering_rmsd[0]))
                print(len(clustering_rmsd[1]))
                
        
        #with open(EXTENSION_PATH%4+"sim_extensions_toutes_aretes_coeff_all1_taille_%d.pickle"%4, 'rb') as fichier_rmsd :
        with open(PATH_MMCIF+"fichiers_rmsd_taille_4_que_carbone_1_normalise.pickle", 'rb') as fichier_rmsd :
            mon_depickler = pickle.Unpickler(fichier_rmsd)
            dico_rmsd = mon_depickler.load()
            
            print(len(dico_rmsd))
            
            compteur = 1
            print(clustering_rmsd)
            couples_dans_groupe = []
            for elt in clustering_rmsd :
                print(len(elt))
                tab_intra_cluster = []
                for i in range(len(elt)) :
                    for j in range(i+1, len(elt)) :
                        cle_1 = "fichier_"+ elt[i] + "_taille_4.pdb"
                        cle_2 = "fichier_"+ elt[j] + "_taille_4.pdb"
                        #cle_1 = elt[i]
                        #cle_2 = elt[j]
                        if (cle_1, cle_2) in dico_rmsd.keys() :
                            tab_intra_cluster.append(dico_rmsd[(cle_1, cle_2)])   
                            couples_dans_groupe.append((cle_1, cle_2))  
                        else :
                            tab_intra_cluster.append(dico_rmsd[(cle_2, cle_1)])
                            couples_dans_groupe.append((cle_2, cle_1))
                
                ax = plt.gca()            
                plt.plot(tab_intra_cluster)
                print(len(tab_intra_cluster))
                #ax.set_ylabel('Valeur de sim')
                ax.set_ylabel('Valeur de RMSD')
                ax.set_yticks(np.arange(0, 1.1, 0.1))
                #plt.title("Distribution des valeurs de sim pour le groupe %d \n du clustering kmeans k=2 sur la RMSD"%compteur)
                plt.title("Distribution des valeurs de RMSD pour le groupe %d \n du clustering kmeans k=2 sur la RMSD"%compteur)
                #plt.savefig(EXTENSION_PATH%4+"distrib_sim_kmeans_rmsd_2groupes_groupe_%d.png"%compteur)
                plt.savefig(EXTENSION_PATH%4+"distrib_rmsd_kmeans_rmsd_2groupes_groupe_%d.png"%compteur)
                #plt.show()
                plt.clf()
                compteur += 1
                
            tab_inter_cluster = []
            for cle in dico_rmsd.keys() :
                if cle not in couples_dans_groupe and cle[0] not in ['fichier_2XD0_V_36_21_taille_4.pdb', 'fichier_3JCS_1_282_1_taille_4.pdb'] and cle[1] not in ['fichier_2XD0_V_36_21_taille_4.pdb', 'fichier_3JCS_1_282_1_taille_4.pdb']:
                    #print(cle)
                    tab_inter_cluster.append(dico_rmsd[cle])
            
            ax = plt.gca()
            plt.plot(tab_inter_cluster)
            ax.set_ylabel('Valeur de sim')
            ax.set_yticks(np.arange(0, 1.1, 0.1))
            print(len(tab_inter_cluster))
            
            #plt.title("Distribution des valeurs de sim entre les deux groupes \n du clustering kmeans k=2 sur la RMSD")
            #plt.savefig(EXTENSION_PATH%4+"distrib_sim_kmeans_rmsd_2groupes_inter_groupe.png")
            plt.title("Distribution des valeurs de RMSD entre les deux groupes \n du clustering kmeans k=2 sur la RMSD")
            plt.savefig(EXTENSION_PATH%4+"distrib_rmsd_kmeans_rmsd_2groupes_inter_groupe_%d.png")
            #plt.show()      
        

''' 23/08/19
Cherche si les elts d'un meme groupe sim 0.7 avec algo PEREZ appartiennent a un meme groupe dans rmsd O.3 algo Perez
'''          
def comp_sim_rmsd_07_03():
    with open("fichier_comp_sim_rmsd_par_seuil_tot_normalise_taille_4.pickle", 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        tab_max_par_seuil = mon_depickler.load()
        
        for elt in tab_max_par_seuil : 
            print(elt)
            if round(elt[0],1) == 0.7 and round(elt[1],1) == 0.3 : 
                clustering_sim = elt[3]
                clustering_rmsd = elt[4]
                
        print(clustering_sim)
        print(clustering_rmsd)
        
        for groupe in clustering_sim :
            pos_groupe = -1
            
            for elt in groupe :
                pos_elt = -1
                compteur = 0
                for groupe2 in clustering_rmsd :
                    if elt in groupe2 :
                        pos_elt = compteur
                    compteur += 1
                
                if pos_groupe == -1  :
                    pos_groupe = pos_elt
                if pos_elt != pos_groupe :
                    print(elt)
                    print(groupe)
                      
                    
                    
             
if __name__ == '__main__':
    
    #distrib_kmeans_sim_rmsd_2paquets()
    comp_sim_rmsd_07_03()
    
    
    #distrib_kmeans()
    #distrib_perez()
    #comp_perez() 
    #comp_kmeans()    
        
        
    
#     with open("script_python_alter_extension.py", 'a') as fichier_py : 
#         fichier_py.write("from aminor.models import Extension\n")
#     tab_max_par_nb_clusters = []
#     for i in range(2, 15) :
#         tab_res = []
#         maxi_jaccard = [(0,0,0)]
#         for _ in range(200) :
#             tab_res.append([comp(4, i)])
#             if tab_res[len(tab_res)-1][0] > maxi_jaccard[0] :
#                 maxi_jaccard = list(tab_res[len(tab_res)-1])
#         print(tab_res)
#         print(maxi_jaccard)
#         tab_max_par_nb_clusters.append(maxi_jaccard)
# # # # #     
#     print(tab_max_par_nb_clusters)
#       
#     with open("fichier_comp_sim_rmsd_par_nb_clusters_normalise.pickle", 'wb') as fichier :
#         mon_pickler = pickle.Pickler(fichier)
#         mon_pickler.dump(tab_max_par_nb_clusters)


    
        #plt.close()
#     
#     with open("fichier_comp_sim_rmsd_par_nb_clusters_tot_normalise.pickle", 'rb') as fichier :
#         mon_depickler = pickle.Unpickler(fichier)
#         tab_max_par_nb_clusters = mon_depickler.load()
#         dico_extension = {}
#         compteur_seuil = 0
#         print(tab_max_par_nb_clusters)
        ## version pour clustering perez
#         for elt in tab_max_par_nb_clusters :
# #             print(elt[3])
# #             print(elt[0])
#             if round(elt[1],1) == 0.5 :
#                 print(elt[1])
#                 print(elt[2])
#                 print(elt[3])
#                 compteur = 1
#                 for groupe in elt[3] :
#                     for elt in groupe :
#                         if elt not in dico_extension.keys() :
#                             dico_extension.update({elt : [[compteur-1]]})
#                         else :
#                             if len(dico_extension[elt]) < compteur_seuil+1 :
#                                 dico_extension[elt].append([])
#                             dico_extension[elt][len(dico_extension[elt])-1].append(compteur-1)     
#                     compteur += 1 
#                     
#                     for cle in dico_extension.keys() :
#                         if len(dico_extension[cle]) < compteur_seuil+1 :
#                             dico_extension[cle].append([])
#                 compteur_seuil +=1
#                 #if compteur_seuil == 2 :
#                     #break
#                     
#         print(dico_extension)
#         ajout_clustering_perez(dico_extension)

## version pour clustering kmeans

#         for elt in tab_max_par_nb_clusters :
# #             print(elt[3])
# #             print(elt[0])
#             if round(elt[0],1) == 2 :
#                 print(elt[1])
#                 print(elt[2])
#                 print(elt[3])
#                 compteur = 1
#                 nb_elts = 0
#                 for groupe in elt[4] :
#                     for elt in groupe :
#                         if elt not in dico_extension.keys() :
#                             dico_extension.update({elt : [compteur]})
#                         else :
#                             dico_extension[elt].append(compteur) 
#                         nb_elts +=1    
#                     compteur += 1 
#                 
#                     for cle in dico_extension.keys() :
#                         if len(dico_extension[cle]) < compteur_seuil :
#                             dico_extension[cle].append("")
#                 print(nb_elts)
#                 compteur_seuil +=1
#                 #if compteur_seuil == 2 :
#                     #break
#                     
#         print(dico_extension)
#         ajout_clustering_perez(dico_extension)
        
# affichage clustering sur Django         
#         print(tab_max_par_nb_clusters[3][0])
#         
#         compteur = 1
#         for groupe in tab_max_par_nb_clusters[3][0][2] :
#             for elt in groupe :
#                 if elt != '' :
#                     alter_base_extension_clustering(elt, compteur,"clustering_rmsd")
#             compteur += 1
#     
#         compteur = 1
#         for groupe in tab_max_par_nb_clusters[3][0][1] :
#             for elt in groupe :
#                 if elt != '' :
#                     alter_base_extension_clustering(elt, compteur,"clustering")
#             compteur += 1
    
    ##test
#     clustering1 = [[1,2,3],[3,4],[5,6]]
#     clustering2 = [[1,2,3], [4,6,5]]
#        
#     print(calcul_sim_jaccard_pour_clustering_recouvrant(clustering2, clustering1,[1,2,3,4,5,6]))


##Normalisation des valeurs de rmsd pour obtenir des valeurs entre 0 et 1
#     with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_1.pickle"%5, 'rb') as fichier_rmsd :
#             mon_depickler = pickle.Unpickler(fichier_rmsd)
#             dico_rmsd = mon_depickler.load()
#             print(len(dico_rmsd))
# #             for x in dico_rmsd.values() :
# #                 if x != None and x < 1.0 :
# #                     print(x)
#  
# 
#             mini = min(list([x for x in dico_rmsd.values() if x != None]))
#             maxi = max(list([x for x in dico_rmsd.values() if x != None]))
#             print(mini)
#             print(maxi)
#                 
#             #print(dico_rmsd.values())
#                 
#             dico_rmsd_normalise = {} 
#             for elt in dico_rmsd.keys() :
#                 if dico_rmsd[elt] != None :
#                     dico_rmsd_normalise.update({elt : (dico_rmsd[elt]-mini)/(maxi-mini)})
#                     cle_1 = elt[0][8:len(elt[0])-13]
#                     cle_2 = elt[1][8:len(elt[1])-13]
#                     alter_base_rmsd(cle_1, cle_2, (dico_rmsd[elt]-mini)/(maxi-mini))
#                     #alter_base_rmsd(cle_1, cle_2, dico_rmsd[elt])
#                 else :
#                     dico_rmsd_normalise.update({elt : None})
#                     print(elt)
#                       
#             #print(dico_rmsd_normalise)
#              
#             with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_1_normalise.pickle"%5, 'wb') as fichier_rmsd_normalise :
#                 mon_pickler = pickle.Pickler(fichier_rmsd_normalise)
#                 mon_pickler.dump(dico_rmsd_normalise)
            
                
            #print(tab_rmsd_normalise)
            
            
            
            