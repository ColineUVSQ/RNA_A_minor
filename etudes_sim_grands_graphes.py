'''
Created on 4 févr. 2019

@author: coline

affichage des composantes connexes associees aux comparaisons des graphes globaux
'''
import pickle
import networkx as nx
import matplotlib.pyplot as plt
from recup_data.calcul_sim import calcul_sim_non_cov_sans_motif
from recup_data.sous_graphe_commun_max_version_grands_graphes import nombre_aretes_arcs

def draw_composantes(graphe_comp) :

    nx.set_node_attributes(graphe_comp, (33,33), "coordonnees")
    pos = nx.circular_layout(graphe_comp)
    
    node_labels=dict([(u, (d["nom"]))for u,d in graphe_comp.nodes(data=True)])

    edge_labels = dict([((u,v), (round(d["sim"][0],2)))for u,v,d in graphe_comp.edges(data=True)])
    
    red_edges = []

    plt.figure(figsize =(9,7))
    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.9)
    
    nx.draw_networkx_nodes(graphe_comp, pos)
    nx.draw_networkx_labels(graphe_comp, pos, labels = node_labels, font_size=8)
    nx.draw_networkx_edges(graphe_comp, pos, edge_color="black")
    nx.draw_networkx_edge_labels(graphe_comp, pos, edge_labels = edge_labels, label_pos=0.3)
    plt.axis('off')
    
    
    plt.show()
    #plt.savefig("composantes_connexes/comp_"+str(round(i,1))+"_"+str(compteur)+".png") # save as png
    plt.close()



with open("fichier_graphes_comp_grands_groupe_1.pickle", 'rb') as fichier_graphe :
        mon_depickler = pickle.Unpickler(fichier_graphe)
        dico_graphe = mon_depickler.load()
        with open("fichiers_pickle/a-minor_test2.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            tab_aminor = mon_depickler.load()
            with open("grands_graphes.pickle", 'rb') as fichier :
                mon_depickler = pickle.Unpickler(fichier)
                dico_graphes = mon_depickler.load()
                
                graphe_complet = nx.Graph()
                compteur = 0
                noeuds_deja_mis = []
                for comp in dico_graphe.keys() :
                    
                    GC = dico_graphe[comp].copy()
                    print(comp[0])
                    print(comp[1])
                    compteur_arc_arete_1, compteur_arc_arete_2 = nombre_aretes_arcs(dico_graphes[comp[0]], dico_graphes[comp[1]])
                    recalcul_sim = (((GC.number_of_nodes() + GC.number_of_edges())*(GC.number_of_nodes() + GC.number_of_edges()))/((dico_graphes[comp[0]].number_of_nodes() + compteur_arc_arete_1)*(dico_graphes[comp[1]].number_of_nodes() + compteur_arc_arete_2)))
                    sim_non_cov_sans_motif = calcul_sim_non_cov_sans_motif(dico_graphes[comp[0]], dico_graphes[comp[1]], GC)
                    
                    premier_noeud = -1
                    deuxieme_noeud = -1
                    for noeud, data in graphe_complet.nodes(data = True) :
                        if data["nom"] == comp[0] :
                            premier_noeud = noeud
                        if data["nom"] == comp[1] :
                            deuxieme_noeud = noeud
                    if premier_noeud == -1 :
                        premier_noeud = compteur
                        compteur += 1
                        graphe_complet.add_node(premier_noeud, nom=comp[0])
                    if deuxieme_noeud == -1 :
                        deuxieme_noeud = compteur
                        graphe_complet.add_node(deuxieme_noeud, nom=comp[1])
                        compteur += 1
                    
                    graphe_complet.add_edge(premier_noeud, deuxieme_noeud, sim=sim_non_cov_sans_motif)

                    
                print(graphe_complet.nodes.data())
                print(graphe_complet.edges.data())
                
                draw_composantes(graphe_complet)