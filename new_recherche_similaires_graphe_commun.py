'''
Created on 18 sept. 2019

@author: coline
'''
from recup_data.new_algo_comparaison import recup_chaines
from recup_data.graphe_commun_clusters import commun_cluster_clique
from recup_data.constantes import CLUSTERING_PEREZ_VERSION_NON_CAN_2,\
    EXTENSION_PATH, EXTENSION_PATH_TAILLE, NEW_EXTENSION_PATH_TAILLE
import pickle

def recherche_similaires(graphe_commun, graphe):
    chaines = recup_chaines(graphe)
    
    chaines_communs = recup_chaines(graphe_commun)
    
    liste_edges_ok = []
    for i in range(4) :
        for elt_commun in chaines_communs[i] :
            for elt in chaines[i] :
                if graphe_commun.nodes[elt_commun]["type"] == graphe.nodes[elt]["type"] :
                    for voisin_commun in graphe_commun[elt_commun] :
                        for edge_commun in graphe_commun[elt_commun][voisin_commun] :
                            for voisin in graphe[elt] : 
                                for edge in graphe[elt][voisin] :
                                    if graphe_commun[elt_commun][voisin_commun][edge_commun]["label"] == graphe[elt][voisin][edge]["label"] and graphe_commun.nodes[voisin_commun]["type"] == graphe.nodes[voisin]["type"] :
                                        chaine_ok = False
                                        for c in graphe_commun.nodes[voisin_commun]["chaine"] :
                                            if c in graphe.nodes[voisin]["chaine"] :
                                                chaine_ok = True
                                                
                                        if chaine_ok :
                                            if (elt_commun, voisin_commun) not in liste_edges_ok :
                                                liste_edges_ok.append((elt_commun, voisin_commun))
    
    print(liste_edges_ok)
if __name__ == '__main__':
    digraphe_commun, liste_cliques = commun_cluster_clique(CLUSTERING_PEREZ_VERSION_NON_CAN_2[11], EXTENSION_PATH%4+"dico_comp_complet_metrique_%s_taille_%s.pickle"%("toutes_aretes_coeff_all1", 4), EXTENSION_PATH_TAILLE%4)
    with open(NEW_EXTENSION_PATH_TAILLE+"fichier_1fjg_6.pickle", 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        graphe = mon_depickler.load()
    recherche_similaires(digraphe_commun, graphe)