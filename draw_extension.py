'''
Created on 10 oct. 2018

@author: Coline Gi

Dessiner les extensions par version matplotlib
'''
import pickle
import networkx as nx
import matplotlib.pyplot as plt
from networkx.classes.function import get_node_attributes
import os
from recup_data.constantes import EXTENSION_PATH_TAILLE,\
    NEW_EXTENSION_PATH_TAILLE

''' (pas utilise normalement) renvoie vrai si le graphe_commun passe en argument contient le sommet passe en argument
renvoie faux sinon 
le parametre num permet de savoir quel index de noeud il faut regarder (premier index pour le premier graphe ou deuxième index pour le deuxième graphe)'''
def contient_sommet(sommet, graphe_commun, num ):
    for noeud in graphe_commun.nodes() :
        if noeud[num] == sommet :
            return True
    return False

''' (pas utilise normalement) renvoie vrai si le graphe_commun passe en argument contient l arete passee en argument
renvoie faux sinon 
le parametre num permet de savoir quel index d'arete il faut regarder (premier index pour le premier graphe ou deuxième index pour le deuxième graphe)'''

def contient_arete(arete, graphe_commun, num):
    for edge in graphe_commun.edges() :
        if edge[0][num] == arete[0] and edge[1][num] == arete[1] or edge[1][num] == arete[0] and edge[0][num] == arete[1]   :
            return True
    return False


def draw_extension(nom_fichier):
# for j in range(1,11) :
#     with open("fichier_affichage_version3_avec_boucles_test_avec_digraph_plus_type_de_liens_avec_liaisons_b53_qui_manquent_cww_non_can_taille_%s.txt"%(j), 'w') as fichier :
# for j in range(6,11) :
    j = 4
    # for element in os.listdir(NEW_EXTENSION_PATH_TAILLE):
    #         if ".pickle" in element and "couples_possibles" not in element and  "graphe" not in element : #and element[:len(element)-7]+".png" not in os.listdir(NEW_EXTENSION_PATH_TAILLE): #and len(element.split("_")) == 5 :
    element = nom_fichier
                #with open(EXTENSION_PATH_TAILLE%j+element, 'rb') as fichier_entree :
    with open(NEW_EXTENSION_PATH_TAILLE+"%s"%(element), 'rb') as fichier_entree :
                            print(element)
                            
                            mon_depickler = pickle.Unpickler(fichier_entree)
                            G = mon_depickler.load()
                            print(G.nodes.data())
                            
                            nx.set_node_attributes(G, (33,33), "coordonnees")
                            G.nodes[1]["coordonnees"] = (0.0,0.5)
                            G.nodes[2]["coordonnees"] = (2.0,0.5)
                            G.nodes[3]["coordonnees"] = (0.0,0.0)
                            G.nodes[4]["coordonnees"] = (2.0,0.0)
                            G.nodes[5]["coordonnees"] = (3.0,0.5)
                             
            #                 fichier.write(str(element)+"\n") 
            #                 fichier.write(str(G.number_of_nodes())+"\n") 
                            print(G.edges.data())
                            
                            for noeud in G.nodes() :
                                if len(G[noeud]) == 0 :
                                    print("ramou")
                                    print(noeud)
                            
                            nodes_list = [u for u,d in G.nodes(data=True) if d["type"] != -1] 
                            print(nodes_list)
                            ordre_noeuds = [1,2,3,4,5]
                            
                            '''on ordonne les noeuds par chaine'''
                            chaines = [[1]]
                            for i in range(1,5) :
                                compteur = i
                                if i != 1 : chaines.append([i])
                                liaison_B53 = True
                                while liaison_B53 :
                                    liaison_B53 = False
                                    temp = compteur
            
                                    for voisin in G.successors(compteur) :
                                        for arc in G[compteur][voisin] :
                                            if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and G[compteur][voisin][arc]["label"] == 'B53' :
                                                liaison_B53 = True
                                                temp = voisin
                                                chaines[len(chaines)-1].append(voisin)
                                                 
                                    for voisin in G.predecessors(compteur) :
                                        for arc in G[voisin][compteur] :
                                            if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and G[voisin][compteur][arc]["label"] == 'B53' :
                                                liaison_B53 = True
                                                temp = voisin
                                                chaines[len(chaines)-1].append(voisin)
                                    compteur = temp
                            
                            for i in range(4) :
                                for elt in chaines[i] :
                                    if elt not in ordre_noeuds :
                                        ordre_noeuds.append(elt)
            #                 print(ordre_noeuds)
                            
                            #ordre_noeuds = list(nodes_list)
                            
                            ''' on place les noeuds les uns apres les autres selon la sequence 
                            en plaçant pour chaque noeud des chaines leur(s) voisin(s) non covalents a l'horizontal '''
                            for noeud in ordre_noeuds :
                                #voisins = G[noeud]
                                coordonnees_noeud = G.nodes[noeud]["coordonnees"]
                                for pred in G.predecessors(noeud) :
                                    if G.nodes[pred]["coordonnees"] == (33,33) :
                                        coordonnees = []
                                        for node in G.nodes() :
                                            coordonnees.append(G.nodes[node]["coordonnees"])
                                        for edge in G[pred][noeud] :
                                            if G[pred][noeud][edge]["label"] == "B53" :
                                                if (coordonnees_noeud[0], coordonnees_noeud[1]-0.5) not in coordonnees :
                                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.5)
                                                elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.5) not in coordonnees :
                                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.5)
                                                elif (coordonnees_noeud[0]-0.5, coordonnees_noeud[1]) not in coordonnees :
                                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.5, coordonnees_noeud[1])
                                                elif (coordonnees_noeud[0]+0.5, coordonnees_noeud[1]) not in coordonnees :
                                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.5, coordonnees_noeud[1])
                                                elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                                elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                                else :
                                                    #fichier.write("probleme\n")
                                                    print("probleme")
                                            #else :
                                            elif G.nodes[pred]["type"] != -1 :
                                                if (coordonnees_noeud[0]-0.75, coordonnees_noeud[1]) not in coordonnees :
                                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.75, coordonnees_noeud[1])
                                                elif (coordonnees_noeud[0]+0.75, coordonnees_noeud[1]) not in coordonnees :
                                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.75, coordonnees_noeud[1])
                                                elif (coordonnees_noeud[0], coordonnees_noeud[1]-0.75) not in coordonnees :
                                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.75)
                                                elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.75) not in coordonnees :
                                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.75)
                                                elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                                elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                                else : 
                                                    #fichier.write("probleme\n")
                                                    print("probleme") 
                                for succ in G.successors(noeud) :
                                    if G.nodes[succ]["coordonnees"] == (33,33) :
                                        coordonnees = []
                                        for node in G.nodes() :
                                            coordonnees.append(G.nodes[node]["coordonnees"])
                                        for edge in G[noeud][succ] :
                                            if G[noeud][succ][edge]["label"] == "B53" :
                                                if (coordonnees_noeud[0], coordonnees_noeud[1]-0.5) not in coordonnees :
                                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.5)
                                                elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.5) not in coordonnees :
                                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.5)
                                                elif (coordonnees_noeud[0]-0.5, coordonnees_noeud[1]) not in coordonnees :
                                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.5, coordonnees_noeud[1])
                                                elif (coordonnees_noeud[0]+0.5, coordonnees_noeud[1]) not in coordonnees :
                                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.5, coordonnees_noeud[1])
                                                elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                                elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                                else :
                                                    #fichier.write("probleme\n")
                                                    print("probleme")
                                            #else : 
                                            elif G.nodes[succ]["type"] != -1 :
                                                if (coordonnees_noeud[0]-0.75, coordonnees_noeud[1]) not in coordonnees :
                                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.75, coordonnees_noeud[1])
                                                elif (coordonnees_noeud[0]+0.75, coordonnees_noeud[1]) not in coordonnees :
                                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.75, coordonnees_noeud[1])
                                                elif (coordonnees_noeud[0], coordonnees_noeud[1]-0.75) not in coordonnees :
                                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.75)
                                                elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.75) not in coordonnees :
                                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.75)
                                                elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                                elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                                else : 
                                                    #fichier.write("probleme\n")
                                                    print("probleme") 
                                                
            #                 
            #                 for noeud in G.nodes() :
            #                     fichier.write(str(noeud) + " " + str(G.nodes[noeud])+"\n")
            #                 for u,v,edata in G.edges(data=True) :
            #                     fichier.write(str(u)+ "'" +str(v) + " " + str(edata["label"])+"\n")
                                
                            
                            '''pour dessiner les graphes en version CaRNAval, on a besoin de stocker les coordonnees '''

                            with open(NEW_EXTENSION_PATH_TAILLE+"%s_avec_coord.pickle"%(element[:len(element)-7]), "wb") as fichier_sortie :
                                mon_pickler = pickle.Pickler(fichier_sortie)
                                mon_pickler.dump(G)
                                        
                            pos = get_node_attributes(G, 'coordonnees')
                            
                                
                            
                            print(pos)
                            for elt in pos :
                                if pos[elt][0] > 30 or pos[elt][1] > 30 :
                                    print(elt)
                            
                            plt.figure(figsize =(5,12))
                            
                            
        #                     with open("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/result_graphes_extension_taille_max/digraphe_commun_0.7_55_taille_10.pickle", 'rb') as fichier : 
        #                         mon_depickler_2 = pickle.Unpickler(fichier)
        #                         commun  = mon_depickler_2.load()
                                       
                            red_edges = [(1,2),(2,1),(1,3),(3,1),(1,5),(5,1),(2,4),(4,2),(3,4),(4,3),(2,5),(5,2)]
                            green_edges = []
                            blue_edges = []
                            black_edges = []
                            weights = []
                            
                            edges_list = [(u,v) for u,v,data in G.edges(data=True) if data["near"] != None]
                            
                            #black_edges = [edge for edge in G.edges() if edge not in red_edges]
                            for u,v,edata in G.edges(data=True) :
                                if (u,v) in edges_list : 
                                    if u in nodes_list and v in nodes_list :
                                        if (u,v) not in red_edges :
                                            if edata["label"] == "B53" :
                                                green_edges.append((u,v))
                                            elif edata["label"] == "CWW" :
                                                blue_edges.append((u,v)) 
                                            else :
                                                black_edges.append((u,v))
                            
        #                     for (u,v) in edges_list :                   
        #                         if contient_arete((u,v), commun, 1) :
        #                             weights.append(3)
        #                         else :
        #                             weights.append(1)
            
                            edge_labels=dict([((u,v,),d["label"])for u,v,d in G.edges(data=True)])
                            #print(edge_labels)
                           # node_labels=dict([(u,(d["nt"], d["type"]))for u,d in G.nodes(data=True)])## if d["type"] != None])
                            
                            node_labels=dict([(u, (d["type"], d["poids"])) for u,d in G.nodes(data=True) if d["type"] != -1 and d["type"] != None]) #else (u, (u)) for u,d in G.nodes(data=True) ])
                            #node_labels=dict([(u, (u)) for u,d in G.nodes(data=True) ])
                            print(node_labels)
                            
        #                     orange_nodes = []
        #                     pink_nodes  = []
        #                     for noeud in nodes_list :
        #                         if contient_sommet(noeud, commun, 1) :
        #                             orange_nodes.append(noeud)
        #                         else :
        #                             pink_nodes.append(noeud)
        #                     node_colors = ['pink' if node in pink_nodes else 'orange' for node in nodes_list]
                            #edge_labels = dict([((u,v,), (d["label"])) for u,v,d in G.edges(data=True) if (u,v) in edges_list]) #else (u, (u)) for u,d in G.nodes(data=True) ])
         #dict([((u,v,),d["label"])for u,v,d in G.edges(data=True) if d["label"] != 'B53' and d["label"] != 'CWW' and ((u,v) not in courbes and (v,u) not in courbes)])
           
                            
                            nx.draw_networkx_nodes(G, pos, nodelist=nodes_list, node_size=150, node_color="pink")
                            nx.draw_networkx_labels(G, pos, labels = node_labels, fontsize=8)
    
                            #edges_list = [(u,v) for u,v,data in G.edges(data=True) if data["long_range"] != None]
                            edge_colors = ['black' if edge in black_edges else 'red' if edge in red_edges else 'blue' if edge in blue_edges else 'green' for edge in edges_list]
                            print(black_edges)
                            print(red_edges)
                            print(edge_colors)
                            #edge_labels = dict([((u,v), (d["long_range"])) for u,v,d in G.edges(data=True)])
                            #nx.draw_networkx_edge_labels(G,pos, edge_labels=edge_labels)
                            print(len(weights))
        #                     print(len(edges_list))
        #                     print(commun.edges.data())
                            #nx.draw_networkx_edge_labels(G,pos, edge_labels = edge_labels, font_size=8)
                            nx.draw_networkx_edges(G,pos, edgelist = edges_list, edge_color=edge_colors)#, width=weights)
                            plt.axis('off')
                            plt.savefig(NEW_EXTENSION_PATH_TAILLE+ "%s.png"%(element[:len(element)-7]), format="png") #transparent=True, dpi=50) # save as png
                            #plt.savefig("Extensions/Metrique_toutes_aretes/graphes_extension_autres_tailles/motif_taille_4_0.8_6/graphe_motif_moyen.png") # save as png
                            
                            #plt.savefig("graphes_extension/fichier_1FJG_A_48_8.svg", format='svg') # save as png
                            plt.show()
                            plt.clf()
                            plt.close()

                
if __name__ == '__main__':
    draw_extension("fichier_4y4o_25_2.pickle")
    